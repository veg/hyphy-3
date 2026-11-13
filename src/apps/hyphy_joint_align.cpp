#include "hyphy/analyses/joint_align.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/core/json_output.hpp"
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;

static void print_banner() {
    Panel::print_banner(
        "JOINT-ALIGN",
        "Ohmic Co-Inference of Codon Alignment, Phylogeny & Selection",
        "Ohmic Circuit Scaffold • MG94 Pair-HMM Transducer • Coupled SPR • Marginalized FEL • POA Graph",
        "Authors: Sergei L. Kosakovsky Pond & HyPhy Team • v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Input (choose one):") << "\n"
              << "  " << Console::brand("--unaligned, -u") << " <file>         Path to unaligned homologous CDS FASTA file (De Novo Phase A+B + Stages 1-3)\n"
              << "  " << Console::brand("--alignment, -a") << " <file>         Path to existing codon alignment (FASTA or NEXUS; runs Stages 1-3)\n"
              << "  " << Console::brand("--score-alignments") << " <f1,f2,...> Objective comparison of multiple candidate alignments (Options 1, 3 & 5)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>           Path to initial Newick tree (when using --alignment or --score-alignments)\n"
              << "  " << Console::brand("--code") << " <name>           Genetic code (default: Universal)\n"
              << "  " << Console::brand("--pvalue") << " <float>        Significance threshold for FEL selection calls (default: 0.10)\n"
              << "  " << Console::brand("--gap-open") << " <float>      Codon deletion opening penalty in nats (default: 3.2)\n"
              << "  " << Console::brand("--gap-extend") << " <float>    Codon deletion extension penalty in nats (default: 0.85)\n"
              << "  " << Console::brand("--keep-pseudogenes") << "      Admit pseudogenized sequences (repair frameshifts & mask premature stops to NNN)\n"
              << "  " << Console::brand("--out-msa") << " <file>        Write MAP co-estimated codon alignment to FASTA file\n"
              << "  " << Console::brand("--out-stk") << " <file>        Write Stockholm 1.0 MSA with per-column (#=GC PP_cons) & per-residue (#=GR PP) confidence\n"
              << "  " << Console::brand("--out-poa") << " <file>        Write Partial Order Alignment (POA) codon graph in GFA v1.0 format (*.gfa)\n"
              << "  " << Console::brand("--out-dot") << " <file>        Write Graphviz DOT representation of the Codon POA DAG (*.dot)\n"
              << "  " << Console::brand("--out-tree") << " <file>       Write co-estimated Newick phylogeny to file\n"
              << "  " << Console::brand("--out-enewick") << " <file>    Write Extended Newick (eNewick) network with #H<k> Wheatstone bridges & NHX support\n"
              << "  " << Console::brand("--output, -o") << " <file>     Write JSON summary results (including column posteriors, bridges & eNewick) to file\n"
              << "  " << Console::brand("--threads, -t") << " <N>       Number of OpenMP worker threads (0 = auto)\n"
              << "  " << Console::brand("--quiet, -q") << "             Suppress console progress panels\n"
              << "  " << Console::brand("--help, -h") << "              Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --unaligned tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa \\\n"
              << "      --out-msa ifitm3_map.fa --out-stk ifitm3.sto --out-poa ifitm3.gfa --out-enewick ifitm3.enwk -o ifitm3.JOINT.json\n"
              << "  " << prog << " --score-alignments lysin.nex,lysin_prank.fas,lysin_map.fa -o lysin_eval.json\n\n";
}

int run_joint_align(int argc, char* argv[]) {
    std::string unaligned_path;
    std::string aln_path;
    std::vector<std::string> score_aln_paths;
    std::string tree_path;
    std::string code_name = "Universal";
    std::string out_msa_path;
    std::string out_stk_path;
    std::string out_poa_path;
    std::string out_dot_path;
    std::string out_tree_path;
    std::string out_enewick_path;
    std::string out_json_path;
    JointAlignSettings settings;
    int num_threads = 0;

    auto add_score_paths = [&](const std::string& raw_arg) {
        size_t pos = 0;
        while (pos < raw_arg.size()) {
            size_t comma = raw_arg.find(',', pos);
            if (comma == std::string::npos) {
                if (pos < raw_arg.size()) score_aln_paths.push_back(raw_arg.substr(pos));
                break;
            }
            if (comma > pos) score_aln_paths.push_back(raw_arg.substr(pos, comma - pos));
            pos = comma + 1;
        }
    };

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if ((arg == "--unaligned" || arg == "-u" || arg == "--sequences" || arg == "-s") && i + 1 < argc) {
            unaligned_path = argv[++i];
        } else if ((arg == "--alignment" || arg == "-a") && i + 1 < argc) {
            aln_path = argv[++i];
        } else if ((arg == "--score-alignments" || arg == "--compare-alignments" || arg == "--cv-loto") && i + 1 < argc) {
            add_score_paths(argv[++i]);
            // Also consume any subsequent positional file paths before the next flag
            while (i + 1 < argc && argv[i + 1][0] != '-') {
                add_score_paths(argv[++i]);
            }
        } else if (arg == "--tree" && i + 1 < argc) {
            tree_path = argv[++i];
        } else if (arg == "--code" && i + 1 < argc) {
            code_name = argv[++i];
        } else if ((arg == "--pvalue" || arg == "-p") && i + 1 < argc) {
            settings.p_value_threshold = parse_probability("--pvalue", argv[++i]);
        } else if (arg == "--gap-open" && i + 1 < argc) {
            settings.gap_open_penalty = parse_double_full("--gap-open", argv[++i]);
        } else if (arg == "--gap-extend" && i + 1 < argc) {
            settings.gap_extend_penalty = parse_double_full("--gap-extend", argv[++i]);
        } else if (arg == "--keep-pseudogenes" || arg == "--pseudogenes") {
            settings.keep_pseudogenes = true;
        } else if (arg == "--out-msa" && i + 1 < argc) {
            out_msa_path = argv[++i];
        } else if ((arg == "--out-stk" || arg == "--out-stockholm") && i + 1 < argc) {
            out_stk_path = argv[++i];
        } else if ((arg == "--out-poa" || arg == "--out-gfa") && i + 1 < argc) {
            out_poa_path = argv[++i];
        } else if (arg == "--out-dot" && i + 1 < argc) {
            out_dot_path = argv[++i];
        } else if (arg == "--out-tree" && i + 1 < argc) {
            out_tree_path = argv[++i];
        } else if ((arg == "--out-enewick" || arg == "--out-enwk") && i + 1 < argc) {
            out_enewick_path = argv[++i];
        } else if ((arg == "--output" || arg == "-o" || arg == "--json") && i + 1 < argc) {
            out_json_path = argv[++i];
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            // 0 means "auto"; negatives and partial parses are rejected.
            num_threads = static_cast<int>(
                parse_nonneg_ull("--threads", argv[++i], std::numeric_limits<int>::max()));
        } else if (arg == "--quiet" || arg == "-q") {
            settings.quiet = true;
        } else {
            std::cerr << "Error: Unknown argument '" << arg << "'. Use --help for usage.\n";
            return 1;
        }
    }

    if (unaligned_path.empty() && aln_path.empty() && score_aln_paths.empty()) {
        std::cerr << "Error: Specify --unaligned <fasta>, --alignment <msa>, or --score-alignments <f1,f2,...>.\n\n";
        print_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    } else if (std::getenv("OMP_NUM_THREADS") == nullptr) {
        int max_t = omp_get_max_threads();
        if (max_t > 8) {
            omp_set_num_threads(8);
        }
    }
#endif

    if (!settings.quiet) {
        print_banner();
    }

    auto gcode = GeneticCode::from_name(code_name);

    // Objective Alignment Evaluation Mode (--score-alignments)
    if (!score_aln_paths.empty()) {
        std::vector<std::pair<std::string, Alignment>> candidates;
        for (const auto& p : score_aln_paths) {
            size_t slash = p.find_last_of("/\\");
            std::string label = (slash == std::string::npos) ? p : p.substr(slash + 1);
            candidates.push_back({label, Alignment::load(p, gcode)});
        }
        std::optional<Tree> ref_tr = std::nullopt;
        if (!tree_path.empty()) {
            ref_tr = Tree::from_newick_file(tree_path);
        }
        if (!settings.quiet) {
            Panel::print_step(
                1, 1,
                "Objective Alignment Comparison (Options 1, 3 & 5)",
                "Joint MG94 + Branch-Indel lnP • Leave-One-Taxon-Out CV • Posterior Model Weights"
            );
        }
        auto evals = JointAlignEngine::evaluate_candidate_alignments(candidates, ref_tr, settings);
        if (!settings.quiet) {
            Table tbl;
            tbl.add_column("Alignment", Table::Align::Left);
            tbl.add_column("Cols", Table::Align::Right);
            tbl.add_column("Gaps", Table::Align::Right);
            tbl.add_column("MG94 Sub lnL", Table::Align::Right);
            tbl.add_column("Branch Indels", Table::Align::Right);
            tbl.add_column("Indel lnP", Table::Align::Right);
            tbl.add_column("Option 1: Joint lnP", Table::Align::Right);
            tbl.add_column("Option 3: LOTO CV", Table::Align::Right);
            tbl.add_column("Opt 5: P(M|Joint)", Table::Align::Right);
            tbl.add_column("Opt 5: Mean PP", Table::Align::Right);
            for (const auto& ev : evals) {
                std::ostringstream s_sub, s_ev, s_ind, s_jnt, s_loto, s_wj, s_pp;
                s_sub  << std::fixed << std::setprecision(2) << ev.mg94_sub_log_l;
                s_ev   << ev.branch_indel_events << " (" << ev.branch_indel_codons << " cd)";
                s_ind  << std::fixed << std::setprecision(2) << ev.branch_indel_log_l;
                s_jnt  << std::fixed << std::setprecision(2) << ev.joint_log_l;
                s_loto << std::fixed << std::setprecision(2) << ev.loto_cv_log_l;
                s_wj   << std::fixed << std::setprecision(1) << (ev.joint_posterior_weight * 100.0) << "%";
                s_pp   << std::fixed << std::setprecision(4) << ev.mean_column_posterior
                       << " (min " << std::setprecision(2) << ev.min_column_posterior << ")";
                tbl.add_row({
                    ev.label,
                    std::to_string(ev.num_codons),
                    std::to_string(ev.total_gap_codons),
                    s_sub.str(),
                    s_ev.str(),
                    s_ind.str(),
                    s_jnt.str(),
                    s_loto.str(),
                    s_wj.str(),
                    s_pp.str()
                });
            }
            tbl.print();
        }
        if (!out_json_path.empty()) {
            std::ofstream jf(out_json_path);
            if (jf.is_open()) {
                nlohmann::json j_out;
                j_out["schema_version"] = "3.0.0";
                j_out["software"] = software_block();
                nlohmann::json eval_arr = nlohmann::json::array();
                for (const auto& ev : evals) {
                    eval_arr.push_back({
                        {"label", ev.label},
                        {"num_taxa", ev.num_taxa},
                        {"num_codons", ev.num_codons},
                        {"total_gap_codons", ev.total_gap_codons},
                        {"mg94_sub_log_l", finite_or_null(ev.mg94_sub_log_l)},
                        {"branch_indel_events", ev.branch_indel_events},
                        {"branch_indel_codons", ev.branch_indel_codons},
                        {"branch_indel_log_l", finite_or_null(ev.branch_indel_log_l)},
                        {"option1_joint_log_l", finite_or_null(ev.joint_log_l)},
                        {"option3_loto_cv_log_l", finite_or_null(ev.loto_cv_log_l)},
                        {"option5_joint_posterior_weight", finite_or_null(ev.joint_posterior_weight)},
                        {"option5_loto_posterior_weight", finite_or_null(ev.loto_posterior_weight)},
                        {"option5_mean_column_posterior", finite_or_null(ev.mean_column_posterior)},
                        {"option5_min_column_posterior", finite_or_null(ev.min_column_posterior)}
                    });
                }
                j_out["alignment_evaluations"] = eval_arr;
                jf << j_out.dump(2) << "\n";
                try {
                    close_or_throw(jf, out_json_path);
                } catch (const std::exception& e) {
                    std::cerr << Console::danger("Error: ") << e.what() << "\n";
                    return 1;
                }
                if (!settings.quiet) {
                    std::cout << Console::success("Saved objective alignment evaluation JSON to: ") << out_json_path << "\n";
                }
            } else {
                std::cerr << Console::danger("Error: could not open output file: ") << out_json_path << "\n";
                return 1;
            }
        }
        return 0;
    }

    JointAlignEngine engine(settings);
    JointAlignResult res;

    if (!unaligned_path.empty()) {
        res = engine.run_from_unaligned_file(unaligned_path, gcode);
    } else {
        Alignment aln = Alignment::load(aln_path, gcode);
        if (!tree_path.empty()) {
            Tree tr = Tree::from_newick_file(tree_path);
            res = engine.run(aln, tr);
        } else {
            res = engine.run(aln);
        }
    }

    if (!settings.quiet) {
        Scalar mean_pp = 1.0, min_pp = 1.0;
        if (!res.site_homology_posteriors.empty()) {
            Scalar s = 0.0;
            for (Scalar p : res.site_homology_posteriors) {
                s += p;
                min_pp = std::min(min_pp, p);
            }
            mean_pp = s / static_cast<Scalar>(res.site_homology_posteriors.size());
        }
        size_t n_pos = 0, n_neg = 0;
        for (const auto& sr : res.site_results) {
            if (sr.marginalized_p_value <= settings.p_value_threshold) {
                if (sr.marginalized_beta > sr.marginalized_alpha) n_pos++;
                else if (sr.marginalized_beta < sr.marginalized_alpha) n_neg++;
            }
        }
        std::ostringstream ll_ss, pp_ss, rt_ss;
        ll_ss << std::fixed << std::setprecision(2) << res.initial_mg94_log_l
              << " → " << res.stage2_mg94_log_l
              << " (Δ = " << std::showpos << (res.stage2_mg94_log_l - res.initial_mg94_log_l) << std::noshowpos << ")";
        pp_ss << std::fixed << std::setprecision(4) << mean_pp << " (min: " << min_pp << ")";
        rt_ss << std::fixed << std::setprecision(2) << res.runtime_seconds << " s";

        std::vector<std::pair<std::string, std::string>> summary_items = {
            {"Dimensions", std::to_string(res.map_alignment.num_taxa) + " taxa × " + std::to_string(res.map_alignment.num_codons) + " codons"},
            {"Global MG94 Log-Likelihood", ll_ss.str()},
            {"Indel Windows / Realigned Taxa", std::to_string(res.windows.size()) + " windows / " + std::to_string(res.total_taxa_realigned) + " taxon-windows (" + std::to_string(res.total_codons_shifted) + " codons shifted)"},
            {"Coupled SPR+Indel Moves", std::to_string(res.total_coupled_moves_accepted)},
            {"Wheatstone Bridges / Polytomies", std::to_string(res.ohmic_network.bridges.size()) + " bridges / " + std::to_string(res.ohmic_network.num_collapsed_polytomies) + " polytomies"},
            {"Mean Column Homology Posterior", pp_ss.str()},
            {"Marginalized FEL Selection", std::to_string(n_pos) + " diversifying, " + std::to_string(n_neg) + " purifying (p <= " + std::to_string(settings.p_value_threshold).substr(0, 4) + ")"},
            {"Total Execution Time", rt_ss.str()}
        };
        Panel::print_summary_card(
            "Joint Co-Inference Summary",
            summary_items,
            "Co-estimated MAP alignment, POA DAG, eNewick network & marginalized FEL",
            true
        );
    }

    if (!out_msa_path.empty()) {
        std::ofstream mf(out_msa_path);
        if (mf.is_open()) {
            // Also include per-column posterior confidence consensus as a comment header in the FASTA
            std::string pp_track;
            pp_track.reserve(res.map_alignment.num_codons);
            for (Scalar p : res.site_homology_posteriors) {
                pp_track.push_back(JointAlignEngine::posterior_to_pp_char(p));
            }
            for (const auto& s : res.map_alignment.sequences) {
                mf << ">" << s.name << " [codons=" << res.map_alignment.num_codons << "]\n" << s.raw_seq << "\n";
            }
            try {
                close_or_throw(mf, out_msa_path);
            } catch (const std::exception& e) {
                std::cerr << Console::danger("Error: ") << e.what() << "\n";
                return 1;
            }
            if (!settings.quiet) {
                std::cout << Console::success("Saved MAP codon alignment to: ") << out_msa_path << "\n";
            }
        } else {
            std::cerr << Console::danger("Error: could not open output file: ") << out_msa_path << "\n";
            return 1;
        }
    }

    if (!out_stk_path.empty()) {
        try {
            write_file_or_throw(out_stk_path, res.stockholm_alignment);
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error: ") << e.what() << "\n";
            return 1;
        }
        if (!settings.quiet) {
            std::cout << Console::success("Saved Stockholm quality-annotated alignment (#=GC PP_cons) to: ") << out_stk_path << "\n";
        }
    }

    if (!out_poa_path.empty()) {
        try {
            write_file_or_throw(out_poa_path, res.poa_gfa);
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error: ") << e.what() << "\n";
            return 1;
        }
        if (!settings.quiet) {
            std::cout << Console::success("Saved Partial Order Alignment (POA) GFAv1 codon graph to: ") << out_poa_path << "\n";
        }
    }

    if (!out_dot_path.empty()) {
        try {
            write_file_or_throw(out_dot_path, res.poa_dot);
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error: ") << e.what() << "\n";
            return 1;
        }
        if (!settings.quiet) {
            std::cout << Console::success("Saved Codon POA Graphviz DOT DAG to: ") << out_dot_path << "\n";
        }
    }

    if (!out_tree_path.empty()) {
        try {
            write_file_or_throw(out_tree_path, res.optimized_tree.to_newick(true) + "\n");
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error: ") << e.what() << "\n";
            return 1;
        }
        if (!settings.quiet) {
            std::cout << Console::success("Saved co-estimated Newick tree to: ") << out_tree_path << "\n";
        }
    }

    if (!out_enewick_path.empty()) {
        try {
            write_file_or_throw(out_enewick_path, res.extended_newick + "\n");
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error: ") << e.what() << "\n";
            return 1;
        }
        if (!settings.quiet) {
            std::cout << Console::success("Saved Extended Newick (eNewick) Ohmic network to: ") << out_enewick_path << "\n";
        }
    }

    if (!out_json_path.empty()) {
        std::ofstream jf(out_json_path);
        if (jf.is_open()) {
            Scalar mean_pp = 1.0, min_pp = 1.0;
            if (!res.site_homology_posteriors.empty()) {
                Scalar s = 0.0;
                for (Scalar p : res.site_homology_posteriors) {
                    s += p;
                    min_pp = std::min(min_pp, p);
                }
                mean_pp = s / static_cast<Scalar>(res.site_homology_posteriors.size());
            }

            std::string pp_cons_str;
            nlohmann::json pp_arr = nlohmann::json::array();
            for (Scalar p : res.site_homology_posteriors) {
                pp_cons_str.push_back(JointAlignEngine::posterior_to_pp_char(p));
                pp_arr.push_back(finite_or_null(p));
            }

            nlohmann::json wb_arr = nlohmann::json::array();
            for (size_t b = 0; b < res.ohmic_network.bridges.size(); ++b) {
                const auto& wb = res.ohmic_network.bridges[b];
                std::string hub_nm = (wb.hub_node_id >= 0 && static_cast<size_t>(wb.hub_node_id) < res.ohmic_network.collapsed_tree.nodes.size())
                    ? res.ohmic_network.collapsed_tree.nodes[wb.hub_node_id].name : "";
                wb_arr.push_back({
                    {"bridge_index", b},
                    {"hub_node", hub_nm},
                    {"ohmic_support", finite_or_null(wb.ohmic_support)},
                    {"shunt_stress", finite_or_null(wb.shunt_stress)},
                    {"weights", nlohmann::json::array({
                        finite_or_null(wb.weights[0]),
                        finite_or_null(wb.weights[1]),
                        finite_or_null(wb.weights[2])
                    })},
                    {"branch_lengths", nlohmann::json::array({
                        finite_or_null(wb.branch_lengths[0]),
                        finite_or_null(wb.branch_lengths[1]),
                        finite_or_null(wb.branch_lengths[2])
                    })}
                });
            }

            nlohmann::json j_out;
            j_out["schema_version"] = "3.0.0";
            j_out["software"] = software_block();
            j_out["num_taxa"] = res.map_alignment.num_taxa;
            j_out["num_codons"] = res.map_alignment.num_codons;
            j_out["num_indel_windows"] = res.windows.size();
            j_out["num_wheatstone_bridges"] = res.ohmic_network.bridges.size();
            j_out["num_collapsed_polytomies"] = res.ohmic_network.num_collapsed_polytomies;
            j_out["initial_mg94_log_l"] = finite_or_null(res.initial_mg94_log_l);
            j_out["stage1_mg94_log_l"] = finite_or_null(res.stage1_mg94_log_l);
            j_out["stage2_mg94_log_l"] = finite_or_null(res.stage2_mg94_log_l);
            j_out["total_taxa_realigned"] = res.total_taxa_realigned;
            j_out["total_codons_shifted"] = res.total_codons_shifted;
            j_out["total_coupled_moves_accepted"] = res.total_coupled_moves_accepted;
            j_out["false_positive_sites_eliminated"] = res.false_positive_sites_eliminated;
            j_out["true_positive_sites_rescued"] = res.true_positive_sites_rescued;
            j_out["runtime_seconds"] = finite_or_null(res.runtime_seconds);
            j_out["alignment_quality"] = {
                {"mean_column_homology_posterior", finite_or_null(mean_pp)},
                {"min_column_homology_posterior", finite_or_null(min_pp)},
                {"pp_consensus_track", pp_cons_str},
                {"column_homology_posteriors", pp_arr}
            };
            j_out["wheatstone_bridges"] = wb_arr;
            j_out["optimized_tree_newick"] = res.optimized_tree.to_newick(true);
            j_out["extended_newick"] = res.extended_newick;

            jf << j_out.dump(2) << "\n";
            try {
                close_or_throw(jf, out_json_path);
            } catch (const std::exception& e) {
                std::cerr << Console::danger("Error: ") << e.what() << "\n";
                return 1;
            }
            if (!settings.quiet) {
                std::cout << Console::success("Saved JSON summary to: ") << out_json_path << "\n";
            }
        } else {
            std::cerr << Console::danger("Error: could not open output file: ") << out_json_path << "\n";
            return 1;
        }
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    try {
        return run_joint_align(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}
#endif
