#include "hyphy/analyses/characterize.hpp"
#include "hyphy/core/console.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;

static void print_characterize_banner() {
    Panel::print_banner(
        "CHARACTERIZE (ONE-BUTTON)",
        "End-to-End Unaligned/Aligned Evolutionary & Selection Characterization",
        "Co-Estimated MSA + POA • Ohmic Network • Compositional Audit • BUSTED-E/aBSREL/RELAX/MEME/dFEL • Model-Free Exchangeability & 1D Domains",
        "Authors: Sergei L. Kosakovsky Pond & HyPhy Team • v3.0.0"
    );
}

static void print_characterize_usage(const char* prog) {
    print_characterize_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Input (choose one):") << "\n"
              << "  " << Console::brand("--unaligned, -u") << " <file>       Path to raw unaligned homologous CDS FASTA file\n"
              << "  " << Console::brand("--alignment, -a") << " <file>       Path to existing codon alignment (FASTA or NEXUS)\n"
              << "  " << Console::brand("--tree, -t") << " <file>            Optional Newick tree file\n\n"
              << Console::bold("Output:") << "\n"
              << "  " << Console::brand("--output, -o") << " <file>          Path for unified JSON characterization report\n"
              << "  " << Console::brand("--html") << " <file>                Path for standalone interactive multi-tab SPA HTML report\n"
              << "  " << Console::brand("--emit-legacy-jsons") << "          Also emit legacy .FEL.json, .MEME.json & .BUSTED.json files\n"
              << "  " << Console::brand("--out-msa") << " <file>             Path to write co-estimated MAP codon alignment (FASTA)\n"
              << "  " << Console::brand("--out-stockholm") << " <file>       Path to write Stockholm 1.0 quality-annotated alignment\n"
              << "  " << Console::brand("--out-poa") << " <file>             Path to write Partial Order Alignment (POA) graph in GFAv1\n\n"
              << Console::bold("Parameters:") << "\n"
              << "  " << Console::brand("--code") << " <name>                Genetic code (default: Universal)\n"
              << "  " << Console::brand("--pvalue, -p") << " <float>         Per-test significance threshold (default: 0.10)\n"
              << "  " << Console::brand("--fdr, -q") << " <float>            Benjamini-Hochberg / Simes FDR q-value threshold (default: 0.10)\n"
              << "  " << Console::brand("--branches, -b") << " <tag>         Optional foreground branch tag for RELAX/Contrast-FEL (default: BIC-gated auto)\n"
              << "  " << Console::brand("--coevolution") << "                Opt-in Tier 3 screen for epistatic co-substituting site pairs\n"
              << "  " << Console::brand("--no-spr") << "                     Disable coupled SPR + indel search\n"
              << "  " << Console::brand("--threads") << " <int>              Number of OpenMP threads\n"
              << "  " << Console::brand("--quiet") << "                      Suppress progress output\n";
}

#ifdef HYPHY3_COMBINED_DRIVER
int run_characterize(int argc, char* argv[])
#else
int main(int argc, char* argv[])
#endif
{
    std::string unaligned_path;
    std::string aln_path;
    std::string tree_path;
    std::string out_json_path;
    std::string out_html_path;
    std::string out_msa_path;
    std::string out_sto_path;
    std::string out_poa_path;
    std::string code_name = "Universal";
    CharacterizeSettings settings;
    int num_threads = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_characterize_usage(argv[0]);
            return 0;
        } else if ((arg == "--unaligned" || arg == "-u") && i + 1 < argc) {
            unaligned_path = argv[++i];
        } else if ((arg == "--alignment" || arg == "-a" || arg == "--msa") && i + 1 < argc) {
            aln_path = argv[++i];
        } else if ((arg == "--tree" || arg == "-t") && i + 1 < argc) {
            tree_path = argv[++i];
        } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
            out_json_path = argv[++i];
        } else if (arg == "--html" && i + 1 < argc) {
            out_html_path = argv[++i];
        } else if (arg == "--emit-legacy-jsons") {
            settings.emit_legacy_jsons = true;
        } else if (arg == "--out-msa" && i + 1 < argc) {
            out_msa_path = argv[++i];
        } else if (arg == "--out-stockholm" && i + 1 < argc) {
            out_sto_path = argv[++i];
        } else if (arg == "--out-poa" && i + 1 < argc) {
            out_poa_path = argv[++i];
        } else if (arg == "--code" && i + 1 < argc) {
            code_name = argv[++i];
        } else if ((arg == "--pvalue" || arg == "-p") && i + 1 < argc) {
            settings.p_value_threshold = std::stod(argv[++i]);
        } else if (arg == "--fdr" && i + 1 < argc) {
            settings.fdr_threshold = std::stod(argv[++i]);
        } else if ((arg == "--branches" || arg == "-b") && i + 1 < argc) {
            settings.foreground_tag = argv[++i];
        } else if (arg == "--coevolution") {
            settings.enable_coevolution = true;
        } else if (arg == "--no-spr") {
            settings.enable_coupled_spr = false;
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--quiet" || arg == "-q") {
            settings.quiet = true;
        } else {
            std::cerr << Console::danger("Error: ") << "Unknown argument '" << arg << "'. Use --help for usage.\n";
            return 1;
        }
    }

    if (unaligned_path.empty() && aln_path.empty()) {
        print_characterize_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
#else
    (void)num_threads;
#endif

    std::string input_label = !unaligned_path.empty() ? unaligned_path : aln_path;
    size_t slash = input_label.find_last_of("/\\");
    std::string short_label = (slash == std::string::npos) ? input_label : input_label.substr(slash + 1);

    if (out_json_path.empty()) {
        out_json_path = input_label + ".CHARACTERIZE.json";
    }
    if (out_html_path.empty()) {
        out_html_path = input_label + ".CHARACTERIZE.html";
    }

    if (!settings.quiet) {
        print_characterize_banner();
        std::vector<std::pair<std::string, std::string>> cfg = {
            {"Input Source", !unaligned_path.empty() ? (unaligned_path + " (Raw Unaligned CDS)") : (aln_path + " (Aligned MSA)")},
            {"Genetic Code", code_name},
            {"Thresholds", "p <= " + std::to_string(settings.p_value_threshold).substr(0, 4) +
                           "  •  Simes/BH FDR q <= " + std::to_string(settings.fdr_threshold).substr(0, 4)},
            {"JSON / SPA HTML Output", out_json_path + "  •  " + out_html_path}
        };
        Panel::print_card("Unified One-Button Configuration", cfg);
    }

    auto gcode = GeneticCode::from_name(code_name);
    CharacterizeEngine engine(settings);
    CharacterizeResult res;

    if (!unaligned_path.empty()) {
        res = engine.run_from_unaligned_file(unaligned_path, gcode);
    } else {
        Alignment aln = Alignment::load(aln_path, gcode);
        if (!tree_path.empty()) {
            Tree tr = Tree::from_newick_file(tree_path);
            res = engine.run_from_alignment(aln, tr);
        } else {
            res = engine.run_from_alignment(aln);
        }
    }

    if (!settings.quiet) {
        // 1. Top Site Selection Atlas Table
        Table site_tbl;
        site_tbl.add_column("Site", Table::Align::Right);
        site_tbl.add_column("Domain", Table::Align::Left);
        site_tbl.add_column("Classification", Table::Align::Left);
        site_tbl.add_column("Simes p (q)", Table::Align::Right);
        site_tbl.add_column("PP (K)", Table::Align::Right);
        site_tbl.add_column("α", Table::Align::Right);
        site_tbl.add_column("β_perv", Table::Align::Right);
        site_tbl.add_column("β+ (Epi)", Table::Align::Right);
        site_tbl.add_column("MEME p", Table::Align::Right);
        site_tbl.add_column("Driving Branches", Table::Align::Left);
        site_tbl.add_column("Exchangeability / Toggle / Convergence", Table::Align::Left);

        std::vector<UnifiedSiteAtlasEntry> ranked_sites = res.site_atlas;
        std::sort(ranked_sites.begin(), ranked_sites.end(), [](const auto& a, const auto& b) {
            bool a_pos = (a.classification.find("Diversifying") != std::string::npos);
            bool b_pos = (b.classification.find("Diversifying") != std::string::npos);
            if (a_pos != b_pos) return a_pos > b_pos;
            return a.simes_omnibus_p < b.simes_omnibus_p;
        });

        size_t shown = 0;
        for (const auto& s : ranked_sites) {
            if (s.classification == "Neutral" || s.classification == "Invariant") continue;
            if (shown >= 18) break;
            shown++;

            std::ostringstream s_simes, s_pp, s_a, s_b, s_bp, s_mp, s_br;
            s_simes << std::fixed << std::setprecision(4) << s.simes_omnibus_p
                    << " (" << std::setprecision(2) << s.simes_omnibus_q << ")";
            s_pp << std::fixed << std::setprecision(3) << s.column_homology_posterior
                 << " (" << s.ensemble_paths_evaluated << ")";
            s_a  << std::fixed << std::setprecision(2) << s.alpha;
            s_b  << std::fixed << std::setprecision(2) << s.beta_pervasive;
            s_bp << std::fixed << std::setprecision(2) << s.beta_plus;
            s_mp << std::fixed << std::setprecision(4) << s.meme_p_value;

            if (s.top_episodic_branches.empty()) {
                s_br << "-";
            } else {
                for (size_t k = 0; k < s.top_episodic_branches.size(); ++k) {
                    if (k > 0) s_br << ",";
                    s_br << s.top_episodic_branches[k];
                }
            }
            std::string exch_summary = s.key_property_driver.empty() ? s.exchangeability_mode : s.key_property_driver;
            if (s.cpg_hypermutable_context) {
                exch_summary += " [CpG]";
            } else if (s.apobec_motif_context) {
                exch_summary += " [APOBEC]";
            }
            site_tbl.add_row({
                std::to_string(s.site_index + 1),
                s.spatial_domain,
                s.classification,
                s_simes.str(),
                s_pp.str(),
                s_a.str(),
                s_b.str(),
                s_bp.str(),
                s_mp.str(),
                s_br.str(),
                exch_summary
            });
        }
        if (site_tbl.num_rows() > 0) {
            std::cout << Console::bold("Top Adaptive & Constrained Sites (Unified Site-Branch-Exchangeability Atlas):") << "\n";
            site_tbl.print();
            std::cout << "\n";
        }

        // 2. Three-Question Executive Summary Card
        std::ostringstream s_q1, s_q2, s_q3_gene, s_q3_dom, s_q3_sites, s_q3_br, s_rt;
        s_q1 << res.joint_align.map_alignment.num_taxa << " taxa × "
             << res.joint_align.map_alignment.num_codons << " codons  │  Mean PP = "
             << std::fixed << std::setprecision(4) << res.alignment_objective.mean_column_posterior
             << "  │  Joint lnP = " << std::setprecision(2) << res.alignment_objective.joint_log_l
             << "  │  LOTO CV = " << res.alignment_objective.loto_cv_log_l;

        s_q2 << (res.composition_test.is_heterogeneous ? "Heterogeneous GC3 Drift (3-Regime Softmax Active)" : "Stationary")
             << " (CP3 p = " << std::fixed << std::setprecision(3) << res.composition_test.p_value_cp3 << ")"
             << "  │  Wheatstone Bridges: " << res.joint_align.ohmic_network.bridges.size()
             << " (" << res.recombination_segments.size() << " localized)";
        if (res.mutational_context_alert.triggered) {
            s_q2 << "  │  [!] " << res.mutational_context_alert.alert_summary;
        }

        s_q3_gene << (res.gene_wide.has_gene_wide_positive_selection ? "POSITIVE SELECTION DETECTED" : "No Gene-Wide Signal")
                  << " (BUSTED-E LRT = " << std::fixed << std::setprecision(2) << res.gene_wide.busted_lrt
                  << ", p = " << std::setprecision(4) << res.gene_wide.busted_p_value
                  << "  │  ω₃ = " << std::setprecision(2) << res.gene_wide.omega_classes[2]
                  << " on " << std::setprecision(1) << (res.gene_wide.omega_weights[2] * 100.0)
                  << "% sites  │  ErrSink = " << std::setprecision(2) << (res.gene_wide.error_sink_weight * 100.0) << "%)";

        s_q3_dom << res.spatial_domains.size() << " 1D domains (";
        for (size_t i = 0; i < res.spatial_domains.size(); ++i) {
            if (i > 0) s_q3_dom << ", ";
            s_q3_dom << (res.spatial_domains[i].start_codon + 1) << "-" << (res.spatial_domains[i].end_codon + 1)
                     << ":" << res.spatial_domains[i].domain_type;
        }
        s_q3_dom << ")  │  SRV Islands: " << res.synonymous_islands.size();
        if (!res.synonymous_islands.empty()) {
            s_q3_dom << " (";
            for (size_t i = 0; i < res.synonymous_islands.size(); ++i) {
                if (i > 0) s_q3_dom << ", ";
                s_q3_dom << "codons " << (res.synonymous_islands[i].start_codon + 1) << "-"
                         << (res.synonymous_islands[i].end_codon + 1) << " [α="
                         << std::fixed << std::setprecision(2) << res.synonymous_islands[i].mean_alpha
                         << ", " << std::setprecision(2) << res.synonymous_islands[i].suppression_ratio << "x]";
            }
            s_q3_dom << ")";
        }

        s_q3_sites << "▲ Episodic (MEME): " << res.sites_episodic_diversifying
                   << "  │  ▲ Pervasive (FEL/dFEL): " << res.sites_pervasive_diversifying
                   << "  │  ▼ Purifying: " << res.sites_pervasive_purifying
                   << "  │  ⇄ Toggle/Convergent: " << res.sites_property_constrained_or_shifting
                   << "  │  Simes q<=" << std::fixed << std::setprecision(2) << settings.fdr_threshold << ": " << res.sites_simes_fdr_significant;

        s_q3_br << "aBSREL Diversifying: " << res.branches_episodic_diversifying
                << "  │  RELAX Intensified: " << res.branches_intensified
                << "  │  RELAX Relaxed: " << res.branches_relaxed
                << "  │  LoF / Core-Invasion: " << res.branches_loss_of_constraint
                << "  │  Clade Split: " << (res.clade_macro_regime.triggered ? res.clade_macro_regime.test_regime_label : "Homogeneous");

        s_rt << std::fixed << std::setprecision(2) << res.runtime_seconds << " s";

        std::vector<std::pair<std::string, std::string>> summary_items = {
            {"1. Is Data Trustworthy? (Homology/LOTO)", s_q1.str()},
            {"2. Confounding Processes (Drift/Network/CpG)", s_q2.str()},
            {"3A. Gene-Wide Adaptation (BUSTED-E)", s_q3_gene.str()},
            {"3B. 1D Domain & Synonymous Architecture", s_q3_dom.str()},
            {"3C. Site-Level Selection & Exchangeability", s_q3_sites.str()},
            {"3D. Lineage Selection & Loss-of-Constraint", s_q3_br.str()},
            {"Total Wall-Clock Time", s_rt.str()}
        };
        if (settings.enable_coevolution) {
            summary_items.insert(summary_items.end() - 1,
                {"Tier 3 Coevolution Pairs (--coevolution)", std::to_string(res.coevolution_pairs.size()) + " co-substituting site pairs"});
        }
        Panel::print_summary_card(
            "One-Button Evolutionary Characterization Summary",
            summary_items,
            "Co-inferred alignment, network, compositional regimes & unified gene/domain/branch/site/exchangeability atlas",
            true
        );
    }

    if (!out_json_path.empty()) {
        std::ofstream jf(out_json_path);
        if (jf.is_open()) {
            jf << CharacterizeEngine::to_json(res) << "\n";
            if (!settings.quiet) {
                std::cout << Console::success("Saved unified JSON characterization to: ") << out_json_path << "\n";
            }
        }
    }
    if (!out_html_path.empty()) {
        std::ofstream hf(out_html_path);
        if (hf.is_open()) {
            hf << CharacterizeEngine::render_standalone_html(res, short_label);
            if (!settings.quiet) {
                std::cout << Console::success("Saved standalone multi-tab SPA HTML report to: ") << out_html_path << "\n";
            }
        }
    }
    if (settings.emit_legacy_jsons) {
        CharacterizeEngine::write_legacy_json_bundle(res, input_label);
        if (!settings.quiet) {
            std::cout << Console::success("Saved legacy JSON bundle (.FEL.json, .MEME.json, .BUSTED.json) for: ") << input_label << "\n";
        }
    }
    if (!out_msa_path.empty()) {
        std::ofstream mf(out_msa_path);
        if (mf.is_open()) {
            for (const auto& s : res.joint_align.map_alignment.sequences) {
                mf << ">" << s.name << " [codons=" << res.joint_align.map_alignment.num_codons << "]\n" << s.raw_seq << "\n";
            }
            if (!settings.quiet) {
                std::cout << Console::success("Saved MAP codon alignment to: ") << out_msa_path << "\n";
            }
        }
    }
    if (!out_sto_path.empty()) {
        std::ofstream sf(out_sto_path);
        if (sf.is_open()) {
            sf << res.joint_align.stockholm_alignment;
            if (!settings.quiet) {
                std::cout << Console::success("Saved Stockholm quality-annotated MSA to: ") << out_sto_path << "\n";
            }
        }
    }
    if (!out_poa_path.empty()) {
        std::ofstream pf(out_poa_path);
        if (pf.is_open()) {
            pf << res.joint_align.poa_gfa;
            if (!settings.quiet) {
                std::cout << Console::success("Saved POA codon graph (GFAv1) to: ") << out_poa_path << "\n";
            }
        }
    }

    return 0;
}
