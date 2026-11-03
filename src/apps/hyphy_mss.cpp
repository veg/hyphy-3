#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/mss.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <algorithm>

using namespace hyphy::core;
using namespace hyphy::analyses;
using namespace hyphy::opt;

static void print_mss_banner() {
    Panel::print_banner(
        "MSS",
        "Multi-class Synonymous Substitution Models & Partition Selector",
        "Codon-pair synonymous exchangeabilities, SynREV / SynREVCodon / SynREV2, & fast DP partition discovery",
        "Citation: Mol Biol Evol. (HyPhy MSS / SynREV) • v3.0.0"
    );
}

static void print_mss_usage(const char* prog) {
    print_mss_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Input Dataset Options (specify one):") << "\n"
              << "  " << Console::brand("--alignment") << " <file>   Path to single-gene codon alignment (FASTA or NEXUS)\n"
              << "  " << Console::brand("--filelist") << " <file>    Path to file listing multiple codon alignments (1 per line, or <aln> <tree>)\n\n"
              << Console::bold("Model & Rate Options:") << "\n"
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file (for single-gene mode; auto-inferred if omitted)\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal)\n"
              << "  " << Console::brand("--model") << " <spec>       Synonymous rate model (alias: --type, --mss):\n"
              << "                       • " << Console::bold("SynREV") << "       : 18 amino-acid synonymous rate classes (default)\n"
              << "                       • " << Console::bold("SynREVCodon") << "  : 68 exchangeable codon-pair synonymous rate classes\n"
              << "                       • " << Console::bold("SynREV2") << "      : Per-AA Transition/Transversion split (~27 classes)\n"
              << "                       • " << Console::bold("SynREV2g") << "     : Per-AA rates + shared global Tv/Ti modifier (19 classes)\n"
              << "                       • " << Console::bold("<path.tsv>") << "   : Custom TSV file (auto-detects 3-col codon, 3-col pair, 5-col GA, or empirical fixed rates)\n"
              << "  " << Console::brand("--neutral") << " <class>    Neutral reference class fixed at rate = 1.0 (e.g. NEUTRAL)\n"
              << "  " << Console::brand("--norm") << " <mode>        Rate normalization constraint (default: flux for --omega Fit, unweighted for --omega Fix, neutral_ref with --neutral):\n"
              << "                       • " << Console::bold("flux") << "         : Equilibrium synonymous flux-weighted mean = 1.0\n"
              << "                       • " << Console::bold("unweighted") << "   : Arithmetic mean of class rates = 1.0\n"
              << "                       • " << Console::bold("neutral_ref") << "  : Reference class fixed at 1.0\n"
              << "  " << Console::brand("--omega") << " <Fix|Fit>    Fix gene omega at MG94 baseline ('Fix') or jointly re-optimize ('Fit', default)\n"
              << "  " << Console::brand("--f81") << "                Zero out nucleotide GTR biases (theta_ij = 1) with coordinate-wise cavity branch optimization\n\n"
              << Console::bold("Fast Partition Selector Options:") << "\n"
              << "  " << Console::brand("--select") << "             Discover optimal synonymous rate partition via 1D profile curvature DP + local search\n"
              << "  " << Console::brand("--classes") << " <K>        Number of synonymous rate classes K for --select (default: 5, max: 12)\n"
              << "  " << Console::brand("--auto-k") << "             Step up K = 2..max_K under --ic criterion (enabled by default with --select)\n"
              << "  " << Console::brand("--fixed-k") << "            Fit only the exact --classes <K> without stepping up from K = 2\n"
              << "  " << Console::brand("--ic") << " <BIC|AIC-c>     Information criterion for partition selection (default: BIC)\n"
              << "  " << Console::brand("--tsv") << " <file>         Export discovered partition & model-averaged rates to TSV (compatible with --model)\n\n"
              << Console::bold("Inference & Output Options:") << "\n"
              << "  " << Console::brand("--no-ci") << "              Skip 95% profile-likelihood confidence intervals\n"
              << "  " << Console::brand("--no-lrt") << "             Skip per-class LRTs (alpha_k == 1) and Holm-Bonferroni p-values\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.MSS.json)\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment tests/data/adh.nex --model SynREV\n"
              << "  " << prog << " --alignment tests/data/adh.nex --model codons.tsv --neutral NEUTRAL\n"
              << "  " << prog << " --alignment tests/data/adh.nex --select --classes 5 --tsv adh.mss.tsv\n"
              << "  " << prog << " --filelist genes.txt --model SynREV2g --norm flux\n\n";
}

int run_mss(int argc, char* argv[]) {
    CommonCLIOptions cli;
    MSSSettings settings;
    settings.select_classes = 5;
    settings.select_auto_k = true;
    std::string filelist_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto status = parse_common_cli_arg(i, argc, argv, cli);
        if (status == CommonArgParseStatus::Handled) continue;
        if (status == CommonArgParseStatus::Error) return 1;

        if ((arg == "--filelist" || arg == "--list" || arg == "--files") && i + 1 < argc) {
            filelist_path = argv[++i];
        } else if ((arg == "--model" || arg == "--type" || arg == "--mss" || arg == "--syn-model") && i + 1 < argc) {
            settings.model_type = argv[++i];
        } else if ((arg == "--neutral" || arg == "--mss-neutral" || arg == "--reference") && i + 1 < argc) {
            settings.neutral_ref = argv[++i];
        } else if ((arg == "--norm" || arg == "--normalization" || arg == "--mss-norm") && i + 1 < argc) {
            std::string n_str = argv[++i];
            settings.normalization_explicit = true;
            if (n_str == "flux" || n_str == "Flux" || n_str == "weighted") {
                settings.normalization = MSSNormalization::Flux;
            } else if (n_str == "unweighted" || n_str == "Unweighted" || n_str == "mean") {
                settings.normalization = MSSNormalization::Unweighted;
            } else if (n_str == "neutral_ref" || n_str == "neutral" || n_str == "ref" || n_str == "NeutralRef") {
                settings.normalization = MSSNormalization::NeutralRef;
            } else {
                std::cerr << Console::danger("Error: ") << "Unknown normalization mode '" << n_str
                          << "'. Choose flux, unweighted, or neutral_ref.\n";
                return 1;
            }
        } else if (arg == "--omega" && i + 1 < argc) {
            std::string o_str = argv[++i];
            settings.fix_omega = (o_str == "Fix" || o_str == "fix" || o_str == "Fixed" || o_str == "fixed");
        } else if (arg == "--fix-omega") {
            settings.fix_omega = true;
        } else if (arg == "--f81" || arg == "--F81") {
            settings.use_f81 = true;
        } else if (arg == "--select" || arg == "--selector" || arg == "--partition-select") {
            settings.select_partition = true;
        } else if ((arg == "--classes" || arg == "--rates" || arg == "-K") && i + 1 < argc) {
            settings.select_classes = std::stoul(argv[++i]);
        } else if (arg == "--auto-k") {
            settings.select_auto_k = true;
        } else if (arg == "--fixed-k" || arg == "--no-auto-k") {
            settings.select_auto_k = false;
        } else if (arg == "--ic" && i + 1 < argc) {
            settings.ic_type = argv[++i];
        } else if ((arg == "--tsv" || arg == "--export-tsv") && i + 1 < argc) {
            settings.output_tsv = argv[++i];
        } else if (arg == "--no-ci") {
            settings.compute_ci = false;
        } else if (arg == "--no-lrt") {
            settings.compute_lrt = false;
        } else if (arg == "--help" || arg == "-h") {
            print_mss_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_mss_usage(argv[0]);
            return 1;
        }
    }

    if (cli.alignment_file.empty() && filelist_path.empty()) {
        std::cerr << "Error: Either --alignment <file> or --filelist <file> is required.\n\n";
        print_mss_usage(argv[0]);
        return 1;
    }

    if (!settings.normalization_explicit) {
        if (!settings.neutral_ref.empty()) {
            settings.normalization = MSSNormalization::NeutralRef;
        } else if (settings.fix_omega) {
            settings.normalization = MSSNormalization::Unweighted;
        }
    }

    finalize_common_cli_options(cli, ".MSS.json");
    if (cli.output_file.empty() && !filelist_path.empty()) {
        cli.output_file = filelist_path + ".MSS.json";
    }

    print_mss_banner();

    auto code = GeneticCode::from_name(cli.code_name);
    MSSAnalyzer analyzer(code);

    if (!cli.alignment_file.empty()) {
        LoadedCodonDataset ds;
        if (!load_codon_dataset(cli, ds, true)) {
            return 1;
        }
        analyzer.add_dataset(cli.alignment_file, ds.tree, ds.aln, settings.use_f81);
    } else {
        std::ifstream fl(filelist_path);
        if (!fl.is_open()) {
            std::cerr << Console::danger("Error: ") << "Cannot open filelist: " << filelist_path << "\n";
            return 1;
        }
        std::string line;
        size_t loaded_count = 0;
        while (std::getline(fl, line)) {
            size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos || line[first] == '#') continue;
            std::istringstream iss(line.substr(first));
            std::string aln_path, tree_path;
            iss >> aln_path;
            if (iss >> tree_path) {
                // Explicit tree path on same line
            }
            CommonCLIOptions item_cli = cli;
            item_cli.alignment_file = aln_path;
            item_cli.tree_file = tree_path;
            LoadedCodonDataset ds;
            if (!load_codon_dataset(item_cli, ds, false)) {
                std::cerr << Console::danger("Error: ") << "Failed to load dataset from filelist: " << aln_path << "\n";
                return 1;
            }
            analyzer.add_dataset(aln_path, ds.tree, ds.aln, settings.use_f81);
            loaded_count++;
        }
        if (loaded_count == 0) {
            std::cerr << Console::danger("Error: ") << "No valid datasets found in filelist: " << filelist_path << "\n";
            return 1;
        }
    }

    size_t total_codons = 0, total_patterns = 0;
    for (const auto& d : analyzer.genes) {
        total_codons += d.aln.num_codons;
        total_patterns += d.aln.patterns.size();
    }

    std::string norm_str = (settings.normalization == MSSNormalization::Flux)
        ? "Flux-Weighted Mean = 1.0"
        : ((settings.normalization == MSSNormalization::Unweighted) ? "Unweighted Mean = 1.0" : "Neutral Reference = 1.0");

    Panel::print_card("MSS Dataset & Model Configuration", {
        {"Input Source", !cli.alignment_file.empty() ? cli.alignment_file : (filelist_path + " (" + std::to_string(analyzer.genes.size()) + " alignments)")},
        {"Total Codon Sites", std::to_string(total_codons) + " (" + std::to_string(total_patterns) + " unique patterns)"},
        {"Synonymous Model", settings.select_partition ? (settings.model_type + " + Fast DP Selector (K <= " + std::to_string(settings.select_classes) + ", " + settings.ic_type + ")") : settings.model_type},
        {"Neutral Reference", settings.neutral_ref.empty() ? "Auto / None" : settings.neutral_ref},
        {"Rate Normalization", norm_str},
        {"Omega Handling", settings.fix_omega ? "Fixed at MG94 baseline" : "Jointly optimized with MSS rates"},
        {"Nucleotide Model", settings.use_f81 ? "F81-like (theta_ij = 1)" : "General Time Reversible (GTR)"},
        {"Genetic Code", code->name}
    });

    Panel::print_step(1, 2, "Fitting Standard MG94 Baseline & Synonymous Flux Weights",
                      "Optimizing nucleotide exchangeabilities, branch lengths & equilibrium codon pair flux");
    Panel::print_step(2, 2, settings.select_partition ? "Running Fast Profile-Curvature DP Partition Selector" : ("Fitting MSS Model (" + settings.model_type + ")"),
                      settings.select_partition ? "68-pair profile screening + exact weighted k-means DP + boundary local search"
                                                : "Coordinate-wise Brent + joint omega/scaler refinement + profile CIs & LRTs");

    MSSResult res = analyzer.run(settings);

    std::cout << "\n  " << Console::muted("Baseline MG94 Log-L : ") << Console::bold(std::to_string(res.baseline_log_l))
              << Console::muted("  |  AIC-c = ") << std::fixed << std::setprecision(2) << res.baseline_aicc
              << Console::muted("  |  BIC = ") << res.baseline_bic << "\n"
              << "  " << Console::muted("MSS Fitted Log-L    : ") << Console::bold(std::to_string(res.mss_log_l))
              << Console::muted("  |  AIC-c = ") << std::fixed << std::setprecision(2) << res.mss_aicc
              << Console::muted("  |  BIC = ") << res.mss_bic << "\n\n";

    if (res.genes.size() > 1) {
        std::cout << Console::bold("  Per-Alignment Summary:") << "\n\n";
        Table ds_table;
        ds_table.add_column("Dataset", Table::Align::Left, 28);
        ds_table.add_column("Taxa", Table::Align::Right, 6);
        ds_table.add_column("Codons", Table::Align::Right, 8);
        ds_table.add_column("MG94 omega", Table::Align::Right, 12);
        ds_table.add_column("MSS omega", Table::Align::Right, 12);
        ds_table.add_column("Tree Scale", Table::Align::Right, 12);
        ds_table.add_column("Syn Length", Table::Align::Right, 12);
        for (const auto& d : res.genes) {
            std::ostringstream o0, o1, ts, sl;
            o0 << std::fixed << std::setprecision(4) << d.baseline_omega;
            o1 << std::fixed << std::setprecision(4) << d.mss_omega;
            ts << std::fixed << std::setprecision(4) << d.tree_scaler;
            sl << std::fixed << std::setprecision(4) << d.synonymous_tree_length;
            ds_table.add_row({
                d.filepath,
                std::to_string(d.num_taxa),
                std::to_string(d.num_codons),
                o0.str(), o1.str(), ts.str(), sl.str()
            });
        }
        ds_table.print();
        std::cout << "\n";
    } else if (res.genes.size() == 1) {
        const auto& d = res.genes[0];
        std::cout << "  " << Console::muted("MG94 omega = ") << std::fixed << std::setprecision(4) << d.baseline_omega
                  << Console::muted("  ->  MSS omega = ") << Console::brand(std::to_string(d.mss_omega))
                  << Console::muted("  |  Synonymous Tree Length = ") << std::fixed << std::setprecision(4) << d.synonymous_tree_length
                  << Console::muted("  |  Non-synonymous Tree Length = ") << d.nonsynonymous_tree_length << "\n\n";
    }

    std::cout << Console::bold("  Inferred Synonymous Substitution Rate Classes (" + res.partition.name + "):") << "\n\n";
    Table rate_table;
    rate_table.add_column("Class", Table::Align::Left, 22);
    rate_table.add_column("Pairs", Table::Align::Right, 6);
    rate_table.add_column("Flux %", Table::Align::Right, 8);
    rate_table.add_column("Rate (alpha)", Table::Align::Right, 13);
    if (settings.compute_ci) {
        rate_table.add_column("95% Profile CI", Table::Align::Center, 22);
    }
    if (settings.compute_lrt) {
        rate_table.add_column("LRT (a=1)", Table::Align::Right, 11);
        rate_table.add_column("p (Holm)", Table::Align::Right, 12);
    }

    size_t sig_classes = 0;
    for (const auto& cr : res.rates) {
        std::ostringstream f_ss, r_ss, ci_ss, lrt_ss, p_ss;
        f_ss << std::fixed << std::setprecision(1) << (cr.equilibrium_flux_share * 100.0) << "%";
        r_ss << std::fixed << std::setprecision(4) << cr.rate;
        if (cr.is_reference) {
            ci_ss << "[fixed reference]";
            lrt_ss << "—";
            p_ss << "—";
        } else if (res.partition.is_empirical_fixed) {
            ci_ss << "[empirical fixed]";
            lrt_ss << "—";
            p_ss << "—";
        } else {
            ci_ss << "[" << std::fixed << std::setprecision(4) << cr.ci_low << ", " << cr.ci_high << "]";
            lrt_ss << std::fixed << std::setprecision(3) << cr.lrt_vs_one;
            p_ss << std::scientific << std::setprecision(3) << cr.p_value_holm;
        }

        std::string row_color = "";
        if (!cr.is_reference && !res.partition.is_empirical_fixed && settings.compute_lrt && cr.p_value_holm <= 0.05) {
            row_color = "\033[38;5;48m";
            sig_classes++;
        }

        std::vector<std::string> row = {
            cr.class_name,
            std::to_string(cr.num_pairs),
            f_ss.str(),
            r_ss.str()
        };
        if (settings.compute_ci) row.push_back(ci_ss.str());
        if (settings.compute_lrt) {
            row.push_back(lrt_ss.str());
            row.push_back(p_ss.str());
        }
        rate_table.add_row(row, row_color);
    }
    rate_table.print();
    std::cout << "\n";

    if (res.selection.enabled) {
        std::cout << Console::bold("  Fast Curvature-Guided DP Partition Selection (" + res.selection.ic_type + "):") << "\n\n";
        std::cout << "    • Null MG94 " << res.selection.ic_type << "       : " << std::fixed << std::setprecision(2) << res.selection.null_mg94_ic << "\n"
                  << "    • SynREVCodon (68-pair) : lnL = " << res.selection.synrev_codon_log_l << ", " << res.selection.ic_type << " = " << res.selection.synrev_codon_ic << "\n"
                  << "    • Selected Partition (K=" << res.selection.num_classes << ") : lnL = " << res.selection.best_partition_log_l
                  << ", " << res.selection.ic_type << " = " << res.selection.best_partition_ic
                  << " (" << res.selection.models_evaluated << " candidate partitions evaluated)\n\n";

        Table sel_table;
        sel_table.add_column("Discovered Class", Table::Align::Left, 20);
        sel_table.add_column("Pairs", Table::Align::Right, 8);
        sel_table.add_column("Class Rate", Table::Align::Right, 14);
        for (size_t k = 0; k < res.selection.class_labels.size(); ++k) {
            size_t count = 0;
            for (const auto& pe : res.selection.pair_entries) {
                if (static_cast<size_t>(pe.assigned_class) == k) count++;
            }
            std::ostringstream cr_ss;
            cr_ss << std::fixed << std::setprecision(4) << res.selection.class_rates[k];
            sel_table.add_row({
                res.selection.class_labels[k],
                std::to_string(count),
                cr_ss.str()
            }, "\033[38;5;48m");
        }
        sel_table.print();
        std::cout << "\n";
        if (!settings.output_tsv.empty()) {
            std::cout << Console::success("✔") << " " << Console::bold("Exported discovered MSS partition TSV to: ")
                      << Console::brand(settings.output_tsv) << "\n\n";
        }
    }

    std::ostringstream lrt_ss, p_ss, time_ss;
    lrt_ss << std::fixed << std::setprecision(4) << res.omnibus_lrt << " (df = " << res.omnibus_df << ")";
    p_ss << std::scientific << std::setprecision(6) << res.omnibus_p_value;
    time_ss << std::fixed << std::setprecision(2) << res.runtime_seconds << " s";

    bool is_sig = (res.omnibus_df > 0 && res.omnibus_p_value < 0.05);
    std::string verdict = is_sig
        ? ("STATISTICALLY SIGNIFICANT synonymous rate heterogeneity across codon classes (Omnibus LRT = " + lrt_ss.str() + ", p = " + p_ss.str() + ")!")
        : "No statistically significant synonymous rate heterogeneity detected against standard MG94 (p >= 0.05).";

    Panel::print_summary_card("MSS Synonymous Substitution Analysis Summary", {
        {"Synonymous Model", res.partition.name + " (" + std::to_string(res.rates.size()) + " classes)"},
        {"Omnibus LRT vs MG94", lrt_ss.str()},
        {"Omnibus p-value", p_ss.str()},
        {"Classes != 1.0 (Holm p <= 0.05)", std::to_string(sig_classes) + " / " + std::to_string(res.rates.size())},
        {"JSON Output File", cli.output_file},
        {"Total Execution Time", time_ss.str()}
    }, verdict, is_sig);

    try {
        MSSAnalyzer::save_json(cli.output_file, res);
        std::cout << Console::success("✔") << " " << Console::bold("Saved MSS JSON report to: ")
                  << Console::brand(cli.output_file) << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Warning: Failed to write JSON output: ") << e.what() << "\n";
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_mss(argc, argv);
}
#endif
