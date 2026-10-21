#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/analyses/relax.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>

using namespace hyphy::core;
using namespace hyphy::analyses;

static void print_relax_banner() {
    Panel::print_banner(
        "RELAX",
        "Test for Selection Relaxation",
        "Detecting relaxed or intensified selection in a phylogenetic framework",
        "Citation: Mol Biol Evol. 32(3): 820-832 (2015) • v3.0.0"
    );
}

static void print_relax_usage(const char* prog) {
    print_relax_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file or 'auto' (optional; auto-infers ML tree via TreeEngine if omitted or 'auto')\n"
              << "  " << Console::brand("--ohmic, --omhic") << "     Infer Ohmic network & Wheatstone bridges (empirical conductances)\n"
              << "  " << Console::brand("--hyphaeon") << "           Opt-in: HyphAeon 384D transformer conductances (needs Python aeon-core)\n"
              << "  " << Console::brand("--embeddings") << " <file>  Path to precomputed [N x D] taxon embeddings TSV/CSV for Ohmic\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal)\n"
              << "  " << Console::brand("--test, --branches") << " <spec> Test set: tag present in tree (e.g. T, TEST, FG), or comma-separated exact\n"
              << "                       branch names, or 'regex:<pattern>' (default: {TEST}/{T}/{FG} tags)\n"
              << "  " << Console::brand("--reference") << " <spec>   Reference set, same syntax (default: {REFERENCE}/{R}/{BG} tags or complement)\n"
              << "  " << Console::brand("--test-regex") << " <re>    Test set by case-insensitive regex on branch names/tags\n"
              << "  " << Console::brand("--reference-regex") << " <re> Reference set by case-insensitive regex on branch names/tags\n"
              << "  " << Console::brand("--pvalue") << " <threshold> Significance threshold for LRT (default: 0.05)\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.RELAX.json)\n"
              << "  " << Console::brand("--html") << " <file>        Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--no-branch-opt") << "      Disable individual branch length refinement (use proportional scaling)\n"
              << "  " << Console::brand("--joint-cycles") << " <N>   Max alternating (distribution <-> branch length) cycles (default: 15)\n"
              << "  " << Console::brand("--relax-normalization") << " <shared|per-class>  Omega-class rate normalization (default: shared =\n"
              << "                       one synonymous rate per branch, branch length = mixture expectation, as in HyPhy 2.5;\n"
              << "                       per-class = legacy hyphy-3 behavior, for comparison)\n"
              << "  " << Console::brand("--relax-mixture") << " <branch-site|site>  Omega mixture level (default: branch-site =\n"
              << "                       P_b = sum_c w_c exp(Q_c t_b) per branch, as in HyPhy 2.5 BS-REL; site = legacy hyphy-3)\n"
              << "  " << Console::brand("--mode") << " <type>        RELAX test mode (supported: 'Classic mode', default)\n"
              << "  " << Console::brand("--rates") << " <N>          Number of omega rate classes (supported: 3, default)\n"
              << "  " << Console::brand("--srv") << " <Yes|No>       Synonymous rate variation (supported: No in standalone hyphy_relax; use hyphy_busted --phenotype --relax)\n"
              << "  " << Console::brand("--multiple-hits") << " <m>  Multi-hit mode (supported: None in standalone hyphy_relax; use hyphy_busted --phenotype --relax)\n"
              << "  " << Console::brand("--error-sink") << " <Yes|No> Error-sink class (supported: No in standalone hyphy_relax; use hyphy_busted --phenotype --relax)\n"
              << "  " << Console::brand("--grid-size") << " <N>       Number of initial candidate points evaluated for RELAX fits (default: 8 null / 11 alt)\n"
              << "  " << Console::brand("--starting-points") << " <N> Number of top-ranked initial seeds refined via Nelder-Mead (default: 1)\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/Fig4E.nex\n"
              << "  " << prog << " --alignment data/cd2.fna --tree auto --output cd2.RELAX.json\n"
              << "  " << prog << " --alignment data/gene.fna --tree data/gene.nwk --test-regex \"PSEUDOGENE\" --threads 8\n\n";
}

int run_relax(int argc, char* argv[]) {
    CommonCLIOptions cli;
    RELAXSettings settings;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto status = parse_common_cli_arg(i, argc, argv, cli);
        if (status == CommonArgParseStatus::Handled) continue;
        if (status == CommonArgParseStatus::Error) return 1;

        if ((arg == "--test" || arg == "--branches") && i + 1 < argc) {
            settings.test_branch_regex = argv[++i];
        } else if (arg == "--reference" && i + 1 < argc) {
            settings.reference_branch_regex = argv[++i];
        } else if ((arg == "--pvalue" || arg == "-p") && i + 1 < argc) {
            settings.p_value_threshold = std::stod(argv[++i]);
        } else if (arg == "--no-branch-opt" || arg == "--no-refine-branches") {
            settings.refine_branch_lengths = false;
        } else if (arg == "--joint-cycles" && i + 1 < argc) {
            settings.max_joint_cycles = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--test-regex" && i + 1 < argc) {
            settings.test_branch_regex = argv[++i];
            settings.test_is_regex = true;
        } else if (arg == "--reference-regex" && i + 1 < argc) {
            settings.reference_branch_regex = argv[++i];
            settings.reference_is_regex = true;
        } else if (arg == "--relax-normalization" && i + 1 < argc) {
            std::string norm = argv[++i];
            if (norm == "per-class" || norm == "perclass") {
                settings.normalization = RELAXNormalization::PerClass;
            } else if (norm == "shared" || norm == "shared-mixture" || norm == "mixture") {
                settings.normalization = RELAXNormalization::SharedMixture;
            } else {
                std::cerr << "Error: --relax-normalization must be 'per-class' or 'shared' (got '" << norm << "').\n";
                return 1;
            }
        } else if (arg == "--relax-mixture" && i + 1 < argc) {
            std::string lvl = argv[++i];
            if (lvl == "site") {
                settings.mixture_level = RELAXMixtureLevel::Site;
            } else if (lvl == "branch-site" || lvl == "branchsite") {
                settings.mixture_level = RELAXMixtureLevel::BranchSite;
            } else {
                std::cerr << "Error: --relax-mixture must be 'site' or 'branch-site' (got '" << lvl << "').\n";
                return 1;
            }
        } else if (arg == "--mode" && i + 1 < argc) {
            std::string m = argv[++i];
            if (m != "Classic mode" && m != "Classic" && m != "classic") {
                std::cerr << "Error: Standalone hyphy_relax implements 'Classic mode' (2-way Test vs Reference); got '" << m << "'.\n";
                return 1;
            }
        } else if (arg == "--rates" && i + 1 < argc) {
            int r = std::stoi(argv[++i]);
            if (r != 3) {
                std::cerr << "Error: Standalone hyphy_relax implements 3 omega classes (--rates 3); got --rates " << r
                          << ". Use 'hyphy_busted --phenotype --relax --rates " << r << "' for arbitrary K.\n";
                return 1;
            }
        } else if (arg == "--srv") {
            std::string val = "Yes";
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                val = argv[++i];
            }
            if (val != "No" && val != "no" && val != "false" && val != "0") {
                std::cerr << "Error: Synonymous rate variation (--srv) is not supported in standalone hyphy_relax. "
                          << "Use 'hyphy_busted --phenotype --relax --srv' instead.\n";
                return 1;
            }
        } else if (arg == "--syn-rates" && i + 1 < argc) {
            int sr = std::stoi(argv[++i]);
            if (sr > 1) {
                std::cerr << "Error: --syn-rates > 1 is not supported in standalone hyphy_relax. "
                          << "Use 'hyphy_busted --phenotype --relax --srv --syn-rates " << sr << "' instead.\n";
                return 1;
            }
        } else if (arg == "--multiple-hits" && i + 1 < argc) {
            std::string mh = argv[++i];
            if (mh != "None" && mh != "none" && mh != "No" && mh != "no") {
                std::cerr << "Error: Multi-hit substitutions (--multiple-hits " << mh << ") are not supported in standalone hyphy_relax. "
                          << "Use 'hyphy_busted --phenotype --relax --multiple-hits " << mh << "' instead.\n";
                return 1;
            }
        } else if (arg == "--error-sink") {
            std::string val = "Yes";
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                val = argv[++i];
            }
            if (val != "No" && val != "no" && val != "false" && val != "0") {
                std::cerr << "Error: Error-sink mixture component (--error-sink) is not supported in standalone hyphy_relax. "
                          << "Use 'hyphy_busted --phenotype --relax --error-sink' instead.\n";
                return 1;
            }
        } else if (arg == "--grid-size" && i + 1 < argc) {
            int n_val = std::stoi(argv[++i]);
            if (n_val <= 0) {
                std::cerr << "Error: --grid-size must be positive.\n";
                return 1;
            }
            settings.grid_size = n_val;
        } else if (arg == "--starting-points" && i + 1 < argc) {
            int n_val = std::stoi(argv[++i]);
            if (n_val <= 0) {
                std::cerr << "Error: --starting-points must be positive.\n";
                return 1;
            }
            settings.starting_points = n_val;
        } else if (arg == "--help" || arg == "-h") {
            print_relax_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_relax_usage(argv[0]);
            return 1;
        }
    }

    if (cli.alignment_file.empty()) {
        std::cerr << "Error: --alignment is required.\n\n";
        print_relax_usage(argv[0]);
        return 1;
    }

    finalize_common_cli_options(cli, ".RELAX.json");
    print_relax_banner();

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string start_iso = Provenance::current_iso8601();
    std::string cli_cmd = join_cli_command(argc, argv);

    LoadedCodonDataset ds;
    if (!load_codon_dataset(cli, ds, false)) {
        return 1;
    }
    auto& code = ds.gcode;
    auto& aln = ds.aln;
    auto& tree = ds.tree;
    const auto& tree_source = ds.tree_source;

    // Initialize and Run RELAX
    auto relax = RELAXAnalyzer::create(tree, aln, settings);

    if (relax.test_branch_names.empty()) {
        std::cerr << Console::danger("Error: No Test branches found! Annotate branches with {TEST}/{T}/{FG} in Newick or pass --test <spec>.\n");
        return 1;
    }

    std::vector<std::pair<std::string, std::string>> card_items = {
        {"Alignment File", cli.alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons) + " (" + std::to_string(aln.num_codons * 3) + " nt)"},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree", tree_source + " (" + std::to_string(tree.num_nodes()) + " nodes, " + std::to_string(tree.num_leaves()) + " leaves)"},
        {"Test Branches", std::to_string(relax.test_branch_names.size()) + " branches"},
        {"Reference Branches", std::to_string(relax.ref_branch_names.size()) + " branches"}
    };
    if (!relax.unclassified_branch_names.empty()) {
        card_items.push_back({"Unclassified Branches", std::to_string(relax.unclassified_branch_names.size()) + " branches"});
    }
    card_items.push_back({"Genetic Code", code->name});
    card_items.push_back({"Omega Mixture",
        std::string(settings.mixture_level == RELAXMixtureLevel::BranchSite ? "branch-site" : "site (legacy)") + ", " +
        (settings.normalization == RELAXNormalization::SharedMixture ? "shared synonymous rate" : "per-class normalization (legacy)")});
    Panel::print_card("Dataset & Phylogeny", card_items);

    Panel::print_step(1, 3, "Fitting Baseline & Separate-Rates Codon Models", "GTR & MG94xREV");

    std::function<void(const std::string&, double)> progress_cb = nullptr;
    if (!cli.show_progress && !cli.force_progress) {
        progress_cb = [](const std::string& stage, double frac) {
            std::cout << "  " << Console::muted("↳ [") << std::setw(3) << static_cast<int>(frac * 100) 
                      << Console::muted("%] ") << stage << std::endl;
        };
    }

    Panel::print_step(2, 3, "Running RELAX Inference Pipeline", 
                      settings.refine_branch_lengths ? "Two-loop L-BFGS branch refinement enabled" : "Proportional branch scaling");

    auto res = relax.run(progress_cb, cli.show_progress, cli.force_progress);

    // Rate distributions
    Table rate_table;
    rate_table.add_column("Class", Table::Align::Center, 5);
    rate_table.add_column("Ref omega", Table::Align::Right, 12);
    rate_table.add_column("Ref Proportion", Table::Align::Right, 14);
    rate_table.add_column("Test omega", Table::Align::Right, 12);
    rate_table.add_column("Test Proportion", Table::Align::Right, 15);

    for (size_t c = 0; c < res.alternative_fit.reference_distribution.omegas.size(); ++c) {
        std::ostringstream w_r_ss, p_r_ss, w_t_ss, p_t_ss;
        w_r_ss << std::fixed << std::setprecision(4) << res.alternative_fit.reference_distribution.omegas[c];
        p_r_ss << std::fixed << std::setprecision(4) << res.alternative_fit.reference_distribution.weights[c];
        w_t_ss << std::fixed << std::setprecision(4) << res.alternative_fit.test_distribution.omegas[c];
        p_t_ss << std::fixed << std::setprecision(4) << res.alternative_fit.test_distribution.weights[c];

        rate_table.add_row({
            std::to_string(c),
            w_r_ss.str(),
            p_r_ss.str(),
            w_t_ss.str(),
            p_t_ss.str()
        });
    }

    if (res.is_failed) {
        std::cerr << Console::danger("Error: RELAX analysis failed to converge or likelihood underflowed: ")
                  << res.error_message << "\n";
        return 1;
    }

    std::ostringstream k_ss, lrt_ss, p_ss, gtr_ss, mg_ss, alt_ss, null_ss, time_ss;
    k_ss << std::fixed << std::setprecision(4) << res.k;
    lrt_ss << std::fixed << std::setprecision(2) << res.lrt;
    p_ss << std::scientific << std::setprecision(4) << res.p_value;
    gtr_ss << std::fixed << std::setprecision(2) << res.gtr_log_likelihood << " (AICc: " << res.gtr_aicc << ")";
    mg_ss << std::fixed << std::setprecision(2) << res.mg94_log_likelihood << " (AICc: " << res.mg94_aicc << ")";
    alt_ss << std::fixed << std::setprecision(2) << res.alternative_fit.log_likelihood << " (AICc: " << res.alternative_fit.aicc << ")";
    null_ss << std::fixed << std::setprecision(2) << res.null_fit.log_likelihood << " (AICc: " << res.null_fit.aicc << ")";
    time_ss << std::fixed << std::setprecision(2) << res.runtime_seconds << " s";

    std::string conclusion;
    if (res.is_significant) {
        if (res.is_relaxed) {
            conclusion = "SIGNIFICANT RELAXATION of selection on Test branches (K = " + k_ss.str() + " < 1, p = " + p_ss.str() + ")";
        } else {
            conclusion = "SIGNIFICANT INTENSIFICATION of selection on Test branches (K = " + k_ss.str() + " > 1, p = " + p_ss.str() + ")";
        }
    } else {
        conclusion = "No significant evidence of relaxation or intensification (p = " + p_ss.str() + " > " + std::to_string(settings.p_value_threshold) + ")";
    }

    std::string json_desc = (cli.json_format == JSONFormat::Legacy) ? " [Legacy]" : " [Modern v3.0]";
    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Relaxation Parameter (K)", k_ss.str()},
        {"Likelihood Ratio Test (LRT)", lrt_ss.str()},
        {"p-value (asymptotic)", p_ss.str() + " [threshold: " + std::to_string(settings.p_value_threshold) + "]"},
        {"RELAX Alternative Log-L", alt_ss.str()},
        {"RELAX Null (K=1) Log-L", null_ss.str()},
        {"Nucleotide GTR Log-L", gtr_ss.str()},
        {"Separate Rates MG94xREV", mg_ss.str()},
        {"JSON Output File", cli.output_file + json_desc},
        {"Total Execution Time", time_ss.str()}
    };

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_runtime = std::chrono::duration<double>(end_time - start_time).count();

    Provenance prov = build_cli_provenance(cli, cli_cmd, start_iso, total_runtime, {
        {"test", settings.test_branch_regex},
        {"pvalue", std::to_string(settings.p_value_threshold)},
        {"refine_branch_lengths", settings.refine_branch_lengths ? "true" : "false"}
    });

    // Save JSON (and optional HTML) before printing final tables/summary
    try {
        res.save_json(cli.output_file, tree, aln, cli.json_format, prov);
        std::string fmt_desc = (cli.json_format == JSONFormat::Legacy) ? "Legacy Datamonkey JSON" : "Modern HyPhy v3.0 JSON";
        std::cout << Console::success("✔") << " " << Console::bold("Saved " + fmt_desc + " report to: ")
                  << Console::brand(cli.output_file) << "\n\n";
        if (!cli.html_file.empty()) {
            hyphy::report::write_html_report(cli.output_file, cli.html_file);
            std::cout << Console::success("✔") << " " << Console::bold("Saved HTML report to: ")
                      << Console::brand(cli.html_file) << "\n\n";
        }
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: Failed to write output: ") << e.what() << "\n";
        return 1;
    }

    std::cout << "\n" << Console::bold("  Inferred Rate Distributions (Alternative Model):") << "\n\n";
    rate_table.print();
    std::cout << "\n";

    Panel::print_summary_card("RELAX Selection Relaxation Analysis Summary", sum_items, conclusion, res.is_significant);
    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_relax(argc, argv);
}
#endif
