#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/fubar.hpp"
#include "hyphy/report/html_report.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;
using namespace hyphy::opt;

static void print_banner() {
    Panel::print_banner(
        "FUBAR",
        "Fast, Unconstrained Bayesian AppRoximation",
        "Site-level detection of positive and purifying selection via bivariate rate grid",
        "Citation: Mol Biol Evol (2013) 30(5):1196–1205 • v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment, -a") << " <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree, -t") << " <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  " << Console::brand("--code") << " <name>            Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  " << Console::brand("--grid, -g") << " <N>           Number of grid points per dimension (default: 20; total D^2 = 400)\n"
              << "  " << Console::brand("--method, -m") << " <name>      Inference method: 'Variational-Bayes' (default), 'Collapsed-Gibbs', or 'Metropolis-Hastings'\n"
              << "  " << Console::brand("--concentration") << " <val>    Dirichlet prior concentration parameter (default: 0.5)\n"
              << "  " << Console::brand("--posterior, -p") << " <float>  Posterior probability cutoff for selection calls (default: 0.90)\n"
              << "  " << Console::brand("--non-zero") << "               Enforce non-zero synonymous rates (min 1e-3)\n"
              << "  " << Console::brand("--tree-scaling") << " <mode>    Tree scaling mode: 'grid' (default, HyPhy 2.5 standard) or 'mg94'\n"
              << "  " << Console::brand("--chains") << " <N>             Number of MCMC chains for MCMC methods (default: 1)\n"
              << "  " << Console::brand("--chain-length") << " <N>       MCMC chain length (default: 2000000)\n"
              << "  " << Console::brand("--burn-in") << " <N>            MCMC burn-in steps (default: 1000000)\n"
              << "  " << Console::brand("--samples") << " <N>            MCMC samples to draw (default: 100)\n"
              << "  " << Console::brand("--threads, -T") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output, -o") << " <file>      Path to output JSON file (default: <alignment>.FUBAR.json)\n"
              << "  " << Console::brand("--html") << " <file>            Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--all-sites") << "              Display all codon sites in console table (default: significant only)\n"
              << "  " << Console::brand("--progress") << "               Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "            Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "               Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " --alignment data/cd2.fna --grid 20 --method Variational-Bayes\n"
              << "  " << prog << " --alignment test.nex --posterior 0.95 --output test.fubar.json\n\n";
}

int run_fubar(int argc, char* argv[]) {
    std::string alignment_file;
    std::string tree_file;
    std::string output_file;
    std::string html_file;
    std::string code_name = "Universal";
    int num_threads = 0;
    
    FUBARSettings settings;
    std::string json_format = "modern";
    bool show_progress = ProgressBar::is_terminal();
    bool all_sites = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--alignment" || arg == "-a") && i + 1 < argc) {
            alignment_file = argv[++i];
        } else if ((arg == "--tree" || arg == "-t") && i + 1 < argc) {
            tree_file = argv[++i];
        } else if (arg == "--code" && i + 1 < argc) {
            code_name = argv[++i];
        } else if ((arg == "--grid" || arg == "-g") && i + 1 < argc) {
            settings.grid_points = std::stoi(argv[++i]);
        } else if ((arg == "--method" || arg == "-m") && i + 1 < argc) {
            settings.method = string_to_fubar_method(argv[++i]);
        } else if (arg == "--concentration" && i + 1 < argc) {
            settings.concentration = std::stod(argv[++i]);
        } else if ((arg == "--posterior" || arg == "-p") && i + 1 < argc) {
            settings.posterior_threshold = std::stod(argv[++i]);
        } else if (arg == "--non-zero") {
            settings.non_zero = true;
        } else if (arg == "--tree-scaling" && i + 1 < argc) {
            settings.tree_scaling = argv[++i];
        } else if (arg == "--chains" && i + 1 < argc) {
            settings.mcmc_chains = std::stoi(argv[++i]);
        } else if (arg == "--chain-length" && i + 1 < argc) {
            settings.mcmc_chain_length = std::stoll(argv[++i]);
        } else if (arg == "--burn-in" && i + 1 < argc) {
            settings.mcmc_burn_in = std::stoll(argv[++i]);
        } else if (arg == "--samples" && i + 1 < argc) {
            settings.mcmc_samples = std::stoll(argv[++i]);
        } else if ((arg == "--threads" || arg == "-T") && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--html" && i + 1 < argc) {
            html_file = argv[++i];
        } else if (arg == "--all-sites") {
            all_sites = true;
        } else if (arg == "--legacy") {
            json_format = "legacy";
        } else if (arg == "--modern") {
            json_format = "modern";
        } else if (arg == "--progress") {
            show_progress = true;
        } else if (arg == "--no-progress") {
            show_progress = false;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        }
    }

    if (alignment_file.empty()) {
        std::cerr << Console::danger("Error: Missing required argument --alignment.\n\n");
        print_usage(argv[0]);
        return 1;
    }

    if (output_file.empty()) {
        output_file = alignment_file + ".FUBAR.json";
    }

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
#endif

    print_banner();

    // 1. Load Alignment & Tree
    std::shared_ptr<const GeneticCode> gcode;
    try {
        gcode = GeneticCode::from_name(code_name);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error loading genetic code: ") << e.what() << "\n";
        return 1;
    }

    Alignment aln;
    try {
        aln = Alignment::load(alignment_file, gcode);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error reading input files: ") << e.what() << "\n";
        return 1;
    }

    Tree tree;
    std::string tree_source;
    if (!tree_file.empty()) {
        try {
            tree = Tree::from_newick_file(tree_file);
            tree_source = tree_file;
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error loading tree: ") << e.what() << "\n";
            return 1;
        }
    } else if (!aln.embedded_tree_newick.empty()) {
        try {
            tree = Tree::from_newick(aln.embedded_tree_newick);
            tree_source = "embedded in NEXUS";
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error parsing embedded tree: ") << e.what() << "\n";
            return 1;
        }
    } else {
        std::cerr << Console::danger("Error: No tree provided via --tree and no embedded tree found in alignment.\n\n");
        print_usage(argv[0]);
        return 1;
    }

    Panel::print_card("Dataset & Configuration", {
        {"Alignment File", alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons)},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree Source", tree_source + " (" + std::to_string(tree.num_nodes()) + " nodes)"},
        {"Genetic Code", gcode->name},
        {"Inference Method", fubar_method_to_string(settings.method)},
        {"Grid Points", std::to_string(settings.grid_points) + " x " + std::to_string(settings.grid_points) + " (" + std::to_string(settings.grid_points * settings.grid_points) + ")"},
        {"Concentration", format_float_val(settings.concentration, 2)},
        {"Posterior Threshold", format_float_val(settings.posterior_threshold, 2)},
#ifdef _OPENMP
        {"OpenMP Worker Threads", std::to_string((num_threads > 0) ? num_threads : omp_get_max_threads())}
#else
        {"OpenMP Worker Threads", "1 (Serial Build)"}
#endif
    });

    // 2. Initialize and Fit
    auto progress_cb = [&](const std::string& msg, double pct) {
        if (!show_progress) return;
        std::cout << Console::brand("==> ") << Console::bold(msg) 
                  << Console::muted(" [" + std::to_string(static_cast<int>(pct * 100)) + "%]") << "\n";
    };

    FUBARAnalyzer fubar = FUBARAnalyzer::create_and_fit(tree, aln, settings, progress_cb);
    FUBARResult result = fubar.run(progress_cb, show_progress);

    // 3. Console Presentation
    FUBARAnalyzer::print_summary(result, settings);

    if (all_sites) {
        Table all_tbl;
        all_tbl.add_column("Codon", Table::Align::Right);
        all_tbl.add_column("alpha", Table::Align::Right);
        all_tbl.add_column("beta", Table::Align::Right);
        all_tbl.add_column("beta - alpha", Table::Align::Right);
        all_tbl.add_column("Prob[alpha>beta]", Table::Align::Right);
        all_tbl.add_column("Prob[alpha<beta]", Table::Align::Right);
        all_tbl.add_column("BayesFactor", Table::Align::Right);
        all_tbl.add_column("Selection Mode");
        for (const auto& sr : result.site_results) {
            bool is_pos = sr.prob_pos >= settings.posterior_threshold;
            bool is_neg = sr.prob_neg >= settings.posterior_threshold;

            std::string mode_str = is_pos ? Console::danger("Diversifying") : (is_neg ? Console::info("Purifying") : Console::muted("Neutral"));
            all_tbl.add_row({
                std::to_string(sr.site),
                format_float_val(sr.alpha, 3),
                format_float_val(sr.beta, 3),
                format_float_val(sr.beta_minus_alpha, 3),
                format_float_val(sr.prob_neg, 4),
                format_float_val(sr.prob_pos, 4),
                format_float_val(sr.bayes_factor, 2),
                mode_str
            });
        }
        std::cout << Console::bold("All Codon Sites:\n");
        all_tbl.print();
        std::cout << "\n";
    } else {
        FUBARAnalyzer::print_results_table(result, settings);
    }

    // 4. Save JSON Output
    std::string cli_cmd = "hyphy3 fubar";
    for (int i = 1; i < argc; ++i) {
        cli_cmd += " ";
        cli_cmd += argv[i];
    }

    try {
        std::string json_str = fubar.to_json(result, json_format, cli_cmd, alignment_file);
        std::ofstream out(output_file);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open output file: " + output_file);
        }
        out << json_str;
        std::cout << Console::bold("Output:\n")
                  << "  " << Console::muted("Analysis results saved to: ") << Console::brand(output_file) << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error writing JSON output: ") << e.what() << "\n";
    }

    if (!html_file.empty()) {
        try {
            hyphy::report::write_html_report(output_file, html_file);
            std::cout << Console::success("✔") << " " << Console::bold("Saved HTML report to: ")
                      << Console::brand(html_file) << "\n\n";
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error writing HTML report: ") << e.what() << "\n";
            return 1;
        }
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_fubar(argc, argv);
}
#endif
