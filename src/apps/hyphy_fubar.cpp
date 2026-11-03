#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/fubar.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <numeric>

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
              << "  " << Console::brand("--tree, -t") << " <file>        Path to Newick tree file or 'auto' (optional; auto-inferred via TreeEngine if omitted or 'auto')\n"
              << "  " << Console::brand("--ohmic, --omhic") << "         Infer Ohmic network & Wheatstone bridges (empirical conductances)\n"
              << "  " << Console::brand("--hyphaeon") << "               Opt-in: HyphAeon 384D transformer conductances (needs Python aeon-core)\n"
              << "  " << Console::brand("--embeddings") << " <file>      Path to precomputed [N x D] taxon embeddings TSV/CSV for Ohmic\n"
              << "  " << Console::brand("--code") << " <name>            Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  " << Console::brand("--grid, -g") << " <N>           Number of grid points per dimension (default: 20; total D^2 = 400)\n"
              << "  " << Console::brand("--method, -m") << " <name>      Inference method: 'Variational-Bayes' (default), 'Collapsed-Gibbs', or 'Metropolis-Hastings'\n"
              << "  " << Console::brand("--concentration") << " <val>    Dirichlet prior concentration parameter (default: 0.5)\n"
              << "  " << Console::brand("--posterior, -p") << " <float>  Posterior probability cutoff for selection calls (default: 0.90)\n"
              << "  " << Console::brand("--non-zero") << "               Enforce non-zero synonymous rates (min 1e-3)\n"
              << "  " << Console::brand("--tree-scaling") << " <mode>    Tree scaling mode: 'grid' (default, HyPhy 2.5 standard) or 'mg94'\n"
              << "  " << Console::brand("--chains") << " <N>             Number of MCMC chains for MCMC methods (default: 5)\n"
              << "  " << Console::brand("--chain-length") << " <N>       MCMC chain length (default: 2000000)\n"
              << "  " << Console::brand("--burn-in") << " <N>            MCMC burn-in steps (default: 1000000)\n"
              << "  " << Console::brand("--samples") << " <N>            MCMC samples to draw (default: 100)\n"
              << "  " << Console::brand("--threads, -T") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output, -o") << " <file>      Path to output JSON file (default: <alignment>.FUBAR.json)\n"
              << "  " << Console::brand("--html") << " <file>            Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--json-format") << " <fmt>      JSON schema format: modern (default) or legacy (HyPhy 2.5)\n"
              << "  " << Console::brand("--legacy-json") << "            Output legacy HyPhy 2.5 JSON schema\n"
              << "  " << Console::brand("--all-sites") << "              Display all codon sites in console table (default: significant only)\n"
              << "  " << Console::brand("--progress") << "               Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "            Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "               Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " --alignment data/cd2.fna --tree auto --output cd2.fubar.json\n"
              << "  " << prog << " --alignment data/cd2.fna --grid 20 --method Variational-Bayes\n"
              << "  " << prog << " --alignment test.nex --posterior 0.95 --output test.fubar.json\n\n";
}

int run_fubar(int argc, char* argv[]) {
    CommonCLIOptions cli;
    FUBARSettings settings;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto status = parse_common_cli_arg(i, argc, argv, cli);
        if (status == CommonArgParseStatus::Handled) continue;
        if (status == CommonArgParseStatus::Error) return 1;

        if ((arg == "--grid" || arg == "-g") && i + 1 < argc) {
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
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << Console::danger("Unknown argument: ") << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (cli.alignment_file.empty()) {
        std::cerr << Console::danger("Error: Missing required argument --alignment.\n\n");
        print_usage(argv[0]);
        return 1;
    }

    finalize_common_cli_options(cli, ".FUBAR.json");
    print_banner();

    LoadedCodonDataset ds;
    if (!load_codon_dataset(cli, ds, false)) {
        return 1;
    }
    auto& gcode = ds.gcode;
    auto& aln = ds.aln;
    auto& tree = ds.tree;
    const auto& tree_source = ds.tree_source;

    Panel::print_card("Dataset & Configuration", {
        {"Alignment File", cli.alignment_file},
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
        {"OpenMP Worker Threads", std::to_string((cli.num_threads > 0) ? cli.num_threads : omp_get_max_threads())}
#else
        {"OpenMP Worker Threads", "1 (Serial Build)"}
#endif
    });

    // 2. Initialize and Fit
    auto progress_cb = [&](const std::string& msg, double pct) {
        if (!cli.show_progress) return;
        std::cout << Console::brand("==> ") << Console::bold(msg) 
                  << Console::muted(" [" + std::to_string(static_cast<int>(pct * 100)) + "%]") << "\n";
    };

    FUBARAnalyzer fubar = FUBARAnalyzer::create_and_fit(tree, aln, settings, progress_cb);
    FUBARResult result = fubar.run(progress_cb, cli.show_progress);

    // 3. Console Presentation
    FUBARAnalyzer::print_summary(result, settings);

    if (cli.all_sites) {
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

    // 4. Save JSON Output (and optional HTML)
    std::string cli_cmd = join_cli_command(argc, argv);
    std::string json_format_str = (cli.json_format == JSONFormat::Legacy) ? "legacy" : "modern";
    std::string json_str = fubar.to_json(result, json_format_str, cli_cmd, cli.alignment_file);
    if (!write_cli_json_and_html(cli, json_str)) {
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_fubar(argc, argv);
}
#endif
