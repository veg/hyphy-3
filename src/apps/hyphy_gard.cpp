#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/analyses/gard.hpp"
#include "hyphy/report/html_report.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;

static void print_gard_banner() {
    Panel::print_banner(
        "GARD",
        "Genetic Algorithm for Recombination Detection",
        "Screening for Recombination Breakpoints & Phylogenetic Incongruence",
        "Citation: Kosakovsky Pond et al. Mol Biol Evol 23(10), 1891–1901 • v3.0.0"
    );
}

static void print_gard_usage(const char* prog) {
    print_gard_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment, -a") << " <file>         Path to sequence alignment (NEXUS or FASTA)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--mode, -m") << " <mode>              Run mode: 'Normal' (default) or 'Faster'\n"
              << "  " << Console::brand("--max-breakpoints, -b") << " <N>     Maximum number of breakpoints (default: 10000)\n"
              << "  " << Console::brand("--threads, -t") << " <N>             Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output, -o") << " <file>           Path to output JSON file (default: <alignment>.GARD.json)\n"
              << "  " << Console::brand("--format") << " <fmt>                Output format: 'modern' (default) or 'legacy'\n"
              << "  " << Console::brand("--legacy-json") << "                  Output HyPhy 2.5 compatible JSON\n"
              << "  " << Console::brand("--html") << " <file>                  Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--no-chart") << "                     Disable dynamic breakpoint track & sparkline\n"
              << "  " << Console::brand("--force-tty") << "                    Force ANSI in-place dynamic chart updates\n"
              << "  " << Console::brand("--quiet, -q") << "                    Suppress progress output\n"
              << "  " << Console::brand("--help, -h") << "                    Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/adh.nex\n"
              << "  " << prog << " --alignment data/HRVI.nex --mode Faster --threads 8\n\n";
}

int run_gard(int argc, char* argv[]) {
    std::string alignment_file;
    std::string output_file;
    std::string html_file;
    std::string run_mode = "Normal";
    std::string output_format = "modern";
    size_t max_breakpoints = 10000;
    int num_threads = 0;
    bool quiet = false;
    bool dynamic_chart = true;
    bool force_tty = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--alignment" || arg == "-a") && i + 1 < argc) {
            alignment_file = argv[++i];
        } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--html" && i + 1 < argc) {
            html_file = argv[++i];
        } else if ((arg == "--mode" || arg == "-m") && i + 1 < argc) {
            run_mode = argv[++i];
        } else if ((arg == "--max-breakpoints" || arg == "-b") && i + 1 < argc) {
            max_breakpoints = std::stoull(argv[++i]);
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--format" && i + 1 < argc) {
            output_format = argv[++i];
        } else if (arg == "--legacy-json") {
            output_format = "legacy";
        } else if (arg == "--no-chart") {
            dynamic_chart = false;
        } else if (arg == "--chart") {
            dynamic_chart = true;
        } else if (arg == "--force-tty") {
            force_tty = true;
        } else if (arg == "--quiet" || arg == "-q") {
            quiet = true;
        } else if (arg == "--help" || arg == "-h") {
            print_gard_usage(argv[0]);
            return 0;
        }
    }

    if (alignment_file.empty()) {
        std::cerr << Console::danger("Error: --alignment argument is required.\n\n");
        print_gard_usage(argv[0]);
        return 1;
    }

    if (num_threads > 0) {
#ifdef _OPENMP
        omp_set_num_threads(num_threads);
#endif
    }

    if (output_file.empty()) {
        output_file = alignment_file + ".GARD.json";
    }

    if (!quiet) {
        print_gard_banner();
    }

    Alignment aln;
    try {
        aln = Alignment::load_nucleotide(alignment_file);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Failed to load alignment: " + std::string(e.what())) << "\n";
        return 1;
    }

    GARDSettings settings;
    settings.input_file = alignment_file;
    settings.output_file = output_file;
    settings.run_mode = run_mode;
    settings.max_breakpoints = max_breakpoints;
    settings.quiet = quiet;
    settings.dynamic_chart = dynamic_chart;
    settings.force_tty = force_tty;
    settings.legacy_json = (output_format == "legacy");

    try {
        GARDAnalyzer analyzer(aln, settings);
        GARDResult result = analyzer.run();

        nlohmann::json json_out;
        if (settings.legacy_json) {
            json_out = result.to_legacy_json();
        } else {
            json_out = result.to_json();
        }

        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cerr << Console::danger("Failed to open output file: " + output_file) << "\n";
            return 1;
        }
        out << json_out.dump(2) << "\n";
        out.close();

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

        if (!quiet) {
            std::string bp_list;
            if (result.best_model.breakpoints.empty()) {
                bp_list = "None (0 breakpoints)";
            } else {
                for (size_t b = 0; b < result.best_model.breakpoints.size(); ++b) {
                    if (b > 0) bp_list += ", ";
                    bp_list += std::to_string(result.best_model.breakpoints[b] + 1);
                }
            }

            std::ostringstream base_ss, best_ss, delta_ss, single_ss;
            base_ss << std::fixed << std::setprecision(2) << result.baseline_aicc;
            best_ss << std::fixed << std::setprecision(2) << result.best_model.aicc;
            Scalar delta_aicc = result.baseline_aicc - result.best_model.aicc;
            delta_ss << std::fixed << std::setprecision(2) << delta_aicc;
            single_ss << std::fixed << std::setprecision(2) << result.single_tree_aicc;

            std::string status_conclusion = result.topological_incongruence
                ? Console::success("✔ Recombination detected: Significant topological incongruence across segments.")
                : Console::muted("ℹ No significant topological incongruence detected across segments.");

            std::vector<std::pair<std::string, std::string>> summary_items = {
                {"Alignment", alignment_file + " (" + std::to_string(result.num_taxa) + " taxa, " + std::to_string(result.num_sites) + " sites)"},
                {"Baseline c-AIC (0 BPs)", base_ss.str()},
                {"Best Model c-AIC", best_ss.str() + " (" + std::to_string(result.best_model.breakpoints.size()) + " breakpoints)"},
                {"Delta c-AIC Improvement", delta_ss.str()},
                {"Breakpoints (1-based)", bp_list},
                {"Single-Tree Constrained c-AIC", single_ss.str()},
                {"Topological Incongruence", result.topological_incongruence ? Console::success("SUPPORTED") : Console::muted("NOT SUPPORTED")},
                {"Total Models Evaluated", std::to_string(result.total_models_evaluated) + " in " + std::to_string(static_cast<int>(result.time_elapsed + 0.5)) + "s"},
                {"Output Results File", Console::brand(output_file)}
            };

            Panel::print_summary_card(
                "GARD Analysis Results",
                summary_items,
                status_conclusion,
                result.topological_incongruence
            );
        }

    } catch (const std::exception& e) {
        std::cerr << Console::danger("GARD analysis error: " + std::string(e.what())) << "\n";
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_gard(argc, argv);
}
#endif
