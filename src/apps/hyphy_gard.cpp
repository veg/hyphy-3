#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/analyses/gard.hpp"

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
              << "  " << Console::brand("--population-size, -p") << " <N>     Population size for CHC GA (default: 32)\n"
              << "  " << Console::brand("--max-stagnant") << " <N>              Max stagnant generations before concluding search (default: 30)\n"
              << "  " << Console::brand("--mutation-rate") << " <float>         Cataclysmic mutation rate (default: 0.35)\n"
              << "  " << Console::brand("--tree-engine") << " <ohmic|nj>        Tree construction method (default: ohmic)\n"
              << "  " << Console::brand("--srv") << " <gamma|none>               Site-to-site rate variation (default: gamma)\n"
              << "  " << Console::brand("--gamma-classes") << " <N>              Number of Gamma rate classes (default: 4)\n"
              << "  " << Console::brand("--bridges") << " <yes|no>               Screen for Wheatstone bridges & seed breakpoints (default: yes)\n"
              << "  " << Console::brand("--topology-test") << " <both|sh|kh|none> Topological incongruence test (default: both)\n"
              << "  " << Console::brand("--test-replicates") << " <N>             Bootstrap replicates for KH/SH tests (default: 1000)\n"
              << "  " << Console::brand("--threads, -t") << " <N>             Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output, -o") << " <file>           Path to output JSON file (default: <alignment>.GARD.json)\n"
              << "  " << Console::brand("--format") << " <fmt>                Output format: 'modern' (default) or 'legacy'\n"
              << "  " << Console::brand("--legacy-json") << "                  Output HyPhy 2.5 compatible JSON\n"
              << "  " << Console::brand("--no-chart") << "                     Disable dynamic breakpoint track & sparkline\n"
              << "  " << Console::brand("--force-tty") << "                    Force ANSI in-place dynamic chart updates\n"
              << "  " << Console::brand("--quiet, -q") << "                    Suppress progress output\n"
              << "  " << Console::brand("--help, -h") << "                    Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/adh.nex\n"
              << "  " << prog << " --alignment data/HRVI.nex --tree-engine ohmic --srv gamma\n"
              << "  " << prog << " --alignment data/HRVI.nex --mode Faster --threads 8\n\n";
}

int run_gard(int argc, char* argv[]) {
    std::string alignment_file;
    std::string output_file;
    std::string run_mode = "Normal";
    std::string output_format = "modern";
    size_t max_breakpoints = 10000;
    size_t population_size = 32;
    bool population_size_set = false;
    size_t max_stagnant = 30;
    bool max_stagnant_set = false;
    Scalar mutation_rate = 0.35;
    bool mutation_rate_set = false;
    int num_threads = 0;
    bool quiet = false;
    bool dynamic_chart = true;
    bool force_tty = false;

    GARDTreeMethod tree_method = GARDTreeMethod::Ohmic;
    GARDSRVMode srv_mode = GARDSRVMode::DiscreteGamma;
    size_t gamma_classes = 4;
    bool allow_bridges = true;
    GARDTopologyTestMode topology_test_mode = GARDTopologyTestMode::Both;
    size_t topology_test_replicates = 1000;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--alignment" || arg == "-a") && i + 1 < argc) {
            alignment_file = argv[++i];
        } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
            output_file = argv[++i];
        } else if ((arg == "--mode" || arg == "-m") && i + 1 < argc) {
            run_mode = argv[++i];
        } else if ((arg == "--max-breakpoints" || arg == "-b") && i + 1 < argc) {
            max_breakpoints = std::stoull(argv[++i]);
        } else if ((arg == "--population-size" || arg == "-p") && i + 1 < argc) {
            population_size = std::stoull(argv[++i]);
            population_size_set = true;
        } else if (arg == "--max-stagnant" && i + 1 < argc) {
            max_stagnant = std::stoull(argv[++i]);
            max_stagnant_set = true;
        } else if (arg == "--mutation-rate" && i + 1 < argc) {
            mutation_rate = std::stod(argv[++i]);
            mutation_rate_set = true;
        } else if (arg == "--tree-engine" && i + 1 < argc) {
            std::string eng = argv[++i];
            if (eng == "nj" || eng == "neighbor_joining") {
                tree_method = GARDTreeMethod::NeighborJoining;
            } else {
                tree_method = GARDTreeMethod::Ohmic;
            }
        } else if (arg == "--srv" && i + 1 < argc) {
            std::string s = argv[++i];
            if (s == "none" || s == "no" || s == "const") {
                srv_mode = GARDSRVMode::None;
            } else {
                srv_mode = GARDSRVMode::DiscreteGamma;
            }
        } else if (arg == "--gamma-classes" && i + 1 < argc) {
            gamma_classes = std::stoull(argv[++i]);
        } else if (arg == "--bridges" && i + 1 < argc) {
            std::string br = argv[++i];
            allow_bridges = (br == "yes" || br == "true" || br == "1");
        } else if (arg == "--no-bridges") {
            allow_bridges = false;
        } else if (arg == "--topology-test" && i + 1 < argc) {
            std::string t = argv[++i];
            if (t == "sh" || t == "SH") {
                topology_test_mode = GARDTopologyTestMode::SH;
            } else if (t == "kh" || t == "KH") {
                topology_test_mode = GARDTopologyTestMode::KH;
            } else if (t == "none" || t == "no" || t == "off") {
                topology_test_mode = GARDTopologyTestMode::None;
            } else {
                topology_test_mode = GARDTopologyTestMode::Both;
            }
        } else if (arg == "--test-replicates" && i + 1 < argc) {
            topology_test_replicates = std::stoull(argv[++i]);
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if ((arg == "--format" || arg == "--json-format") && i + 1 < argc) {
            std::string fmt = argv[++i];
            JSONFormat jf;
            if (!parse_json_format(fmt, jf)) {
                std::cerr << Console::danger("Error: unknown --json-format '") << fmt
                          << "'. Expected 'modern' or 'legacy'.\n";
                return 1;
            }
            output_format = (jf == JSONFormat::Legacy) ? "legacy" : "modern";
        } else if (arg == "--legacy-json") {
            output_format = "legacy";
        } else if (arg == "--modern-json") {
            output_format = "modern";
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
    if (population_size_set) settings.population_size = population_size;
    if (max_stagnant_set) settings.max_stagnant_generations = max_stagnant;
    if (mutation_rate_set) settings.mutation_rate = mutation_rate;
    settings.tree_method = tree_method;
    settings.srv_mode = srv_mode;
    settings.gamma_rate_classes = gamma_classes;
    settings.allow_bridges = allow_bridges;
    settings.topology_test_mode = topology_test_mode;
    settings.topology_test_replicates = topology_test_replicates;
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

            std::string status_conclusion;
            bool rec_detected = false;
            if (result.topology_test_result.run) {
                if (result.topology_test_result.any_incongruence_significant) {
                    status_conclusion = Console::success("✔ Recombination detected: Significant topological incongruence across segments (KH/SH p < 0.05).");
                    rec_detected = true;
                } else if (result.topological_incongruence) {
                    status_conclusion = Console::warning("⚠ Breakpoints improve c-AIC, but topological incongruence is not statistically significant (may reflect rate variation).");
                } else {
                    status_conclusion = Console::muted("ℹ No significant topological incongruence detected across segments.");
                }
            } else {
                if (result.topological_incongruence) {
                    status_conclusion = Console::success("✔ Recombination detected: Multiple trees improve c-AIC over single tree.");
                    rec_detected = true;
                } else if (result.ohmic_bridges_detected && result.ohmic_delta_aicc > 2.0) {
                    status_conclusion = Console::success("✔ Reticulation parsimoniously captured by Wheatstone bridge without data fragmentation!");
                } else {
                    status_conclusion = Console::muted("ℹ No significant topological incongruence detected across segments.");
                }
            }

            std::string tree_engine_str = (result.tree_method == GARDTreeMethod::Ohmic)
                ? "Ohmic Spectral Scaffolding (O(N log N))"
                : "Neighbor-Joining (O(N^3))";

            std::string srv_str = (result.srv_mode == GARDSRVMode::DiscreteGamma)
                ? ("Discrete Gamma (" + std::to_string(result.gamma_rates.size()) + " classes, alpha = " + [&]{
                    std::ostringstream ss;
                    ss << std::fixed << std::setprecision(3) << result.gamma_alpha;
                    return ss.str();
                }() + ")")
                : "None (Constant rates)";

            std::string topo_test_str = (settings.topology_test_mode == GARDTopologyTestMode::Both) ? "Kishino-Hasegawa (KH) & Shimodaira-Hasegawa (SH)"
                : (settings.topology_test_mode == GARDTopologyTestMode::SH) ? "Shimodaira-Hasegawa (SH)"
                : (settings.topology_test_mode == GARDTopologyTestMode::KH) ? "Kishino-Hasegawa (KH)" : "None";
            if (settings.topology_test_mode != GARDTopologyTestMode::None) {
                topo_test_str += " (" + std::to_string(settings.topology_test_replicates) + " RELL reps)";
            }

            std::vector<std::pair<std::string, std::string>> summary_items = {
                {"Alignment", alignment_file + " (" + std::to_string(result.num_taxa) + " taxa, " + std::to_string(result.num_sites) + " sites)"},
                {"Tree Search Engine", tree_engine_str},
                {"Rate Variation", srv_str},
                {"Topology Incongruence Test", topo_test_str},
                {"Baseline c-AIC (0 BPs)", base_ss.str()},
                {"Best Model c-AIC", best_ss.str() + " (" + std::to_string(result.best_model.breakpoints.size()) + " breakpoints)"},
                {"Delta c-AIC Improvement", delta_ss.str()},
                {"Breakpoints (1-based)", bp_list},
                {"Single-Tree Constrained c-AIC", single_ss.str()},
                {"Topological Incongruence", result.topological_incongruence ? Console::success("SUPPORTED") : Console::muted("NOT SUPPORTED")},
                {"Total Models Evaluated", std::to_string(result.total_models_evaluated) + " in " + std::to_string(static_cast<int>(result.time_elapsed + 0.5)) + "s"},
                {"Output Results File", Console::brand(output_file)}
            };

            if (result.ohmic_bridges_detected) {
                std::ostringstream ohm_ss, ohm_d_ss;
                ohm_ss << std::fixed << std::setprecision(2) << result.ohmic_aicc;
                ohm_d_ss << std::fixed << std::setprecision(2) << result.ohmic_delta_aicc;
                summary_items.insert(summary_items.begin() + 4, {
                    "Ohmic Network (Bridges)", ohm_ss.str() + " (Δ: " + ohm_d_ss.str() + ", " + std::to_string(result.num_bridges) + " bridge(s))"
                });
                if (!result.bridge_candidate_breakpoints.empty()) {
                    std::string seeds_str;
                    for (size_t k = 0; k < result.bridge_candidate_breakpoints.size(); ++k) {
                        if (k > 0) seeds_str += ", ";
                        seeds_str += std::to_string(result.bridge_candidate_breakpoints[k] + 1);
                    }
                    summary_items.insert(summary_items.begin() + 5, {
                        "Bridge Seeded Breakpoints", seeds_str
                    });
                }
            }

            Panel::print_summary_card(
                "GARD Analysis Results",
                summary_items,
                status_conclusion,
                rec_detected || (result.ohmic_bridges_detected && result.ohmic_delta_aicc > 2.0)
            );

            if (result.topology_test_result.run && !result.topology_test_result.breakpoint_tests.empty()) {
                Table table;
                table.add_column("Breakpoint", Table::Align::Right, 10);
                table.add_column("Left Segment", Table::Align::Center, 14);
                table.add_column("Right Segment", Table::Align::Center, 14);
                table.add_column("ΔlnL (L)", Table::Align::Right, 10);
                table.add_column("ΔlnL (R)", Table::Align::Right, 10);
                table.add_column("KH p (L / R)", Table::Align::Center, 16);
                table.add_column("SH p (L / R)", Table::Align::Center, 16);
                table.add_column("Significant", Table::Align::Center, 14);

                for (const auto& bp : result.topology_test_result.breakpoint_tests) {
                    std::ostringstream dl_ss, dr_ss, kh_ss, sh_ss;
                    dl_ss << std::fixed << std::setprecision(2) << bp.delta_lnl_left;
                    dr_ss << std::fixed << std::setprecision(2) << bp.delta_lnl_right;
                    kh_ss << std::fixed << std::setprecision(4) << bp.kh_p_value_left << " / " << bp.kh_p_value_right;
                    sh_ss << std::fixed << std::setprecision(4) << bp.sh_p_value_left << " / " << bp.sh_p_value_right;

                    std::string left_seg = "[" + std::to_string(bp.left_interval.first) + ".." + std::to_string(bp.left_interval.second) + "]";
                    std::string right_seg = "[" + std::to_string(bp.right_interval.first) + ".." + std::to_string(bp.right_interval.second) + "]";

                    std::string sig_str = bp.statistically_significant 
                        ? Console::success("Yes (p < 0.05)") 
                        : Console::muted("No");

                    table.add_row({
                        std::to_string(bp.breakpoint_site),
                        left_seg,
                        right_seg,
                        dl_ss.str(),
                        dr_ss.str(),
                        kh_ss.str(),
                        sh_ss.str(),
                        sig_str
                    });
                }

                std::cout << Console::bold(" Post-Hoc Topological Incongruence Tests (RELL Bootstrap):\n");
                table.print();
                std::cout << "\n";
            }
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
