#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/core/preferences.hpp"
#include "hyphy/analyses/expcm_optimizer.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/analyses/persite_optimizer.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;

static void print_banner() {
    Panel::print_banner(
        "ExpCM",
        "Experimentally-informed Codon Model",
        "Site-specific dN/dS (ω) analysis with LRT and FDR correction",
        "v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment FILE") << "     Path to codon alignment (FASTA or NEXUS)\n"
              << "  " << Console::brand("--tree FILE") << "          Path to Newick tree file\n"
              << "  " << Console::brand("--preferences FILE") << "   Path to site-specific preferences (CSV)\n\n"
              << Console::bold("Parameter options:") << "\n"
              << "  " << Console::brand("--beta VALUE") << "         Stringency parameter [default: auto-optimize]\n"
              << "  " << Console::brand("--omega VALUE") << "        dN/dS ratio [default: auto-optimize]\n"
              << "  " << Console::brand("--kappa VALUE") << "        Transition/transversion ratio [default: auto-optimize]\n"
              << "  " << Console::brand("--fix-beta") << "           Fix β to provided value\n"
              << "  " << Console::brand("--fix-omega") << "          Fix ω to provided value\n"
              << "  " << Console::brand("--fix-kappa") << "          Fix κ to provided value\n"
              << "  " << Console::brand("--optimize-branches") << "  Optimize branch lengths [default: fixed]\n\n"
              << Console::bold("Output options:") << "\n"
              << "  " << Console::brand("--output FILE") << "        Output JSON file [default: expcm.json]\n"
              << "  " << Console::brand("--omegabysite FILE") << "   Per-site ω results (TSV format) [optional]\n"
              << "  " << Console::brand("--fdr-level VALUE") << "    FDR threshold for significance [default: 0.05]\n\n"
              << Console::bold("Convergence options:") << "\n"
              << "  " << Console::brand("--gradient-tol VALUE") << "  Gradient tolerance [default: 1e-6]\n"
              << "  " << Console::brand("--relative-tol VALUE") << "  Relative tolerance [default: 1e-10]\n"
              << "  " << Console::brand("--max-iterations N") << "    Maximum iterations [default: 1000]\n"
              << "  " << Console::brand("--threads N") << "           Number of OpenMP threads (0 = auto)\n"
              << "  " << Console::brand("--verbose") << "            Verbose output\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data.fasta --tree data.nwk --preferences data.csv --output results.json\n"
              << "  " << prog << " --alignment data.fasta --tree data.nwk --preferences data.csv"
              << " --beta 3.367 --omega 0.979 --kappa 8.829 --fix-beta --fix-omega --fix-kappa\n"
              << "  " << prog << " --alignment data.fasta --tree data.nwk --preferences data.csv"
              << " --omegabysite results_omegabysite.txt --threads 4\n\n";
}

int run_expcm(int argc, char* argv[]) {
    std::string alignment_file;
    std::string tree_file;
    std::string preferences_file;
    std::string code_name = "Universal";
    std::string output_file = "expcm.json";
    std::string omegabysite_file;

    double beta_val = 0.0;
    double omega_val = 0.0;
    double kappa_val = 0.0;
    bool fix_beta = false;
    bool fix_omega = false;
    bool fix_kappa = false;
    bool optimize_branch_lengths = false;

    double fdr_level = 0.05;
    double gradient_tol = 1e-6;
    double relative_tol = 1e-10;
    int max_iterations = 1000;
    int num_threads = 0;
    bool verbose = false;

    try {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if ((arg == "--alignment" || arg == "-a") && i + 1 < argc) {
                alignment_file = argv[++i];
            } else if ((arg == "--tree" || arg == "-t") && i + 1 < argc) {
                tree_file = argv[++i];
            } else if (arg == "--preferences" && i + 1 < argc) {
                preferences_file = argv[++i];
            } else if ((arg == "--code" || arg == "-c") && i + 1 < argc) {
                code_name = argv[++i];
            } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
                output_file = argv[++i];
            } else if (arg == "--omegabysite" && i + 1 < argc) {
                omegabysite_file = argv[++i];
            } else if (arg == "--beta" && i + 1 < argc) {
                beta_val = parse_double_full("--beta", argv[++i]);
            } else if (arg == "--omega" && i + 1 < argc) {
                omega_val = parse_double_full("--omega", argv[++i]);
            } else if (arg == "--kappa" && i + 1 < argc) {
                kappa_val = parse_double_full("--kappa", argv[++i]);
            } else if (arg == "--fix-beta") {
                fix_beta = true;
            } else if (arg == "--fix-omega") {
                fix_omega = true;
            } else if (arg == "--fix-kappa") {
                fix_kappa = true;
            } else if (arg == "--optimize-branches") {
                optimize_branch_lengths = true;
            } else if (arg == "--fdr-level" && i + 1 < argc) {
                fdr_level = parse_probability("--fdr-level", argv[++i]);
            } else if (arg == "--gradient-tol" && i + 1 < argc) {
                gradient_tol = parse_double_full("--gradient-tol", argv[++i]);
            } else if (arg == "--relative-tol" && i + 1 < argc) {
                relative_tol = parse_double_full("--relative-tol", argv[++i]);
            } else if (arg == "--max-iterations" && i + 1 < argc) {
                max_iterations = static_cast<int>(parse_nonneg_ull(
                    "--max-iterations", argv[++i], std::numeric_limits<int>::max(), 1));
            } else if ((arg == "--threads" || arg == "-T") && i + 1 < argc) {
                num_threads = static_cast<int>(
                    parse_nonneg_ull("--threads", argv[++i], std::numeric_limits<int>::max()));
            } else if (arg == "--verbose") {
                verbose = true;
            } else if (arg == "--help" || arg == "-h") {
                print_usage(argv[0]);
                return 0;
            } else {
                std::cerr << Console::danger("Error: unknown or incomplete argument '") << arg
                          << "'. Use --help for usage.\n";
                return 1;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }

    if (alignment_file.empty() || preferences_file.empty()) {
        std::cerr << Console::danger("Error: ")
                  << "Required arguments missing (--alignment, --tree, --preferences)\n\n";
        print_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
        LikelihoodEngine::set_user_requested_threads(num_threads);
    } else if (const char* env_t = std::getenv("OMP_NUM_THREADS")) {
        int env_threads = std::atoi(env_t);
        if (env_threads > 0) {
            LikelihoodEngine::set_user_requested_threads(env_threads);
        }
    } else {
        int hw = omp_get_max_threads();
        omp_set_num_threads(std::min(hw, 8));
        LikelihoodEngine::set_user_requested_threads(0);
    }
#else
    if (num_threads > 0) {
        LikelihoodEngine::set_user_requested_threads(num_threads);
    }
#endif

    print_banner();

    auto start_time = std::chrono::high_resolution_clock::now();

    std::shared_ptr<const GeneticCode> gcode;
    try {
        gcode = GeneticCode::from_name(code_name);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error resolving genetic code: ") << e.what() << "\n";
        return 1;
    }

    Alignment aln;
    try {
        aln = Alignment::load(alignment_file, gcode);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error loading alignment: ") << e.what() << "\n";
        return 1;
    }

    if (verbose) {
        std::cout << "✓ Loaded alignment: " << alignment_file << "\n";
    }

    Tree tree;
    try {
        if (!tree_file.empty()) {
            tree = Tree::from_newick_file(tree_file);
        } else if (!aln.embedded_tree_newick.empty()) {
            tree = Tree::from_newick(aln.embedded_tree_newick);
        } else {
            std::cerr << Console::danger("Error: ")
                      << "Required arguments missing (--tree or embedded tree in alignment)\n\n";
            print_usage(argv[0]);
            return 1;
        }
        LikelihoodEngine::reconcile_tree_leaf_names(tree, aln);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error loading tree: ") << e.what() << "\n";
        return 1;
    }

    if (verbose) {
        std::cout << "✓ Loaded tree: " << (tree_file.empty() ? "<embedded>" : tree_file) << "\n";
    }

    Preferences prefs;
    try {
        prefs = Preferences::from_csv(preferences_file);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error loading preferences: ") << e.what() << "\n";
        return 1;
    }

    if (verbose) {
        std::cout << "✓ Loaded preferences: " << preferences_file << "\n";
    }

    // Use empirical nucleotide frequencies from the alignment if valid and non-zero
    Vector4 empirical_phi = {0.35, 0.20, 0.25, 0.20};
    if ((aln.nuc_frequencies.array() > 0.0).all() && std::fabs(aln.nuc_frequencies.sum() - 1.0) < 1e-6) {
        empirical_phi = aln.nuc_frequencies;
    }

    Panel::print_card("Dataset & Phylogeny", {
        {"Alignment File", alignment_file},
        {"Sequences / Taxa", std::to_string(aln.sequences.size())},
        {"Codon Sites", std::to_string(prefs.num_sites())},
        {"Tree File", (tree_file.empty() ? std::string("<embedded>") : tree_file) + " (" + std::to_string(tree.num_nodes()) + " nodes)"},
        {"Genetic Code", gcode->name}
    });

    double beta = beta_val;
    double omega = omega_val;
    double kappa = kappa_val;
    double global_lnl = 0.0;

    if (!fix_beta || !fix_omega || !fix_kappa) {
        Panel::print_step(1, 2, "Global Parameter Optimization",
                         "Optimizing β, ω, κ with L-BFGS");

        try {
            ExpCMOptimizer optimizer(tree, aln, *gcode, prefs);

            ExpCMParameters initial_params;
            initial_params.beta = beta_val > 0 ? beta_val : 3.0;
            initial_params.omega = omega_val > 0 ? omega_val : 1.0;
            initial_params.kappa = kappa_val > 0 ? kappa_val : 2.0;
            initial_params.phi = empirical_phi;

            ExpCMOptimizer::ConvergenceCriteria criteria;
            criteria.gradient_tolerance = gradient_tol;
            criteria.relative_tol = relative_tol;
            criteria.max_iterations = max_iterations;

            auto result = optimizer.optimize(initial_params, ExpCMOptimizer::Bounds(), criteria);

            Eigen::Vector3d real_params = result.get_real_params();
            beta = fix_beta && beta_val > 0 ? beta_val : real_params(0);
            omega = fix_omega && omega_val > 0 ? omega_val : real_params(1);
            kappa = fix_kappa && kappa_val > 0 ? kappa_val : real_params(2);
            global_lnl = result.optimal_likelihood;
            if (fix_beta || fix_omega || fix_kappa) {
                ExpCMParameters fixed_p(beta, omega, kappa, empirical_phi);
                global_lnl = LikelihoodEngine::compute_expcm_log_likelihood(
                    tree, aln, *gcode, prefs, fixed_p
                );
            }

            std::cout << "  " << Console::muted("Log-Likelihood: ") << Console::bold(std::to_string(global_lnl)) << "\n"
                      << "  " << Console::muted("Parameters: β=") << std::fixed << std::setprecision(4) << beta
                      << Console::muted(" ω=") << omega
                      << Console::muted(" κ=") << kappa << "\n"
                      << "  " << Console::muted("Iterations: ") << result.iterations << "\n\n";

            if (verbose) {
                std::cout << "✓ Global optimization completed\n";
            }
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error in global optimization: ") << e.what() << "\n";
            return 1;
        }
    } else {
        Panel::print_step(1, 2, "Computing Likelihood",
                         "With fixed parameters");

        try {
            ExpCMParameters params;
            params.beta = beta;
            params.omega = omega;
            params.kappa = kappa;
            params.phi = empirical_phi;

            global_lnl = LikelihoodEngine::compute_expcm_log_likelihood(
                tree, aln, *gcode, prefs, params
            );

            std::cout << "  " << Console::muted("Log-Likelihood: ") << Console::bold(std::to_string(global_lnl)) << "\n"
                      << "  " << Console::muted("Parameters: β=") << std::fixed << std::setprecision(4) << beta
                      << Console::muted(" ω=") << omega
                      << Console::muted(" κ=") << kappa << "\n\n";

            if (verbose) {
                std::cout << "✓ Likelihood computation completed\n";
            }
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error computing likelihood: ") << e.what() << "\n";
            return 1;
        }
    }

    if (optimize_branch_lengths) {
        try {
            ExpCMParameters params(beta, omega, kappa, empirical_phi);
            const double inv_phi = (std::sqrt(5.0) - 1.0) * 0.5;
            for (auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                double orig_bl = std::max(node.branch_length, 1e-6);
                double a = std::log(std::max(1e-8, orig_bl * 0.1));
                double b = std::log(std::min(10.0, std::max(orig_bl * 10.0, 0.5)));
                double c = b - inv_phi * (b - a);
                double d = a + inv_phi * (b - a);
                node.branch_length = std::exp(c);
                double fc = LikelihoodEngine::compute_expcm_log_likelihood(tree, aln, *gcode, prefs, params);
                node.branch_length = std::exp(d);
                double fd = LikelihoodEngine::compute_expcm_log_likelihood(tree, aln, *gcode, prefs, params);
                for (int iter = 0; iter < 12; ++iter) {
                    if (fc > fd) {
                        b = d;
                        d = c;
                        fd = fc;
                        c = b - inv_phi * (b - a);
                        node.branch_length = std::exp(c);
                        fc = LikelihoodEngine::compute_expcm_log_likelihood(tree, aln, *gcode, prefs, params);
                    } else {
                        a = c;
                        c = d;
                        fc = fd;
                        d = a + inv_phi * (b - a);
                        node.branch_length = std::exp(d);
                        fd = LikelihoodEngine::compute_expcm_log_likelihood(tree, aln, *gcode, prefs, params);
                    }
                }
                double cand_bl = std::exp(0.5 * (a + b));
                node.branch_length = cand_bl;
                double cand_ll = LikelihoodEngine::compute_expcm_log_likelihood(tree, aln, *gcode, prefs, params);
                if (cand_ll >= global_lnl) {
                    global_lnl = cand_ll;
                } else {
                    node.branch_length = orig_bl;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error optimizing branch lengths: ") << e.what() << "\n";
            return 1;
        }
    }

    Panel::print_step(2, 2, "Per-Site ω Analysis",
                     "Computing site-specific selection statistics");

    std::vector<PersiteOptimizer::SiteResult> site_results;
    try {
        site_results = PersiteOptimizer::optimize_all_sites(
            tree, aln, *gcode, prefs, beta, kappa, empirical_phi
        );

        PersiteOptimizer::apply_fdr_correction(site_results, fdr_level);

        int sig_sites = 0;
        for (const auto& result : site_results) {
            if (result.q_value < fdr_level) sig_sites++;
        }

        std::cout << "  " << Console::muted("Sites analyzed: ") << site_results.size() << "\n"
                  << "  " << Console::muted("Significant sites (q<") << fdr_level << "): " << sig_sites << "\n"
                  << "  " << Console::muted("ω range: [") << std::fixed << std::setprecision(4);

        double min_omega = site_results.empty() ? 0.0 : 1e10;
        double max_omega = 0.0;
        for (const auto& r : site_results) {
            min_omega = std::min(min_omega, r.omega_mle);
            max_omega = std::max(max_omega, r.omega_mle);
        }

        std::cout << min_omega << ", " << max_omega << "]\n\n";

        if (verbose) {
            std::cout << "✓ Per-site analysis completed\n";
        }
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error in per-site analysis: ") << e.what() << "\n";
        return 1;
    }

    auto ends_with_json = [](const std::string& path) -> bool {
        if (path.size() < 5) return false;
        std::string ext = path.substr(path.size() - 5);
        for (char& ch : ext) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        return ext == ".json";
    };

    try {
        if (ends_with_json(output_file)) {
            nlohmann::json j;
            j["method"] = "ExpCM";
            j["input"] = {
                {"alignment_file", alignment_file},
                {"tree_file", tree_file},
                {"preferences_file", preferences_file},
                {"num_sequences", aln.sequences.size()},
                {"num_sites", prefs.num_sites()},
                {"genetic_code", gcode->name}
            };
            j["fits"] = {
                {"log_likelihood", finite_or_null(global_lnl)},
                {"beta", finite_or_null(beta)},
                {"omega", finite_or_null(omega)},
                {"kappa", finite_or_null(kappa)},
                {"phi", {
                    finite_or_null(empirical_phi(0)),
                    finite_or_null(empirical_phi(1)),
                    finite_or_null(empirical_phi(2)),
                    finite_or_null(empirical_phi(3))
                }}
            };
            j["fdr_level"] = finite_or_null(fdr_level);

            nlohmann::json sites_arr = nlohmann::json::array();
            for (const auto& r : site_results) {
                double dlnl = r.lnl_alt - r.lnl_null;
                sites_arr.push_back({
                    {"site", r.site_id + 1},
                    {"omega", finite_or_null(r.omega_mle)},
                    {"log_likelihood_null", finite_or_null(r.lnl_null)},
                    {"log_likelihood_alt", finite_or_null(r.lnl_alt)},
                    {"delta_lnl", finite_or_null(dlnl)},
                    {"lrt", finite_or_null(r.test_statistic)},
                    {"p_value", finite_or_null(r.p_value)},
                    {"q_value", finite_or_null(r.q_value)}
                });
            }
            j["sites"] = std::move(sites_arr);
            write_file_or_throw(output_file, j.dump(2) + "\n");
        } else {
            PersiteOptimizer::write_omegabysite_file(site_results, output_file, true);
        }
        std::cout << Console::success("✓ ") << "Results written to: " << output_file << "\n";
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error writing output: ") << e.what() << "\n";
        return 1;
    }

    if (!omegabysite_file.empty()) {
        try {
            PersiteOptimizer::write_omegabysite_file(site_results, omegabysite_file, true);
            std::cout << Console::success("✓ ") << "Per-site results written to: " << omegabysite_file << "\n";
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error writing per-site file: ") << e.what() << "\n";
            return 1;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    Panel::print_card("Execution Summary", {
        {"Total Time", std::to_string(duration.count()) + " seconds"},
        {"Global ln L", std::to_string(global_lnl)},
        {"β", std::to_string(beta)},
        {"ω", std::to_string(omega)},
        {"κ", std::to_string(kappa)}
    });

    return 0;
}
