#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/busted.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;
using namespace hyphy::opt;

static void print_banner() {
    Panel::print_banner(
        "BUSTED",
        "Branch-site Unrestricted Statistical Test",
        "Gene-wide identification of episodic diversifying positive selection",
        "Citation: Mol Biol Evol. 32: 1365-1371 (2015) • v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal)\n"
              << "  " << Console::brand("--rates") << " <N>          Number of omega rate categories (default: 3)\n"
              << "  " << Console::brand("--srv") << "                Enable synonymous rate variation across sites (BUSTED-S)\n"
              << "  " << Console::brand("--syn-rates") << " <N>      Number of synonymous rate categories (default: 3)\n"
              << "  " << Console::brand("--auto-k") << "             Automatically select optimal K via AICc step-up\n"
              << "  " << Console::brand("--multiple-hits") << " <M>  Multi-nucleotide substitutions: None (default), Double, Double+Triple\n"
              << "  " << Console::brand("--no-branch-opt") << "      Disable individual branch length refinement (use proportional scaling)\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.BUSTED.json)\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " --alignment tests/data/adh.nex --auto-k --threads 8\n\n";
}

int run_busted(int argc, char* argv[]) {
    std::string alignment_file;
    std::string tree_file;
    std::string output_file;
    std::string code_name = "Universal";
    int num_threads = 0;
    JSONFormat json_format = JSONFormat::ModernV3;
    BUSTEDSettings settings;
    bool show_progress = ProgressBar::is_terminal();
    bool force_progress = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--alignment" && i + 1 < argc) {
            alignment_file = argv[++i];
        } else if (arg == "--tree" && i + 1 < argc) {
            tree_file = argv[++i];
        } else if (arg == "--code" && i + 1 < argc) {
            code_name = argv[++i];
        } else if (arg == "--rates" && i + 1 < argc) {
            settings.num_rate_classes = std::stoul(argv[++i]);
        } else if (arg == "--srv") {
            settings.srv = true;
        } else if (arg == "--syn-rates" && i + 1 < argc) {
            settings.srv = true;
            settings.num_syn_rate_classes = std::stoul(argv[++i]);
        } else if (arg == "--auto-k") {
            settings.auto_select_k = true;
        } else if (arg == "--multiple-hits" && i + 1 < argc) {
            settings.multiple_hits = argv[++i];
        } else if (arg == "--no-branch-opt") {
            settings.refine_branch_lengths = false;
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--json-format" && i + 1 < argc) {
            std::string fmt = argv[++i];
            if (fmt == "legacy" || fmt == "datamonkey") {
                json_format = JSONFormat::Legacy;
            } else {
                json_format = JSONFormat::ModernV3;
            }
        } else if (arg == "--legacy-json") {
            json_format = JSONFormat::Legacy;
        } else if (arg == "--modern-json") {
            json_format = JSONFormat::ModernV3;
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--progress") {
            show_progress = true;
            force_progress = true;
        } else if (arg == "--no-progress") {
            show_progress = false;
            force_progress = false;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (alignment_file.empty()) {
        std::cerr << "Error: --alignment is required.\n\n";
        print_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
#endif

    if (output_file.empty()) {
        output_file = alignment_file + ".BUSTED.json";
    }

    print_banner();

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string start_iso = Provenance::current_iso8601();
    std::string cli_cmd = argv[0];
    for (int i = 1; i < argc; ++i) {
        cli_cmd += " ";
        cli_cmd += argv[i];
    }

    // 1. Load Genetic Code
    std::shared_ptr<const GeneticCode> code;
    try {
        code = GeneticCode::from_name(code_name);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error loading genetic code '") << code_name << "': " << e.what() << "\n";
        return 1;
    }

    // 2. Load Alignment
    Alignment aln;
    try {
        aln = Alignment::load(alignment_file, code);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error loading alignment: ") << e.what() << "\n";
        return 1;
    }

    // 3. Load Tree
    Tree tree;
    if (!tree_file.empty()) {
        try {
            tree = Tree::from_newick_file(tree_file);
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error loading tree: ") << e.what() << "\n";
            return 1;
        }
    } else if (!aln.embedded_tree_newick.empty()) {
        try {
            tree = Tree::from_newick(aln.embedded_tree_newick);
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error parsing embedded tree: ") << e.what() << "\n";
            return 1;
        }
    } else {
        std::cerr << Console::danger("Error: No tree specified and no embedded tree found in alignment.\n");
        return 1;
    }

    Panel::print_card("Dataset & Model Settings", {
        {"Alignment File", alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons) + " (" + std::to_string(aln.num_codons * 3) + " nt)"},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree", std::to_string(tree.num_nodes()) + " nodes (" + std::to_string(tree.num_leaves()) + " leaves)"},
        {"Rate Classes (K)", std::to_string(settings.num_rate_classes) + (settings.auto_select_k ? " (Auto-K enabled)" : "")},
        {"Synonymous Variation", settings.srv ? "Enabled (" + std::to_string(settings.num_syn_rate_classes) + " classes)" : "Disabled"},
        {"Multiple Hits", settings.multiple_hits},
        {"Genetic Code", code->name}
    });

    // 4. Initialize and run BUSTED
    Panel::print_step(1, 2, "Fitting Baseline GTR & Global MG94 Models", "Exchangeability rates & branch lengths");
    BUSTEDAnalyzer analyzer = BUSTEDAnalyzer::create_and_fit(tree, aln);

    std::cout << "  " << Console::muted("Nucleotide GTR Log-L : ") << Console::bold(std::to_string(analyzer.gtr_log_l)) << "\n"
              << "  " << Console::muted("Global MG94 Log-L    : ") << Console::bold(std::to_string(analyzer.mg94_log_l))
              << Console::muted("  |  omega = ") << Console::brand(std::to_string(analyzer.mg94_omega)) << "\n\n";

    Panel::print_step(2, 2, "Running BUSTED Mixture Models", "ECM + SQUAREM acceleration");
    if (settings.auto_select_k) {
        std::cout << "  " << Console::muted("↳ Automatic model selection enabled (testing K = 1.." + std::to_string(settings.max_k) + ")...\n");
    }
    BUSTEDResult res = analyzer.run(settings, show_progress, force_progress);

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_runtime = std::chrono::duration<double>(end_time - start_time).count();

    // Summary table
    Table rate_table;
    rate_table.add_column("Category", Table::Align::Center, 8);
    rate_table.add_column("Unconstrained omega", Table::Align::Right, 20);
    rate_table.add_column("Unconstrained Weight", Table::Align::Right, 20);
    rate_table.add_column("Constrained omega", Table::Align::Right, 18);
    rate_table.add_column("Constrained Weight", Table::Align::Right, 18);

    size_t num_k = res.unconstrained.test_distribution.omegas.size();
    for (size_t k = 0; k < num_k; ++k) {
        std::ostringstream u_w_ss, u_p_ss, c_w_ss, c_p_ss;
        u_w_ss << std::fixed << std::setprecision(4) << res.unconstrained.test_distribution.omegas[k];
        u_p_ss << std::fixed << std::setprecision(4) << res.unconstrained.test_distribution.weights[k];
        c_w_ss << std::fixed << std::setprecision(4) << res.constrained.test_distribution.omegas[k];
        c_p_ss << std::fixed << std::setprecision(4) << res.constrained.test_distribution.weights[k];

        std::string row_color = "";
        if (k + 1 == num_k && res.unconstrained.test_distribution.omegas[k] > 1.0) {
            row_color = "\033[38;5;48m"; // Highlight positive selection class in emerald
        }

        rate_table.add_row({
            "omega_" + std::to_string(k + 1),
            u_w_ss.str(),
            u_p_ss.str(),
            c_w_ss.str(),
            c_p_ss.str()
        }, row_color);
    }

    std::cout << "\n" << Console::bold("  Inferred Rate Distributions:") << "\n\n";
    rate_table.print();
    std::cout << "\n";

    if (!res.unconstrained.test_distribution.syn_rates.empty() && res.unconstrained.test_distribution.syn_rates.size() > 1) {
        std::cout << Console::bold("  Synonymous Site-to-Site Rates (SRV):") << "\n";
        for (size_t m = 0; m < res.unconstrained.test_distribution.syn_rates.size(); ++m) {
            std::cout << "    alpha_" << (m + 1) << " = " << std::fixed << std::setprecision(4)
                      << res.unconstrained.test_distribution.syn_rates[m]
                      << " (weight = " << res.unconstrained.test_distribution.syn_weights[m] << ")\n";
        }
        std::cout << "\n";
    }

    if (res.settings.multiple_hits != "None") {
        std::cout << Console::bold("  Multi-hit Substitution Rates:") << "\n";
        std::cout << "    delta (double-hit rate) : " << std::fixed << std::setprecision(4) << res.unconstrained.delta
                  << " (fraction: " << res.unconstrained.frac_delta * 100.0 << "%)\n";
        if (res.settings.multiple_hits == "Double+Triple") {
            std::cout << "    psi (triple-hit rate)   : " << std::fixed << std::setprecision(4) << res.unconstrained.psi
                      << " (fraction: " << res.unconstrained.frac_psi * 100.0 << "%)\n";
        }
        std::cout << "\n";
    }

    std::ostringstream lrt_ss, p_ss, u_ss, c_ss, time_ss;
    lrt_ss << std::fixed << std::setprecision(4) << res.lrt;
    p_ss << std::scientific << std::setprecision(6) << res.p_value;
    u_ss << std::fixed << std::setprecision(2) << res.unconstrained.log_likelihood << " (AICc: " << res.unconstrained.aicc << ", K = " << res.optimal_k << ")";
    c_ss << std::fixed << std::setprecision(2) << res.constrained.log_likelihood << " (AICc: " << res.constrained.aicc << ", omega_" << res.optimal_k << " = 1.0)";
    time_ss << std::fixed << std::setprecision(2) << total_runtime << " s";

    std::string conclusion;
    bool is_sig = (res.p_value < 0.05);
    if (is_sig) {
        conclusion = "STATISTICALLY SIGNIFICANT evidence of episodic diversifying selection (LRT = " + lrt_ss.str() + ", p = " + p_ss.str() + " < 0.05)!";
    } else {
        conclusion = "No statistically significant evidence of episodic diversifying selection (p = " + p_ss.str() + " >= 0.05).";
    }

    std::string json_desc = (json_format == JSONFormat::Legacy) ? " [Legacy]" : " [Modern v3.0]";
    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Likelihood Ratio Test (LRT)", lrt_ss.str()},
        {"p-value (mixture distribution)", p_ss.str() + " [threshold: 0.05]"},
        {"Unconstrained Model Fit", u_ss.str()},
        {"Constrained Null Model Fit", c_ss.str()},
        {"Optimal Model Selection", "K = " + std::to_string(res.optimal_k) + " rate classes"},
        {"JSON Output File", output_file + json_desc},
        {"Total Execution Time", time_ss.str()}
    };

    Panel::print_summary_card("BUSTED Episodic Selection Analysis Summary", sum_items, conclusion, is_sig);

    // Build Provenance
    Provenance prov;
    prov.invocation.cli_command = cli_cmd;
    prov.invocation.working_directory = Provenance::get_cwd();
    prov.invocation.arguments["alignment"] = alignment_file;
    if (!tree_file.empty()) prov.invocation.arguments["tree"] = tree_file;
    prov.invocation.arguments["code"] = code_name;
    prov.invocation.arguments["rates"] = std::to_string(settings.num_rate_classes);
    prov.invocation.arguments["srv"] = settings.srv ? "true" : "false";
    prov.invocation.arguments["syn_rates"] = std::to_string(settings.num_syn_rate_classes);
    prov.invocation.arguments["auto_k"] = settings.auto_select_k ? "true" : "false";
    prov.invocation.arguments["multiple_hits"] = settings.multiple_hits;
    prov.invocation.arguments["refine_branch_lengths"] = settings.refine_branch_lengths ? "true" : "false";
    prov.invocation.arguments["json_format"] = (json_format == JSONFormat::Legacy) ? "legacy" : "modern_v3";

    size_t aln_sz = 0;
    std::string aln_hash = crypto::SHA256::hash_file(alignment_file, &aln_sz);
    prov.inputs["alignment"] = {alignment_file, "FASTA/NEXUS", aln_hash, aln_sz};

    if (!tree_file.empty()) {
        size_t tree_sz = 0;
        std::string tree_hash = crypto::SHA256::hash_file(tree_file, &tree_sz);
        prov.inputs["tree"] = {tree_file, "Newick", tree_hash, tree_sz};
    }

    prov.execution.start_time = start_iso;
    prov.execution.end_time = Provenance::current_iso8601();
    prov.execution.wall_time_seconds = total_runtime;
#ifdef _OPENMP
    prov.execution.cpu_threads = (num_threads > 0) ? num_threads : omp_get_max_threads();
#else
    prov.execution.cpu_threads = 1;
#endif
    prov.execution.hostname = Provenance::get_hostname();
    prov.execution.os = Provenance::detect_os();
    prov.execution.compiler = Provenance::detect_compiler();

    // Save JSON output
    try {
        analyzer.save_json(output_file, res, json_format, prov);
        std::string fmt_desc = (json_format == JSONFormat::Legacy) ? "Legacy Datamonkey JSON" : "Modern HyPhy v3.0 JSON";
        std::cout << Console::success("✔") << " " << Console::bold("Saved " + fmt_desc + " report to: ")
                  << Console::brand(output_file) << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Warning: Failed to write JSON output: ") << e.what() << "\n";
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_busted(argc, argv);
}
#endif
