#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
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
    std::cout << "\n=======================================================\n"
              << "   HYPHY 3: Branch-site Unrestricted Statistical       \n"
              << "            Test for Episodic Diversification          \n"
              << "            (BUSTED)                                   \n"
              << "=======================================================\n"
              << " Citation: Gene-wide identification of episodic         \n"
              << "           selection (2015), Mol Biol Evol. 32: 1365-71\n"
              << " Version:  3.0.0 (Modern C++20 Core)                   \n"
              << "=======================================================\n\n";
}

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [OPTIONS]\n\n"
              << "Required arguments:\n"
              << "  --alignment <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << "Optional arguments:\n"
              << "  --tree <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  --code <name>        Genetic code (default: Universal)\n"
              << "  --rates <N>          Number of omega rate categories (default: 3)\n"
              << "  --srv                Enable synonymous rate variation across sites (BUSTED-S)\n"
              << "  --syn-rates <N>      Number of synonymous rate categories (default: 3)\n"
              << "  --auto-k             Automatically select optimal K via AICc step-up\n"
              << "  --multiple-hits <M>  Multi-nucleotide substitutions: None (default), Double, Double+Triple\n"
              << "  --no-branch-opt      Disable individual branch length refinement (use proportional scaling)\n"
              << "  --threads <N>        Number of OpenMP worker threads\n"
              << "  --output <file>      Path to output JSON file (default: <alignment>.BUSTED.json)\n"
              << "  --help, -h           Show this help message\n\n"
              << "Examples:\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " --alignment tests/data/adh.nex --auto-k --threads 8\n\n";
}

int run_busted(int argc, char* argv[]) {
    std::string alignment_file;
    std::string tree_file;
    std::string output_file;
    std::string code_name = "Universal";
    int num_threads = 0;
    BUSTEDSettings settings;

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
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
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

    // 1. Load Genetic Code
    std::shared_ptr<const GeneticCode> code;
    try {
        code = GeneticCode::from_name(code_name);
    } catch (const std::exception& e) {
        std::cerr << "Error loading genetic code '" << code_name << "': " << e.what() << "\n";
        return 1;
    }

    // 2. Load Alignment
    std::cout << "[1/4] Loading alignment from: " << alignment_file << "\n";
    Alignment aln;
    try {
        aln = Alignment::load(alignment_file, code);
    } catch (const std::exception& e) {
        std::cerr << "Error loading alignment: " << e.what() << "\n";
        return 1;
    }
    std::cout << "      Sequences: " << aln.num_taxa
              << " | Sites: " << aln.num_codons
              << " | Unique Patterns: " << aln.patterns.size() << "\n";

    // 3. Load Tree
    Tree tree;
    if (!tree_file.empty()) {
        std::cout << "[2/4] Loading tree from: " << tree_file << "\n";
        try {
            tree = Tree::from_newick_file(tree_file);
        } catch (const std::exception& e) {
            std::cerr << "Error loading tree: " << e.what() << "\n";
            return 1;
        }
    } else if (!aln.embedded_tree_newick.empty()) {
        std::cout << "[2/4] Using embedded tree from alignment file\n";
        try {
            tree = Tree::from_newick(aln.embedded_tree_newick);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing embedded tree: " << e.what() << "\n";
            return 1;
        }
    } else {
        std::cerr << "Error: No tree specified and no embedded tree found in alignment.\n";
        return 1;
    }
    std::cout << "      Taxa: " << tree.num_leaves()
              << " | Internal nodes: " << (tree.num_nodes() - tree.num_leaves()) << "\n";

    // 4. Initialize and run BUSTED
    std::cout << "[3/4] Fitting baseline Nucleotide GTR and Global MG94...\n";
    BUSTEDAnalyzer analyzer = BUSTEDAnalyzer::create_and_fit(tree, aln);

    std::cout << "      GTR Log-Likelihood : " << std::fixed << std::setprecision(2) << analyzer.gtr_log_l << "\n";
    std::cout << "      MG94 Log-Likelihood: " << std::fixed << std::setprecision(2) << analyzer.mg94_log_l
              << " (omega = " << std::setprecision(4) << analyzer.mg94_omega << ")\n";

    std::cout << "[4/4] Running BUSTED mixture models (ECM + SQUAREM)...\n";
    if (settings.auto_select_k) {
        std::cout << "      Automatic model selection enabled (testing K = 1.." << settings.max_k << ")...\n";
    }
    BUSTEDResult res = analyzer.run(settings);

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_runtime = std::chrono::duration<double>(end_time - start_time).count();

    // Summary output
    std::cout << "\n=======================================================\n"
              << "                   BUSTED RESULTS                      \n"
              << "=======================================================\n";
    if (settings.auto_select_k) {
        std::cout << "Selected Model: K = " << res.optimal_k << " rate categories (via AICc step-up)\n\n";
    }
    std::cout << "Unconstrained Model (K = " << res.optimal_k << "):\n"
              << "  Log-Likelihood : " << std::setprecision(2) << res.unconstrained.log_likelihood << "\n"
              << "  AIC-c          : " << res.unconstrained.aicc << "\n"
              << "  Tree Scale     : " << std::setprecision(4) << res.unconstrained.tree_scale << "\n";
    for (size_t k = 0; k < res.unconstrained.test_distribution.omegas.size(); ++k) {
        std::cout << "  omega_" << (k + 1) << " = " << std::setprecision(4) << res.unconstrained.test_distribution.omegas[k]
                  << " (p = " << std::setprecision(4) << res.unconstrained.test_distribution.weights[k] << ")";
        if (k < res.unconstrained.test_distribution.annotations.size() && !res.unconstrained.test_distribution.annotations[k].empty()) {
            std::cout << " [" << res.unconstrained.test_distribution.annotations[k] << "]";
        }
        std::cout << "\n";
    }
    if (!res.unconstrained.test_distribution.syn_rates.empty() && res.unconstrained.test_distribution.syn_rates.size() > 1) {
        std::cout << "  Synonymous Site-to-Site Rates (SRV):\n";
        for (size_t m = 0; m < res.unconstrained.test_distribution.syn_rates.size(); ++m) {
            std::cout << "    alpha_" << (m + 1) << " = " << std::setprecision(4) << res.unconstrained.test_distribution.syn_rates[m]
                      << " (p = " << std::setprecision(4) << res.unconstrained.test_distribution.syn_weights[m] << ")\n";
        }
    }
    if (res.settings.multiple_hits != "None") {
        std::cout << "  Multi-hit Substitutions:\n";
        std::cout << "    delta (double-hit rate) : " << std::setprecision(4) << res.unconstrained.delta
                  << " (fraction: " << std::setprecision(4) << res.unconstrained.frac_delta * 100.0 << "%)\n";
        if (res.settings.multiple_hits == "Double+Triple") {
            std::cout << "    psi (triple-hit rate)   : " << std::setprecision(4) << res.unconstrained.psi
                      << " (fraction: " << std::setprecision(4) << res.unconstrained.frac_psi * 100.0 << "%)\n";
        }
    }
    std::cout << "\n";

    std::cout << "Constrained Null Model (omega_" << res.optimal_k << " = 1.0):\n"
              << "  Log-Likelihood : " << std::setprecision(2) << res.constrained.log_likelihood << "\n"
              << "  AIC-c          : " << res.constrained.aicc << "\n";
    for (size_t k = 0; k < res.constrained.test_distribution.omegas.size(); ++k) {
        std::cout << "  omega_" << (k + 1) << " = " << std::setprecision(4) << res.constrained.test_distribution.omegas[k]
                  << " (p = " << std::setprecision(4) << res.constrained.test_distribution.weights[k] << ")";
        if (k < res.constrained.test_distribution.annotations.size() && !res.constrained.test_distribution.annotations[k].empty()) {
            std::cout << " [" << res.constrained.test_distribution.annotations[k] << "]";
        }
        std::cout << "\n";
    }
    if (res.settings.multiple_hits != "None") {
        std::cout << "  Multi-hit Substitutions:\n";
        std::cout << "    delta (double-hit rate) : " << std::setprecision(4) << res.constrained.delta
                  << " (fraction: " << std::setprecision(4) << res.constrained.frac_delta * 100.0 << "%)\n";
        if (res.settings.multiple_hits == "Double+Triple") {
            std::cout << "    psi (triple-hit rate)   : " << std::setprecision(4) << res.constrained.psi
                      << " (fraction: " << std::setprecision(4) << res.constrained.frac_psi * 100.0 << "%)\n";
        }
    }
    std::cout << "\n";

    std::cout << "Hypothesis Test for Episodic Diversifying Positive Selection:\n"
              << "  Likelihood Ratio Test (LRT) = " << std::setprecision(4) << res.lrt << "\n"
              << "  p-value                     = " << std::setprecision(6) << res.p_value << "\n";

    if (res.p_value < 0.05) {
        std::cout << "  Conclusion: Statistically significant evidence of episodic selection (p < 0.05)!\n";
    } else {
        std::cout << "  Conclusion: No statistically significant evidence of episodic selection.\n";
    }

    std::cout << "\nTotal runtime: " << std::setprecision(2) << total_runtime << " s\n"
              << "=======================================================\n\n";

    // Save JSON output
    try {
        auto j = analyzer.to_json(res);
        std::ofstream out(output_file);
        if (out.is_open()) {
            out << j.dump(2) << "\n";
            std::cout << "Saved Datamonkey JSON to: " << output_file << "\n";
        } else {
            std::cerr << "Warning: Could not open output file " << output_file << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to write JSON output: " << e.what() << "\n";
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_busted(argc, argv);
}
#endif
