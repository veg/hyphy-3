#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/analyses/relax.hpp"

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

static void print_relax_banner() {
    std::cout << "\n=======================================================\n"
              << "   HYPHY 3: RELAX (Test for Selection Relaxation)       \n"
              << "=======================================================\n"
              << " Citation: RELAX: Detecting Relaxed Selection in a      \n"
              << "           Phylogenetic Framework                       \n"
              << "           (2015), Mol Biol Evol. 32(3): 820-832        \n"
              << " Version:  3.0.0 (Modern C++20 Core)                   \n"
              << "=======================================================\n\n";
}

static void print_relax_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [OPTIONS]\n\n"
              << "Required arguments:\n"
              << "  --alignment <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << "Optional arguments:\n"
              << "  --tree <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  --code <name>        Genetic code (default: Universal)\n"
              << "  --test <regex|names> Branch names or regex to mark as Test set (default: use {T} tags)\n"
              << "  --pvalue <threshold> Significance threshold for LRT (default: 0.05)\n"
              << "  --threads <N>        Number of OpenMP worker threads\n"
              << "  --output <file>      Path to output JSON file (default: <alignment>.RELAX.json)\n"
              << "  --progress           Force interactive progress bar\n"
              << "  --no-progress        Disable progress bar\n"
              << "  --help, -h           Show this help message\n\n"
              << "Examples:\n"
              << "  " << prog << " --alignment data/Fig4E.nex\n"
              << "  " << prog << " --alignment data/gene.fna --tree data/gene.nwk --test \".*PSEUDOGENE.*\" --threads 8\n\n";
}

int run_relax(int argc, char* argv[]) {
    std::string alignment_file;
    std::string tree_file;
    std::string output_file;
    std::string code_name = "Universal";
    int num_threads = 0;
    RELAXSettings settings;
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
        } else if (arg == "--test" && i + 1 < argc) {
            settings.test_branch_regex = argv[++i];
        } else if (arg == "--pvalue" && i + 1 < argc) {
            settings.p_value_threshold = std::stod(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--progress") {
            show_progress = true;
            force_progress = true;
        } else if (arg == "--no-progress") {
            show_progress = false;
            force_progress = false;
        } else if (arg == "--help" || arg == "-h") {
            print_relax_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_relax_usage(argv[0]);
            return 1;
        }
    }

    if (alignment_file.empty()) {
        std::cerr << "Error: --alignment is required.\n\n";
        print_relax_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
#endif

    if (output_file.empty()) {
        output_file = alignment_file + ".RELAX.json";
    }

    print_relax_banner();

    // 1. Load Genetic Code
    std::shared_ptr<const GeneticCode> code;
    try {
        code = GeneticCode::from_name(code_name);
    } catch (const std::exception& e) {
        std::cerr << "Error: Unknown genetic code '" << code_name << "'.\n";
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
    std::cout << "      Taxa: " << aln.num_taxa << ", Codons: " << aln.num_codons
              << ", Unique Patterns: " << aln.patterns.size() << "\n";

    // 3. Load Tree
    std::cout << "[2/4] Resolving phylogeny and branch annotations...\n";
    Tree tree;
    if (!tree_file.empty()) {
        try {
            tree = Tree::from_newick_file(tree_file);
        } catch (const std::exception& e) {
            std::cerr << "Error loading tree from file: " << e.what() << "\n";
            return 1;
        }
    } else if (!aln.embedded_tree_newick.empty()) {
        try {
            tree = Tree::from_newick(aln.embedded_tree_newick);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing embedded tree: " << e.what() << "\n";
            return 1;
        }
    } else {
        std::cerr << "Error: No tree provided and no tree embedded in alignment.\n";
        return 1;
    }

    // 4. Initialize and Run RELAX
    auto relax = RELAXAnalyzer::create(tree, aln, settings);
    std::cout << "      Nodes: " << tree.num_nodes() << ", Leaves: " << tree.num_leaves() << "\n";
    std::cout << "      Test branches: " << relax.test_branch_names.size()
              << ", Reference branches: " << relax.ref_branch_names.size() << "\n";

    if (relax.test_branch_names.empty()) {
        std::cerr << "Error: No Test branches found! Annotate branches with {T} in Newick or pass --test <regex>.\n";
        return 1;
    }

    std::cout << "[3/4] Running RELAX model inference pipeline...\n";

    std::function<void(const std::string&, double)> progress_cb = nullptr;
    if (!show_progress && !force_progress) {
        progress_cb = [](const std::string& stage, double frac) {
            std::cout << "      [" << std::setw(3) << static_cast<int>(frac * 100) << "%] " << stage << std::endl;
        };
    }

    auto res = relax.run(progress_cb, show_progress, force_progress);

    // 5. Display Summary
    std::cout << "\n=======================================================\n"
              << "   RELAX Model Fitting Summary                         \n"
              << "=======================================================\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Nucleotide GTR Log-L    : " << std::setw(9) << res.gtr_log_likelihood
              << " (AICc: " << res.gtr_aicc << ", params: " << res.gtr_parameters << ")\n";
    std::cout << "  MG94xREV Sep Rates Log-L: " << std::setw(9) << res.mg94_log_likelihood
              << " (AICc: " << res.mg94_aicc << ", params: " << res.mg94_parameters << ")\n";
    std::cout << "    omega (Reference)     : " << std::setprecision(4) << res.mg94_omega_R << "\n";
    std::cout << "    omega (Test)          : " << std::setprecision(4) << res.mg94_omega_T << "\n";
    std::cout << "  RELAX Alternative Log-L : " << std::setprecision(2) << std::setw(9) << res.alternative_fit.log_likelihood
              << " (AICc: " << res.alternative_fit.aicc << ", params: " << res.alternative_fit.parameters << ")\n";
    std::cout << "  RELAX Null (K=1) Log-L  : " << std::setprecision(2) << std::setw(9) << res.null_fit.log_likelihood
              << " (AICc: " << res.null_fit.aicc << ", params: " << res.null_fit.parameters << ")\n";
    std::cout << "=======================================================\n\n";

    // Rate distributions
    std::cout << "### Inferred Rate Distributions (Alternative Model)\n\n";
    std::cout << "| Class | Reference omega | Reference Prop | Test omega | Test Prop |\n";
    std::cout << "| :---: | :---: | :---: | :---: | :---: |\n";
    for (size_t c = 0; c < res.alternative_fit.reference_distribution.omegas.size(); ++c) {
        std::cout << "| " << c
                  << " | " << std::fixed << std::setprecision(4) << std::setw(15) << res.alternative_fit.reference_distribution.omegas[c]
                  << " | " << std::fixed << std::setprecision(4) << std::setw(14) << res.alternative_fit.reference_distribution.weights[c]
                  << " | " << std::fixed << std::setprecision(4) << std::setw(10) << res.alternative_fit.test_distribution.omegas[c]
                  << " | " << std::fixed << std::setprecision(4) << std::setw(9) << res.alternative_fit.test_distribution.weights[c]
                  << " |\n";
    }

    std::cout << "\n### Hypothesis Test for Selection Relaxation\n\n";
    std::cout << "  Relaxation parameter (K) : " << std::fixed << std::setprecision(4) << res.k << "\n";
    std::cout << "  Likelihood Ratio Test    : " << std::fixed << std::setprecision(2) << res.lrt << "\n";
    std::cout << "  Asymptotic p-value       : " << std::scientific << std::setprecision(4) << res.p_value << "\n";
    std::cout << "  Significance threshold   : " << std::fixed << std::setprecision(2) << settings.p_value_threshold << "\n\n";

    if (res.is_significant) {
        if (res.is_relaxed) {
            std::cout << "  **Result: Significant RELAXATION of selection on Test branches** (K = "
                      << std::fixed << std::setprecision(4) << res.k << " < 1, p = "
                      << std::scientific << std::setprecision(3) << res.p_value << " <= "
                      << std::fixed << std::setprecision(2) << settings.p_value_threshold << ").\n";
        } else {
            std::cout << "  **Result: Significant INTENSIFICATION of selection on Test branches** (K = "
                      << std::fixed << std::setprecision(4) << res.k << " > 1, p = "
                      << std::scientific << std::setprecision(3) << res.p_value << " <= "
                      << std::fixed << std::setprecision(2) << settings.p_value_threshold << ").\n";
        }
    } else {
        std::cout << "  **Result: No significant evidence of relaxation or intensification** (p = "
                  << std::scientific << std::setprecision(3) << res.p_value << " > "
                  << std::fixed << std::setprecision(2) << settings.p_value_threshold << ").\n";
    }

    std::cout << "  Total execution time     : " << std::fixed << std::setprecision(2) << res.runtime_seconds << "s\n\n";

    // 6. Write JSON
    std::cout << "[4/4] Writing JSON results to: " << output_file << "\n";
    std::ofstream out(output_file);
    if (!out) {
        std::cerr << "Warning: Could not open output file for writing: " << output_file << "\n";
    } else {
        out << res.to_json(tree, aln).dump(2) << "\n";
        std::cout << "      Saved successfully.\n";
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_relax(argc, argv);
}
#endif
