#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/analyses/absrel.hpp"

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

static void print_absrel_banner() {
    std::cout << "\n=======================================================\n"
              << "   HYPHY 3: Adaptive Branch-Site Random Effects        \n"
              << "            Likelihood (aBSREL)                        \n"
              << "=======================================================\n"
              << " Citation: Less Is More: An Adaptive Branch-Site        \n"
              << "           Random Effects Model for Efficient          \n"
              << "           Detection of Episodic Diversifying Selection\n"
              << "           (2015), Mol Biol Evol. 32(5): 1342-1353     \n"
              << " Version:  3.0.0 (Modern C++20 Core)                   \n"
              << "=======================================================\n\n";
}

static void print_absrel_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [OPTIONS]\n\n"
              << "Required arguments:\n"
              << "  --alignment <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << "Optional arguments:\n"
              << "  --tree <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  --code <name>        Genetic code (default: Universal)\n"
              << "  --branches <set>     Branches to test: All, Internal, Leaves (default: All)\n"
              << "  --pvalue <threshold> Holm-Bonferroni p-value threshold (default: 0.05)\n"
              << "  --max-rates <N>      Maximum rate classes per branch (default: 3)\n"
              << "  --threads <N>        Number of OpenMP worker threads\n"
              << "  --output <file>      Path to output JSON file (default: <alignment>.ABSREL.json)\n"
              << "  --help, -h           Show this help message\n\n"
              << "Examples:\n"
              << "  " << prog << " --alignment data/bglobin.nex\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk --branches Internal --threads 8\n\n";
}

int run_absrel(int argc, char* argv[]) {
    std::string alignment_file;
    std::string tree_file;
    std::string output_file;
    std::string code_name = "Universal";
    int num_threads = 0;
    ABSRELSettings settings;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--alignment" && i + 1 < argc) {
            alignment_file = argv[++i];
        } else if (arg == "--tree" && i + 1 < argc) {
            tree_file = argv[++i];
        } else if (arg == "--code" && i + 1 < argc) {
            code_name = argv[++i];
        } else if (arg == "--branches" && i + 1 < argc) {
            settings.test_branches = argv[++i];
        } else if (arg == "--pvalue" && i + 1 < argc) {
            settings.p_threshold = std::stod(argv[++i]);
        } else if (arg == "--max-rates" && i + 1 < argc) {
            settings.max_rate_classes = std::stoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            print_absrel_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_absrel_usage(argv[0]);
            return 1;
        }
    }

    if (alignment_file.empty()) {
        std::cerr << "Error: --alignment is required.\n\n";
        print_absrel_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
#endif

    if (output_file.empty()) {
        output_file = alignment_file + ".ABSREL.json";
    }

    print_absrel_banner();

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
    std::cout << "[2/4] Resolving phylogeny...\n";
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
    std::cout << "      Nodes: " << tree.num_nodes() << ", Leaves: " << tree.num_leaves() << "\n";

    // 4. Initialize and Run aBSREL
    std::cout << "[3/4] Running aBSREL model inference...\n";
    auto absrel = ABSRELAnalyzer::create(tree, aln, settings);

    auto progress_cb = [](const std::string& stage, double frac) {
        std::cout << "      [" << std::setw(3) << static_cast<int>(frac * 100) << "%] " << stage << std::endl;
    };

    auto res = absrel.run(progress_cb);

    // 5. Display Summary
    std::cout << "\n=======================================================\n"
              << "   aBSREL Model Fitting Summary                        \n"
              << "=======================================================\n"
              << "  Nucleotide GTR Log-L    : " << std::fixed << std::setprecision(2) << res.gtr_fit.log_likelihood
              << " (AICc: " << res.gtr_fit.aicc << ")\n"
              << "  Baseline MG94xREV Log-L : " << res.baseline_fit.log_likelihood
              << " (AICc: " << res.baseline_fit.aicc << ")\n"
              << "  Full Adaptive Log-L     : " << res.full_adaptive_fit.log_likelihood
              << " (AICc: " << res.full_adaptive_fit.aicc << ")\n"
              << "  Tested Branches         : " << res.tested_branches.size() << "\n"
              << "  Positive Branches       : " << res.positive_branches.size() << "\n"
              << "  Execution Time          : " << std::setprecision(2) << res.runtime_seconds << "s\n"
              << "=======================================================\n\n";

    // Branch table
    std::cout << "### Tested Branches Summary\n\n";
    std::cout << "| Branch | Rates | Full Length | Test LRT | Uncorrected p | Corrected p | Sites @ EBF>=100 | Selected? |\n";
    std::cout << "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n";

    for (const auto& bname : res.tested_branches) {
        const auto& bres = res.branches.at(bname);
        std::cout << "| " << std::left << std::setw(15) << bname
                  << " | " << std::setw(5) << bres.rate_classes
                  << " | " << std::fixed << std::setprecision(4) << std::setw(11) << bres.full_branch_length
                  << " | " << std::setprecision(2) << std::setw(8) << bres.lrt
                  << " | " << std::scientific << std::setprecision(2) << std::setw(13) << bres.uncorrected_p_value
                  << " | " << std::setw(11) << bres.corrected_p_value
                  << " | " << std::defaultfloat << std::setw(16) << bres.sites_ebf_100
                  << " | " << (bres.is_positive ? "**Yes**" : "No") << " |\n";
    }

    if (!res.positive_branches.empty()) {
        std::cout << "\n**Positive selection detected** on " << res.positive_branches.size()
                  << " branches at Holm-Bonferroni p <= " << settings.p_threshold << ":\n";
        for (const auto& bname : res.positive_branches) {
            const auto& bres = res.branches.at(bname);
            std::cout << "  - **" << bname << "**: p_adj = " << std::scientific << std::setprecision(3)
                      << bres.corrected_p_value << " (LRT = " << std::fixed << std::setprecision(2)
                      << bres.lrt << ", " << bres.rate_classes << " rate classes, "
                      << bres.sites_ebf_100 << " sites with EBF >= 100)\n";
        }
    } else {
        std::cout << "\nNo branches were detected to be under positive selection at p <= "
                  << settings.p_threshold << ".\n";
    }

    // 6. Write JSON
    std::cout << "\n[4/4] Writing JSON results to: " << output_file << "\n";
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
    return run_absrel(argc, argv);
}
#endif
