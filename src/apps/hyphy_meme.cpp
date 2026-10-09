#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/meme.hpp"

#include <iostream>
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
              << "   HYPHY 3: Mixed Effects Model of Evolution (MEME)    \n"
              << "=======================================================\n"
              << " Citation: Detecting Individual Sites Subject to       \n"
              << "           Episodic Diversifying Selection (2012)      \n"
              << "           PLoS Genet 8(7): e1002764                   \n"
              << " Version:  3.0.0 (Modern C++20 Core)                   \n"
              << "=======================================================\n\n";
}

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [OPTIONS]\n\n"
              << "Required arguments:\n"
              << "  --alignment <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << "Optional arguments:\n"
              << "  --tree <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  --code <name>        Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  --threads <N>        Number of OpenMP worker threads\n"
              << "  --output <file>      Path to output JSON file (default: <alignment>.MEME.json)\n"
              << "  --pvalue <float>     P-value significance threshold (default: 0.1)\n"
              << "  --help, -h           Show this help message\n\n"
              << "Examples:\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk --output cd2.MEME.json\n"
              << "  " << prog << " --alignment tests/data/adh.nex --threads 8\n\n";
}

int run_meme(int argc, char* argv[]) {
    std::string alignment_file;
    std::string tree_file;
    std::string output_file;
    std::string code_name = "Universal";
    int num_threads = 0;
    Scalar pvalue_threshold = 0.1;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--alignment" && i + 1 < argc) {
            alignment_file = argv[++i];
        } else if (arg == "--tree" && i + 1 < argc) {
            tree_file = argv[++i];
        } else if (arg == "--code" && i + 1 < argc) {
            code_name = argv[++i];
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--pvalue" && i + 1 < argc) {
            pvalue_threshold = std::stod(argv[++i]);
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
        output_file = alignment_file + ".MEME.json";
    }

    print_banner();

    auto start_time = std::chrono::high_resolution_clock::now();

    // 0. Genetic Code
    std::shared_ptr<const GeneticCode> gcode;
    try {
        gcode = GeneticCode::from_name(code_name);
        std::cout << "> Using Genetic Code: " << gcode->name 
                  << " (" << gcode->num_sense_codons << " sense codons)\n";
    } catch (const std::exception& e) {
        std::cerr << "Error resolving genetic code: " << e.what() << "\n";
        return 1;
    }

    // 1. Load Alignment
    std::cout << "> Loading alignment from '" << alignment_file << "'...\n";
    Alignment aln;
    try {
        aln = Alignment::load(alignment_file, gcode);
    } catch (const std::exception& e) {
        std::cerr << "Error loading alignment: " << e.what() << "\n";
        return 1;
    }
    std::cout << "  Sequences: " << aln.num_taxa << "\n"
              << "  Codon Sites: " << aln.num_codons << "\n"
              << "  Unique Site Patterns: " << aln.patterns.size() << "\n\n";

    // 2. Load Tree
    Tree tree;
    if (!tree_file.empty()) {
        std::cout << "> Loading tree from '" << tree_file << "'...\n";
        try {
            tree = Tree::from_newick_file(tree_file);
        } catch (const std::exception& e) {
            std::cerr << "Error loading tree: " << e.what() << "\n";
            return 1;
        }
    } else if (!aln.embedded_tree_newick.empty()) {
        std::cout << "> Using embedded tree from NEXUS alignment...\n";
        tree = Tree::from_newick(aln.embedded_tree_newick);
    } else {
        std::cerr << "Error: No tree provided and no embedded tree found in alignment.\n";
        return 1;
    }
    std::cout << "  Tree leaves: " << tree.num_leaves() << "\n"
              << "  Total nodes: " << tree.num_nodes() << "\n\n";

    // 3. Multi-phase Global Fitting
    auto progress_cb = [](const std::string& msg, double pct) {
        std::cout << "> " << msg << "... (" << std::fixed << std::setprecision(0) << pct * 100 << "%)\n";
    };

    auto fit_start = std::chrono::high_resolution_clock::now();
    auto analyzer = MEMEAnalyzer::create_and_fit(tree, std::move(aln), pvalue_threshold, progress_cb);
    auto fit_end = std::chrono::high_resolution_clock::now();
    double fit_duration = std::chrono::duration<double, std::milli>(fit_end - fit_start).count();

    std::cout << "\n  --- Global Model Fit Summary ---\n"
              << "  GTR Log-Likelihood:  " << std::fixed << std::setprecision(2) << analyzer.gtr_log_l << "\n"
              << "  MG94 Log-Likelihood: " << std::fixed << std::setprecision(2) << analyzer.global_log_l << "\n"
              << "  Global omega (dN/dS): " << std::fixed << std::setprecision(4) 
              << (analyzer.base_params.beta / analyzer.base_params.alpha) << "\n"
              << "  Fitting time: " << std::setprecision(1) << fit_duration << " ms\n\n";

    // 4. Site-level MEME testing
    std::cout << "> Running site-by-site MEME testing across " 
              << analyzer.aln.patterns.size() << " unique patterns...\n";
    auto run_start = std::chrono::high_resolution_clock::now();
    auto results = analyzer.run();
    auto run_end = std::chrono::high_resolution_clock::now();
    double run_duration = std::chrono::duration<double, std::milli>(run_end - run_start).count();

    std::cout << "  Site testing completed in: " << std::fixed << std::setprecision(1) 
              << run_duration << " ms\n\n";

    // 5. Report Significant Sites
    std::vector<const MEMESiteResult*> positive_sites;
    for (const auto& r : results) {
        if (r.p_value <= pvalue_threshold && r.beta_plus > r.alpha) {
            positive_sites.push_back(&r);
        }
    }

    std::cout << "=========================================================================================\n"
              << " Codon |   alpha   |   beta1 (p1)   |    beta+ (p+)    |   LRT   |  p-value | # Branches \n"
              << "=========================================================================================\n";

    for (const auto* r : positive_sites) {
        std::cout << " " << std::setw(5) << (r->site_index + 1) << " | "
                  << std::fixed << std::setw(9) << std::setprecision(3) << r->alpha << " | "
                  << std::setw(6) << std::setprecision(2) << r->beta1 << " (" 
                  << std::setw(4) << std::setprecision(2) << r->p1 << ") | "
                  << std::setw(8) << std::setprecision(2) << r->beta_plus << " ("
                  << std::setw(4) << std::setprecision(2) << r->p_plus << ") | "
                  << std::setw(7) << std::setprecision(3) << r->lrt << " | "
                  << std::setw(8) << std::setprecision(4) << r->p_value << " | "
                  << std::setw(10) << r->branches_under_selection << "\n";
    }

    std::cout << "=========================================================================================\n"
              << "Found " << positive_sites.size() 
              << " site(s) under episodic diversifying positive selection at p <= " 
              << pvalue_threshold << "\n\n";

    // 6. Save JSON
    std::cout << "> Writing results to '" << output_file << "'...\n";
    try {
        analyzer.save_json(output_file);
        std::cout << "  Successfully saved HyPhy JSON report.\n\n";
    } catch (const std::exception& e) {
        std::cerr << "Error writing JSON: " << e.what() << "\n";
        return 1;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_duration = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "Total Analysis Runtime: " << std::fixed << std::setprecision(2) 
              << (total_duration / 1000.0) << " seconds.\n\n";

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_meme(argc, argv);
}
#endif
