#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/fel.hpp"

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
              << "       HYPHY 3: Fixed Effects Likelihood (FEL)         \n"
              << "=======================================================\n"
              << " Citation: Not So Different After All (2005) MBE 22:1208\n"
              << " Version:  3.0.0 (Modern C++20 Core)\n"
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
              << "  --output <file>      Path to output JSON file (default: <alignment>.FEL.json)\n"
              << "  --pvalue <float>     P-value significance threshold (default: 0.1)\n"
              << "  --full-model         Perform branch length re-optimization under full codon model (default)\n"
              << "  --quick              Disable full branch re-optimization (proportional branch scaling)\n"
              << "  --progress           Force interactive progress bar\n"
              << "  --no-progress        Disable progress bar\n"
              << "  --help, -h           Show this help message\n\n"
              << "Examples:\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk --output cd2.FEL.json\n"
              << "  " << prog << " --alignment tests/data/COXI.nex --code Vertebrate-mtDNA\n\n";
}

int run_fel(int argc, char* argv[]) {
    std::string alignment_file;
    std::string tree_file;
    std::string output_file;
    std::string code_name = "Universal";
    int num_threads = 0;
    Scalar pvalue_threshold = 0.1;
    bool show_progress = ProgressBar::is_terminal();
    bool force_progress = false;
    bool full_model = true;

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
        } else if (arg == "--full-model") {
            full_model = true;
        } else if (arg == "--no-full-model" || arg == "--quick") {
            full_model = false;
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
        output_file = alignment_file + ".FEL.json";
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
              << "  Codon sites: " << aln.num_codons << "\n"
              << "  Unique site patterns: " << aln.patterns.size() << "\n\n";

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
        try {
            tree = Tree::from_newick(aln.embedded_tree_newick);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing embedded tree: " << e.what() << "\n";
            return 1;
        }
    } else {
        std::cerr << "Error: No tree provided via --tree and no embedded tree found in alignment.\n\n";
        print_usage(argv[0]);
        return 1;
    }
    std::cout << "  Tree nodes: " << tree.num_nodes() << " (" << aln.num_taxa << " leaves)\n\n";

    // 3. Global Model Phase
    // Phase 3A: Nucleotide GTR Model
    std::cout << "### Phase 1: Fitting Nucleotide GTR Model...\n";
    GTRFitter gtr_fitter(tree, aln);
    auto gtr_res = gtr_fitter.fit();
    std::cout << "  GTR Log-Likelihood: " << std::fixed << std::setprecision(4) << gtr_res.log_likelihood
              << " (" << gtr_res.iterations << " iterations)\n"
              << "  Substitution biases:\n"
              << "    AC: " << gtr_res.params.theta_AC
              << "    AT: " << gtr_res.params.theta_AT
              << "    CG: " << gtr_res.params.theta_CG
              << "    CT: " << gtr_res.params.theta_CT
              << "    GT: " << gtr_res.params.theta_GT << "\n\n";

    // Phase 3B: Global MG94xREV Model
    std::cout << "### Phase 2: Refining under Global MG94xREV Model"
              << (full_model ? " (Full branch length re-optimization)...\n" : " (Proportional branch scaling)...\n");
    MG94Parameters base_p;
    base_p.theta_AC = gtr_res.params.theta_AC;
    base_p.theta_AT = gtr_res.params.theta_AT;
    base_p.theta_CG = gtr_res.params.theta_CG;
    base_p.theta_CT = gtr_res.params.theta_CT;
    base_p.theta_GT = gtr_res.params.theta_GT;

    MG94Fitter mg_fitter(gtr_res.tree, aln);
    FitResult mg_res;
    if (full_model) {
        std::function<void(const std::string&, double)> mg_cb = nullptr;
        if (!show_progress && !force_progress) {
            mg_cb = [](const std::string& msg, double) {
                std::cout << "  > " << msg << "...\n";
            };
        }
        mg_res = mg_fitter.fit_full_model(1.0, base_p, mg_cb);
        base_p = mg_res.params;
    } else {
        mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
        base_p.alpha = 1.0;
        base_p.beta = mg_res.x_opt(0);
    }

    std::cout << "  MG94 Log-Likelihood: " << std::fixed << std::setprecision(4) << mg_res.log_likelihood << "\n"
              << "  Global omega (dN/dS): " << base_p.beta << "\n\n";

    // 4. Site-by-Site Testing Phase
    std::cout << "### Phase 3: Testing " << aln.num_codons << " codon sites for selection (OpenMP accelerated)...\n";
    FELAnalyzer fel(mg_res.tree, aln, base_p);
    fel.p_value_threshold = pvalue_threshold;
    fel.global_log_l = mg_res.log_likelihood;
    size_t K_mg = 5 + 1 + (mg_res.tree.num_nodes() - 1);
    size_t N_codons = aln.num_codons;
    fel.global_aicc = -2.0 * mg_res.log_likelihood + 2.0 * K_mg * N_codons / (N_codons - K_mg - 1.0);

    fel.has_gtr_fit = true;
    fel.gtr_log_l = gtr_res.log_likelihood;
    fel.gtr_aicc = gtr_res.aicc;
    fel.gtr_rates = gtr_res.params;
    for (const auto& node : gtr_res.tree.nodes) {
        if (node.id != gtr_res.tree.root_id) {
            fel.gtr_branch_lengths[node.name] = node.branch_length;
        }
    }

    auto site_results = fel.run(show_progress, force_progress);

    std::cout << "\n--------------------------------------------------------------------------------\n"
              << " Codon |   alpha   |   beta    | alpha=beta |    LRT    |  p-value  | Selection \n"
              << "--------------------------------------------------------------------------------\n";

    for (size_t c = 0; c < site_results.size(); ++c) {
        const auto& r = site_results[c];
        std::string sel_label = "-";
        if (r.p_value <= pvalue_threshold) {
            if (r.beta > r.alpha) {
                sel_label = "Positive";
            } else if (r.alpha > r.beta) {
                sel_label = "Negative";
            }
        }

        std::cout << std::setw(6) << (c + 1) << " | "
                  << std::setw(9) << std::fixed << std::setprecision(4) << r.alpha << " | "
                  << std::setw(9) << std::fixed << std::setprecision(4) << r.beta << " | "
                  << std::setw(10) << std::fixed << std::setprecision(4) << r.alpha_null << " | "
                  << std::setw(9) << std::fixed << std::setprecision(4) << r.lrt << " | "
                  << std::setw(9) << std::fixed << std::setprecision(4) << r.p_value << " | "
                  << sel_label << "\n";
    }
    std::cout << "--------------------------------------------------------------------------------\n\n";

    // 5. Summary
    auto summary = fel.get_summary();
    std::cout << "### Selection Summary (p-value <= " << pvalue_threshold << "):\n"
              << "  Tested sites: " << summary.tested_sites << "\n"
              << "  Diversifying positive selection: " << summary.positively_selected << "\n"
              << "  Purifying negative selection:   " << summary.negatively_selected << "\n\n";

    // 6. Output JSON
    try {
        fel.save_json(output_file, alignment_file);
        std::cout << "> Saved Datamonkey-compatible JSON report to: '" << output_file << "'\n";
    } catch (const std::exception& e) {
        std::cerr << "Error writing JSON: " << e.what() << "\n";
        return 1;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "> Total runtime: " << std::fixed << std::setprecision(2) << elapsed_ms << " ms\n\n";

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_fel(argc, argv);
}
#endif
