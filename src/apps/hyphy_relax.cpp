#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/console.hpp"
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
    Panel::print_banner(
        "RELAX",
        "Test for Selection Relaxation",
        "Detecting relaxed or intensified selection in a phylogenetic framework",
        "Citation: Mol Biol Evol. 32(3): 820-832 (2015) • v3.0.0"
    );
}

static void print_relax_usage(const char* prog) {
    print_relax_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal)\n"
              << "  " << Console::brand("--test") << " <regex|names> Branch names or regex to mark as Test set (default: use {T} tags)\n"
              << "  " << Console::brand("--pvalue") << " <threshold> Significance threshold for LRT (default: 0.05)\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.RELAX.json)\n"
              << "  " << Console::brand("--no-branch-opt") << "      Disable individual branch length refinement (use proportional scaling)\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
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
        } else if (arg == "--no-branch-opt" || arg == "--no-refine-branches") {
            settings.refine_branch_lengths = false;
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
        std::cerr << Console::danger("Error: Unknown genetic code '") << code_name << "'.\n";
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
            std::cerr << Console::danger("Error loading tree from file: ") << e.what() << "\n";
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
        std::cerr << Console::danger("Error: No tree provided and no tree embedded in alignment.\n");
        return 1;
    }

    // 4. Initialize and Run RELAX
    auto relax = RELAXAnalyzer::create(tree, aln, settings);

    if (relax.test_branch_names.empty()) {
        std::cerr << Console::danger("Error: No Test branches found! Annotate branches with {T} in Newick or pass --test <regex>.\n");
        return 1;
    }

    Panel::print_card("Dataset & Phylogeny", {
        {"Alignment File", alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons) + " (" + std::to_string(aln.num_codons * 3) + " nt)"},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree", std::to_string(tree.num_nodes()) + " nodes (" + std::to_string(tree.num_leaves()) + " leaves)"},
        {"Test Branches", std::to_string(relax.test_branch_names.size()) + " branches"},
        {"Reference Branches", std::to_string(relax.ref_branch_names.size()) + " branches"},
        {"Genetic Code", code->name}
    });

    Panel::print_step(1, 3, "Fitting Baseline & Separate-Rates Codon Models", "GTR & MG94xREV");

    std::function<void(const std::string&, double)> progress_cb = nullptr;
    if (!show_progress && !force_progress) {
        progress_cb = [](const std::string& stage, double frac) {
            std::cout << "  " << Console::muted("↳ [") << std::setw(3) << static_cast<int>(frac * 100) 
                      << Console::muted("%] ") << stage << std::endl;
        };
    }

    Panel::print_step(2, 3, "Running RELAX Inference Pipeline", 
                      settings.refine_branch_lengths ? "Two-loop L-BFGS branch refinement enabled" : "Proportional branch scaling");

    auto res = relax.run(progress_cb, show_progress, force_progress);

    // Rate distributions
    Table rate_table;
    rate_table.add_column("Class", Table::Align::Center, 5);
    rate_table.add_column("Ref omega", Table::Align::Right, 12);
    rate_table.add_column("Ref Proportion", Table::Align::Right, 14);
    rate_table.add_column("Test omega", Table::Align::Right, 12);
    rate_table.add_column("Test Proportion", Table::Align::Right, 15);

    for (size_t c = 0; c < res.alternative_fit.reference_distribution.omegas.size(); ++c) {
        std::ostringstream w_r_ss, p_r_ss, w_t_ss, p_t_ss;
        w_r_ss << std::fixed << std::setprecision(4) << res.alternative_fit.reference_distribution.omegas[c];
        p_r_ss << std::fixed << std::setprecision(4) << res.alternative_fit.reference_distribution.weights[c];
        w_t_ss << std::fixed << std::setprecision(4) << res.alternative_fit.test_distribution.omegas[c];
        p_t_ss << std::fixed << std::setprecision(4) << res.alternative_fit.test_distribution.weights[c];

        rate_table.add_row({
            std::to_string(c),
            w_r_ss.str(),
            p_r_ss.str(),
            w_t_ss.str(),
            p_t_ss.str()
        });
    }

    std::cout << "\n" << Console::bold("  Inferred Rate Distributions (Alternative Model):") << "\n\n";
    rate_table.print();
    std::cout << "\n";

    std::ostringstream k_ss, lrt_ss, p_ss, gtr_ss, mg_ss, alt_ss, null_ss, time_ss;
    k_ss << std::fixed << std::setprecision(4) << res.k;
    lrt_ss << std::fixed << std::setprecision(2) << res.lrt;
    p_ss << std::scientific << std::setprecision(4) << res.p_value;
    gtr_ss << std::fixed << std::setprecision(2) << res.gtr_log_likelihood << " (AICc: " << res.gtr_aicc << ")";
    mg_ss << std::fixed << std::setprecision(2) << res.mg94_log_likelihood << " (AICc: " << res.mg94_aicc << ")";
    alt_ss << std::fixed << std::setprecision(2) << res.alternative_fit.log_likelihood << " (AICc: " << res.alternative_fit.aicc << ")";
    null_ss << std::fixed << std::setprecision(2) << res.null_fit.log_likelihood << " (AICc: " << res.null_fit.aicc << ")";
    time_ss << std::fixed << std::setprecision(2) << res.runtime_seconds << " s";

    std::string conclusion;
    if (res.is_significant) {
        if (res.is_relaxed) {
            conclusion = "SIGNIFICANT RELAXATION of selection on Test branches (K = " + k_ss.str() + " < 1, p = " + p_ss.str() + ")";
        } else {
            conclusion = "SIGNIFICANT INTENSIFICATION of selection on Test branches (K = " + k_ss.str() + " > 1, p = " + p_ss.str() + ")";
        }
    } else {
        conclusion = "No significant evidence of relaxation or intensification (p = " + p_ss.str() + " > " + std::to_string(settings.p_value_threshold) + ")";
    }

    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Relaxation Parameter (K)", k_ss.str()},
        {"Likelihood Ratio Test (LRT)", lrt_ss.str()},
        {"p-value (asymptotic)", p_ss.str() + " [threshold: " + std::to_string(settings.p_value_threshold) + "]"},
        {"RELAX Alternative Log-L", alt_ss.str()},
        {"RELAX Null (K=1) Log-L", null_ss.str()},
        {"Nucleotide GTR Log-L", gtr_ss.str()},
        {"Separate Rates MG94xREV", mg_ss.str()},
        {"Total Execution Time", time_ss.str()}
    };

    Panel::print_summary_card("RELAX Selection Relaxation Analysis Summary", sum_items, conclusion, res.is_significant);

    // Save JSON
    std::ofstream out(output_file);
    if (!out) {
        std::cerr << Console::danger("Warning: Could not open output file for writing: ") << output_file << "\n";
    } else {
        out << res.to_json(tree, aln).dump(2) << "\n";
        std::cout << Console::success("✔") << " " << Console::bold("Saved Datamonkey-compatible JSON report to: ")
                  << Console::brand(output_file) << "\n\n";
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_relax(argc, argv);
}
#endif
