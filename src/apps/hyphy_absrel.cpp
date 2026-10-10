#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/console.hpp"
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
    Panel::print_banner(
        "aBSREL",
        "Adaptive Branch-Site Random Effects Likelihood",
        "Lineage-specific test for episodic diversifying selection",
        "Citation: Mol Biol Evol. 32(5): 1342-1353 (2015) • v3.0.0"
    );
}

static void print_absrel_usage(const char* prog) {
    print_absrel_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal)\n"
              << "  " << Console::brand("--branches") << " <set>     Branches to test: All, Internal, Leaves (default: All)\n"
              << "  " << Console::brand("--pvalue") << " <threshold> Holm-Bonferroni p-value threshold (default: 0.05)\n"
              << "  " << Console::brand("--max-rates") << " <N>      Maximum rate classes per branch (default: 3)\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.ABSREL.json)\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
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
        } else if (arg == "--progress") {
            show_progress = true;
            force_progress = true;
        } else if (arg == "--no-progress") {
            show_progress = false;
            force_progress = false;
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

    // 4. Initialize and Run aBSREL
    auto absrel = ABSRELAnalyzer::create(tree, aln, settings);

    Panel::print_card("Dataset & Lineage Testing", {
        {"Alignment File", alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons) + " (" + std::to_string(aln.num_codons * 3) + " nt)"},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree", std::to_string(tree.num_nodes()) + " nodes (" + std::to_string(tree.num_leaves()) + " leaves)"},
        {"Tested Branches", settings.test_branches},
        {"Max Rate Classes", std::to_string(settings.max_rate_classes)},
        {"Genetic Code", code->name}
    });

    Panel::print_step(1, 2, "Fitting Baseline GTR & MG94xREV Models", "Exchangeability rates & branch lengths");

    std::function<void(const std::string&, double)> progress_cb = nullptr;
    if (!show_progress && !force_progress) {
        progress_cb = [](const std::string& stage, double frac) {
            std::cout << "  " << Console::muted("↳ [") << std::setw(3) << static_cast<int>(frac * 100) 
                      << Console::muted("%] ") << stage << std::endl;
        };
    }

    Panel::print_step(2, 2, "Running Adaptive Model Selection & Lineage Testing", "OpenMP parallel");
    auto res = absrel.run(progress_cb, show_progress, force_progress);

    // Branch table
    Table branch_table;
    branch_table.add_column("Branch", Table::Align::Left, 16);
    branch_table.add_column("Rates", Table::Align::Center, 5);
    branch_table.add_column("Full Length", Table::Align::Right, 11);
    branch_table.add_column("Test LRT", Table::Align::Right, 8);
    branch_table.add_column("Uncorrected p", Table::Align::Right, 13);
    branch_table.add_column("Corrected p", Table::Align::Right, 12);
    branch_table.add_column("EBF>=100", Table::Align::Right, 9);
    branch_table.add_column("Selection", Table::Align::Left, 14);

    for (const auto& bname : res.tested_branches) {
        const auto& bres = res.branches.at(bname);
        std::ostringstream len_ss, lrt_ss, un_p_ss, corr_p_ss;
        len_ss << std::fixed << std::setprecision(4) << bres.full_branch_length;
        lrt_ss << std::fixed << std::setprecision(2) << bres.lrt;
        un_p_ss << std::scientific << std::setprecision(2) << bres.uncorrected_p_value;
        corr_p_ss << std::scientific << std::setprecision(2) << bres.corrected_p_value;

        std::string sel_label = Console::muted("-");
        std::string row_color = "";
        if (bres.is_positive) {
            sel_label = Console::success("▲ Positive");
            row_color = "\033[38;5;48m";
        }

        branch_table.add_row({
            bname,
            std::to_string(bres.rate_classes),
            len_ss.str(),
            lrt_ss.str(),
            un_p_ss.str(),
            corr_p_ss.str(),
            std::to_string(bres.sites_ebf_100),
            sel_label
        }, row_color);
    }

    std::cout << "\n" << Console::bold("  Tested Branches Summary:") << "\n\n";
    branch_table.print();
    std::cout << "\n";

    std::ostringstream gtr_ss, base_ss, full_ss, time_ss;
    gtr_ss << std::fixed << std::setprecision(2) << res.gtr_fit.log_likelihood << " (AICc: " << res.gtr_fit.aicc << ")";
    base_ss << std::fixed << std::setprecision(2) << res.baseline_fit.log_likelihood << " (AICc: " << res.baseline_fit.aicc << ")";
    full_ss << std::fixed << std::setprecision(2) << res.full_adaptive_fit.log_likelihood << " (AICc: " << res.full_adaptive_fit.aicc << ")";
    time_ss << std::fixed << std::setprecision(2) << res.runtime_seconds << " s";

    std::string conclusion;
    bool is_sig = !res.positive_branches.empty();
    if (is_sig) {
        conclusion = "EPISODIC DIVERSIFYING SELECTION detected on " + std::to_string(res.positive_branches.size()) + 
                     " branch(es) at Holm-Bonferroni p <= " + std::to_string(settings.p_threshold) + "!";
    } else {
        conclusion = "No branches detected to be under positive selection at Holm-Bonferroni p <= " + std::to_string(settings.p_threshold) + ".";
    }

    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Tested Branches", std::to_string(res.tested_branches.size()) + " (" + settings.test_branches + ")"},
        {"Positively Selected Branches", 
         (is_sig ? Console::success(std::to_string(res.positive_branches.size())) : "0") + 
         " [Holm-Bonferroni p <= " + std::to_string(settings.p_threshold) + "]"},
        {"Baseline MG94xREV Log-L", base_ss.str()},
        {"Full Adaptive aBSREL Log-L", full_ss.str()},
        {"Nucleotide GTR Log-L", gtr_ss.str()},
        {"JSON Output File", output_file},
        {"Total Execution Time", time_ss.str()}
    };

    Panel::print_summary_card("aBSREL Lineage Selection Analysis Summary", sum_items, conclusion, is_sig);

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
    return run_absrel(argc, argv);
}
#endif
