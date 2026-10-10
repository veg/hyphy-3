#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/meme.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;
using namespace hyphy::opt;

static void print_banner() {
    Panel::print_banner(
        "MEME",
        "Mixed Effects Model of Evolution",
        "Detecting individual sites subject to episodic diversifying selection",
        "Citation: PLoS Genet 8(7): e1002764 (2012) • v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file (optional if embedded in NEXUS)\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.MEME.json)\n"
              << "  " << Console::brand("--pvalue") << " <float>     P-value significance threshold (default: 0.1)\n"
              << "  " << Console::brand("--full-model") << "         Perform branch length re-optimization under full codon model (default)\n"
              << "  " << Console::brand("--quick") << "              Disable full branch re-optimization (proportional branch scaling)\n"
              << "  " << Console::brand("--all-sites") << "          Display all codon sites in console table (default: significant only)\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
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
    bool show_progress = ProgressBar::is_terminal();
    bool force_progress = false;
    bool full_model = true;
    bool all_sites = false;

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
        } else if (arg == "--all-sites") {
            all_sites = true;
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
            std::cerr << Console::danger("Unknown argument: ") << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (alignment_file.empty()) {
        std::cerr << Console::danger("Error: --alignment is required.\n\n");
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
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error resolving genetic code: ") << e.what() << "\n";
        return 1;
    }

    // 1. Load Alignment
    Alignment aln;
    try {
        aln = Alignment::load(alignment_file, gcode);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error loading alignment: ") << e.what() << "\n";
        return 1;
    }

    // 2. Load Tree
    Tree tree;
    if (!tree_file.empty()) {
        try {
            tree = Tree::from_newick_file(tree_file);
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error loading tree: ") << e.what() << "\n";
            return 1;
        }
    } else if (!aln.embedded_tree_newick.empty()) {
        tree = Tree::from_newick(aln.embedded_tree_newick);
    } else {
        std::cerr << Console::danger("Error: No tree provided and no embedded tree found in alignment.\n");
        return 1;
    }

    Panel::print_card("Dataset & Phylogeny", {
        {"Alignment File", alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons) + " (" + std::to_string(aln.num_codons * 3) + " nt)"},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree", std::to_string(tree.num_nodes()) + " nodes (" + std::to_string(tree.num_leaves()) + " leaves)"},
        {"Genetic Code", gcode->name + " (" + std::to_string(gcode->num_sense_codons) + " sense codons)"}
    });

    // 3. Multi-phase Global Fitting
    Panel::print_step(1, 2, "Fitting Baseline GTR & Global MG94 Codon Models", 
                      full_model ? "Full branch re-optimization" : "Proportional branch scaling");
    auto progress_cb = [](const std::string& msg, double pct) {
        std::cout << "  " << Console::muted("↳ [") << std::setw(3) << static_cast<int>(pct * 100)
                  << Console::muted("%] ") << msg << "...\n";
    };

    auto fit_start = std::chrono::high_resolution_clock::now();
    auto analyzer = MEMEAnalyzer::create_and_fit(tree, std::move(aln), pvalue_threshold, progress_cb, full_model);
    auto fit_end = std::chrono::high_resolution_clock::now();
    double fit_duration = std::chrono::duration<double, std::milli>(fit_end - fit_start).count();

    std::cout << "  " << Console::muted("GTR Log-Likelihood  : ") << Console::bold(std::to_string(analyzer.gtr_log_l)) << "\n"
              << "  " << Console::muted("MG94 Log-Likelihood : ") << Console::bold(std::to_string(analyzer.global_log_l))
              << Console::muted("  |  omega = ") << Console::brand(std::to_string(analyzer.base_params.beta / analyzer.base_params.alpha))
              << Console::muted("  |  Fit time: ") << std::fixed << std::setprecision(1) << fit_duration << " ms\n\n";

    // 4. Site-level MEME testing
    Panel::print_step(2, 2, "Testing Codon Sites for Episodic Selection", "OpenMP parallel");
    auto results = analyzer.run(show_progress, force_progress);

    // 5. Report Sites
    Table table;
    table.add_column("Codon", Table::Align::Right, 6);
    table.add_column("alpha (dS)", Table::Align::Right, 10);
    table.add_column("beta1 (p1)", Table::Align::Right, 14);
    table.add_column("beta+ (p+)", Table::Align::Right, 14);
    table.add_column("LRT", Table::Align::Right, 8);
    table.add_column("p-value", Table::Align::Right, 9);
    table.add_column("# Branches", Table::Align::Right, 10);

    std::vector<const MEMESiteResult*> positive_sites;
    for (const auto& r : results) {
        if (r.p_value <= pvalue_threshold && r.beta_plus > r.alpha) {
            positive_sites.push_back(&r);
        }
    }

    std::vector<const MEMESiteResult*> sites_to_display;
    if (all_sites) {
        for (const auto& r : results) sites_to_display.push_back(&r);
    } else if (!positive_sites.empty()) {
        sites_to_display = positive_sites;
    } else {
        std::vector<size_t> sorted_by_lrt(results.size());
        std::iota(sorted_by_lrt.begin(), sorted_by_lrt.end(), 0);
        std::sort(sorted_by_lrt.begin(), sorted_by_lrt.end(), [&](size_t i, size_t j) {
            return results[i].lrt > results[j].lrt;
        });
        size_t n_top = std::min(sorted_by_lrt.size(), size_t(5));
        for (size_t i = 0; i < n_top; ++i) sites_to_display.push_back(&results[sorted_by_lrt[i]]);
    }

    for (const auto* r : sites_to_display) {
        std::ostringstream a_ss, b1_ss, bp_ss, lrt_ss, p_ss;
        a_ss << std::fixed << std::setprecision(3) << r->alpha;
        b1_ss << std::fixed << std::setprecision(2) << r->beta1 << " (" << std::setprecision(2) << r->p1 << ")";
        bp_ss << std::fixed << std::setprecision(2) << r->beta_plus << " (" << std::setprecision(2) << r->p_plus << ")";
        lrt_ss << std::fixed << std::setprecision(3) << r->lrt;
        p_ss << std::fixed << std::setprecision(4) << r->p_value;

        std::string row_color = "";
        if (r->p_value <= pvalue_threshold && r->beta_plus > r->alpha) {
            row_color = "\033[38;5;48m"; // Emerald highlight
        }

        table.add_row({
            std::to_string(r->site_index + 1),
            a_ss.str(),
            b1_ss.str(),
            bp_ss.str(),
            lrt_ss.str(),
            p_ss.str(),
            std::to_string(r->branches_under_selection)
        }, row_color);
    }

    std::cout << "\n";
    table.print();

    if (!all_sites) {
        if (!positive_sites.empty()) {
            std::cout << Console::muted("  Showing " + std::to_string(positive_sites.size()) + 
                                       " site(s) under episodic diversifying selection (p <= " + 
                                       std::to_string(pvalue_threshold) + 
                                       ") out of " + std::to_string(results.size()) + 
                                       ". (Use --all-sites to display all)\n\n");
        } else {
            std::cout << Console::muted("  No sites reached significance (p <= " + 
                                       std::to_string(pvalue_threshold) + 
                                       "). Displaying top 5 sites by LRT. (Use --all-sites to display all)\n\n");
        }
    } else {
        std::cout << "\n";
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(end_time - start_time).count();

    auto fmt_pct = [](size_t n, size_t tot) {
        if (tot == 0) return std::string("0.0%");
        double p = (100.0 * n) / tot;
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << p << "%";
        return ss.str();
    };

    std::ostringstream gtr_ss, mg_ss, time_ss, omega_ss;
    gtr_ss << std::fixed << std::setprecision(2) << analyzer.gtr_log_l;
    mg_ss << std::fixed << std::setprecision(2) << analyzer.global_log_l;
    time_ss << std::fixed << std::setprecision(2) << total_sec << " s";
    omega_ss << std::fixed << std::setprecision(4) << (analyzer.base_params.beta / analyzer.base_params.alpha);

    std::string conclusion = std::to_string(positive_sites.size()) + 
        " site(s) under episodic diversifying selection identified (p <= " + std::to_string(pvalue_threshold) + ")";

    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Tested Codon Sites", std::to_string(results.size())},
        {"Episodic Diversifying Sites", 
         Console::success(std::to_string(positive_sites.size())) + " sites (" + 
         fmt_pct(positive_sites.size(), results.size()) + ") [p <= " + std::to_string(pvalue_threshold) + "]"},
        {"Global dN/dS (beta/alpha)", omega_ss.str()},
        {"Nucleotide GTR Log-L", gtr_ss.str()},
        {"Global MG94xREV Log-L", mg_ss.str()},
        {"JSON Output File", output_file},
        {"Total Execution Time", time_ss.str()}
    };

    Panel::print_summary_card("MEME Episodic Selection Analysis Summary", sum_items, conclusion, !positive_sites.empty());

    // Save JSON
    try {
        analyzer.save_json(output_file);
        std::cout << Console::success("✔") << " " << Console::bold("Saved Datamonkey-compatible JSON report to: ")
                  << Console::brand(output_file) << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error writing JSON: ") << e.what() << "\n";
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_meme(argc, argv);
}
#endif
