#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/fel.hpp"

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
        "FEL",
        "Fixed Effects Likelihood",
        "Site-level test for diversifying and purifying selection",
        "Citation: Not So Different After All (2005) MBE 22:1208 • v3.0.0"
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
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.FEL.json)\n"
              << "  " << Console::brand("--pvalue") << " <float>     P-value significance threshold (default: 0.1)\n"
              << "  " << Console::brand("--json-format") << " <type> Output JSON schema format: 'modern' (v3.0, default) or 'legacy' (Datamonkey 2.5)\n"
              << "  " << Console::brand("--legacy-json") << "        Emit legacy Datamonkey/HyPhy 2.5 JSON schema\n"
              << "  " << Console::brand("--modern-json") << "        Emit modern HyPhy 3.0 JSON schema with provenance and columnar data\n"
              << "  " << Console::brand("--full-model") << "         Perform branch length re-optimization under full codon model (default)\n"
              << "  " << Console::brand("--quick") << "              Disable full branch re-optimization (proportional branch scaling)\n"
              << "  " << Console::brand("--all-sites") << "          Display all codon sites in console table (default: significant only)\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk --output cd2.FEL.json\n"
              << "  " << prog << " --alignment tests/data/COXI.nex --code Vertebrate-mtDNA --legacy-json\n\n";
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
    bool all_sites = false;
    JSONFormat json_format = JSONFormat::ModernV3;

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
        } else if (arg == "--json-format" && i + 1 < argc) {
            json_format = parse_json_format(argv[++i]);
        } else if (arg == "--legacy-json") {
            json_format = JSONFormat::Legacy;
        } else if (arg == "--modern-json") {
            json_format = JSONFormat::ModernV3;
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
        output_file = alignment_file + ".FEL.json";
    }

    print_banner();

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string start_iso = Provenance::current_iso8601();
    std::string cli_cmd = argv[0];
    for (int i = 1; i < argc; ++i) {
        cli_cmd += " ";
        cli_cmd += argv[i];
    }

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
    std::string tree_source;
    if (!tree_file.empty()) {
        try {
            tree = Tree::from_newick_file(tree_file);
            tree_source = tree_file;
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error loading tree: ") << e.what() << "\n";
            return 1;
        }
    } else if (!aln.embedded_tree_newick.empty()) {
        try {
            tree = Tree::from_newick(aln.embedded_tree_newick);
            tree_source = "embedded in NEXUS";
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error parsing embedded tree: ") << e.what() << "\n";
            return 1;
        }
    } else {
        std::cerr << Console::danger("Error: No tree provided via --tree and no embedded tree found in alignment.\n\n");
        print_usage(argv[0]);
        return 1;
    }

    // Display Dataset Card
    Panel::print_card("Dataset & Phylogeny", {
        {"Alignment File", alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons) + " (" + std::to_string(aln.num_codons * 3) + " nt)"},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree Source", tree_source + " (" + std::to_string(tree.num_nodes()) + " nodes, " + std::to_string(aln.num_taxa) + " leaves)"},
        {"Genetic Code", gcode->name + " (" + std::to_string(gcode->num_sense_codons) + " sense codons)"}
    });

    // 3. Global Model Phase
    // Step 1: Nucleotide GTR Model
    Panel::print_step(1, 3, "Fitting Nucleotide GTR Model", "Branch lengths & exchangeability rates");
    GTRFitter gtr_fitter(tree, aln);
    auto gtr_res = gtr_fitter.fit();
    std::cout << "  " << Console::muted("Log-Likelihood : ") << Console::bold(std::to_string(gtr_res.log_likelihood))
              << Console::muted("  |  AICc: ") << std::fixed << std::setprecision(2) << gtr_res.aicc
              << Console::muted("  |  Iter: ") << gtr_res.iterations << "\n"
              << "  " << Console::muted("Biases         : ")
              << "AC=" << std::setprecision(4) << gtr_res.params.theta_AC << "  "
              << "AT=" << gtr_res.params.theta_AT << "  "
              << "CG=" << gtr_res.params.theta_CG << "  "
              << "CT=" << gtr_res.params.theta_CT << "  "
              << "GT=" << gtr_res.params.theta_GT << "\n\n";

    // Step 2: Global MG94xREV Model
    Panel::print_step(2, 3, "Refining Global MG94xREV Codon Model", 
                      full_model ? "Full branch re-optimization" : "Proportional branch scaling");
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
                std::cout << "  " << Console::muted("↳ ") << msg << "...\n";
            };
        }
        mg_res = mg_fitter.fit_full_model(1.0, base_p, mg_cb);
        base_p = mg_res.params;
    } else {
        mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
        base_p.alpha = 1.0;
        base_p.beta = mg_res.x_opt(0);
    }

    std::cout << "  " << Console::muted("Log-Likelihood : ") << Console::bold(std::to_string(mg_res.log_likelihood))
              << Console::muted("  |  Global dN/dS (omega): ") << Console::brand(std::to_string(base_p.beta)) << "\n\n";

    // Step 3: Site-by-Site Testing Phase
    Panel::print_step(3, 3, "Testing Codon Sites for Selection", "OpenMP parallel");
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

    // Site results table
    Table table;
    table.add_column("Codon", Table::Align::Right, 6);
    table.add_column("alpha (dS)", Table::Align::Right, 10);
    table.add_column("beta (dN)", Table::Align::Right, 10);
    table.add_column("alpha=beta", Table::Align::Right, 10);
    table.add_column("LRT", Table::Align::Right, 8);
    table.add_column("p-value", Table::Align::Right, 9);
    table.add_column("Selection", Table::Align::Left, 14);

    std::vector<size_t> sig_indices;
    for (size_t c = 0; c < site_results.size(); ++c) {
        if (site_results[c].p_value <= pvalue_threshold) {
            sig_indices.push_back(c);
        }
    }

    std::vector<size_t> sites_to_display;
    if (all_sites) {
        sites_to_display.resize(site_results.size());
        std::iota(sites_to_display.begin(), sites_to_display.end(), 0);
    } else if (!sig_indices.empty()) {
        sites_to_display = sig_indices;
    } else {
        // Top 5 by LRT
        std::vector<size_t> sorted_by_lrt(site_results.size());
        std::iota(sorted_by_lrt.begin(), sorted_by_lrt.end(), 0);
        std::sort(sorted_by_lrt.begin(), sorted_by_lrt.end(), [&](size_t i, size_t j) {
            return site_results[i].lrt > site_results[j].lrt;
        });
        size_t n_top = std::min(sorted_by_lrt.size(), size_t(5));
        sites_to_display.assign(sorted_by_lrt.begin(), sorted_by_lrt.begin() + n_top);
    }

    for (size_t c : sites_to_display) {
        const auto& r = site_results[c];
        std::string sel_label = Console::muted("-");
        std::string row_color = "";
        if (r.p_value <= pvalue_threshold) {
            if (r.beta > r.alpha) {
                sel_label = Console::success("▲ Diversifying");
                row_color = "\033[38;5;48m";
            } else if (r.alpha > r.beta) {
                sel_label = Console::info("▼ Purifying");
                row_color = "\033[38;5;75m";
            }
        }

        std::ostringstream a_ss, b_ss, n_ss, lrt_ss, p_ss;
        a_ss << std::fixed << std::setprecision(4) << r.alpha;
        b_ss << std::fixed << std::setprecision(4) << r.beta;
        n_ss << std::fixed << std::setprecision(4) << r.alpha_null;
        lrt_ss << std::fixed << std::setprecision(4) << r.lrt;
        p_ss << std::fixed << std::setprecision(4) << r.p_value;

        table.add_row({
            std::to_string(c + 1),
            a_ss.str(),
            b_ss.str(),
            n_ss.str(),
            lrt_ss.str(),
            p_ss.str(),
            sel_label
        }, row_color);
    }

    std::cout << "\n";
    table.print();

    if (!all_sites) {
        if (!sig_indices.empty()) {
            std::cout << Console::muted("  Showing " + std::to_string(sig_indices.size()) + 
                                       " significant site(s) at p <= " + std::to_string(pvalue_threshold) + 
                                       " out of " + std::to_string(aln.num_codons) + 
                                       ". (Use --all-sites to display every site)\n\n");
        } else {
            std::cout << Console::muted("  No sites reached significance threshold (p <= " + 
                                       std::to_string(pvalue_threshold) + 
                                       "). Displaying top 5 sites by LRT. (Use --all-sites to display all)\n\n");
        }
    } else {
        std::cout << "\n";
    }

    auto summary = fel.get_summary();
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
    gtr_ss << std::fixed << std::setprecision(2) << gtr_res.log_likelihood << " (AICc: " << gtr_res.aicc << ")";
    mg_ss << std::fixed << std::setprecision(2) << mg_res.log_likelihood << " (AICc: " << fel.global_aicc << ")";
    time_ss << std::fixed << std::setprecision(2) << total_sec << " s";
    omega_ss << std::fixed << std::setprecision(4) << base_p.beta;

    std::string json_desc = (json_format == JSONFormat::Legacy) ? " [Legacy]" : " [Modern v3.0]";
    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Tested Codon Sites", std::to_string(summary.tested_sites)},
        {"Diversifying Selection (Positive)", 
         Console::success(std::to_string(summary.positively_selected)) + " sites (" + 
         fmt_pct(summary.positively_selected, summary.tested_sites) + ") [p <= " + std::to_string(pvalue_threshold) + "]"},
        {"Purifying Selection (Negative)", 
         Console::info(std::to_string(summary.negatively_selected)) + " sites (" + 
         fmt_pct(summary.negatively_selected, summary.tested_sites) + ") [p <= " + std::to_string(pvalue_threshold) + "]"},
        {"Global dN/dS (beta/alpha)", omega_ss.str()},
        {"Nucleotide GTR Fit", gtr_ss.str()},
        {"Global MG94xREV Fit", mg_ss.str()},
        {"JSON Output File", output_file + json_desc},
        {"Total Execution Time", time_ss.str()}
    };

    std::string conclusion = std::to_string(summary.positively_selected) + 
        " diversifying site(s) and " + std::to_string(summary.negatively_selected) + 
        " purifying site(s) identified (p <= " + std::to_string(pvalue_threshold) + ")";

    Panel::print_summary_card("FEL Selection Analysis Summary", sum_items, conclusion, true);

    // Build Provenance
    Provenance prov;
    prov.invocation.cli_command = cli_cmd;
    prov.invocation.working_directory = Provenance::get_cwd();
    prov.invocation.arguments["alignment"] = alignment_file;
    if (!tree_file.empty()) prov.invocation.arguments["tree"] = tree_file;
    prov.invocation.arguments["code"] = code_name;
    prov.invocation.arguments["pvalue_threshold"] = std::to_string(pvalue_threshold);
    prov.invocation.arguments["full_model"] = full_model ? "true" : "false";
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
    prov.execution.wall_time_seconds = total_sec;
#ifdef _OPENMP
    prov.execution.cpu_threads = (num_threads > 0) ? num_threads : omp_get_max_threads();
#else
    prov.execution.cpu_threads = 1;
#endif
    prov.execution.hostname = Provenance::get_hostname();
    prov.execution.os = Provenance::detect_os();
    prov.execution.compiler = Provenance::detect_compiler();

    // Save JSON
    try {
        fel.save_json(output_file, alignment_file, "", json_format, prov);
        std::string fmt_desc = (json_format == JSONFormat::Legacy) ? "Legacy Datamonkey JSON" : "Modern HyPhy v3.0 JSON";
        std::cout << Console::success("✔") << " " << Console::bold("Saved " + fmt_desc + " report to: ")
                  << Console::brand(output_file) << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error writing JSON: ") << e.what() << "\n";
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_fel(argc, argv);
}
#endif
