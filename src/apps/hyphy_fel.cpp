#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/fel.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <numeric>

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
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file or 'auto' (optional; auto-infers ML tree via TreeEngine if omitted or 'auto')\n"
              << "  " << Console::brand("--ohmic") << "              Build Ohmic Phylogenetic Network with Wheatstone bridges\n"
              << "  " << Console::brand("--hyphaeon") << "           Opt-in: HyphAeon 384D transformer conductances (needs Python aeon-core)\n"
              << "  " << Console::brand("--no-hyphaeon") << "        Use empirical codon conductances only (default)\n"
              << "  " << Console::brand("--embeddings") << " <file>  Path to precomputed [N x 384] HyphAeon embeddings (.bin/.csv/.tsv)\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  " << Console::brand("--branches") << " <set>     Branches to test: 'All' (default), 'Internal', 'Leaves', 'Unlabeled', tag name, or comma-separated names\n"
              << "  " << Console::brand("--branch-regex") << " <reg> Regular expression to match branch names to test\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.FEL.json)\n"
              << "  " << Console::brand("--html") << " <file>        Also write a self-contained interactive HTML report\n"
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
              << "  " << prog << " --alignment data/cd2.fna --tree auto --output cd2.FEL.json\n"
              << "  " << prog << " --alignment data/cd2.fna --ohmic --output cd2.FEL.json\n"
              << "  " << prog << " --alignment tests/data/COXI.nex --code Vertebrate-mtDNA --legacy-json\n\n";
}

int run_fel(int argc, char* argv[]) {
    CommonCLIOptions opts;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (parse_common_cli_arg(i, argc, argv, opts)) {
            if (opts.has_error) {
                return 1;
            }
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << Console::danger("Unknown argument: ") << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (opts.alignment_file.empty()) {
        std::cerr << Console::danger("Error: --alignment is required.\n\n");
        print_usage(argv[0]);
        return 1;
    }

    finalize_common_cli_options(opts, ".FEL.json");
    print_banner();

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string start_iso = Provenance::current_iso8601();
    std::string cli_cmd = join_cli_command(argc, argv);

    LoadedCodonDataset data;
    if (!load_codon_dataset(opts, data)) {
        return 1;
    }
    auto& [gcode, aln, tree, ohmic_net, tree_source] = data;
    const auto& alignment_file = opts.alignment_file;
    const auto& output_file = opts.output_file;
    const auto& test_branches = opts.branches;
    const auto& test_branch_regex = opts.branch_regex;
    const Scalar pvalue_threshold = opts.pvalue_threshold;
    const bool full_model = opts.full_model;
    const bool all_sites = opts.all_sites;
    const bool show_progress = opts.show_progress;
    const bool force_progress = opts.force_progress;
    const JSONFormat json_format = opts.json_format;

    // Validate branch selection upfront before global model fitting
    try {
        (void)tree.select_branches(test_branches, test_branch_regex);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error in branch selection: ") << e.what() << "\n";
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
    auto sel = gtr_res.tree.select_branches(test_branches, test_branch_regex, {}, true);
    const std::vector<bool>* p_sel = (sel.num_selected > 0 && sel.num_unselected > 0) ? &sel.is_selected : nullptr;
    FitResult mg_res;
    if (full_model) {
        std::function<void(const std::string&, double)> mg_cb = nullptr;
        if (!show_progress && !force_progress) {
            mg_cb = [](const std::string& msg, double) {
                std::cout << "  " << Console::muted("↳ ") << msg << "...\n";
            };
        }
        mg_res = mg_fitter.fit_full_model(1.0, base_p, mg_cb, 250, 1e-5, p_sel);
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
    FELSettings fel_settings;
    fel_settings.p_value_threshold = pvalue_threshold;
    fel_settings.full_model = full_model;
    fel_settings.test_branches = test_branches;
    fel_settings.test_branch_regex = test_branch_regex;

    FELAnalyzer fel(mg_res.tree, aln, base_p, fel_settings);
    fel.global_log_l = mg_res.log_likelihood;
    size_t K_mg = 5 + (p_sel ? 2 : 1) + (mg_res.tree.num_nodes() - 1);
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

    if (fel.num_background > 0) {
        std::cout << "  " << Console::muted("Tested branches : ") 
                  << Console::bold(std::to_string(fel.num_tested)) << Console::muted(" test, ") 
                  << Console::bold(std::to_string(fel.num_background)) << Console::muted(" background\n\n");
    }

    auto site_results = fel.run(show_progress, force_progress);

    auto render_rate_spark = [](double alpha, double beta) -> std::string {
        auto segs = [](double v) -> int {
            if (v <= 1e-4) return 0;
            if (v < 0.5) return 1;
            if (v < 3.0) return 2;
            return 3;
        };
        int na = segs(alpha);
        int nb = segs(beta);

        std::string sa_filled, sa_empty, sb_filled, sb_empty;
        for (int i = 0; i < na; ++i) sa_filled += "▰";
        for (int i = 0; i < 3 - na; ++i) sa_empty += "▱";
        for (int i = 0; i < nb; ++i) sb_filled += "▰";
        for (int i = 0; i < 3 - nb; ++i) sb_empty += "▱";

        if (Console::color_enabled()) {
            return "\033[38;5;75m" + sa_filled + "\033[0m" +
                   "\033[2m" + sa_empty + "\033[0m" +
                   "\033[38;5;244m:\033[0m" +
                   "\033[38;5;48m" + sb_filled + "\033[0m" +
                   "\033[2m" + sb_empty + "\033[0m";
        } else {
            return sa_filled + sa_empty + ":" + sb_filled + sb_empty;
        }
    };

    // Site results table
    Table table;
    table.add_column("Codon", Table::Align::Right, 6);
    if (fel.num_background > 0) {
        table.add_column("alpha : beta (test rates)", Table::Align::Right, 26);
        table.add_column("beta (bg)", Table::Align::Right, 10);
    } else {
        table.add_column("alpha : beta (rates)", Table::Align::Right, 24);
    }
    table.add_column("p-value", Table::Align::Right, 9);
    table.add_column("Selection", Table::Align::Left, 14);
    table.add_column("Substitutions", Table::Align::Left, 40);

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
        // Top 5 by p-value
        std::vector<size_t> sorted_by_lrt(site_results.size());
        std::iota(sorted_by_lrt.begin(), sorted_by_lrt.end(), 0);
        std::sort(sorted_by_lrt.begin(), sorted_by_lrt.end(), [&](size_t i, size_t j) {
            return site_results[i].p_value < site_results[j].p_value;
        });
        size_t n_top = std::min(sorted_by_lrt.size(), size_t(5));
        sites_to_display.assign(sorted_by_lrt.begin(), sorted_by_lrt.begin() + n_top);
    }

    for (size_t c : sites_to_display) {
        const auto& r = site_results[c];
        std::string sel_label = Console::muted("-");
        if (r.p_value <= pvalue_threshold) {
            if (r.beta > r.alpha) {
                sel_label = Console::success("▲ Diversifying");
            } else if (r.alpha > r.beta) {
                sel_label = Console::info("▼ Purifying");
            }
        }

        std::ostringstream ab_ss, p_ss;
        ab_ss << std::fixed << std::setprecision(2) << r.alpha << " : " << std::setprecision(2) << r.beta
              << "   " << render_rate_spark(r.alpha, r.beta);

        if (r.p_value < 0.0001) {
            p_ss << "<0.0001";
        } else {
            p_ss << std::fixed << std::setprecision(4) << r.p_value;
        }

        std::string subs_str = (fel.num_background > 0) 
            ? r.subs.summary_string_partitioned 
            : r.subs.summary_string;

        if (fel.num_background > 0) {
            std::ostringstream bg_ss;
            bg_ss << std::fixed << std::setprecision(2) << r.beta_nuisance;
            table.add_row({
                std::to_string(c + 1),
                ab_ss.str(),
                bg_ss.str(),
                p_ss.str(),
                sel_label,
                subs_str
            });
        } else {
            table.add_row({
                std::to_string(c + 1),
                ab_ss.str(),
                p_ss.str(),
                sel_label,
                subs_str
            });
        }
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
                                       "). Displaying top 5 sites by p-value. (Use --all-sites to display all)\n\n");
        }
    }

    std::string rate_desc = (fel.num_background > 0)
        ? "alpha : beta (test rates) : Synonymous (alpha) and non-synonymous (beta) rates on tested branches\n"
          "    beta (bg)                : Non-synonymous rate on background branches"
        : "alpha : beta (rates)      : Synonymous (alpha) and non-synonymous (beta) substitution rates";

    std::cout << Console::bold("  Table Legend:\n")
              << Console::muted("  • " + rate_desc + "\n"
                                "  • Rates Sparkline     : α (cyan) : β (green), scaled by rate magnitude:\n"
                                "                          ▱▱▱ zero (0.0) | ▰▱▱ low (<0.5) | ▰▰▱ moderate (0.5..3.0) | ▰▰▰ elevated (>=3.0)\n"
                                "  • Selection           : ▲ Diversifying (beta > alpha, p <= " + std::to_string(pvalue_threshold) + ") | "
                                "▼ Purifying (beta < alpha, p <= " + std::to_string(pvalue_threshold) + ")\n"
                                "  • Substitutions       : A ↔ B (undirected exchange), (k×) recurrence, [2nt/3nt] multiple hits\n"
                                "                          ⚡ high perplexity substitution (Shannon surprise >= 10.0 bits; rare jump / short branch)\n\n");

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
    std::string branch_summary = (fel.num_background > 0) 
        ? (std::to_string(fel.num_tested) + " test, " + std::to_string(fel.num_background) + " background")
        : ("All (" + std::to_string(fel.num_tested) + " branches)");

    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Branch Set", branch_summary},
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

    Provenance prov = build_cli_provenance(opts, cli_cmd, start_iso, total_sec);
    nlohmann::json out_j = (json_format == JSONFormat::Legacy)
        ? fel.to_legacy_json(alignment_file, "")
        : fel.to_modern_json(alignment_file, "", prov);
    if (!write_cli_json_and_html(opts, out_j)) {
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_fel(argc, argv);
}
#endif
