#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/meme.hpp"

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
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file or 'auto' (optional; auto-infers ML tree via TreeEngine if omitted or 'auto')\n"
              << "  " << Console::brand("--ohmic") << "              Build Ohmic Phylogenetic Network with Wheatstone bridges\n"
              << "  " << Console::brand("--hyphaeon") << "           Opt-in: HyphAeon 384D transformer conductances (needs Python aeon-core)\n"
              << "  " << Console::brand("--no-hyphaeon") << "        Use empirical codon conductances only (default)\n"
              << "  " << Console::brand("--embeddings") << " <file>  Path to precomputed [N x 384] HyphAeon embeddings (.bin/.csv/.tsv)\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  " << Console::brand("--branches") << " <spec>    Branches to test: All (default), Internal, Leaves, Unlabeled, {TAG}, or names\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.MEME.json)\n"
              << "  " << Console::brand("--html") << " <file>        Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--pvalue") << " <float>     P-value significance threshold (default: 0.1)\n"
              << "  " << Console::brand("--json-format") << " <type> Output JSON schema format: 'modern' (v3.0, default) or 'legacy' (Datamonkey 2.5)\n"
              << "  " << Console::brand("--legacy-json") << "        Emit legacy Datamonkey/HyPhy 2.5 JSON schema\n"
              << "  " << Console::brand("--modern-json") << "        Emit modern HyPhy 3.0 JSON schema with provenance and columnar data\n"
              << "  " << Console::brand("--full-model") << "         Perform branch length re-optimization under full codon model (default)\n"
              << "  " << Console::brand("--quick") << "              Disable full branch re-optimization (proportional branch scaling)\n"
              << "  " << Console::brand("--rates") << " <N>          Number of site omega classes (supported: 2, default)\n"
              << "  " << Console::brand("--multiple-hits") << " <m>  Multi-hit mode (supported: None; use hyphy_busted for multi-hit mixtures)\n"
              << "  " << Console::brand("--resample") << " <N>       Parametric bootstrap replicates (supported: 0)\n"
              << "  " << Console::brand("--impute-states") << " <b > Ancestral state imputation mode (supported: No)\n"
              << "  " << Console::brand("--precision") << " <mode>   Optimization precision ('standard' or 'fast'/'reduced' [= --quick])\n"
              << "  " << Console::brand("--all-sites") << "          Display all codon sites in console table (default: significant only)\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk --output cd2.MEME.json\n"
              << "  " << prog << " --alignment data/cd2.fna --tree auto --output cd2.MEME.json\n"
              << "  " << prog << " --alignment data/cd2.fna --ohmic --output cd2.MEME.json\n"
              << "  " << prog << " --alignment tests/data/adh.nex --threads 8 --legacy-json\n\n";
}

int run_meme(int argc, char* argv[]) {
    CommonCLIOptions opts;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (parse_common_cli_arg(i, argc, argv, opts)) {
            if (opts.has_error) {
                return 1;
            }
        } else if (arg == "--test" && i + 1 < argc) {
            opts.branches = argv[++i];
        } else if (arg == "--rates" && i + 1 < argc) {
            int r = std::stoi(argv[++i]);
            if (r != 2) {
                std::cerr << Console::danger("Error: ") << "hyphy_meme implements the canonical 2-rate MEME mixture (--rates 2); got --rates " << r << ".\n";
                return 1;
            }
        } else if ((arg == "--multiple-hits" || arg == "--site-multihit") && i + 1 < argc) {
            std::string mh = argv[++i];
            if (mh != "None" && mh != "none" && mh != "No" && mh != "no") {
                std::cerr << Console::danger("Error: ") << "hyphy_meme supports single-hit MG94xREV only (--multiple-hits None); got '" << mh
                          << "'. Use hyphy_busted --multiple-hits " << mh << " for multi-hit episodic selection testing.\n";
                return 1;
            }
        } else if (arg == "--resample" && i + 1 < argc) {
            int n_res = std::stoi(argv[++i]);
            if (n_res > 0) {
                std::cerr << Console::danger("Error: ") << "Parametric bootstrap (--resample > 0) is not supported in hyphy_meme.\n";
                return 1;
            }
        } else if (arg == "--impute-states") {
            std::string val = "Yes";
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                val = argv[++i];
            }
            if (val != "No" && val != "no" && val != "false" && val != "0") {
                std::cerr << Console::danger("Error: ") << "--impute-states is not supported as a standalone output block in hyphy_meme.\n";
                return 1;
            }
        } else if (arg == "--precision" && i + 1 < argc) {
            std::string prec = argv[++i];
            if (prec == "fast" || prec == "faster" || prec == "reduced") {
                opts.full_model = false;
            } else if (prec == "standard" || prec == "default") {
                opts.full_model = true;
            } else {
                std::cerr << Console::danger("Error: ") << "--precision must be 'standard', 'fast', or 'reduced' (got '" << prec << "').\n";
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

    finalize_common_cli_options(opts, ".MEME.json");
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
    const auto& branches_spec = opts.branches;
    const Scalar pvalue_threshold = opts.pvalue_threshold;
    const bool full_model = opts.full_model;
    const bool all_sites = opts.all_sites;
    const bool show_progress = opts.show_progress;
    const bool force_progress = opts.force_progress;
    const JSONFormat json_format = opts.json_format;

    Tree::BranchSelection branch_sel;
    try {
        branch_sel = tree.select_branches(branches_spec, "", {}, false);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }

    Panel::print_card("Dataset & Phylogeny", {
        {"Alignment File", alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons) + " (" + std::to_string(aln.num_codons * 3) + " nt)"},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree Source", tree_source + " (" + std::to_string(tree.num_nodes()) + " nodes, " + std::to_string(tree.num_leaves()) + " leaves)"},
        {"Tested Branches", branches_spec + " (" + std::to_string(branch_sel.num_selected) + " test, " + std::to_string(branch_sel.num_unselected) + " background)"},
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
    auto analyzer = MEMEAnalyzer::create_and_fit(tree, std::move(aln), pvalue_threshold, progress_cb, full_model, branches_spec);
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
    auto render_rate_spark = [](double alpha, double beta_plus) -> std::string {
        auto segs = [](double v) -> int {
            if (v <= 1e-4) return 0;
            if (v < 0.5) return 1;
            if (v < 3.0) return 2;
            return 3;
        };
        int na = segs(alpha);
        int nb = segs(beta_plus);

        std::string sa_filled, sa_empty, sb_filled, sb_empty;
        for (int i = 0; i < na; ++i) sa_filled += "▰";
        for (int i = 0; i < 3 - na; ++i) sa_empty += "▱";
        for (int i = 0; i < nb; ++i) sb_filled += "▰";
        for (int i = 0; i < 3 - nb; ++i) sb_empty += "▱";

        if (Console::color_enabled()) {
            return "\033[38;5;75m" + sa_filled + "\033[0m" +
                   "\033[2m" + sa_empty + "\033[0m" +
                   "\033[38;5;244m:\033[0m" +
                   "\033[38;5;208m" + sb_filled + "\033[0m" +
                   "\033[2m" + sb_empty + "\033[0m";
        } else {
            return sa_filled + sa_empty + ":" + sb_filled + sb_empty;
        }
    };

    Table table;
    table.add_column("Codon", Table::Align::Right, 6);
    table.add_column("alpha", Table::Align::Right, 8);
    table.add_column("beta1 (p1)", Table::Align::Right, 15);
    table.add_column("beta+ (p+)", Table::Align::Right, 16);
    table.add_column("Rates", Table::Align::Right, 9);
    table.add_column("LRT", Table::Align::Right, 7);
    table.add_column("p-value", Table::Align::Right, 9);
    table.add_column("# Branches", Table::Align::Right, 11);
    table.add_column("Selection", Table::Align::Left, 14);
    table.add_column("Substitutions", Table::Align::Left, 40);

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
        a_ss << std::fixed << std::setprecision(2) << r->alpha;
        b1_ss << std::fixed << std::setprecision(2) << r->beta1 << " (" << std::setprecision(2) << r->p1 << ")";
        bp_ss << std::fixed << std::setprecision(2) << r->beta_plus << " (" << std::setprecision(2) << r->p_plus << ")";
        lrt_ss << std::fixed << std::setprecision(2) << r->lrt;

        if (r->p_value < 0.0001) {
            p_ss << "<0.0001";
        } else {
            p_ss << std::fixed << std::setprecision(4) << r->p_value;
        }

        bool is_episodic = (r->p_value <= pvalue_threshold && r->beta_plus > r->alpha);
        std::string sel_label = is_episodic ? Console::brand("▲ Episodic") : Console::muted("-");
        std::string row_color = is_episodic ? "\033[38;5;208m" : ""; // Braun Signal Orange highlight

        std::string subs_str = r->subs.summary_string.empty() ? "Conserved" : r->subs.summary_string;

        table.add_row({
            std::to_string(r->site_index + 1),
            a_ss.str(),
            b1_ss.str(),
            bp_ss.str(),
            render_rate_spark(r->alpha, r->beta_plus),
            lrt_ss.str(),
            p_ss.str(),
            std::to_string(r->branches_under_selection),
            sel_label,
            subs_str
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

    std::cout << Console::bold("  Table Legend:\n")
              << Console::muted("  • alpha : beta1 (p1) : beta+ (p+) : Synonymous (alpha), constrained/neutral non-synonymous (beta1 with weight p1),\n"
                                "                                      and episodic non-synonymous (beta+ with weight p+ = 1 - p1)\n"
                                "  • Rates Sparkline                 : α (cyan) : β+ (orange), scaled by rate magnitude:\n"
                                "                                      ▱▱▱ zero (0.0) | ▰▱▱ low (<0.5) | ▰▰▱ moderate (0.5..3.0) | ▰▰▰ elevated (>=3.0)\n"
                                "  • Selection                       : ▲ Episodic diversifying selection (beta+ > alpha, p <= " + std::to_string(pvalue_threshold) + ")\n"
                                "  • # Branches                      : Number of branches with Empirical Bayes Factor (EBF) >= 100 for beta+\n"
                                "  • Substitutions                   : A ↔ B (undirected exchange), (k×) recurrence, [2nt/3nt] multiple hits\n"
                                "                                      ⚡ high perplexity substitution (Shannon surprise >= 10.0 bits; rare jump / short branch)\n\n");

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

    std::string json_desc = (json_format == JSONFormat::Legacy) ? " [Legacy]" : " [Modern v3.0]";
    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Tested Codon Sites", std::to_string(results.size())},
        {"Episodic Diversifying Sites", 
         Console::success(std::to_string(positive_sites.size())) + " sites (" + 
         fmt_pct(positive_sites.size(), results.size()) + ") [p <= " + std::to_string(pvalue_threshold) + "]"},
        {"Global dN/dS (beta/alpha)", omega_ss.str()},
        {"Nucleotide GTR Log-L", gtr_ss.str()},
        {"Global MG94xREV Log-L", mg_ss.str()},
        {"JSON Output File", output_file + json_desc},
        {"Total Execution Time", time_ss.str()}
    };

    Panel::print_summary_card("MEME Episodic Selection Analysis Summary", sum_items, conclusion, !positive_sites.empty());

    analyzer.input_filepath = alignment_file;
    Provenance prov = build_cli_provenance(opts, cli_cmd, start_iso, total_sec);
    nlohmann::json out_j = (json_format == JSONFormat::Legacy)
        ? analyzer.to_legacy_json()
        : analyzer.to_modern_json(prov);
    if (!write_cli_json_and_html(opts, out_j)) {
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_meme(argc, argv);
}
#endif
