#include "hyphy/analyses/directional_fel.hpp"
#include "hyphy/core/cli_common.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace hyphy::core;
using namespace hyphy::analyses;

void print_dfel_banner() {
    Panel::print_banner(
        "DIRECTIONAL FEL (dFEL)",
        "Non-Stationary Fixed Effects Likelihood",
        "Site-level selection testing with auto-rooting & compositional heterogeneity screening",
        "HyPhy 3 • Sergei L. Kosakovsky Pond & HyPhy Team • v3.0.0"
    );
}

void print_dfel_usage(const char* prog) {
    print_dfel_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " " 
              << Console::accent("--alignment") << " <file> "
              << Console::muted("[OPTIONS]") << "\n\n";

    std::cout << Console::bold("Required Arguments:") << "\n"
              << "  " << Console::brand("--alignment, -a  ") << "  Path to in-frame codon alignment (FASTA or NEXUS)\n\n";

    std::cout << Console::bold("Model & Rooting Options:") << "\n"
              << "  " << Console::brand("--tree, -t       ") << "  Path to Newick phylogenetic tree or 'auto' (optional; auto-inferred via TreeEngine if omitted or 'auto')\n"
              << "  " << Console::brand("--ohmic, --omhic ") << "  Infer Ohmic network & Wheatstone bridges (empirical conductances)\n"
              << "  " << Console::brand("--hyphaeon    ") << "     Opt-in: HyphAeon 384D transformer conductances (needs Python aeon-core)\n"
              << "  " << Console::brand("--embeddings <f> ") << "  Path to precomputed [N x D] taxon embeddings TSV/CSV for Ohmic\n"
              << "  " << Console::brand("--root, -r       ") << "  Rooting strategy: 'auto' (ML auto-rooting, default), 'keep' (input root), or <taxon>\n"
              << "  " << Console::brand("--heterogeneity  ") << "  P-value threshold to trigger directional regimes (default: 0.05)\n"
              << "  " << Console::brand("--regimes, -k    ") << "  Compositional regimes: 'auto' (default: 3 if heterogeneous, 1 if homogeneous), 1, or 3\n"
              << "  " << Console::brand("--code, -c       ") << "  Genetic code name (default: Universal)\n"
              << "  " << Console::brand("--pvalue, -p     ") << "  Significance threshold for selection testing (default: 0.1)\n"
              << "  " << Console::brand("--quick          ") << "  Quick mode: estimate global tree scale and omega (skip full branch refit)\n\n";

    std::cout << Console::bold("Output & Performance Options:") << "\n"
              << "  " << Console::brand("--output, -o     ") << "  Output JSON filepath (default: <alignment>.dFEL.json)\n"
              << "  " << Console::brand("--html <file>    ") << "  Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--modern-json    ") << "  Export HyPhy 3.0 Modern JSON schema (default)\n"
              << "  " << Console::brand("--legacy-json    ") << "  Export backward-compatible HyPhy 2.5 legacy JSON format\n"
              << "  " << Console::brand("--threads, -T    ") << "  Number of OpenMP worker threads (default: system auto)\n"
              << "  " << Console::brand("--all-sites      ") << "  Print results table for all sites, not just significant\n"
              << "  " << Console::brand("--no-progress    ") << "  Disable interactive terminal visualizer\n"
              << "  " << Console::brand("--help, -h       ") << "  Show this help message\n\n";

    std::cout << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " -a tests/data/adh.nex --root auto --regimes 3\n"
              << "  " << prog << " -a tests/data/adh.nex --root Drosophila_melanogaster\n\n";
}

static std::string render_rate_spark(Scalar alpha, Scalar beta) {
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
}

int run_dfel(int argc, char* argv[]) {
    CommonCLIOptions cli;
    std::string root_strategy = "auto";
    std::string regimes_opt = "auto";
    Scalar pvalue_threshold = 0.1;
    Scalar het_threshold = 0.05;
    bool full_model = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto status = parse_common_cli_arg(i, argc, argv, cli);
        if (status == CommonArgParseStatus::Handled) continue;
        if (status == CommonArgParseStatus::Error) return 1;

        if ((arg == "--root" || arg == "-r") && i + 1 < argc) {
            root_strategy = argv[++i];
        } else if ((arg == "--regimes" || arg == "-k") && i + 1 < argc) {
            regimes_opt = argv[++i];
        } else if ((arg == "--heterogeneity" || arg == "--heterogeneity-threshold") && i + 1 < argc) {
            het_threshold = std::stod(argv[++i]);
        } else if ((arg == "--pvalue" || arg == "-p") && i + 1 < argc) {
            pvalue_threshold = std::stod(argv[++i]);
        } else if (arg == "--quick" || arg == "--no-full-model") {
            full_model = false;
        } else if (arg == "--full-model") {
            full_model = true;
        } else if (arg == "--help" || arg == "-h") {
            print_dfel_usage(argv[0]);
            return 0;
        } else {
            std::cerr << Console::danger("Unknown argument: ") << arg << "\n\n";
            print_dfel_usage(argv[0]);
            return 1;
        }
    }

    if (cli.alignment_file.empty()) {
        std::cerr << Console::danger("Error: --alignment is required.\n\n");
        print_dfel_usage(argv[0]);
        return 1;
    }

    finalize_common_cli_options(cli, ".dFEL.json");
    print_dfel_banner();

    auto start_time = std::chrono::high_resolution_clock::now();

    LoadedCodonDataset ds;
    if (!load_codon_dataset(cli, ds, false)) {
        return 1;
    }
    auto& gcode = ds.gcode;
    auto& aln = ds.aln;
    auto& tree = ds.tree;
    const auto& tree_source = ds.tree_source;

    Panel::print_card("Dataset & Phylogeny", {
        {"Alignment File", cli.alignment_file},
        {"Sequences / Taxa", std::to_string(aln.num_taxa)},
        {"Codon Sites", std::to_string(aln.num_codons) + " (" + std::to_string(aln.num_codons * 3) + " nt)"},
        {"Unique Patterns", std::to_string(aln.patterns.size())},
        {"Tree Source", tree_source + " (" + std::to_string(tree.num_nodes()) + " nodes)"},
        {"Root Strategy", root_strategy},
        {"Genetic Code", gcode->name + " (" + std::to_string(gcode->num_sense_codons) + " sense codons)"}
    });

    DirectionalFELSettings settings;
    settings.p_value_threshold = pvalue_threshold;
    settings.heterogeneity_threshold = het_threshold;
    settings.root_strategy = root_strategy;
    settings.full_model = full_model;
    if (regimes_opt == "1") settings.num_regimes = 1;
    else if (regimes_opt == "3") settings.num_regimes = 3;
    else settings.num_regimes = 0; // auto

    // Step 1: Upfront Compositional Heterogeneity Test
    Panel::print_step(1, 4, "Compositional Heterogeneity Screening", "Contingency Chi-Square Tests Across Taxa");
    auto comp_res = DirectionalFELAnalyzer::run_composition_test(aln, het_threshold);

    Table comp_tbl;
    comp_tbl.add_column("Partition", Table::Align::Left, 25);
    comp_tbl.add_column("Chi-Square", Table::Align::Right, 12);
    comp_tbl.add_column("df", Table::Align::Right, 6);
    comp_tbl.add_column("p-value", Table::Align::Right, 10);
    comp_tbl.add_column("Composition Status", Table::Align::Left, 24);
    auto fmt_p = [](Scalar p) {
        if (p < 0.0001) return std::string("<0.0001");
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4) << p;
        return ss.str();
    };

    std::string overall_status = (comp_res.p_value_overall < het_threshold)
        ? Console::danger("Heterogeneous (p < 0.05)")
        : Console::success("Homogeneous (p >= 0.05)");

    std::string cp3_status = (comp_res.p_value_cp3 < het_threshold)
        ? Console::danger("Heterogeneous (p < 0.05)")
        : Console::success("Homogeneous (p >= 0.05)");

    comp_tbl.add_row({"Overall (All Codon Pos)", SiteVisualizer::format_scalar(comp_res.chi2_overall, 2), std::to_string(comp_res.df_overall), fmt_p(comp_res.p_value_overall), overall_status});
    comp_tbl.add_row({"Codon Position 1 (CP1)", SiteVisualizer::format_scalar(comp_res.chi2_cp1, 2), std::to_string(comp_res.df_cp1), fmt_p(comp_res.p_value_cp1), (comp_res.p_value_cp1 < het_threshold) ? Console::danger("Heterogeneous") : Console::muted("Homogeneous")});
    comp_tbl.add_row({"Codon Position 2 (CP2)", SiteVisualizer::format_scalar(comp_res.chi2_cp2, 2), std::to_string(comp_res.df_cp2), fmt_p(comp_res.p_value_cp2), (comp_res.p_value_cp2 < het_threshold) ? Console::danger("Heterogeneous") : Console::muted("Homogeneous")});
    comp_tbl.add_row({"Codon Position 3 (CP3)", SiteVisualizer::format_scalar(comp_res.chi2_cp3, 2), std::to_string(comp_res.df_cp3), fmt_p(comp_res.p_value_cp3), cp3_status});
    comp_tbl.print();

    // Show GC extremes
    if (!comp_res.taxa.empty()) {
        auto sorted_taxa = comp_res.taxa;
        std::sort(sorted_taxa.begin(), sorted_taxa.end(), [](const TaxonComposition& a, const TaxonComposition& b) {
            return a.gc_cp3 < b.gc_cp3;
        });
        std::cout << "  " << Console::muted("CP3 GC Range : ")
                  << Console::info(sorted_taxa.front().name + " (" + SiteVisualizer::format_scalar(sorted_taxa.front().gc_cp3 * 100.0, 1) + "%)")
                  << Console::muted(" to ")
                  << Console::success(sorted_taxa.back().name + " (" + SiteVisualizer::format_scalar(sorted_taxa.back().gc_cp3 * 100.0, 1) + "%)")
                  << "\n\n";
    }

    // Step 2: Nucleotide GTR Model & Auto-Rooting
    Panel::print_step(2, 4, "Nucleotide GTR & Directional Auto-Rooting", "Evaluating candidate root branches");
    GTRFitter gtr_fitter(tree, aln);
    auto gtr_res = gtr_fitter.fit();

    auto root_res = DirectionalFELAnalyzer::auto_select_root(gtr_res.tree, aln, gtr_res.params, settings);

    Table root_tbl;
    root_tbl.add_column("Rank", Table::Align::Right, 5);
    root_tbl.add_column("Candidate Branch", Table::Align::Left, 24);
    root_tbl.add_column("Type", Table::Align::Left, 10);
    root_tbl.add_column("Log-L", Table::Align::Right, 12);
    root_tbl.add_column("Δ ln L", Table::Align::Right, 10);
    root_tbl.add_column("Root Weight", Table::Align::Right, 12);
    root_tbl.add_column("Status", Table::Align::Left, 16);

    size_t display_max = std::min(size_t(5), root_res.candidates.size());
    for (size_t i = 0; i < display_max; ++i) {
        const auto& c = root_res.candidates[i];
        bool is_selected = (i == root_res.chosen_candidate_idx);
        std::string status_str = is_selected ? Console::success("★ Selected Root") : "";
        std::string weight_str = SiteVisualizer::format_scalar(c.weight * 100.0, 1) + "%";

        root_tbl.add_row({
            std::to_string(c.rank),
            c.branch_name,
            c.is_terminal ? "Terminal" : "Internal",
            SiteVisualizer::format_scalar(c.log_l, 2),
            SiteVisualizer::format_scalar(c.delta_log_l, 2),
            weight_str,
            status_str
        });
    }
    root_tbl.print();
    std::cout << "  " << Console::muted("Auto-Rooting Decision: ") 
              << Console::bold(root_res.chosen_branch_name) 
              << " (Weight: " << SiteVisualizer::format_scalar(root_res.chosen_weight * 100.0, 1) << "%)\n\n";

    // Step 3: Directional MG94 Model Refinement
    Panel::print_step(3, 4, "Refining Global Directional MG94 Model", "Lineage-specific composition regimes");
    auto dfel = DirectionalFELAnalyzer::create_and_fit(
        root_res.best_tree,
        std::move(aln),
        settings,
        nullptr
    );

    Table reg_tbl;
    reg_tbl.add_column("Regime ID", Table::Align::Right, 10);
    reg_tbl.add_column("Regime Name", Table::Align::Left, 14);
    reg_tbl.add_column("Directional Drift (ψ)", Table::Align::Right, 22);
    reg_tbl.add_column("Mean GC", Table::Align::Right, 10);
    reg_tbl.add_column("Branches Assigned", Table::Align::Right, 18);
    for (const auto& reg : dfel.regimes) {
        reg_tbl.add_row({
            std::to_string(reg.id),
            reg.name,
            SiteVisualizer::format_scalar(reg.psi, 3),
            SiteVisualizer::format_scalar(reg.mean_gc * 100.0, 1) + "%",
            std::to_string(reg.branch_count)
        });
    }
    reg_tbl.print();

    std::cout << "  " << Console::muted("Global Log-L : ") << Console::bold(std::to_string(dfel.global_log_l))
              << Console::muted("  |  AICc: ") << std::fixed << std::setprecision(2) << dfel.global_aicc
              << Console::muted("  |  Tree Length: ") << std::fixed << std::setprecision(3) << dfel.baseline_tree_length << "\n\n";

    // Step 4: Site-Level Testing
    Panel::print_step(4, 4, "Site-by-Site Directional Selection Testing", "Testing alpha vs beta with directed substitutions");
    auto results = dfel.run(cli.show_progress, cli.force_progress);

    // Candidates Table
    std::vector<size_t> sig_indices;
    for (size_t s = 0; s < results.size(); ++s) {
        if (results[s].p_value <= pvalue_threshold) {
            sig_indices.push_back(s);
        }
    }

    std::vector<size_t> display_indices;
    if (cli.all_sites) {
        display_indices.resize(results.size());
        for (size_t s = 0; s < results.size(); ++s) display_indices[s] = s;
    } else if (!sig_indices.empty()) {
        display_indices = sig_indices;
    } else {
        // Show top 5 lowest p-values
        std::vector<size_t> all_idx(results.size());
        for (size_t s = 0; s < results.size(); ++s) all_idx[s] = s;
        std::sort(all_idx.begin(), all_idx.end(), [&](size_t a, size_t b) {
            return results[a].p_value < results[b].p_value;
        });
        size_t n_top = std::min(size_t(5), all_idx.size());
        display_indices.assign(all_idx.begin(), all_idx.begin() + n_top);
    }

    Table sites_tbl;
    sites_tbl.add_column("Site", Table::Align::Right, 6);
    sites_tbl.add_column("Rates (α : β)", Table::Align::Right, 24);
    sites_tbl.add_column("LRT", Table::Align::Right, 8);
    sites_tbl.add_column("p-value", Table::Align::Right, 10);
    sites_tbl.add_column("Selection", Table::Align::Left, 16);
    sites_tbl.add_column("Directed Substitutions", Table::Align::Left, 36);

    for (size_t s : display_indices) {
        const auto& r = results[s];
        std::string sel_label = Console::muted("Neutral");
        if (r.p_value <= pvalue_threshold) {
            if (r.beta > r.alpha) {
                sel_label = Console::success("▲ Diversifying");
            } else if (r.alpha > r.beta) {
                sel_label = Console::info("▼ Purifying");
            }
        }

        std::ostringstream ab_ss;
        ab_ss << std::fixed << std::setprecision(2) << r.alpha << " : " << std::setprecision(2) << r.beta
              << "  " << render_rate_spark(r.alpha, r.beta);

        sites_tbl.add_row({
            std::to_string(s + 1),
            ab_ss.str(),
            SiteVisualizer::format_scalar(r.lrt, 2),
            fmt_p(r.p_value),
            sel_label,
            r.subs.summary_string
        });
    }

    std::cout << "\n";
    sites_tbl.print();

    std::cout << Console::bold("  Table Legend:\n")
              << Console::muted("  • Rates (α : β)       : Synonymous (alpha) and non-synonymous (beta) substitution rates\n"
                                "  • Selection           : ▲ Diversifying (beta > alpha) | ▼ Purifying (beta < alpha) at p <= " + std::to_string(pvalue_threshold) + "\n"
                                "  • Directed Subs       : A → B directed ancestral transition along rooted phylogeny\n"
                                "                          [GC↑] = mutation towards G/C | [AT↑] = mutation towards A/T\n"
                                "                          ⚡ high perplexity substitution (Shannon surprise >= 10.0 bits)\n\n");

    // Summary Card
    auto summary = dfel.get_summary();
    auto end_time = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(end_time - start_time).count();

    Panel::print_card("Directional FEL Execution Summary", {
        {"Compositional Heterogeneity", comp_res.is_heterogeneous ? "Detected (Regimes Activated)" : "Homogeneous"},
        {"Selected Root Branch", root_res.chosen_branch_name + " (Support: " + SiteVisualizer::format_scalar(root_res.chosen_weight * 100.0, 1) + "%)"},
        {"Directional Regimes", std::to_string(dfel.regimes.size())},
        {"Tested Codon Sites", std::to_string(summary.tested_sites)},
        {"Diversifying Selection (Positive)", Console::success(std::to_string(summary.positively_selected)) + " sites (" + format_percentage(summary.positively_selected, summary.tested_sites) + ")"},
        {"Purifying Selection (Negative)", Console::info(std::to_string(summary.negatively_selected)) + " sites (" + format_percentage(summary.negatively_selected, summary.tested_sites) + ")"},
        {"Neutral Sites", std::to_string(summary.neutral_sites) + " sites (" + format_percentage(summary.neutral_sites, summary.tested_sites) + ")"},
        {"Global MG94 Fit", SiteVisualizer::format_scalar(dfel.global_log_l, 2) + " (AICc: " + SiteVisualizer::format_scalar(dfel.global_aicc, 2) + ")"},
        {"Output JSON", cli.output_file},
        {"Total Runtime", SiteVisualizer::format_scalar(total_sec, 2) + " s"}
    });

    // Export JSON
    try {
        dfel.export_json(cli.output_file, cli.json_format, cli.alignment_file, dfel.tree.to_newick(true));
        std::cout << Console::success("✓ ") << "Analysis results written to " << Console::bold(cli.output_file) << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error exporting JSON: ") << e.what() << "\n";
        return 1;
    }

    if (!cli.html_file.empty()) {
        try {
            hyphy::report::write_html_report(cli.output_file, cli.html_file);
            std::cout << Console::success("✔") << " " << Console::bold("Saved HTML report to: ")
                      << Console::brand(cli.html_file) << "\n\n";
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error writing HTML report: ") << e.what() << "\n";
            return 1;
        }
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    try {
        return run_dfel(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}
#endif
