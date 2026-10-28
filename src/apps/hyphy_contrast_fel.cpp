#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/contrast_fel.hpp"

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
        "CONTRAST-FEL",
        "Differential Selection Between Branch Sets",
        "Empirical Bayes Ridge Shrinkage & Hypothesis Testing Across Lineages",
        "Citation: Mol Biol Evol (2021) 38(3):1184–1198 • v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file or 'auto' (optional; auto-infers ML tree via TreeEngine if omitted or 'auto')\n"
              << "  " << Console::brand("--ohmic, --omhic") << "     Infer Ohmic network & Wheatstone bridges (empirical conductances)\n"
              << "  " << Console::brand("--hyphaeon") << "           Opt-in: HyphAeon 384D transformer conductances (needs Python aeon-core)\n"
              << "  " << Console::brand("--embeddings") << " <file>  Path to precomputed [N x D] taxon embeddings TSV/CSV for Ohmic\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  " << Console::brand("--branches") << " <spec>    Branch sets to contrast: 'auto' (default, uses tree tags), 'Internal,Leaves', or tag names\n"
              << "  " << Console::brand("--shrinkage") << " <val>    Empirical Bayes Ridge shrinkage lambda (default: 1.0; 0 to disable)\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.CONTRAST-FEL.json)\n"
              << "  " << Console::brand("--html") << " <file>        Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--json-format") << " <fmt>  JSON schema format: modern (default) or legacy (HyPhy 2.5)\n"
              << "  " << Console::brand("--legacy-json") << "        Output legacy HyPhy 2.5 JSON schema\n"
              << "  " << Console::brand("--pvalue") << " <float>     Nominal p-value threshold for site display (default: 0.05)\n"
              << "  " << Console::brand("--qvalue") << " <float>     FDR q-value threshold for reporting (default: 0.20)\n"
              << "  " << Console::brand("--all-sites") << "          Display all codon sites in console table (default: significant only)\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " --alignment data/cd2.fna --tree auto --branches Internal,Leaves\n"
              << "  " << prog << " --alignment data/cd2.fna --branches Internal,Leaves --shrinkage 1.0\n"
              << "  " << prog << " --alignment test.nex --pvalue 0.05 --qvalue 0.10\n\n";
}

int run_contrast_fel(int argc, char* argv[]) {
    CommonCLIOptions cli;
    std::string branch_specs = "auto";
    Scalar pvalue_threshold = 0.05;
    Scalar qvalue_threshold = 0.20;
    Scalar shrinkage_lambda = 1.0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--branches" || arg == "--branch-set" || arg == "--branch-sets") && i + 1 < argc) {
            branch_specs = argv[++i];
            continue;
        } else if (arg == "--shrinkage" && i + 1 < argc) {
            shrinkage_lambda = std::stod(argv[++i]);
            continue;
        } else if (arg == "--no-shrinkage") {
            shrinkage_lambda = 0.0;
            continue;
        } else if ((arg == "--pvalue" || arg == "--p-value" || arg == "-p") && i + 1 < argc) {
            pvalue_threshold = std::stod(argv[++i]);
            continue;
        } else if ((arg == "--qvalue" || arg == "--q-value" || arg == "-q") && i + 1 < argc) {
            qvalue_threshold = std::stod(argv[++i]);
            continue;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        }

        auto status = parse_common_cli_arg(i, argc, argv, cli);
        if (status == CommonArgParseStatus::Handled) continue;
        if (status == CommonArgParseStatus::Error) return 1;

        std::cerr << Console::danger("Unknown argument: ") << arg << "\n";
        print_usage(argv[0]);
        return 1;
    }
    if (branch_specs == "auto" && cli.branches != "All" && !cli.branches.empty()) {
        branch_specs = cli.branches;
    }

    if (cli.alignment_file.empty()) {
        std::cerr << Console::danger("Error: --alignment is required.\n\n");
        print_usage(argv[0]);
        return 1;
    }

    finalize_common_cli_options(cli, ".CONTRAST-FEL.json");
    print_banner();

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string start_iso = Provenance::current_iso8601();
    std::string cli_cmd = join_cli_command(argc, argv);

    LoadedCodonDataset ds;
    if (!load_codon_dataset(cli, ds, false, /*allow_fallback_code=*/false)) {
        return 1;
    }
    auto& gcode = ds.gcode;
    auto& aln = ds.aln;
    auto& tree = ds.tree;
    const auto& tree_source = ds.tree_source;
    size_t num_codons = aln.num_codons;

    std::cout << Console::bold("Input Data:\n")
              << "  " << Console::muted("Alignment file : ") << Console::bold(cli.alignment_file) << "\n";
    if (!cli.tree_file.empty()) {
        std::cout << "  " << Console::muted("Tree file      : ") << Console::bold(cli.tree_file) << "\n";
    }
    std::cout << "  " << Console::muted("Genetic code   : ") << Console::bold(gcode->name) << "\n"
              << "  " << Console::muted("Shrinkage λ    : ") << Console::bold(std::to_string(shrinkage_lambda)) 
              << (shrinkage_lambda > 0.0 ? Console::success(" [Empirical Bayes Ridge Enabled]") : Console::warning(" [Unpenalized MLE]")) << "\n"
              << "  " << Console::muted("Tree source    : ") << Console::bold(tree_source) << "\n"
              << "  " << Console::muted("Sequences      : ") << Console::bold(std::to_string(aln.num_taxa)) << "\n"
              << "  " << Console::muted("Codon sites    : ") << Console::bold(std::to_string(aln.num_codons)) << "\n"
              << "  " << Console::muted("Unique patterns: ") << Console::bold(std::to_string(aln.patterns.size())) << "\n\n";

    ContrastFELSettings settings;
    settings.p_value_threshold = pvalue_threshold;
    settings.q_value_threshold = qvalue_threshold;
    settings.shrinkage_lambda = shrinkage_lambda;
    settings.branch_specs = branch_specs;

    auto progress_cb = [&](const std::string& task, [[maybe_unused]] double p) {
        std::cout << "  " << Console::muted("• ") << task << "...\n";
    };

    auto contrast_fel = ContrastFELAnalyzer::create_and_fit(std::move(tree), std::move(aln), settings, progress_cb);

    size_t K = contrast_fel.testable_groups.size();
    std::cout << "\n" << Console::bold("Branch Sets Partition:\n");
    for (size_t k = 0; k < K; ++k) {
        const auto& grp = contrast_fel.testable_groups[k];
        std::cout << "  " << Console::brand("Group " + std::to_string(k + 1) + " [" + grp.name + "]") 
                  << " : " << grp.branch_count << " branches (total length: " 
                  << std::fixed << std::setprecision(3) << grp.total_branch_length << ")\n";
    }
    if (contrast_fel.has_background) {
        std::cout << "  " << Console::muted("Background [Nuisance]") 
                  << " : " << contrast_fel.background_group.branch_count << " branches (total length: "
                  << std::fixed << std::setprecision(3) << contrast_fel.background_group.total_branch_length << ")\n";
    }
    std::cout << "\n";

    auto site_results = contrast_fel.run(cli.show_progress, cli.force_progress);

    // Sparkline renderer for 2 groups
    auto render_contrast_spark = [](double alpha, double b0, double b1) -> std::string {
        auto segs = [](double v) -> int {
            if (v <= 1e-4) return 0;
            if (v < 0.5) return 1;
            if (v < 3.0) return 2;
            return 3;
        };
        int na = segs(alpha);
        int nb0 = segs(b0);
        int nb1 = segs(b1);

        std::string sa_filled, sa_empty, sb0_filled, sb0_empty, sb1_filled, sb1_empty;
        for (int i = 0; i < na; ++i) sa_filled += "▰";
        for (int i = 0; i < 3 - na; ++i) sa_empty += "▱";
        for (int i = 0; i < nb0; ++i) sb0_filled += "▰";
        for (int i = 0; i < 3 - nb0; ++i) sb0_empty += "▱";
        for (int i = 0; i < nb1; ++i) sb1_filled += "▰";
        for (int i = 0; i < 3 - nb1; ++i) sb1_empty += "▱";

        std::string arrow = (b0 < b1) ? " ◄ " : ((b0 > b1) ? " ► " : " = ");

        if (Console::color_enabled()) {
            return "\033[38;5;75m" + sa_filled + "\033[0m\033[2m" + sa_empty + "\033[0m : " +
                   "\033[38;5;39m" + sb0_filled + "\033[0m\033[2m" + sb0_empty + "\033[0m" +
                   "\033[38;5;244m" + arrow + "\033[0m" +
                   "\033[38;5;208m" + sb1_filled + "\033[0m\033[2m" + sb1_empty + "\033[0m";
        } else {
            return sa_filled + sa_empty + " : " + sb0_filled + sb0_empty + arrow + sb1_filled + sb1_empty;
        }
    };

    // Console table formatting
    Table table;
    table.add_column("Codon", Table::Align::Right, 6);
    table.add_column("alpha", Table::Align::Right, 7);

    if (K == 2) {
        std::string col_rates = "beta [" + contrast_fel.testable_groups[0].name + " : " + contrast_fel.testable_groups[1].name + "]";
        table.add_column(col_rates, Table::Align::Right, 24);
        table.add_column("Contrast Spark", Table::Align::Center, 18);
    } else {
        table.add_column("beta rates (by group)", Table::Align::Right, 28);
    }

    table.add_column("PLRT", Table::Align::Right, 8);
    table.add_column("p-value", Table::Align::Right, 9);
    table.add_column("q-value", Table::Align::Right, 9);
    table.add_column("Contrast Shift", Table::Align::Left, 22);
    table.add_column("Substitutions by Group", Table::Align::Left, 38);

    std::vector<size_t> sig_indices;
    for (size_t c = 0; c < site_results.size(); ++c) {
        if (site_results[c].p_value <= pvalue_threshold) {
            sig_indices.push_back(c);
        }
    }

    std::vector<size_t> sites_to_display;
    if (cli.all_sites) {
        sites_to_display.resize(site_results.size());
        std::iota(sites_to_display.begin(), sites_to_display.end(), 0);
    } else if (!sig_indices.empty()) {
        sites_to_display = sig_indices;
    } else {
        std::vector<size_t> sorted_by_pval(site_results.size());
        std::iota(sorted_by_pval.begin(), sorted_by_pval.end(), 0);
        std::sort(sorted_by_pval.begin(), sorted_by_pval.end(), [&](size_t i, size_t j) {
            return site_results[i].p_value < site_results[j].p_value;
        });
        size_t n_top = std::min(sorted_by_pval.size(), size_t(5));
        sites_to_display.assign(sorted_by_pval.begin(), sorted_by_pval.begin() + n_top);
    }

    for (size_t c : sites_to_display) {
        const auto& r = site_results[c];
        std::ostringstream a_ss, p_ss, q_ss, lrt_ss;
        a_ss << std::fixed << std::setprecision(2) << r.alpha;
        lrt_ss << std::fixed << std::setprecision(2) << r.penalized_lrt;

        if (r.p_value < 0.0001) p_ss << "<0.0001";
        else p_ss << std::fixed << std::setprecision(4) << r.p_value;

        if (r.q_value < 0.0001) q_ss << "<0.0001";
        else q_ss << std::fixed << std::setprecision(4) << r.q_value;

        std::string tag_str = r.contrast_tag;
        if (r.p_value <= pvalue_threshold) {
            if (r.q_value <= qvalue_threshold) {
                tag_str = Console::success(tag_str);
            } else {
                tag_str = Console::accent(tag_str);
            }
        } else {
            tag_str = Console::muted("-");
        }

        std::string subs_str = r.subs.summary_string_partitioned;

        if (K == 2) {
            std::ostringstream b_ss;
            b_ss << std::fixed << std::setprecision(2) << r.beta[0] << " : " << std::setprecision(2) << r.beta[1];
            std::string spark_str = render_contrast_spark(r.alpha, r.beta[0], r.beta[1]);

            table.add_row({
                std::to_string(c + 1),
                a_ss.str(),
                b_ss.str(),
                spark_str,
                lrt_ss.str(),
                p_ss.str(),
                q_ss.str(),
                tag_str,
                subs_str
            });
        } else {
            std::ostringstream b_ss;
            for (size_t k = 0; k < K; ++k) {
                if (k > 0) b_ss << " / ";
                b_ss << std::fixed << std::setprecision(2) << r.beta[k];
            }

            table.add_row({
                std::to_string(c + 1),
                a_ss.str(),
                b_ss.str(),
                lrt_ss.str(),
                p_ss.str(),
                q_ss.str(),
                tag_str,
                subs_str
            });
        }
    }

    std::cout << "\n";
    table.print();

    if (!cli.all_sites) {
        if (!sig_indices.empty()) {
            std::cout << Console::muted("  Showing " + std::to_string(sig_indices.size()) + 
                                       " site(s) with differential selection at p <= " + std::to_string(pvalue_threshold) + 
                                       " out of " + std::to_string(num_codons) + 
                                       ". (Use --all-sites to display all)\n\n");
        } else {
            std::cout << Console::muted("  No sites reached nominal significance (p <= " + 
                                       std::to_string(pvalue_threshold) + 
                                       "). Displaying top 5 sites by p-value. (Use --all-sites to display all)\n\n");
        }
    }

    // Legend
    std::cout << Console::bold("  Table Legend:\n")
              << Console::muted("  • Rates Sparkline  : α (cyan) : β₁ (blue) ◄/► β₂ (orange), arrow indicates elevated non-synonymous rate\n"
                                "                       ▱▱▱ zero (0.0) | ▰▱▱ low (<0.5) | ▰▰▱ moderate (0.5..3.0) | ▰▰▰ elevated (>=3.0)\n"
                                "  • Regularization   : Empirical Bayes Ridge Shrinkage (λ = " + std::to_string(shrinkage_lambda) + 
                                ") penalizes spurious rate differences on sparse branches\n"
                                "  • Contrast Shift   : G1 ≪ G2 (Adaptive / Relaxed / Accelerated) identifies the biological direction of selection\n"
                                "  • Substitutions    : Group-partitioned codon changes, [2nt/3nt] multiple hits, ⚡ Shannon surprise >= 10.0 bits\n\n");

    // FDR Summary Box
    auto summary = contrast_fel.get_summary();
    std::cout << Console::bold("  False Discovery Rate (FDR) Summary:\n");
    if (summary.significant_q_sites > 0) {
        std::cout << "  " << Console::success("Found " + std::to_string(summary.significant_q_sites) + 
                                            " site(s) passing alignment-wide FDR threshold (q <= " + 
                                            std::to_string(qvalue_threshold) + "):\n\n");

        Table fdr_table;
        fdr_table.add_column("Codon", Table::Align::Right, 8);
        fdr_table.add_column("PLRT", Table::Align::Right, 10);
        fdr_table.add_column("p-value", Table::Align::Right, 12);
        fdr_table.add_column("q-value (FDR)", Table::Align::Right, 14);
        fdr_table.add_column("Directional Contrast", Table::Align::Left, 26);

        for (size_t c = 0; c < site_results.size(); ++c) {
            if (site_results[c].q_value <= qvalue_threshold) {
                const auto& r = site_results[c];
                std::ostringstream p_ss, q_ss, lrt_ss;
                lrt_ss << std::fixed << std::setprecision(2) << r.penalized_lrt;
                p_ss << std::scientific << std::setprecision(3) << r.p_value;
                q_ss << std::scientific << std::setprecision(3) << r.q_value;
                fdr_table.add_row({
                    std::to_string(c + 1),
                    lrt_ss.str(),
                    p_ss.str(),
                    q_ss.str(),
                    Console::success(r.contrast_tag)
                });
            }
        }
        fdr_table.print();
        std::cout << "\n";
    } else {
        std::cout << "  " << Console::muted("No sites pass alignment-wide FDR threshold of q <= " + 
                                           std::to_string(qvalue_threshold) + ".\n\n");
    }

    if (K > 2 && !site_results.empty() && !site_results[0].pairwise.empty()) {
        std::cout << Console::bold("  Pairwise Group Contrasts Summary (p <= " + std::to_string(pvalue_threshold) + "):\n\n");
        Table pw_table;
        pw_table.add_column("Group Pair", Table::Align::Left, 24);
        pw_table.add_column("Uncorrected p <= " + std::to_string(pvalue_threshold), Table::Align::Right, 22);
        pw_table.add_column("Holm-Bonferroni p <= " + std::to_string(pvalue_threshold), Table::Align::Right, 26);
        for (size_t pw_idx = 0; pw_idx < site_results[0].pairwise.size(); ++pw_idx) {
            const auto& pw0 = site_results[0].pairwise[pw_idx];
            size_t cnt_unadj = 0, cnt_adj = 0;
            for (const auto& r : site_results) {
                if (pw_idx < r.pairwise.size()) {
                    if (r.pairwise[pw_idx].p_value <= pvalue_threshold) cnt_unadj++;
                    if (r.pairwise[pw_idx].p_value_adjusted <= pvalue_threshold) cnt_adj++;
                }
            }
            pw_table.add_row({
                pw0.name_a + " vs " + pw0.name_b,
                std::to_string(cnt_unadj),
                std::to_string(cnt_adj)
            });
        }
        pw_table.print();
        std::cout << "\n";
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(end_time - start_time).count();

    Provenance prov = build_cli_provenance(cli, cli_cmd, start_iso, total_sec, {
        {"branches", branch_specs},
        {"shrinkage", std::to_string(shrinkage_lambda)},
        {"pvalue", std::to_string(pvalue_threshold)},
        {"qvalue", std::to_string(qvalue_threshold)}
    });

    try {
        contrast_fel.save_json(cli.output_file, cli.alignment_file, contrast_fel.tree.to_newick(), prov, cli.json_format);
        std::cout << Console::bold("Output:\n")
                  << "  " << Console::muted("Analysis results saved to: ") << Console::brand(cli.output_file) << "\n\n";
        if (!cli.html_file.empty()) {
            hyphy::report::write_html_report(cli.output_file, cli.html_file);
            std::cout << Console::success("✔") << " " << Console::bold("Saved HTML report to: ")
                      << Console::brand(cli.html_file) << "\n\n";
        }
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error writing output: ") << e.what() << "\n";
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    try {
        return run_contrast_fel(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}
#endif
