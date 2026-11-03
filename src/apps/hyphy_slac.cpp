#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/slac.hpp"

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
        "SLAC",
        "Single-Likelihood Ancestor Counting",
        "Counting-based selection test using joint maximum likelihood ancestral reconstruction",
        "Citation: Kosakovsky Pond & Frost (2005) MBE 22(5): 1208-1222 • v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>   Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>        Path to Newick tree file or 'auto' (optional; auto-inferred via TreeEngine if omitted or 'auto')\n"
              << "  " << Console::brand("--ohmic, --omhic") << "     Infer Ohmic network & Wheatstone bridges (empirical conductances)\n"
              << "  " << Console::brand("--hyphaeon") << "           Opt-in: HyphAeon 384D transformer conductances (needs Python aeon-core)\n"
              << "  " << Console::brand("--embeddings") << " <file>  Path to precomputed [N x D] taxon embeddings TSV/CSV for Ohmic\n"
              << "  " << Console::brand("--code") << " <name>        Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  " << Console::brand("--branches") << " <set>     Branches to test: 'All' (default), 'Internal', 'Leaves', 'Unlabeled', tag name, or comma-separated names\n"
              << "  " << Console::brand("--branch-regex") << " <reg> Regular expression to match branch names to test\n"
              << "  " << Console::brand("--samples") << " <N>        Number of ancestral state samples for confidence intervals (default: 100, 0 to skip)\n"
              << "  " << Console::brand("--pvalue") << " <float>     P-value significance threshold (default: 0.1)\n"
              << "  " << Console::brand("--full-model") << "         Perform branch length re-optimization under full codon model (default)\n"
              << "  " << Console::brand("--quick") << "              Disable full branch re-optimization (proportional branch scaling)\n"
              << "  " << Console::brand("--all-sites") << "          Display all codon sites in console table (default: significant only)\n"
              << "  " << Console::brand("--threads") << " <N>        Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>      Path to output JSON file (default: <alignment>.SLAC.json)\n"
              << "  " << Console::brand("--html") << " <file>        Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--json-format") << " <type> Output JSON schema format: 'modern' (v3.0, default) or 'legacy' (Datamonkey 2.5)\n"
              << "  " << Console::brand("--legacy-json") << "        Emit legacy Datamonkey/HyPhy 2.5 JSON schema\n"
              << "  " << Console::brand("--modern-json") << "        Emit modern HyPhy 3.0 JSON schema with provenance and columnar data\n"
              << "  " << Console::brand("--progress") << "           Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "        Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "           Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk --output cd2.SLAC.json\n"
              << "  " << prog << " --alignment data/cd2.fna --tree auto --output cd2.SLAC.json\n"
              << "  " << prog << " --alignment tests/data/adh.nex --samples 100 --pvalue 0.05\n"
              << "  " << prog << " --alignment data/cd2.fna --branches Internal --legacy-json\n\n";
}

int run_slac(int argc, char* argv[]) {
    CommonCLIOptions cli;
    std::string test_branches = "All";
    std::string test_branch_regex = "";
    size_t samples = 100;
    Scalar pvalue_threshold = 0.1;
    bool full_model = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto status = parse_common_cli_arg(i, argc, argv, cli);
        if (status == CommonArgParseStatus::Handled) continue;
        if (status == CommonArgParseStatus::Error) return 1;

        if ((arg == "--branches" || arg == "--test-branches") && i + 1 < argc) {
            test_branches = argv[++i];
        } else if ((arg == "--branch-regex" || arg == "--test-branch-regex") && i + 1 < argc) {
            test_branch_regex = argv[++i];
        } else if (arg == "--samples" && i + 1 < argc) {
            samples = static_cast<size_t>(std::stoul(argv[++i]));
        } else if ((arg == "--pvalue" || arg == "-p") && i + 1 < argc) {
            pvalue_threshold = std::stod(argv[++i]);
        } else if (arg == "--full-model") {
            full_model = true;
        } else if (arg == "--quick" || arg == "--no-full-model") {
            full_model = false;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (cli.alignment_file.empty()) {
        std::cerr << "Error: --alignment argument is required.\n\n";
        print_usage(argv[0]);
        return 1;
    }

    finalize_common_cli_options(cli, ".SLAC.json");
    print_banner();

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string start_iso = Provenance::current_iso8601();
    std::string cli_cmd = join_cli_command(argc, argv);

    std::cout << Console::bold("Loading alignment:") << " " << cli.alignment_file << "\n";
    LoadedCodonDataset ds;
    if (!load_codon_dataset(cli, ds, false)) {
        return 1;
    }
    auto& gcode = ds.gcode;
    auto& aln = ds.aln;
    auto& tree = ds.tree;
    std::cout << Console::bold("Phylogeny:") << " " << ds.tree_source << "\n";

    try {
        (void)tree.select_branches(test_branches, test_branch_regex);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error in branch selection: ") << e.what() << "\n";
        return 1;
    }

    std::cout << Console::muted("  • Sequences : ") << Console::accent(std::to_string(aln.num_taxa)) << "\n"
              << Console::muted("  • Codon Sites: ") << Console::accent(std::to_string(aln.num_codons)) << "\n"
              << Console::muted("  • Tree Leaves: ") << Console::accent(std::to_string(tree.num_leaves())) << "\n"
              << Console::muted("  • Genetic Code: ") << Console::accent(gcode->name) << "\n"
              << Console::muted("  • Tested Set: ") << Console::accent(test_branches) << "\n"
              << Console::muted("  • Ancestral Samples: ") << Console::accent(std::to_string(samples)) << "\n\n";

    SLACSettings settings;
    settings.p_value_threshold = pvalue_threshold;
    settings.samples = samples;
    settings.test_branches = test_branches;
    settings.test_branch_regex = test_branch_regex;
    settings.full_model = full_model;
    settings.all_sites = cli.all_sites;

    std::function<void(const std::string&, double)> progress_cb = nullptr;
    if (!cli.show_progress && !cli.force_progress) {
        progress_cb = [](const std::string& msg, double) {
            std::cout << "  " << Console::muted("↳ ") << msg << "...\n";
        };
    }

    SLACAnalyzer analyzer = SLACAnalyzer::create_and_fit(tree, aln, settings, progress_cb);
    SLACResult result = analyzer.run(progress_cb, cli.show_progress, cli.force_progress);

    // Results table
    Table table;
    table.add_column("Site", Table::Align::Right, 6);

    if (samples > 0) {
        table.add_column("S [med, IQR]", Table::Align::Center, 18);
        table.add_column("N [med, IQR]", Table::Align::Center, 18);
        table.add_column("dS [med, IQR]", Table::Align::Center, 18);
        table.add_column("dN [med, IQR]", Table::Align::Center, 18);
        table.add_column("dN - dS", Table::Align::Right, 10);
        table.add_column("P [dN/dS > 1]", Table::Align::Right, 13);
        table.add_column("P [dN/dS < 1]", Table::Align::Right, 13);
        table.add_column("Selection", Table::Align::Left, 16);
    } else {
        table.add_column("S", Table::Align::Right, 6);
        table.add_column("N", Table::Align::Right, 6);
        table.add_column("E[S]", Table::Align::Right, 8);
        table.add_column("E[N]", Table::Align::Right, 8);
        table.add_column("dS", Table::Align::Right, 8);
        table.add_column("dN", Table::Align::Right, 8);
        table.add_column("dN - dS", Table::Align::Right, 10);
        table.add_column("P [dN/dS > 1]", Table::Align::Right, 13);
        table.add_column("P [dN/dS < 1]", Table::Align::Right, 13);
        table.add_column("Selection", Table::Align::Left, 16);
    }

    std::vector<size_t> sig_indices;
    for (size_t s = 0; s < result.total_sites; ++s) {
        Scalar p_pos = result.resolved_sites[s].p_positive;
        Scalar p_neg = result.resolved_sites[s].p_negative;
        if (p_pos <= pvalue_threshold || p_neg <= pvalue_threshold) {
            sig_indices.push_back(s);
        }
    }

    std::vector<size_t> display_indices;
    if (cli.all_sites) {
        display_indices.resize(result.total_sites);
        std::iota(display_indices.begin(), display_indices.end(), 0);
    } else if (!sig_indices.empty()) {
        display_indices = sig_indices;
    } else {
        std::vector<size_t> sorted_indices(result.total_sites);
        std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
        std::sort(sorted_indices.begin(), sorted_indices.end(), [&](size_t a, size_t b) {
            Scalar min_p_a = std::min(result.resolved_sites[a].p_positive, result.resolved_sites[a].p_negative);
            Scalar min_p_b = std::min(result.resolved_sites[b].p_positive, result.resolved_sites[b].p_negative);
            return min_p_a < min_p_b;
        });
        size_t show_n = std::min(sorted_indices.size(), static_cast<size_t>(5));
        display_indices.assign(sorted_indices.begin(), sorted_indices.begin() + show_n);
    }

    auto fmt_val = [](Scalar v, int prec = 3) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(prec) << v;
        return ss.str();
    };

    auto fmt_iqr = [&](Scalar med, Scalar lo, Scalar hi) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << med << " [" << lo << "-" << hi << "]";
        return ss.str();
    };

    for (size_t s : display_indices) {
        const auto& r = result.resolved_sites[s];
        std::string sel_label = Console::muted("Neutral");
        if (r.p_positive <= pvalue_threshold) {
            sel_label = Console::success("▲ Diversifying");
        } else if (r.p_negative <= pvalue_threshold) {
            sel_label = Console::info("▼ Purifying");
        }

        std::string ds_str = r.has_dS ? fmt_val(r.dS) : "-";
        std::string dn_str = r.has_dN ? fmt_val(r.dN) : "-";
        std::string diff_str = r.has_dN_minus_dS ? fmt_val(r.dN_minus_dS) : "-";

        if (samples > 0 && !r.quantiles_median.empty()) {
            std::string s_iqr = fmt_iqr(r.quantiles_median[2], r.quantiles_2_5[2], r.quantiles_97_5[2]);
            std::string n_iqr = fmt_iqr(r.quantiles_median[3], r.quantiles_2_5[3], r.quantiles_97_5[3]);
            std::string ds_iqr = r.has_dS ? fmt_iqr(r.quantiles_median[5], r.quantiles_2_5[5], r.quantiles_97_5[5]) : "-";
            std::string dn_iqr = r.has_dN ? fmt_iqr(r.quantiles_median[6], r.quantiles_2_5[6], r.quantiles_97_5[6]) : "-";

            table.add_row({
                std::to_string(s + 1),
                s_iqr,
                n_iqr,
                ds_iqr,
                dn_iqr,
                diff_str,
                fmt_val(r.p_positive, 4),
                fmt_val(r.p_negative, 4),
                sel_label
            });
        } else {
            table.add_row({
                std::to_string(s + 1),
                fmt_val(r.S, 1),
                fmt_val(r.N, 1),
                fmt_val(r.ES, 2),
                fmt_val(r.EN, 2),
                ds_str,
                dn_str,
                diff_str,
                fmt_val(r.p_positive, 4),
                fmt_val(r.p_negative, 4),
                sel_label
            });
        }
    }

    std::cout << "\n";
    table.print();

    if (!cli.all_sites) {
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

    // Summary Panel
    auto end_time = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(end_time - start_time).count();

    std::string branch_summary = (result.num_background_branches > 0)
        ? (std::to_string(result.num_tested_branches) + " test, " + std::to_string(result.num_background_branches) + " background")
        : ("All (" + std::to_string(result.num_tested_branches) + " branches)");

    std::vector<std::pair<std::string, std::string>> sum_items = {
        {"Branch Set", branch_summary},
        {"Tested Codon Sites", std::to_string(result.total_sites)},
        {"Diversifying Selection (Positive)", 
         Console::success(std::to_string(result.diversifying_sites)) + " sites (" + 
         format_percentage(result.diversifying_sites, result.total_sites) + ") [p <= " + std::to_string(pvalue_threshold) + "]"},
        {"Purifying Selection (Negative)", 
         Console::info(std::to_string(result.purifying_sites)) + " sites (" + 
         format_percentage(result.purifying_sites, result.total_sites) + ") [p <= " + std::to_string(pvalue_threshold) + "]"},
        {"Neutral Sites", std::to_string(result.neutral_sites) + " sites (" + format_percentage(result.neutral_sites, result.total_sites) + ")"},
        {"Global dN/dS (omega)", fmt_val(result.global_omega, 4)},
        {"Nucleotide GTR Log-L", fmt_val(result.gtr_log_l, 2) + " (AICc: " + fmt_val(result.gtr_aicc, 2) + ")"},
        {"Global MG94xREV Log-L", fmt_val(result.baseline_log_l, 2) + " (AICc: " + fmt_val(result.baseline_aicc, 2) + ")"},
        {"Total Tree Length", fmt_val(result.total_tree_length, 4)},
        {"Analysis Runtime", fmt_val(total_sec, 2) + " s"}
    };

    Panel::print_summary_card("SLAC Analysis Summary", sum_items);

    nlohmann::json out_json;
    if (cli.json_format == JSONFormat::Legacy) {
        out_json = result.to_legacy_json(cli.alignment_file, tree.to_newick());
    } else {
        Provenance prov = build_cli_provenance(cli, cli_cmd, start_iso, total_sec, {
            {"branches", test_branches},
            {"samples", std::to_string(samples)},
            {"pvalue", std::to_string(pvalue_threshold)},
            {"full_model", full_model ? "true" : "false"}
        });
        out_json = result.to_modern_json(cli.alignment_file, tree.to_newick(), prov);
    }

    if (!write_cli_json_and_html(cli, out_json.dump(2))) {
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_slac(argc, argv);
}
#endif
