#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/analyses/forecast.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace hyphy::core;
using namespace hyphy::analyses;

static void print_forecast_banner() {
    Panel::print_banner(
        "FORECAST",
        "Mutational Fitness & Evolutionary Trajectory Engine",
        "Reversible Site-Specific Codon Landscapes via FEL-SRV + G+S-PRIME + SynREVCodon + AncAA",
        "HyPhy 3 Evolutionary Forecasting & DMS Prediction • v3.0.0"
    );
}

static void print_forecast_usage(const char* prog) {
    print_forecast_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>          Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>               Path to Newick tree file or 'auto' (default: embedded or auto-inferred)\n"
              << "  " << Console::brand("--code") << " <name>               Genetic code (default: Universal)\n"
              << "  " << Console::brand("--branches") << " <spec>           Lineage partition: 'Internal' (default), 'All', 'Leaves', or tag\n"
              << "  " << Console::brand("--property-set") << " <set>        Biophysical set: 3PROP (default), 5PROP, LCAP, Atchley, 4PROP, 2PROP\n"
              << "  " << Console::brand("--property-file") << " <file>      Custom amino acid properties JSON file\n"
              << "  " << Console::brand("--shrinkage") << " <val>           Empirical Bayes Ridge penalty on site lambdas (default: 0.5)\n"
              << "  " << Console::brand("--prior") << " <mode>              Shrinkage center: 'global' (default, G-PRIME EB) or 'zero'\n"
              << "  " << Console::brand("--syn-model") << " <model>         Synonymous model: SynREVCodon (default), SynREV, SynREV2g, MG94, or TSV\n"
              << "  " << Console::brand("--syn-shrinkage") << " <val>       Log-normal EB regularization on synonymous pair rates (default: 0.5)\n"
              << "  " << Console::brand("--ancestral-exponent") << " <val>  Exponent on internal-node AA posterior tilt (default: 0.5; 0 disables)\n"
              << "  " << Console::brand("--horizon") << " <float>           Forecast horizon in expected codon substitutions/site (default: 0.05)\n"
              << "  " << Console::brand("--reference") << " <taxon>         Focal/wildtype taxon name (default: ML trunk consensus codon)\n"
              << "  " << Console::brand("--dms-tsv") << " <file>            Write tidy [L x 20] amino-acid DMS & stationary preference TSV\n"
              << "  " << Console::brand("--codon-tsv") << " <file>          Write tidy [L x 61] codon instantaneous rate & horizon log-prob TSV\n"
              << "  " << Console::brand("--output") << " <file>             Path to output JSON file (default: <alignment>.FORECAST.json)\n"
              << "  " << Console::brand("--threads") << " <N>               Number of OpenMP worker threads (0 = auto)\n"
              << "  " << Console::brand("--help, -h") << "                  Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment benchmarks/data/cd2.fna --tree benchmarks/data/cd2.nwk --dms-tsv cd2.dms.tsv\n"
              << "  " << prog << " --alignment H3N2_HA.fa --branches Internal --property-set 3PROP --syn-model SynREVCodon --reference A/Perth/16/2009\n\n";
}

int run_forecast(int argc, char* argv[]) {
    CommonCLIOptions cli;
    cli.branches = "Internal";

    ForecastSettings settings;
    std::string dms_tsv_file;
    std::string codon_tsv_file;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--branches" || arg == "--test-branches") && i + 1 < argc) {
            settings.branches = argv[++i];
            cli.branches = settings.branches;
            continue;
        }
        if ((arg == "--pvalue" || arg == "-p") && i + 1 < argc) {
            settings.p_value_threshold = parse_probability("--pvalue", argv[++i]);
            cli.pvalue_threshold = settings.p_value_threshold;
            continue;
        }

        auto status = parse_common_cli_arg(i, argc, argv, cli);
        if (status == CommonArgParseStatus::Handled) continue;
        if (status == CommonArgParseStatus::Error) return 1;

        if ((arg == "--property-set" || arg == "--properties") && i + 1 < argc) {
            settings.property_set = argv[++i];
        } else if (arg == "--property-file" && i + 1 < argc) {
            settings.custom_property_file = argv[++i];
        } else if (arg == "--shrinkage" && i + 1 < argc) {
            settings.shrinkage_lambda = parse_double_full("--shrinkage", argv[++i]);
        } else if ((arg == "--prior" || arg == "--shrinkage-prior") && i + 1 < argc) {
            settings.shrinkage_prior = argv[++i];
            std::string lp = settings.shrinkage_prior;
            std::transform(lp.begin(), lp.end(), lp.begin(), ::tolower);
            if (lp != "global" && lp != "g-prime" && lp != "gprime" && lp != "zero") {
                std::cerr << Console::danger("Error: unknown --prior '") << settings.shrinkage_prior
                          << "'. Expected 'global' or 'zero'.\n";
                return 1;
            }
        } else if ((arg == "--syn-model" || arg == "--mss") && i + 1 < argc) {
            settings.syn_model = argv[++i];
        } else if (arg == "--syn-shrinkage" && i + 1 < argc) {
            settings.syn_shrinkage = parse_double_full("--syn-shrinkage", argv[++i]);
        } else if ((arg == "--ancestral-exponent" || arg == "--anc-exp") && i + 1 < argc) {
            settings.ancestral_exponent = parse_double_full("--ancestral-exponent", argv[++i]);
        } else if (arg == "--ancestral-pseudocount" && i + 1 < argc) {
            settings.ancestral_pseudocount = parse_double_full("--ancestral-pseudocount", argv[++i]);
        } else if (arg == "--horizon" && i + 1 < argc) {
            settings.horizon = parse_double_full("--horizon", argv[++i]);
            if (settings.horizon <= 0.0) {
                std::cerr << Console::danger("Error: --horizon must be > 0.\n");
                return 1;
            }
        } else if ((arg == "--reference" || arg == "--ref") && i + 1 < argc) {
            settings.reference_taxon = argv[++i];
        } else if (arg == "--dms-tsv" && i + 1 < argc) {
            dms_tsv_file = argv[++i];
        } else if (arg == "--codon-tsv" && i + 1 < argc) {
            codon_tsv_file = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            print_forecast_usage(argv[0]);
            return 0;
        } else {
            std::cerr << Console::danger("Unknown argument: ") << arg << "\n";
            print_forecast_usage(argv[0]);
            return 1;
        }
    }

    if (cli.alignment_file.empty()) {
        std::cerr << Console::danger("Error: Missing required --alignment argument.\n\n");
        print_forecast_usage(argv[0]);
        return 1;
    }

    finalize_common_cli_options(cli, ".FORECAST.json");
    print_forecast_banner();

    LoadedCodonDataset ds;
    if (!load_codon_dataset(cli, ds, false)) {
        return 1;
    }

    std::cout << "\n" << Console::bold("Configuration:") << "\n"
              << "  " << Console::muted("Alignment file : ") << Console::bold(cli.alignment_file) << "\n"
              << "  " << Console::muted("Tree source    : ") << Console::bold(ds.tree_source) << "\n"
              << "  " << Console::muted("Sequences/Sites: ") << Console::bold(std::to_string(ds.aln.num_taxa) + " taxa × " + std::to_string(ds.aln.num_codons) + " codons") << "\n"
              << "  " << Console::muted("Branches       : ") << Console::bold(settings.branches) << "\n"
              << "  " << Console::muted("Property model : ") << Console::bold(settings.custom_property_file.empty() ? settings.property_set : settings.custom_property_file)
              << Console::success(" [prior=" + settings.shrinkage_prior + ", λ=" + std::to_string(settings.shrinkage_lambda) + "]") << "\n"
              << "  " << Console::muted("Synonymous MSS : ") << Console::bold(settings.syn_model) << "\n"
              << "  " << Console::muted("AncAA exponent : ") << Console::bold(std::to_string(settings.ancestral_exponent)) << "\n"
              << "  " << Console::muted("Horizon Δt     : ") << Console::bold(std::to_string(settings.horizon) + " expected codon subs/site") << "\n\n";

    auto start_iso = Provenance::current_iso8601();
    ForecastAnalyzer analyzer(ds.tree, ds.aln, settings);
    ForecastResult res = analyzer.run(cli.show_progress, cli.force_progress);

    std::cout << Console::bold("Global Evolutionary Landscape Summary:") << "\n"
              << "  " << Console::muted("• Effective branches : ") << Console::bold(res.effective_branches) << "\n"
              << "  " << Console::muted("• Global ω (G-PRIME) : ") << std::fixed << std::setprecision(3) << res.global_omega
              << " (lnL = " << std::setprecision(2) << res.prime_global_log_l << ")\n"
              << "  " << Console::muted("• Global λ prior     : [");
    for (size_t k = 0; k < res.global_lambdas.size(); ++k) {
        std::cout << (k > 0 ? ", " : "") << res.property_names[k] << "=" << std::fixed << std::setprecision(2) << res.global_lambdas[k];
    }
    std::cout << "]\n";

    if (res.effective_syn_model != "MG94") {
        std::cout << "  " << Console::muted("• Synonymous ") << std::setw(8) << std::left << res.effective_syn_model << Console::muted(": ")
                  << "ΔlnL = +" << std::fixed << std::setprecision(2) << std::max<Scalar>(0.0, res.mss_log_l - res.mss_baseline_log_l)
                  << " (" << res.mss_num_classes << " classes, LRT=" << std::setprecision(2) << res.mss_omnibus_lrt
                  << ", p=" << std::scientific << std::setprecision(2) << res.mss_omnibus_p_value << std::fixed << ")\n";
    }

    // Top most mutable & most constrained sites summary table
    std::vector<size_t> order(res.site_results.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return res.site_results[a].expected_substitutions_horizon > res.site_results[b].expected_substitutions_horizon;
    });

    Table tab;
    tab.add_column("Site", Table::Align::Right);
    tab.add_column("Focal", Table::Align::Center);
    tab.add_column("α (SRV)", Table::Align::Right);
    tab.add_column("β (Test)", Table::Align::Right);
    tab.add_column("ω_eff", Table::Align::Right);
    tab.add_column("N_eff(AA)", Table::Align::Right);
    tab.add_column("P(mut|Δt)", Table::Align::Right);
    tab.add_column("Top Predicted AA Shift", Table::Align::Left);

    auto fmt = [](Scalar v, int p) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(p) << v;
        return ss.str();
    };

    size_t top_show = std::min<size_t>(10, order.size());
    for (size_t idx = 0; idx < top_show; ++idx) {
        size_t s = order[idx];
        const auto& sr = res.site_results[s];
        int focal_a = aa_char_to_index(sr.focal_aa);
        int best_alt_a = -1;
        Scalar best_dms = -1e30;
        for (int a = 0; a < 20; ++a) {
            if (a == focal_a) continue;
            if (res.aa_dms_matrix(s, a) > best_dms) {
                best_dms = res.aa_dms_matrix(s, a);
                best_alt_a = a;
            }
        }
        std::string top_mut = (best_alt_a >= 0)
            ? (std::string(1, sr.focal_aa) + "->" + res.aa_alphabet[best_alt_a] + " (" + fmt(best_dms, 2) + ")")
            : "-";

        tab.add_row({
            std::to_string(sr.site_index + 1),
            sr.focal_codon + "(" + std::string(1, sr.focal_aa) + ")",
            fmt(sr.alpha, 2),
            fmt(sr.beta, 2),
            fmt(sr.effective_omega, 2),
            fmt(sr.effective_aa_count, 2),
            fmt(sr.expected_substitutions_horizon, 4),
            top_mut
        });
    }

    std::cout << "\n" << Console::bold("Top Forecasted Mutable Sites (Horizon Δt = " + fmt(settings.horizon, 3) + "):") << "\n";
    tab.print();
    std::cout << "\n";

    Provenance prov = build_cli_provenance(cli, join_cli_command(argc, argv), start_iso, res.runtime_seconds);
    std::string json_str = res.to_modern_json(analyzer.tree, ds.aln, prov, cli.alignment_file).dump(2);
    if (!write_cli_json_and_html(cli, json_str)) {
        return 1;
    }

    if (!dms_tsv_file.empty()) {
        res.write_dms_tsv(dms_tsv_file);
        std::cout << Console::success("✓") << " Wrote [L × 20] DMS & stationary AA matrix TSV to " << Console::bold(dms_tsv_file) << "\n";
    }
    if (!codon_tsv_file.empty()) {
        res.write_codon_tsv(codon_tsv_file);
        std::cout << Console::success("✓") << " Wrote [L × " << res.num_sense_codons << "] codon forecast TSV to " << Console::bold(codon_tsv_file) << "\n";
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    try {
        return run_forecast(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}
#endif
