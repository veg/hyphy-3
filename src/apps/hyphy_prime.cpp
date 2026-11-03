#include "hyphy/core/types.hpp"
#include "hyphy/core/cli_common.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/prime.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <fstream>

using namespace hyphy::core;
using namespace hyphy::analyses;
using namespace hyphy::opt;

static void print_prime_banner() {
    Panel::print_banner(
        "PRIME",
        "PRoperty Informed Models of Evolution",
        "Inference of Amino Acid Biophysical Property Selection & Evolutionary Constraints",
        "Citation: PRoperty Informed Models of Evolution (PRIME) • v3.0.0"
    );
}

static void print_prime_usage(const char* prog) {
    print_prime_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment") << " <file>       Path to codon alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--tree") << " <file>            Path to Newick tree file or 'auto' (optional; auto-inferred via TreeEngine if omitted or 'auto')\n"
              << "  " << Console::brand("--ohmic, --omhic") << "         Infer Ohmic network & Wheatstone bridges (empirical conductances)\n"
              << "  " << Console::brand("--hyphaeon") << "               Opt-in: HyphAeon 384D transformer conductances (needs Python aeon-core)\n"
              << "  " << Console::brand("--embeddings") << " <file>      Path to precomputed [N x D] taxon embeddings TSV/CSV for Ohmic\n"
              << "  " << Console::brand("--code") << " <name>            Genetic code (default: Universal; e.g. Vertebrate-mtDNA)\n"
              << "  " << Console::brand("--property-set") << " <set>     Property set: 5PROP (default), Atchley, 3PROP, 4PROP, 2PROP\n"
              << "  " << Console::brand("--property-file") << " <file>   Custom amino acid properties JSON file\n"
              << "  " << Console::brand("--shrinkage") << " <val>        Empirical Bayes Ridge shrinkage lambda (default: 0.5; 0 to disable)\n"
              << "  " << Console::brand("--branches") << " <spec>        Branches to test: 'All' (default), 'Internal', 'Leaves', or regex\n"
              << "  " << Console::brand("--threads") << " <N>            Number of OpenMP worker threads\n"
              << "  " << Console::brand("--output") << " <file>          Path to output JSON file (default: <alignment>.PRIME.json)\n"
              << "  " << Console::brand("--html") << " <file>            Also write a self-contained interactive HTML report\n"
              << "  " << Console::brand("--pvalue") << " <float>         Omnibus FDR threshold for site reporting (default: 0.10)\n"
              << "  " << Console::brand("--pvalue2") << " <float>        Property-level adjusted p-value threshold (default: 0.10)\n"
              << "  " << Console::brand("--format") << " <fmt>           Output format: 'modern' (default), 'legacy', or 'both'\n"
              << "  " << Console::brand("--all-sites") << "              Display all codon sites in console table (default: significant only)\n"
              << "  " << Console::brand("--progress") << "               Force interactive progress bar\n"
              << "  " << Console::brand("--no-progress") << "            Disable progress bar\n"
              << "  " << Console::brand("--help, -h") << "               Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " --alignment data/cd2.fna --tree auto --output cd2.PRIME.json\n"
              << "  " << prog << " --alignment tests/data/bglobin.nex --property-set 3PROP\n"
              << "  " << prog << " --alignment data/yokoyama.nex --property-set Atchley --shrinkage 1.0\n\n";
}

int run_prime(int argc, char* argv[]) {
    CommonCLIOptions cli;
    std::string property_set = "5PROP";
    std::string property_file = "";
    std::string branch_spec = "All";
    std::string output_format = "modern";
    Scalar pvalue_threshold = 0.10;
    Scalar stage2_pvalue = 0.10;
    Scalar shrinkage_lambda = 0.5;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--format" || arg == "--json-format") && i + 1 < argc) {
            std::string fmt = argv[++i];
            if (fmt == "both") {
                output_format = "both";
                cli.json_format = JSONFormat::ModernV3;
            } else {
                JSONFormat jf;
                if (!parse_json_format(fmt, jf)) {
                    std::cerr << Console::danger("Error: unknown --json-format '") << fmt
                              << "'. Expected 'modern' or 'legacy'.\n";
                    return 1;
                }
                cli.json_format = jf;
                output_format = (jf == JSONFormat::Legacy) ? "legacy" : "modern";
            }
            continue;
        }

        auto status = parse_common_cli_arg(i, argc, argv, cli);
        if (status == CommonArgParseStatus::Handled) {
            output_format = (cli.json_format == JSONFormat::Legacy) ? "legacy" : "modern";
            continue;
        }
        if (status == CommonArgParseStatus::Error) return 1;

        if (arg == "--property-set" && i + 1 < argc) {
            property_set = argv[++i];
        } else if (arg == "--property-file" && i + 1 < argc) {
            property_file = argv[++i];
        } else if (arg == "--shrinkage" && i + 1 < argc) {
            shrinkage_lambda = std::stod(argv[++i]);
        } else if (arg == "--no-shrinkage") {
            shrinkage_lambda = 0.0;
        } else if ((arg == "--branches" || arg == "--test-branches") && i + 1 < argc) {
            branch_spec = argv[++i];
        } else if ((arg == "--pvalue" || arg == "-p") && i + 1 < argc) {
            pvalue_threshold = std::stod(argv[++i]);
        } else if (arg == "--pvalue2" && i + 1 < argc) {
            stage2_pvalue = std::stod(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            print_prime_usage(argv[0]);
            return 0;
        } else {
            std::cerr << Console::danger("Unknown argument: ") << arg << "\n";
            print_prime_usage(argv[0]);
            return 1;
        }
    }

    if (cli.alignment_file.empty()) {
        std::cerr << Console::danger("Error: Missing required --alignment argument.\n\n");
        print_prime_usage(argv[0]);
        return 1;
    }

    finalize_common_cli_options(cli, ".PRIME.json");
    print_prime_banner();

    std::cout << "\n" << Console::bold("Configuration:") << "\n"
              << "  " << Console::muted("Alignment file : ") << Console::bold(cli.alignment_file) << "\n";
    if (!cli.tree_file.empty()) {
        std::cout << "  " << Console::muted("Tree file      : ") << Console::bold(cli.tree_file) << "\n";
    }
    std::cout << "  " << Console::muted("Genetic code   : ") << Console::bold(cli.code_name) << "\n"
              << "  " << Console::muted("Property set   : ") << Console::bold(property_file.empty() ? property_set : ("Custom: " + property_file)) << "\n"
              << "  " << Console::muted("Shrinkage λ    : ") << Console::bold(std::to_string(shrinkage_lambda))
              << (shrinkage_lambda > 0.0 ? Console::success(" [Empirical Bayes Ridge Enabled]") : Console::warning(" [Unpenalized MLE]")) << "\n"
              << "  " << Console::muted("Branches       : ") << Console::bold(branch_spec) << "\n\n";

    LoadedCodonDataset ds;
    if (!load_codon_dataset(cli, ds, false)) {
        return 1;
    }
    auto& aln = ds.aln;
    auto& tree = ds.tree;
    const auto& tree_source = ds.tree_source;

    std::cout << "  " << Console::muted("Tree source    : ") << Console::bold(tree_source) << "\n"
              << "  " << Console::muted("Sequences      : ") << Console::bold(std::to_string(aln.num_taxa)) << "\n"
              << "  " << Console::muted("Codon sites    : ") << Console::bold(std::to_string(aln.num_codons)) << "\n"
              << "  " << Console::muted("Unique patterns: ") << Console::bold(std::to_string(aln.patterns.size())) << "\n\n";

    PRIMESettings settings;
    settings.property_set = property_set;
    settings.custom_property_file = property_file;
    settings.shrinkage_lambda = shrinkage_lambda;
    settings.p_value_threshold = pvalue_threshold;
    settings.stage2_p_value = stage2_pvalue;
    settings.test_branches = branch_spec;

    std::cout << Console::bold("Step 1:") << " Nucleotide GTR & Codon Baseline Fitting...\n";
    auto prime_analyzer = PRIMEAnalyzer::create_and_fit(tree, aln, settings);

    std::cout << "  " << Console::muted("• GTR Log-L    : ") << std::fixed << std::setprecision(2) << prime_analyzer.gtr_log_l << "\n"
              << "  " << Console::muted("• Baseline Tree: ") << std::fixed << std::setprecision(4) << prime_analyzer.baseline_tree_length << " substitutions/site\n";

    if (settings.fit_global) {
        std::cout << "  " << Console::muted("• G-PRIME Log-L: ") << std::fixed << std::setprecision(2) << prime_analyzer.global_fit.log_likelihood
                  << " (beta=" << std::setprecision(3) << prime_analyzer.global_fit.beta << ")\n";
        std::cout << "  " << Console::muted("• Global λ     : [");
        for (size_t k = 0; k < prime_analyzer.global_fit.lambdas.size(); ++k) {
            std::cout << (k > 0 ? ", " : "") << std::fixed << std::setprecision(2) << prime_analyzer.global_fit.lambdas[k];
        }
        std::cout << "]\n";
    }

    std::cout << "\n" << Console::bold("Step 2:") << " Running Site-by-Site S-PRIME Model Fitting...\n";
    auto results = prime_analyzer.run(nullptr, cli.show_progress, cli.force_progress);

    auto format_float = [](Scalar val, int prec) -> std::string {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(prec) << val;
        return ss.str();
    };
    auto format_pvalue = [](Scalar p) -> std::string {
        std::ostringstream ss;
        if (p < 0.0001) ss << "<0.0001";
        else if (p < 0.001) ss << std::scientific << std::setprecision(2) << p;
        else ss << std::fixed << std::setprecision(4) << p;
        return ss.str();
    };

    // Format Console Results Table
    Table tab;
    tab.add_column("Site", Table::Align::Right);
    tab.add_column("alpha", Table::Align::Right);
    tab.add_column("beta", Table::Align::Right);
    tab.add_column("Omnibus LRT", Table::Align::Right);
    tab.add_column("p-value", Table::Align::Right);
    tab.add_column("q-value", Table::Align::Right);
    tab.add_column("Significant Properties (Effect / Adj p)", Table::Align::Left);
    tab.add_column("Subs", Table::Align::Right);
    tab.add_column("AA", Table::Align::Right);

    for (const auto& sr : results.site_results) {
        bool is_sig = sr.is_significant(pvalue_threshold);
        if (!cli.all_sites && !is_sig) continue;

        std::string site_str = std::to_string(sr.site_index + 1);
        std::string a_str = format_float(sr.alpha, 2);
        std::string b_str = format_float(sr.beta, 2);
        std::string lrt_str = format_float(sr.omnibus_lrt, 2);
        std::string p_str = format_pvalue(sr.omnibus_p_value);
        std::string q_str = format_pvalue(sr.omnibus_q_value);

        std::string prop_summary = "None";
        std::vector<std::string> sig_props;
        for (const auto& p_inf : sr.properties) {
            if (p_inf.classification != "Neutral") {
                std::string short_name = p_inf.property_name;
                if (short_name.find("Hydrophobicity") != std::string::npos) short_name = "Hydro";
                else if (short_name.find("Volume") != std::string::npos) short_name = "Vol";
                else if (short_name.find("Isoelectric") != std::string::npos) short_name = "pI";
                else if (short_name.find("AlphaHelix") != std::string::npos) short_name = "Helix";
                else if (short_name.find("BetaSheet") != std::string::npos) short_name = "Sheet";

                std::string entry = short_name + " (" +
                    (p_inf.classification == "Conserved" ? Console::accent("Cons") : Console::danger("Chg")) +
                    ", p=" + format_pvalue(p_inf.adjusted_p_value) + ")";
                sig_props.push_back(entry);
            }
        }
        if (!sig_props.empty()) {
            prop_summary = "";
            for (size_t s = 0; s < sig_props.size(); ++s) {
                if (s > 0) prop_summary += "; ";
                prop_summary += sig_props[s];
            }
        }

        std::string subs_str = std::to_string(sr.nonsyn_subs);
        std::string aa_str = std::to_string(sr.unique_aa);

        if (is_sig) {
            site_str = Console::bold(site_str);
            q_str = Console::success(q_str);
        }

        tab.add_row({site_str, a_str, b_str, lrt_str, p_str, q_str, prop_summary, subs_str, aa_str});
    }

    std::cout << "\n";
    tab.print();

    // Summary Box
    std::cout << "\n" << Console::bold("PRIME Analysis Summary:") << "\n";
    std::cout << "  " << Console::muted("Total sites analyzed      : ") << Console::bold(std::to_string(results.total_sites)) << "\n"
              << "  " << Console::muted("Invariant sites           : ") << Console::bold(std::to_string(results.invariant_sites)) << "\n"
              << "  " << Console::muted("Omnibus significant (q<=") << pvalue_threshold << "): "
              << Console::bold(std::to_string(results.significant_sites)) << " sites\n";

    std::cout << "  " << Console::muted("Property Breakdown        :\n");
    for (const auto& p : results.properties) {
        size_t c_cons = results.conserved_counts.at(p.name);
        size_t c_chg = results.changing_counts.at(p.name);
        std::cout << "    • " << Console::brand(p.name) << ": "
                  << Console::accent(std::to_string(c_cons) + " conserved") << ", "
                  << Console::danger(std::to_string(c_chg) + " changing") << "\n";
    }
    std::cout << "  " << Console::muted("Total runtime             : ") << std::fixed << std::setprecision(2) << results.runtime_seconds << "s\n\n";

    // Write Output JSON File (and optional HTML)
    std::string cli_cmd = join_cli_command(argc, argv);
    std::string final_json_str;
    if (output_format == "legacy") {
        final_json_str = results.to_legacy_json(tree, aln, cli.alignment_file).dump(2);
    } else if (output_format == "both") {
        auto mj = results.to_modern_json(tree, aln, cli_cmd, cli.alignment_file);
        mj["legacy_mle"] = results.to_legacy_json(tree, aln, cli.alignment_file)["MLE"];
        final_json_str = mj.dump(2);
    } else {
        final_json_str = results.to_modern_json(tree, aln, cli_cmd, cli.alignment_file).dump(2);
    }

    if (!write_cli_json_and_html(cli, final_json_str)) {
        return 1;
    }

    return 0;
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    return run_prime(argc, argv);
}
#endif
