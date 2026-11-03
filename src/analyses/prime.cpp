// Out-of-line JSON emitter definitions for prime.
// Moved from include/hyphy/analyses/prime.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/prime.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json PRIMEResult::to_modern_json(
    const Tree& tree,
    const Alignment& aln,
    const std::string& cli_cmd,
    const std::string& input_filepath
) const {
    nlohmann::json j;
    j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/prime.v3.schema.json";
    j["schema_version"] = "3.0.0";

    // PRIME's schema describes this exact analysis-block shape
    // (type/authors/citation/contact rather than the shared category/
    // citations contract), so it is kept verbatim.
    j["analysis"] = {
        {"id", "prime"},
        {"name", "PRoperty Informed Models of Evolution (PRIME)"},
        {"type", "selection_property"},
        {"version", "3.0.0"},
        {"authors", "Sergei L. Kosakovsky Pond"},
        {"citation", "Pond, SLK et al. PRoperty Informed Models of Evolution (PRIME)"},
        {"contact", "spond@temple.edu"}
    };

    j["software"] = {
        {"name", "hyphy"},
        {"version", "3.0.0"},
        {"homepage", "https://hyphy.org"}
    };

    Provenance auto_prov;
    auto_prov.invocation.cli_command = cli_cmd;
    auto_prov.invocation.working_directory = Provenance::get_cwd();
    if (!input_filepath.empty()) {
        size_t sz = 0;
        std::string hash = crypto::SHA256::hash_file(input_filepath, &sz);
        if (!hash.empty()) {
            auto_prov.inputs["alignment"] = {input_filepath, "FASTA/NEXUS", hash, sz};
        }
    }
    auto_prov.execution.start_time = Provenance::current_iso8601();
    auto_prov.execution.end_time = auto_prov.execution.start_time;
    auto_prov.execution.wall_time_seconds = runtime_seconds;
    auto_prov.execution.cpu_threads = 1;
    auto_prov.execution.hostname = Provenance::get_hostname();
    auto_prov.execution.os = Provenance::detect_os();
    auto_prov.execution.compiler = Provenance::detect_compiler();
    j["provenance"] = auto_prov.to_json();

    // Dataset — PRIME's schema pins genetic_code.id as a string of the
    // numeric code id, so this stays a custom block.
    const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
    j["dataset"] = {
        {"taxa_count", aln.num_taxa},
        {"codon_sites", aln.num_codons},
        {"nucleotide_sites", aln.num_nucleotides},
        {"unique_patterns", aln.patterns.size()},
        {"taxa", aln.taxon_names},
        {"genetic_code", {
            {"id", std::to_string(gcode.id)},
            {"name", gcode.name},
            {"sense_codons", gcode.num_sense_codons},
            {"stop_codons", gcode.stop_codons}
        }}
    };

    // Properties metadata
    nlohmann::json prop_meta = nlohmann::json::array();
    for (size_t k = 0; k < properties.size(); ++k) {
        nlohmann::json p_entry;
        p_entry["index"] = k + 1;
        p_entry["name"] = properties[k].name;
        p_entry["description"] = properties[k].description;
        p_entry["mean"] = properties[k].mean;
        p_entry["std_dev"] = properties[k].std_dev;
        prop_meta.push_back(p_entry);
    }
    j["properties"] = prop_meta;

    // Phylogeny — shared block; the "MLE" key preserves upstream's
    // per-node branch-length label.
    j["phylogeny"] = phylogeny_block(
        tree, "",
        [](const TreeNode& node) {
            return std::map<std::string, Scalar>{{"MLE", node.branch_length}};
        });
    j["phylogeny"]["tree_lengths"]["baseline"] =
        finite_or_null(global_fit.baseline_tree_length);

    // Model fits
    nlohmann::json gpr = fit_entry(
        "Global property-informed codon model fit across all sites",
        global_fit.log_likelihood, global_fit.aicc, 0);
    gpr.erase("parameters_count"); // PRIME's schema does not carry it
    gpr["alpha"] = finite_or_null(global_fit.alpha);
    gpr["beta"] = finite_or_null(global_fit.beta);
    gpr["lambdas"] = global_fit.lambdas;
    j["model_fits"] = {{"G-PRIME", gpr}};

    // Statistical tests summary
    nlohmann::json summary_j = {
        {"method", "Site-Level PRoperty Informed Models of Evolution (S-PRIME)"},
        {"test_statistic", "Likelihood Ratio Test (LRT)"},
        {"threshold", settings.p_value_threshold},
        {"stage2_threshold", settings.stage2_p_value},
        {"shrinkage_lambda", settings.shrinkage_lambda},
        {"total_sites", static_cast<int>(total_sites)},
        {"invariant_sites", static_cast<int>(invariant_sites)},
        {"significant_sites", static_cast<int>(significant_sites)},
        {"properties_count", static_cast<int>(properties.size())}
    };
    for (const auto& p : properties) {
        summary_j["conserved_" + p.name] = conserved_counts.at(p.name);
        summary_j["changing_" + p.name] = changing_counts.at(p.name);
    }
    j["statistical_tests"] = {{"site_level_summary", summary_j}};

    // Site results (columnar; ColumnarResults enforces equal-length)
    ColumnarResults table;
    std::vector<int> col_site, col_ns, col_ua;
    std::vector<double> col_a, col_b, col_fa, col_fb, col_llp, col_llf,
                        col_olrt, col_op, col_oq, col_tbl, col_rr;
    std::map<std::string, std::vector<double>> col_lams, col_pv, col_adj_pv;
    std::map<std::string, std::vector<std::string>> col_cls;

    for (const auto& p : properties) {
        col_lams[p.name] = {};
        col_pv[p.name] = {};
        col_adj_pv[p.name] = {};
        col_cls[p.name] = {};
    }

    for (const auto& sr : site_results) {
        col_site.push_back(static_cast<int>(sr.site_index + 1));
        col_a.push_back(sr.alpha);
        col_b.push_back(sr.beta);
        col_fa.push_back(sr.fel_alpha);
        col_fb.push_back(sr.fel_beta);
        col_llp.push_back(sr.log_l_alt);
        col_llf.push_back(sr.log_l_fel);
        col_olrt.push_back(sr.omnibus_lrt);
        col_op.push_back(sr.omnibus_p_value);
        col_oq.push_back(sr.omnibus_q_value);
        col_tbl.push_back(sr.total_branch_length);
        col_ns.push_back(static_cast<int>(sr.nonsyn_subs));
        col_ua.push_back(static_cast<int>(sr.unique_aa));
        col_rr.push_back(sr.redundancy_ratio);

        for (const auto& p_inf : sr.properties) {
            col_lams[p_inf.property_name].push_back(p_inf.lambda);
            col_pv[p_inf.property_name].push_back(p_inf.p_value);
            col_adj_pv[p_inf.property_name].push_back(p_inf.adjusted_p_value);
            col_cls[p_inf.property_name].push_back(p_inf.classification);
        }
    }

    table.add_column("site", {"integer", "1-based codon site index"}, col_site);
    table.add_column("alpha", {"number", "Synonymous rate at site under PRIME"}, col_a);
    table.add_column("beta", {"number", "Baseline non-synonymous rate under PRIME"}, col_b);
    table.add_column("fel_alpha", {"number", "Synonymous rate under neutral FEL"}, col_fa);
    table.add_column("fel_beta", {"number", "Non-synonymous rate under neutral FEL"}, col_fb);
    table.add_column("log_l_prime", {"number", "Site log-likelihood under PRIME"}, col_llp);
    table.add_column("log_l_fel", {"number", "Site log-likelihood under FEL"}, col_llf);
    table.add_column("omnibus_lrt", {"number", "LRT statistic for any property effect"}, col_olrt);
    table.add_column("omnibus_p_value", {"number", "Omnibus p-value"}, col_op);
    table.add_column("omnibus_q_value", {"number", "Benjamini-Hochberg FDR q-value"}, col_oq);
    table.add_column("total_branch_length", {"number", "Total branch length contributing to site"}, col_tbl);
    table.add_column("nonsynonymous_substitutions", {"integer", "Inferred non-synonymous substitution count"}, col_ns);
    table.add_column("unique_amino_acids", {"integer", "Number of unique observed amino acids"}, col_ua);
    table.add_column("redundancy_ratio", {"number", "Substitutions to unique amino acids ratio"}, col_rr);

    for (const auto& p : properties) {
        table.add_column("lambda_" + p.name, {"number", "Constraint parameter for " + p.name}, col_lams[p.name]);
        table.add_column("p_" + p.name, {"number", "p-value for " + p.name}, col_pv[p.name]);
        table.add_column("adj_p_" + p.name, {"number", "Holm-Bonferroni adjusted p-value for " + p.name}, col_adj_pv[p.name]);
        table.add_column("class_" + p.name, {"string", "Classification for " + p.name + " (Conserved, Changing, Neutral)"}, col_cls[p.name]);
    }
    j["site_results"] = table.to_json();

    return j;
}

nlohmann::json PRIMEResult::to_legacy_json(
    const Tree& tree,
    const Alignment& aln,
    const std::string& input_filepath
) const {
    nlohmann::json j;

    j["analysis"] = {
        {"authors", "Sergei L. Kosakovsky Pond"},
        {"citation", "Pond, SLK et al. PRoperty Informed Models of Evolution (PRIME)"},
        {"contact", "spond@temple.edu"},
        {"info", "PRIME (PRoperty Informed Model of Evolution)."},
        {"requirements", "in-frame codon alignment and a phylogenetic tree"},
        {"residue_properties", nullptr},
        {"version", "0.2.6"}
    };

    j["input"] = {
        {"file name", input_filepath},
        {"number of sequences", aln.num_taxa},
        {"number of sites", aln.num_codons},
        {"partition count", 1},
        {"trees", {{"0", tree.to_newick()}}}
    };

    // Legacy MLE Table Headers
    // Index 0..11: alpha, beta, FEL alpha, FEL beta, Total branch length, # subs, # aa, PRIME LogL, FEL LogL, p-value, q-value, R
    // Followed by 3 columns per property: &lambda;k, pk, LogLk
    nlohmann::json headers = nlohmann::json::array();
    headers.push_back({"&alpha;", "Synonymous substitution rate at a site"});
    headers.push_back({"&beta;", "Property independent non-synonymous rate a site"});
    headers.push_back({"FEL &alpha;", "Synonymous substitution rate at a site (FEL model)"});
    headers.push_back({"FEL &beta;", "Log of property independent non-synonymous rate a site"});
    headers.push_back({"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"});
    headers.push_back({"# subs", "The number of amino-acid substitutions that occurred at this site"});
    headers.push_back({"# aa", "The number of unique amino-acids occurring at that site"});
    headers.push_back({"PRIME LogL", "Site Log-likelihood under the PRIME model"});
    headers.push_back({"FEL LogL", "Site Log-likelihood under the FEL model"});
    headers.push_back({"p-value", "Omnibus p-value (any property is important)"});
    headers.push_back({"q-value", "Omnibus q-value (any property is important)"});
    headers.push_back({"R", "Redundancy ratio (subs / unique aa)"});

    for (size_t k = 0; k < properties.size(); ++k) {
        std::string k_str = std::to_string(k + 1);
        headers.push_back({"&lambda;" + k_str, "Importance for " + properties[k].name});
        headers.push_back({"p" + k_str, "p-value for non-zero effect of " + properties[k].name});
        headers.push_back({"LogL" + k_str, "Log likelihood when there is no effect of " + properties[k].name});
    }

    nlohmann::json content_rows = nlohmann::json::array();
    for (const auto& sr : site_results) {
        nlohmann::json row = nlohmann::json::array();
        row.push_back(sr.alpha);
        row.push_back(sr.beta);
        row.push_back(sr.fel_alpha);
        row.push_back(sr.fel_beta);
        row.push_back(sr.total_branch_length);
        row.push_back(static_cast<int>(sr.nonsyn_subs));
        row.push_back(static_cast<int>(sr.unique_aa));
        row.push_back(sr.log_l_alt);
        row.push_back(sr.log_l_fel);
        row.push_back(sr.omnibus_p_value);
        row.push_back(sr.omnibus_q_value);
        row.push_back(sr.redundancy_ratio);

        for (const auto& p_inf : sr.properties) {
            row.push_back(p_inf.lambda);
            row.push_back(p_inf.adjusted_p_value); // Stage 2 adjusted p-value
            row.push_back(p_inf.log_l_null);
        }
        content_rows.push_back(row);
    }

    j["MLE"] = {
        {"headers", headers},
        {"content", {{"0", content_rows}}}
    };

    // Legacy "test results"
    nlohmann::json sig_sites_map = nlohmann::json::object();
    size_t sig_idx = 0;
    for (const auto& sr : site_results) {
        if (sr.is_significant(settings.p_value_threshold)) {
            nlohmann::json site_entry;
            site_entry["codon"] = static_cast<int>(sr.site_index);
            site_entry["partiton"] = 0;
            site_entry["q-value"] = sr.omnibus_q_value;

            nlohmann::json prop_obj = nlohmann::json::object();
            for (const auto& p_inf : sr.properties) {
                if (p_inf.classification != "Neutral") {
                    prop_obj[p_inf.property_name] = p_inf.adjusted_p_value;
                }
            }
            site_entry["properties"] = prop_obj;
            sig_sites_map[std::to_string(sig_idx++)] = site_entry;
        }
    }

    j["test results"] = {
        {"terms.qvalue", settings.p_value_threshold},
        {"Global test significant sites", sig_sites_map}
    };

    j["fits"] = {
        {"Global MG94xREV", {{"Log Likelihood", global_fit.log_likelihood}}},
        {"Nucleotide GTR", {{"Log Likelihood", global_fit.gtr_log_l}}}
    };

    j["runtime"] = "3.0.0";
    j["timers"] = {{"Total time", runtime_seconds}};

    return j;
}

} // namespace hyphy::analyses
