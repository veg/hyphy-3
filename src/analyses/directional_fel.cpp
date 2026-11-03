// Out-of-line JSON emitter definitions for directional_fel.
// Moved from include/hyphy/analyses/directional_fel.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/directional_fel.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json DirectionalFELAnalyzer::to_modern_json(const std::string& input_filepath, const std::string& tree_string) const {
        AnalysisMetadata meta{
            "dfel",
            "Directional Fixed Effects Likelihood",
            "3.0.0",
            "selection",
            "Directional FEL detects site-level selection while accounting for lineage-specific compositional heterogeneity and auto-rooting via non-reversible dynamics",
            {},
            {
                {"p_value_threshold", p_value_threshold},
                {"heterogeneity_threshold", settings.heterogeneity_threshold},
                {"root_strategy", settings.root_strategy},
                {"num_regimes", regimes.size()}
            }
        };

        nlohmann::json j;
        j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/dfel.v3.schema.json";
        j["schema_version"] = meta.version;
        j["analysis"] = meta.to_json();
        j["software"] = software_block();
        j["software"]["citation"] = "Directional Evolution & Fixed Effects Likelihood under Non-Stationary Composition Dynamics (2026)";

        // Provenance — dFEL's signature carries no cli_cmd, so build the
        // block inline with a fixed invocation label.
        {
            Provenance prov;
            prov.invocation.cli_command = "hyphy3 dfel";
            prov.invocation.working_directory = Provenance::get_cwd();
            if (!input_filepath.empty()) {
                size_t sz = 0;
                std::string hash = crypto::SHA256::hash_file(input_filepath, &sz);
                if (!hash.empty()) {
                    prov.inputs["alignment"] = {input_filepath, "FASTA/NEXUS", hash, sz};
                }
            }
            prov.execution.start_time = Provenance::current_iso8601();
            prov.execution.end_time = prov.execution.start_time;
            prov.execution.hostname = Provenance::get_hostname();
            prov.execution.os = Provenance::detect_os();
            prov.execution.compiler = Provenance::detect_compiler();
            j["provenance"] = prov.to_json();
        }
        j["provenance"]["created_at"] = Provenance::current_iso8601();

        // Dataset — dfel.v3.schema.json pins genetic_code as a plain string
        // and patterns_count (not the shared unique_patterns key), so this
        // block stays custom.
        j["dataset"]["file"] = input_filepath;
        j["dataset"]["taxa_count"] = aln.num_taxa;
        j["dataset"]["codon_sites"] = aln.num_codons;
        j["dataset"]["patterns_count"] = aln.patterns.size();
        j["dataset"]["genetic_code"] = aln.code ? aln.code->name : "Universal";

        // Phylogeny — schema-pinned custom shape (rooted auto-root details).
        j["phylogeny"]["rooted"] = true;
        j["phylogeny"]["root_branch"] = root_result.chosen_branch_name;
        j["phylogeny"]["root_log_likelihood"] = finite_or_null(root_result.chosen_log_l);
        j["phylogeny"]["root_weight"] = finite_or_null(root_result.chosen_weight);
        j["phylogeny"]["tree_length"] = finite_or_null(baseline_tree_length);
        if (!tree_string.empty()) {
            j["phylogeny"]["newick"] = tree_string;
        }

        // Compositional heterogeneity test
        nlohmann::json comp_j;
        comp_j["is_heterogeneous"] = comp_test.is_heterogeneous;
        comp_j["overall"] = { {"chi2", comp_test.chi2_overall}, {"df", comp_test.df_overall}, {"p_value", comp_test.p_value_overall} };
        comp_j["codon_position_1"] = { {"chi2", comp_test.chi2_cp1}, {"df", comp_test.df_cp1}, {"p_value", comp_test.p_value_cp1} };
        comp_j["codon_position_2"] = { {"chi2", comp_test.chi2_cp2}, {"df", comp_test.df_cp2}, {"p_value", comp_test.p_value_cp2} };
        comp_j["codon_position_3"] = { {"chi2", comp_test.chi2_cp3}, {"df", comp_test.df_cp3}, {"p_value", comp_test.p_value_cp3} };

        nlohmann::json taxa_arr = nlohmann::json::array();
        for (const auto& tc : comp_test.taxa) {
            taxa_arr.push_back({
                {"taxon", tc.name},
                {"gc_overall", tc.gc_overall},
                {"gc_cp3", tc.gc_cp3}
            });
        }
        comp_j["taxa"] = taxa_arr;
        j["compositional_heterogeneity"] = comp_j;

        // Auto-rooting candidates
        nlohmann::json root_cands = nlohmann::json::array();
        for (const auto& c : root_result.candidates) {
            root_cands.push_back({
                {"rank", c.rank},
                {"branch", c.branch_name},
                {"is_terminal", c.is_terminal},
                {"log_likelihood", c.log_l},
                {"delta_log_l", c.delta_log_l},
                {"weight", c.weight}
            });
        }
        j["root_selection"]["candidates"] = root_cands;

        // Regimes
        nlohmann::json reg_arr = nlohmann::json::array();
        for (const auto& reg : regimes) {
            reg_arr.push_back({
                {"id", reg.id},
                {"name", reg.name},
                {"psi", reg.psi},
                {"mean_gc", reg.mean_gc},
                {"branch_count", reg.branch_count}
            });
        }
        j["directional_model"]["regimes"] = reg_arr;

        nlohmann::json branch_regimes_j = nlohmann::json::object();
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                const auto& w = branch_regime_weights[node.id];
                branch_regimes_j[node.name] = {
                    {"regime_id", branch_regime_ids[node.id]},
                    {"weights", nlohmann::json::array({w[0], w[1], w[2]})}
                };
            }
        }
        j["directional_model"]["branch_regimes"] = branch_regimes_j;

        // Model fits — schema pins these two fit names
        j["model_fits"]["gtr"]["log_likelihood"] = finite_or_null(gtr_log_l);
        j["model_fits"]["gtr"]["aicc"] = finite_or_null(gtr_aicc);
        j["model_fits"]["directional_mg94"]["log_likelihood"] = finite_or_null(global_log_l);
        j["model_fits"]["directional_mg94"]["aicc"] = finite_or_null(global_aicc);

        // Site results (columnar; ColumnarResults enforces equal-length)
        ColumnarResults table;
        std::vector<int> sites_vec;
        std::vector<double> alpha_vec, beta_vec, alpha_null_vec, lrt_vec, pval_vec;
        std::vector<std::string> class_vec, subs_vec;

        for (const auto& r : site_results) {
            sites_vec.push_back(static_cast<int>(r.site_index + 1));
            alpha_vec.push_back(r.alpha);
            beta_vec.push_back(r.beta);
            alpha_null_vec.push_back(r.alpha_null);
            lrt_vec.push_back(r.lrt);
            pval_vec.push_back(r.p_value);

            std::string c_type = "neutral";
            if (r.p_value <= p_value_threshold) {
                if (r.beta > r.alpha) c_type = "diversifying";
                else if (r.alpha > r.beta) c_type = "purifying";
            }
            class_vec.push_back(c_type);
            subs_vec.push_back(r.subs.summary_string);
        }

        table.add_column("site", {"integer", "1-based codon site index"}, sites_vec);
        table.add_column("alpha", {"number", "Synonymous rate under the directional model"}, alpha_vec);
        table.add_column("beta", {"number", "Non-synonymous rate under the directional model"}, beta_vec);
        table.add_column("alpha_null", {"number", "Synonymous rate under the null model"}, alpha_null_vec);
        table.add_column("lrt", {"number", "Likelihood ratio test statistic"}, lrt_vec);
        table.add_column("p_value", {"number", "Site-level p-value"}, pval_vec);
        table.add_column("classification", {"string", "diversifying, purifying, or neutral"}, class_vec);
        table.add_column("substitutions", {"string", "Summary of inferred substitutions at the site"}, subs_vec);
        j["site_results"] = table.to_json();

        return j;
    }

nlohmann::json DirectionalFELAnalyzer::to_legacy_json(const std::string& input_filepath, const std::string& tree_string) const {
        nlohmann::json j;
        j["analysis"]["authors"] = "Sergei L Kosakovsky Pond and HyPhy Team";
        j["analysis"]["citation"] = "Directional Fixed Effects Likelihood under Non-Stationary Composition Dynamics (2026)";
        j["analysis"]["contact"] = "spond@temple.edu";
        j["analysis"]["info"] = "Directional FEL detects site-level selection while accounting for lineage-specific compositional heterogeneity and auto-rooting via non-reversible dynamics";
        j["analysis"]["requirements"] = "in-frame codon alignment and a phylogenetic tree";
        j["analysis"]["version"] = "3.0.0";

        j["MLE"]["headers"] = nlohmann::json::array({
            nlohmann::json::array({"alpha", "Synonymous substitution rate at a site"}),
            nlohmann::json::array({"beta", "Non-synonymous substitution rate at a site"}),
            nlohmann::json::array({"alpha=beta", "The rate estimate under the neutral model"}),
            nlohmann::json::array({"LRT", "Likelihood ratio test statistic"}),
            nlohmann::json::array({"p-value", "Asymptotic p-value for evidence of selection"}),
            nlohmann::json::array({"Total branch length", "The total length of branches"})
        });

        nlohmann::json content = nlohmann::json::array();
        for (const auto& r : site_results) {
            content.push_back({
                r.alpha,
                r.beta,
                r.alpha_null,
                r.lrt,
                r.p_value,
                r.total_branch_length
            });
        }
        j["MLE"]["content"]["0"] = content;

        j["input"]["file name"] = input_filepath;
        j["input"]["genetic code"] = aln.code ? aln.code->name : "Universal";
        j["input"]["number of sequences"] = aln.num_taxa;
        j["input"]["number of sites"] = aln.num_codons;
        if (!tree_string.empty()) {
            j["input"]["trees"]["0"] = tree_string;
        }

        return j;
    }

} // namespace hyphy::analyses
