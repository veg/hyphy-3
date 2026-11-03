// Out-of-line JSON emitter definitions for absrel.
// Moved from include/hyphy/analyses/absrel.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/absrel.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json ABSRELResult::to_legacy_json(const Tree& tree, const Alignment& aln) const {
        (void)tree;
        nlohmann::json j;

        // Analysis metadata
        j["analysis"] = {
            {"info", "aBSREL (Adaptive branch-site random effects likelihood) uses an adaptive random effects branch-site model framework to test whether each branch has evolved under positive selection, inferring an optimal number of rate categories per branch."},
            {"reference", "Less Is More: An Adaptive Branch-Site Random Effects Model for Efficient Detection of Episodic Diversifying Selection (2015). Mol Biol Evol 32 (5): 1342-1353."},
            {"version", "3.0.0"}
        };

        // Input
        j["input"] = {
            {"number of sequences", aln.num_taxa},
            {"number of sites", aln.num_codons},
            {"partition count", 1},
            {"trees", {{"0", aln.embedded_tree_newick}}}
        };

        // Fits
        j["fits"] = {
            {"Nucleotide GTR", {
                {"Log Likelihood", gtr_fit.log_likelihood},
                {"estimated parameters", gtr_fit.parameters},
                {"AIC-c", gtr_fit.aicc},
                {"display order", 0}
            }},
            {"Baseline MG94xREV", {
                {"Log Likelihood", baseline_fit.log_likelihood},
                {"estimated parameters", baseline_fit.parameters},
                {"AIC-c", baseline_fit.aicc},
                {"display order", 1}
            }},
            {"Full adaptive model", {
                {"Log Likelihood", full_adaptive_fit.log_likelihood},
                {"estimated parameters", full_adaptive_fit.parameters},
                {"AIC-c", full_adaptive_fit.aicc},
                {"display order", 2}
            }}
        };

        // Branch Attributes
        nlohmann::json ba = nlohmann::json::object();
        for (const auto& [name, bres] : branches) {
            nlohmann::json b_dict;
            b_dict["Baseline MG94xREV"] = bres.baseline_branch_length;
            b_dict["Baseline MG94xREV omega ratio"] = bres.baseline_omega;
            b_dict["Rate classes"] = bres.rate_classes;
            b_dict["Full adaptive model"] = bres.full_branch_length;
            b_dict["Full adaptive model (synonymous subs/site)"] = bres.full_es;
            b_dict["Full adaptive model (non-synonymous subs/site)"] = bres.full_en;

            // Rate distributions: [[omega, weight], ...]
            nlohmann::json rdist = nlohmann::json::array();
            for (size_t k = 0; k < bres.rate_distribution.rates.size(); ++k) {
                rdist.push_back({bres.rate_distribution.rates[k], bres.rate_distribution.weights[k]});
            }
            b_dict["Rate Distributions"] = rdist;

            if (bres.is_tested) {
                b_dict["LRT"] = bres.lrt;
                b_dict["Uncorrected P-value"] = bres.uncorrected_p_value;
                b_dict["Corrected P-value"] = bres.corrected_p_value;
                b_dict["Sites @ EBF>=100"] = bres.sites_ebf_100;
            }
            ba[name] = b_dict;
        }
        j["branch attributes"] = {{"0", ba}};

        // Test Results
        j["test results"] = {
            {"P-value threshold", p_threshold},
            {"tested", tested_branches.size()},
            {"positive test results", positive_branches.size()}
        };

        // Tested map
        nlohmann::json tested_map = nlohmann::json::object();
        for (const auto& [name, bres] : branches) {
            tested_map[name] = bres.is_tested ? "test" : "background";
        }
        j["tested"] = {{"0", tested_map}};

        // Timers
        j["timers"] = {
            {"Overall", {{"timer", static_cast<int>(runtime_seconds)}, {"order", 0}}}
        };

        return j;
    }

nlohmann::json ABSRELResult::to_modern_json(const Tree& tree, const Alignment& aln, const Provenance& prov) const {
        AnalysisMetadata meta{
            "absrel",
            "Adaptive Branch-Site Random Effects Likelihood",
            "3.0.0",
            "branch_site_selection",
            "aBSREL (Adaptive branch-site random effects likelihood) uses an adaptive random effects branch-site model framework to test whether each branch has evolved under positive selection, inferring an optimal number of rate categories per branch.",
            {citation("Smith MD, Wertheim JO, Weaver S, Murrell B, Scheffler K, Kosakovsky Pond SL (2015). Less Is More: An Adaptive Branch-Site Random Effects Model for Efficient Detection of Episodic Diversifying Selection. Mol Biol Evol 32(5): 1342-1353.",
                      "10.1093/molbev/msv022", "25697341")},
            {
                {"p_threshold", p_threshold},
                {"tested_branches", tested_branches.size()}
            }
        };

        nlohmann::json j = new_analysis_json(meta, prov);
        // Preserve measured runtime when the caller supplied no real provenance.
        if (prov.invocation.cli_command.empty() && prov.inputs.empty() &&
            j["provenance"].contains("execution")) {
            j["provenance"]["execution"]["wall_time_seconds"] = runtime_seconds;
        }
        j["dataset"] = dataset_block(aln);
        j["phylogeny"] = phylogeny_block(
            tree, "",
            [this](const TreeNode& node) {
                std::map<std::string, Scalar> lengths;
                auto it = branches.find(node.name);
                if (it != branches.end()) {
                    lengths["baseline_mg94"] = it->second.baseline_branch_length;
                    lengths["full_adaptive"] = it->second.full_branch_length;
                }
                return lengths;
            },
            [this](const TreeNode& node) -> nlohmann::json {
                auto it = branches.find(node.name);
                if (it == branches.end()) return nlohmann::json::object();
                const auto& bres = it->second;
                nlohmann::json entry;
                entry["rate_classes"] = bres.rate_classes;
                entry["baseline_omega"] = finite_or_null(bres.baseline_omega);
                entry["is_tested"] = bres.is_tested;
                entry["is_positive"] = bres.is_positive;
                const auto& rd = bres.rate_distribution;
                if (rd.rates.size() != rd.weights.size()) {
                    throw std::runtime_error(
                        "absrel: rate_distribution length mismatch on branch '" +
                        node.name + "' (rates=" + std::to_string(rd.rates.size()) +
                        ", weights=" + std::to_string(rd.weights.size()) + ")");
                }
                nlohmann::json rdist = nlohmann::json::array();
                for (size_t k = 0; k < rd.rates.size(); ++k) {
                    rdist.push_back({
                        {"omega", finite_or_null(rd.rates[k])},
                        {"proportion", finite_or_null(rd.weights[k])}
                    });
                }
                entry["rate_distribution"] = rdist;
                if (bres.is_tested) {
                    entry["lrt"] = finite_or_null(bres.lrt);
                    entry["uncorrected_p_value"] = finite_or_null(bres.uncorrected_p_value);
                    entry["corrected_p_value"] = finite_or_null(bres.corrected_p_value);
                    entry["sites_ebf_100"] = bres.sites_ebf_100;
                }
                return entry;
            });

        // Model fits
        if (gtr_fit.parameters > 0 || gtr_fit.log_likelihood != 0.0) {
            j["model_fits"]["nucleotide_gtr"] = nucleotide_gtr_fit(
                aln, gtr_params, gtr_fit.log_likelihood, gtr_fit.aicc,
                static_cast<int>(gtr_fit.parameters));
        }
        j["model_fits"]["baseline_mg94"] = fit_entry(
            "Baseline MG94xREV model with individual branch omegas and branch lengths",
            baseline_fit.log_likelihood, baseline_fit.aicc,
            static_cast<int>(baseline_fit.parameters));
        j["model_fits"]["full_adaptive"] = fit_entry(
            "Full adaptive aBSREL mixture model with branch-specific rate categories",
            full_adaptive_fit.log_likelihood, full_adaptive_fit.aicc,
            static_cast<int>(full_adaptive_fit.parameters));

        // Statistical tests
        j["statistical_tests"]["branch_level_summary"] = {
            {"method", "Adaptive Branch-Site Random Effects Likelihood (aBSREL)"},
            {"test_statistic", "Likelihood Ratio Test (LRT)"},
            {"distribution", "Asymptotic mixture distribution"},
            {"multiple_testing_correction", "Holm-Bonferroni step-down procedure"},
            {"threshold", p_threshold},
            {"tested_branches_count", static_cast<int>(tested_branches.size())},
            {"positive_branches_count", static_cast<int>(positive_branches.size())},
            {"positive_branches", positive_branches}
        };

        // Branch results (columnar; ColumnarResults enforces equal-length)
        ColumnarResults table;
        std::vector<std::string> col_branch;
        std::vector<int> col_k, col_ebf;
        std::vector<double> col_base_omega, col_len, col_lrt, col_raw_p, col_adj_p;
        std::vector<bool> col_pos;

        for (const auto& bname : tested_branches) {
            auto it = branches.find(bname);
            if (it != branches.end()) {
                const auto& bres = it->second;
                col_branch.push_back(bname);
                col_k.push_back(bres.rate_classes);
                col_base_omega.push_back(bres.baseline_omega);
                col_len.push_back(bres.full_branch_length);
                col_lrt.push_back(bres.lrt);
                col_raw_p.push_back(bres.uncorrected_p_value);
                col_adj_p.push_back(bres.corrected_p_value);
                col_ebf.push_back(bres.sites_ebf_100);
                col_pos.push_back(bres.is_positive);
            }
        }

        table.add_column("branch", {"string", "Branch name in phylogeny"}, col_branch);
        table.add_column("rate_classes", {"integer", "Number of inferred omega rate classes"}, col_k);
        table.add_column("baseline_omega", {"number", "Baseline single omega estimate"}, col_base_omega);
        table.add_column("full_branch_length", {"number", "Estimated branch length under full adaptive model", "substitutions/site"}, col_len);
        table.add_column("lrt", {"number", "Likelihood ratio test statistic for episodic positive selection"}, col_lrt);
        table.add_column("uncorrected_p_value", {"number", "Asymptotic uncorrected p-value"}, col_raw_p);
        table.add_column("corrected_p_value", {"number", "Holm-Bonferroni corrected p-value"}, col_adj_p);
        table.add_column("sites_ebf_100", {"integer", "Number of sites with Empirical Bayes Factor (EBF) >= 100 favoring selection"}, col_ebf);
        table.add_column("positive", {"boolean", "Flag indicating episodic diversifying selection at threshold"}, col_pos);
        j["branch_results"] = table.to_json();

        return j;
    }

} // namespace hyphy::analyses
