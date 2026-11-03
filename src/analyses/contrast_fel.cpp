// Out-of-line JSON emitter definitions for contrast_fel.
#include "hyphy/analyses/contrast_fel.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json ContrastFELAnalyzer::to_modern_json(
        const std::string& input_filepath,
        const std::string& tree_string,
        const Provenance& prov) const {
        nlohmann::json j;

        j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/contrast_fel.v3.schema.json";
        j["schema_version"] = "3.0.0";
        j["analysis"] = {
            {"id", "contrast-fel"},
            {"name", "Contrast-FEL (Fixed Effects Likelihood)"},
            {"version", "3.0.0"},
            {"description", "Test for differential selective pressures between branch sets using Empirical Bayes Ridge Shrinkage"},
            {"reference", "Mol Biol Evol (2021) 38(3):1184–1198"},
            {"shrinkage_lambda", shrinkage_lambda}
        };
        j["software"] = software_block();

        j["provenance"] = prov.to_json();

        // Dataset
        j["dataset"]["file"] = input_filepath;
        j["dataset"]["sequences"] = aln.num_taxa;
        j["dataset"]["codon_sites"] = aln.num_codons;
        j["dataset"]["nucleotide_sites"] = aln.num_nucleotides;
        j["dataset"]["taxa"] = aln.taxon_names;

        // Phylogeny — contrast-FEL's node entries carry group membership and a
        // scalar (rescaled) branch length rather than the shared
        // branch_lengths map, so this block stays custom.
        j["phylogeny"]["newick"] = tree_string.empty() ? tree.to_newick() : tree_string;
        j["phylogeny"]["tree_length"] = finite_or_null(baseline_tree_length);
        nlohmann::json nodes_json = nlohmann::json::object();
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                int gid = node_group_id[node.id];
                std::string grp_name = (gid >= 0 && gid < static_cast<int>(testable_groups.size()))
                    ? testable_groups[gid].name
                    : (has_background ? "Background" : "Reference");
                nodes_json[node.name] = {
                    {"type", node.is_leaf ? "leaf" : "internal"},
                    {"branch_length", (conversion_factor > 0.0) ? (node.branch_length / conversion_factor) : node.branch_length},
                    {"group", grp_name},
                    {"is_testable", (gid >= 0)}
                };
            }
        }
        j["phylogeny"]["nodes"] = nodes_json;

        // Model fits
        if (has_gtr_fit) {
            j["model_fits"]["nucleotide_gtr"]["log_likelihood"] = finite_or_null(gtr_log_l);
            j["model_fits"]["nucleotide_gtr"]["aicc"] = finite_or_null(gtr_aicc);
            j["model_fits"]["nucleotide_gtr"]["rates"] = {
                {"AC", finite_or_null(gtr_rates.theta_AC)},
                {"AT", finite_or_null(gtr_rates.theta_AT)},
                {"CG", finite_or_null(gtr_rates.theta_CG)},
                {"CT", finite_or_null(gtr_rates.theta_CT)},
                {"GT", finite_or_null(gtr_rates.theta_GT)}
            };
        }

        j["model_fits"]["codon_mg94"]["log_likelihood"] = finite_or_null(global_log_l);
        j["model_fits"]["codon_mg94"]["aicc"] = finite_or_null(global_aicc);
        j["model_fits"]["codon_mg94"]["omega"] = finite_or_null(
            (base_params.alpha > 0.0) ? (base_params.beta / base_params.alpha) : 1.0);

        // Summary
        auto sum = get_summary();
        j["statistical_tests"]["summary"] = {
            {"total_sites", static_cast<int>(sum.total_sites)},
            {"degrees_of_freedom", static_cast<int>(testable_groups.size() > 1 ? testable_groups.size() - 1 : 1)},
            {"significant_p_threshold", p_value_threshold},
            {"significant_p_sites", static_cast<int>(sum.significant_p_sites)},
            {"significant_q_threshold", q_value_threshold},
            {"significant_q_sites", static_cast<int>(sum.significant_q_sites)},
            {"shrinkage_lambda", shrinkage_lambda}
        };

        // Pairwise summary if K >= 2
        if (!site_results.empty() && !site_results[0].pairwise.empty()) {
            nlohmann::json pw_summary = nlohmann::json::array();
            size_t num_pairs = site_results[0].pairwise.size();
            for (size_t pw_idx = 0; pw_idx < num_pairs; ++pw_idx) {
                const auto& pw0 = site_results[0].pairwise[pw_idx];
                int sig_p_cnt = 0;
                int sig_adj_cnt = 0;
                for (const auto& r : site_results) {
                    if (pw_idx < r.pairwise.size()) {
                        if (r.pairwise[pw_idx].p_value <= p_value_threshold) sig_p_cnt++;
                        if (r.pairwise[pw_idx].p_value_adjusted <= p_value_threshold) sig_adj_cnt++;
                    }
                }
                pw_summary.push_back({
                    {"group_a", pw0.name_a},
                    {"group_b", pw0.name_b},
                    {"significant_p_sites", sig_p_cnt},
                    {"significant_adjusted_p_sites", sig_adj_cnt}
                });
            }
            j["statistical_tests"]["pairwise_contrasts"] = pw_summary;
        }

        // Site results (columnar; ColumnarResults enforces equal-length)
        ColumnarResults table;
        const size_t num_pairs = (!site_results.empty()) ? site_results[0].pairwise.size() : 0;
        std::vector<int> col_site;
        std::vector<double> col_a, col_a_null, col_b_null, col_raw_lrt, col_plrt, col_pval, col_qval;
        std::vector<std::string> col_tag, col_subs;
        std::vector<std::vector<double>> col_betas(testable_groups.size());
        std::vector<double> col_beta_bg;
        std::vector<std::vector<double>> col_pw_lrt(num_pairs), col_pw_pval(num_pairs), col_pw_padj(num_pairs);

        for (size_t c = 0; c < site_results.size(); ++c) {
            const auto& r = site_results[c];
            col_site.push_back(static_cast<int>(c + 1));
            col_a.push_back(r.alpha);
            col_a_null.push_back(r.alpha_null);
            col_b_null.push_back(r.beta_null);
            col_raw_lrt.push_back(r.raw_lrt);
            col_plrt.push_back(r.penalized_lrt);
            col_pval.push_back(r.p_value);
            col_qval.push_back(r.q_value);
            col_tag.push_back(r.contrast_tag);
            col_subs.push_back(r.subs.summary_string_partitioned);

            for (size_t k = 0; k < testable_groups.size(); ++k) {
                col_betas[k].push_back(r.beta[k]);
            }
            if (has_background) col_beta_bg.push_back(r.beta_background);
            if (testable_groups.size() > 2 && num_pairs > 0) {
                for (size_t pw_idx = 0; pw_idx < num_pairs; ++pw_idx) {
                    // Regularized/invariant sites carry no pairwise inference;
                    // keep columns aligned with neutral placeholders (0 LRT, 1.0 p),
                    // matching the legacy emitter's convention.
                    bool has_pw = pw_idx < r.pairwise.size();
                    col_pw_lrt[pw_idx].push_back(has_pw ? r.pairwise[pw_idx].lrt : 0.0);
                    col_pw_pval[pw_idx].push_back(has_pw ? r.pairwise[pw_idx].p_value : 1.0);
                    col_pw_padj[pw_idx].push_back(has_pw ? r.pairwise[pw_idx].p_value_adjusted : 1.0);
                }
            }
        }

        table.add_column("site", {"integer", "1-based codon site index"}, col_site);
        table.add_column("alpha", {"number", "Site synonymous substitution rate"}, col_a);
        table.add_column("alpha_null", {"number", "Synonymous rate under H0"}, col_a_null);
        table.add_column("beta_null", {"number", "Shared non-synonymous rate under H0"}, col_b_null);
        table.add_column("raw_lrt", {"number", "Unpenalized likelihood ratio test statistic"}, col_raw_lrt);
        table.add_column("penalized_lrt", {"number", "Penalized likelihood ratio test statistic"}, col_plrt);
        table.add_column("p_value", {"number", "Asymptotic p-value from chi^2_{K-1}"}, col_pval);
        table.add_column("q_value", {"number", "Benjamini-Hochberg false discovery rate"}, col_qval);
        table.add_column("contrast_tag", {"string", "Directional contrast classification"}, col_tag);
        table.add_column("substitutions", {"string", "Partitioned ancestral substitutions"}, col_subs);

        for (size_t k = 0; k < testable_groups.size(); ++k) {
            table.add_column("beta_" + testable_groups[k].name,
                             {"number", "Non-synonymous rate for group " + testable_groups[k].name},
                             col_betas[k]);
        }
        if (has_background) {
            table.add_column("beta_background",
                             {"number", "Nuisance non-synonymous rate on background branches"},
                             col_beta_bg);
        }
        if (testable_groups.size() > 2 && num_pairs > 0) {
            for (size_t pw_idx = 0; pw_idx < num_pairs; ++pw_idx) {
                const auto& pw0 = site_results[0].pairwise[pw_idx];
                std::string suffix = pw0.name_a + "_vs_" + pw0.name_b;
                table.add_column("lrt_" + suffix,
                                 {"number", "Pairwise penalized LRT for " + pw0.name_a + " vs " + pw0.name_b},
                                 col_pw_lrt[pw_idx]);
                table.add_column("p_value_" + suffix,
                                 {"number", "Uncorrected pairwise p-value for " + pw0.name_a + " vs " + pw0.name_b},
                                 col_pw_pval[pw_idx]);
                table.add_column("p_adjusted_" + suffix,
                                 {"number", "Holm-Bonferroni adjusted pairwise p-value for " + pw0.name_a + " vs " + pw0.name_b},
                                 col_pw_padj[pw_idx]);
            }
        }
        j["site_results"] = table.to_json();

        return j;
    }

nlohmann::json ContrastFELAnalyzer::to_legacy_json(
        const std::string& input_filepath,
        const std::string& tree_string) const {
        nlohmann::json j;

        j["analysis"] = {
            {"authors", "Sergei L Kosakovsky Pond, Blair Perry, Velislava N Petrova, Stephen D Shank, Ben Murrell, Todd A Castoe"},
            {"citation", "Contrast-FEL — A Test for Differences in Selective Pressures at Individual Sites among Clades and Sets of Branches (2021). _Mol Biol Evol_ 38(3): 1184–1198"},
            {"contact", "spond@temple.edu"},
            {"info", "Contrast-FEL (Fixed Effects Likelihood) tests for differences in non-synonymous substitution rates at individual sites between two or more sets of branches."},
            {"requirements", "in-frame codon alignment and a phylogenetic tree with branch sets"},
            {"version", "3.0.0"},
            {"settings", {
                {"pvalue", p_value_threshold},
                {"qvalue", q_value_threshold},
                {"srv", 1},
                {"permutation", 0}
            }}
        };

        j["input"] = {
            {"file name", input_filepath},
            {"genetic code", aln.code ? aln.code->name : "Universal"},
            {"number of sequences", aln.num_taxa},
            {"number of sites", aln.num_codons},
            {"partition count", 1}
        };
        j["input"]["trees"]["0"] = tree_string.empty() ? tree.to_newick() : tree_string;
        // Tested branches mapping
        nlohmann::json tested_map = nlohmann::json::object();
        nlohmann::json branch_attr = nlohmann::json::object();
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                int gid = node_group_id[node.id];
                std::string grp_name = (gid >= 0 && gid < static_cast<int>(testable_groups.size()))
                    ? testable_groups[gid].name
                    : "background";
                tested_map[node.name] = grp_name;

                Scalar orig_bl = (conversion_factor > 0.0) ? (node.branch_length / conversion_factor) : node.branch_length;
                branch_attr[node.name]["Global MG94xREV"] = orig_bl;
                if (has_gtr_fit) {
                    branch_attr[node.name]["Nucleotide GTR"] = gtr_branch_lengths.count(node.name) ? gtr_branch_lengths.at(node.name) : orig_bl;
                }
                branch_attr[node.name]["original name"] = node.name;
            }
        }
        j["tested"]["0"] = tested_map;
        j["branch attributes"]["0"] = branch_attr;
        j["branch attributes"]["attributes"]["Global MG94xREV"] = {{"attribute type", "branch length"}, {"display order", 1}};
        if (has_gtr_fit) {
            j["branch attributes"]["attributes"]["Nucleotide GTR"] = {{"attribute type", "branch length"}, {"display order", 0}};
        }

        // Fits
        if (has_gtr_fit) {
            nlohmann::json gtr_ef = nlohmann::json::array();
            for (int i = 0; i < 4; ++i) gtr_ef.push_back({aln.nuc_frequencies(i)});
            j["fits"]["Nucleotide GTR"] = {
                {"Log Likelihood", gtr_log_l},
                {"AIC-c", gtr_aicc},
                {"display order", 0},
                {"Equilibrium frequencies", gtr_ef},
                {"Rate Distributions", {
                    {"Substitution rate from nucleotide A to nucleotide C", gtr_rates.theta_AC},
                    {"Substitution rate from nucleotide A to nucleotide G", 1.0},
                    {"Substitution rate from nucleotide A to nucleotide T", gtr_rates.theta_AT},
                    {"Substitution rate from nucleotide C to nucleotide G", gtr_rates.theta_CG},
                    {"Substitution rate from nucleotide C to nucleotide T", gtr_rates.theta_CT},
                    {"Substitution rate from nucleotide G to nucleotide T", gtr_rates.theta_GT}
                }}
            };
        }

        nlohmann::json mg_ef = nlohmann::json::array();
        for (int i = 0; i < aln.codon_frequencies_f3x4.size(); ++i) mg_ef.push_back({aln.codon_frequencies_f3x4(i)});
        Scalar omega_mle = (base_params.alpha > 0.0) ? (base_params.beta / base_params.alpha) : 1.0;
        nlohmann::json mg_rates = nlohmann::json::object();
        for (const auto& grp : testable_groups) {
            mg_rates["non-synonymous/synonymous rate ratio for *" + grp.name + "*"] = nlohmann::json::array({nlohmann::json::array({omega_mle, 1.0})});
        }
        if (has_background) {
            mg_rates["non-synonymous/synonymous rate ratio for *background*"] = nlohmann::json::array({nlohmann::json::array({omega_mle, 1.0})});
        }
        j["fits"]["Global MG94xREV"] = {
            {"Log Likelihood", global_log_l},
            {"AIC-c", global_aicc},
            {"display order", 1},
            {"Equilibrium frequencies", mg_ef},
            {"Rate Distributions", mg_rates}
        };

        // MLE headers & content
        nlohmann::json headers = nlohmann::json::array();
        headers.push_back(nlohmann::json::array({"alpha", "Synonymous substitution rate at a site"}));
        for (const auto& grp : testable_groups) {
            headers.push_back(nlohmann::json::array({"beta (" + grp.name + ")", "Non-synonymous substitution rate at a site for " + grp.name + " branches"}));
        }
        if (has_background) {
            headers.push_back(nlohmann::json::array({"beta (background)", "Non-synonymous substitution rate at a site for background branches"}));
        }
        headers.push_back(nlohmann::json::array({"beta (overall)", "Non-synonymous substitution rate at a site under the null model"}));
        headers.push_back(nlohmann::json::array({"P-value (overall)", "P-value for the overall test of difference in selective pressure among tested branch groups"}));
        headers.push_back(nlohmann::json::array({"Q-value (overall)", "False discovery rate q-value for the overall test"}));

        size_t num_pairs = (!site_results.empty()) ? site_results[0].pairwise.size() : 0;
        for (size_t pw_idx = 0; pw_idx < num_pairs; ++pw_idx) {
            const auto& pw0 = site_results[0].pairwise[pw_idx];
            headers.push_back(nlohmann::json::array({
                "P-value for " + pw0.name_a + " vs " + pw0.name_b,
                "P-value for difference between " + pw0.name_a + " and " + pw0.name_b + " branches"
            }));
        }
        headers.push_back(nlohmann::json::array({"Permutation P-value", "Permutation p-value for difference in selective pressure"}));
        headers.push_back(nlohmann::json::array({"Total branch length", "The total length of branches contributing to inference at this site"}));

        j["MLE"]["headers"] = headers;

        nlohmann::json content = nlohmann::json::array();
        for (const auto& r : site_results) {
            nlohmann::json row = nlohmann::json::array();
            row.push_back(r.alpha);
            for (size_t k = 0; k < testable_groups.size(); ++k) {
                row.push_back(r.beta[k]);
            }
            if (has_background) {
                row.push_back(r.beta_background);
            }
            row.push_back(r.beta_null);
            row.push_back(r.p_value);
            row.push_back(r.q_value);
            for (size_t pw_idx = 0; pw_idx < num_pairs; ++pw_idx) {
                row.push_back((pw_idx < r.pairwise.size()) ? r.pairwise[pw_idx].p_value : 1.0);
            }
            row.push_back(1.0); // Permutation P-value placeholder when permutation=0
            row.push_back(r.total_tested_branch_length);
            content.push_back(row);
        }
        j["MLE"]["content"]["0"] = content;

        // Substitutions block
        nlohmann::json subs_json = nlohmann::json::object();
        for (size_t s = 0; s < site_results.size(); ++s) {
            const auto& r = site_results[s];
            std::string site_str = std::to_string(s);
            if (r.subs.events.empty()) {
                subs_json[site_str] = nullptr;
            } else {
                nlohmann::json site_obj = nlohmann::json::object();
                site_obj["root"] = r.subs.root_codon;
                for (const auto& ev : r.subs.events) {
                    site_obj[ev.branch_name] = ev.child_codon;
                }
                subs_json[site_str] = site_obj;
            }
        }
        j["substitutions"]["0"] = subs_json;

        return j;
    }

} // namespace hyphy::analyses
