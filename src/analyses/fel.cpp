// Out-of-line JSON emitter definitions for fel.
// Moved from include/hyphy/analyses/fel.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/fel.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json FELAnalyzer::to_legacy_json(const std::string& input_filepath, const std::string& tree_string) const {
        nlohmann::json j;

        // Analysis block
        j["analysis"]["authors"] = "Sergei L Kosakovsky Pond and Simon DW Frost";
        j["analysis"]["citation"] = "Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection (2005). _Mol Biol Evol_ 22 (5): 1208-1222";
        j["analysis"]["contact"] = "spond@temple.edu";
        j["analysis"]["info"] = "FEL (Fixed Effects Likelihood) estimates site-wise synonymous (alpha) and non-synonymous (beta) rates, and uses a likelihood ratio test to determine if beta != alpha at a site.";
        j["analysis"]["requirements"] = "in-frame codon alignment and a phylogenetic tree";
        j["analysis"]["settings"]["ci"] = 0;
        j["analysis"]["settings"]["multihit"] = "None";
        j["analysis"]["settings"]["pvalue"] = p_value_threshold;
        j["analysis"]["settings"]["resample"] = 0;
        j["analysis"]["settings"]["srv"] = 1;
        j["analysis"]["version"] = "3.0-next";

        // MLE block
        j["MLE"]["headers"] = nlohmann::json::array({
            nlohmann::json::array({"alpha", "Synonymous substitution rate at a site"}),
            nlohmann::json::array({"beta", "Non-synonymous substitution rate at a site"}),
            nlohmann::json::array({"alpha=beta", "The rate estimate under the neutral model"}),
            nlohmann::json::array({"LRT", "Likelihood ratio test statistic for beta = alpha, versus beta &neq; alpha"}),
            nlohmann::json::array({"p-value", "Asymptotic p-value for evidence of selection, i.e. beta &neq; alpha"}),
            nlohmann::json::array({"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"})
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

        // Input block
        j["input"]["file name"] = input_filepath;
        j["input"]["genetic code"] = aln.code ? aln.code->name : "Universal";
        j["input"]["number of sequences"] = aln.num_taxa;
        j["input"]["number of sites"] = aln.num_codons;
        j["input"]["partition count"] = 1;
        if (!tree_string.empty()) {
            j["input"]["trees"]["0"] = tree_string;
        }

        // Branch attributes
        nlohmann::json branch_attr;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar orig_bl = (conversion_factor > 0.0) ? (node.branch_length / conversion_factor) : node.branch_length;
                branch_attr[node.name]["Global MG94xREV"] = orig_bl;
                if (has_gtr_fit) {
                    branch_attr[node.name]["Nucleotide GTR"] = gtr_branch_lengths.count(node.name) ? gtr_branch_lengths.at(node.name) : orig_bl;
                }
                branch_attr[node.name]["original name"] = node.name;
            }
        }
        j["branch attributes"]["0"] = branch_attr;
        j["branch attributes"]["attributes"]["Global MG94xREV"]["attribute type"] = "branch length";
        j["branch attributes"]["attributes"]["Global MG94xREV"]["display order"] = 1;
        if (has_gtr_fit) {
            j["branch attributes"]["attributes"]["Nucleotide GTR"]["attribute type"] = "branch length";
            j["branch attributes"]["attributes"]["Nucleotide GTR"]["display order"] = 0;
        }

        // Fits
        if (has_gtr_fit) {
            j["fits"]["Nucleotide GTR"]["Log Likelihood"] = gtr_log_l;
            j["fits"]["Nucleotide GTR"]["AIC-c"] = gtr_aicc;
            j["fits"]["Nucleotide GTR"]["display order"] = 0;
            nlohmann::json gtr_ef = nlohmann::json::array();
            for (int i = 0; i < 4; ++i) {
                gtr_ef.push_back({aln.nuc_frequencies(i)});
            }
            j["fits"]["Nucleotide GTR"]["Equilibrium frequencies"] = gtr_ef;
            j["fits"]["Nucleotide GTR"]["Rate Distributions"] = {
                {"Substitution rate from nucleotide A to nucleotide C", gtr_rates.theta_AC},
                {"Substitution rate from nucleotide A to nucleotide G", 1.0},
                {"Substitution rate from nucleotide A to nucleotide T", gtr_rates.theta_AT},
                {"Substitution rate from nucleotide C to nucleotide G", gtr_rates.theta_CG},
                {"Substitution rate from nucleotide C to nucleotide T", gtr_rates.theta_CT},
                {"Substitution rate from nucleotide G to nucleotide T", gtr_rates.theta_GT}
            };
        }

        j["fits"]["Global MG94xREV"]["Log Likelihood"] = global_log_l;
        j["fits"]["Global MG94xREV"]["AIC-c"] = global_aicc;
        j["fits"]["Global MG94xREV"]["display order"] = 1;
        nlohmann::json mg_ef = nlohmann::json::array();
        for (int i = 0; i < aln.codon_frequencies_f3x4.size(); ++i) {
            mg_ef.push_back({aln.codon_frequencies_f3x4(i)});
        }
        j["fits"]["Global MG94xREV"]["Equilibrium frequencies"] = mg_ef;
        Scalar omega_mle = (base_params.alpha > 0.0) ? (base_params.beta / base_params.alpha) : 1.0;
        if (num_background > 0) {
            j["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"] = {
                {omega_mle, 1.0}
            };
            j["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *background*"] = {
                {omega_mle, 1.0}
            };
        } else {
            j["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"] = {
                {omega_mle, 1.0}
            };
        }

        // Tested branches (HyPhy 2.5 datamonkey compatibility)
        nlohmann::json tested_branches = nlohmann::json::object();
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                tested_branches[node.name] = is_test_branch[node.id] ? "test" : "background";
            }
        }
        j["tested"]["0"] = tested_branches;

        // Substitutions block (HyPhy 2.5 datamonkey / hyphy-vision compatibility)
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

nlohmann::json FELAnalyzer::to_modern_json(const std::string& input_filepath,
                                  const std::string& tree_string,
                                  const Provenance& prov) const {
        auto json_settings = nlohmann::json{
            {"pvalue_threshold", p_value_threshold},
            {"srv", true},
            {"ci", false},
            {"full_model", full_model},
            {"test_branches", settings.test_branches}
        };
        if (!settings.test_branch_regex.empty()) {
            json_settings["test_branch_regex"] = settings.test_branch_regex;
        }

        AnalysisMetadata meta{
            "fel",
            "Fixed Effects Likelihood",
            "3.0.0",
            "site_selection",
            "FEL (Fixed Effects Likelihood) estimates site-wise synonymous (alpha) and non-synonymous (beta) rates, and uses a likelihood ratio test to determine if beta != alpha at a site.",
            {citation("Kosakovsky Pond SL, Frost SDW (2005). Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection. Mol Biol Evol 22(5): 1208-1222.",
                      "10.1093/molbev/msi105", "15703277")},
            std::move(json_settings)
        };

        nlohmann::json j = new_analysis_json(meta, prov, input_filepath);
        j["dataset"] = dataset_block(aln);
        j["phylogeny"] = phylogeny_block(
            tree, tree_string,
            mg94_branch_lengths_fn(conversion_factor, has_gtr_fit, gtr_branch_lengths),
            [&](const TreeNode& node) -> nlohmann::json {
                if (node.id == tree.root_id) return nlohmann::json::object();
                return {{"tested", is_test_branch[node.id]}};
            });

        // Model fits
        if (has_gtr_fit) {
            j["model_fits"]["nucleotide_gtr"] = nucleotide_gtr_fit(
                aln, gtr_rates, gtr_log_l, gtr_aicc,
                5 + static_cast<int>(tree.nodes.size()) - 1);
        }
        j["model_fits"]["codon_mg94"] = fit_entry(
            "Muse-Gaut 1994 x GTR codon model with shared global omega",
            global_log_l, global_aicc,
            6 + static_cast<int>(tree.nodes.size()) - 1);
        j["model_fits"]["codon_mg94"]["global_parameters"] = {
            {"omega", (base_params.alpha > 0.0) ? (base_params.beta / base_params.alpha) : 1.0},
            {"tree_scale", 1.0}
        };

        // Statistical tests
        auto summary = get_summary();
        j["statistical_tests"]["site_level_summary"] = {
            {"method", "Fixed Effects Likelihood (FEL)"},
            {"test_statistic", "Likelihood Ratio Test (LRT)"},
            {"distribution", "Asymptotic Chi-squared"},
            {"degrees_of_freedom", 1},
            {"null_hypothesis", "beta = alpha (Neutral evolution)"},
            {"alternative_hypothesis", "beta != alpha (Positive or negative selection)"},
            {"threshold", p_value_threshold},
            {"total_sites", static_cast<int>(summary.tested_sites)},
            {"diversifying_count", static_cast<int>(summary.positively_selected)},
            {"purifying_count", static_cast<int>(summary.negatively_selected)},
            {"neutral_count", static_cast<int>(summary.tested_sites - summary.positively_selected - summary.negatively_selected)},
            {"tested_branches_count", static_cast<int>(num_tested)},
            {"background_branches_count", static_cast<int>(num_background)},
            {"total_tested_tree_length", (conversion_factor > 0.0) ? (baseline_tested_tree_length / conversion_factor) : baseline_tested_tree_length}
        };

        // Site results (columnar; ColumnarResults enforces equal-length)
        ColumnarResults table;
        const size_t n = site_results.size();
        const bool has_bg = num_background > 0;
        std::vector<int> col_site, col_syn, col_nonsyn, col_syn_bg, col_nonsyn_bg,
                         col_1nt, col_2nt, col_3nt;
        std::vector<double> col_alpha, col_beta, col_beta_bg, col_alpha_null,
                            col_lrt, col_pval, col_tbl, col_surprise, col_mut_bl;
        std::vector<std::string> col_class, col_history;

        col_site.reserve(n); col_alpha.reserve(n); col_beta.reserve(n);
        if (has_bg) col_beta_bg.reserve(n);
        col_alpha_null.reserve(n); col_lrt.reserve(n); col_pval.reserve(n);
        col_tbl.reserve(n); col_class.reserve(n);
        col_syn.reserve(n); col_nonsyn.reserve(n);
        if (has_bg) { col_syn_bg.reserve(n); col_nonsyn_bg.reserve(n); }
        col_1nt.reserve(n); col_2nt.reserve(n); col_3nt.reserve(n);
        col_surprise.reserve(n); col_mut_bl.reserve(n); col_history.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            const auto& r = site_results[i];
            col_site.push_back(static_cast<int>(i + 1));
            col_alpha.push_back(r.alpha);
            col_beta.push_back(r.beta);
            if (has_bg) col_beta_bg.push_back(r.beta_nuisance);
            col_alpha_null.push_back(r.alpha_null);
            col_lrt.push_back(r.lrt);
            col_pval.push_back(r.p_value);
            col_tbl.push_back(r.total_branch_length);

            if (r.p_value <= p_value_threshold) {
                if (r.beta > r.alpha) {
                    col_class.push_back("diversifying");
                } else if (r.alpha > r.beta) {
                    col_class.push_back("purifying");
                } else {
                    col_class.push_back("neutral");
                }
            } else {
                col_class.push_back("neutral");
            }

            col_syn.push_back(static_cast<int>(r.subs.syn_test));
            col_nonsyn.push_back(static_cast<int>(r.subs.nonsyn_test));
            if (has_bg) {
                col_syn_bg.push_back(static_cast<int>(r.subs.syn_bg));
                col_nonsyn_bg.push_back(static_cast<int>(r.subs.nonsyn_bg));
            }
            col_1nt.push_back(static_cast<int>(r.subs.subs_1nt));
            col_2nt.push_back(static_cast<int>(r.subs.subs_2nt));
            col_3nt.push_back(static_cast<int>(r.subs.subs_3nt));
            col_surprise.push_back(r.subs.mean_surprise_bits);
            col_mut_bl.push_back(r.subs.mean_branch_length);
            col_history.push_back(has_bg ? r.subs.summary_string_partitioned : r.subs.summary_string);
        }

        table.add_column("site", {"integer", "Alignment codon site index (1-based)", "1-based codon position"}, col_site);
        table.add_column("alpha", {"number", "Synonymous substitution rate estimate (alpha)", "substitutions/site"}, col_alpha);
        table.add_column("beta", {"number", "Non-synonymous substitution rate estimate on tested branches (beta)", "substitutions/site"}, col_beta);
        if (has_bg) {
            table.add_column("beta_background", {"number", "Non-synonymous substitution rate estimate on background branches (beta_background)", "substitutions/site"}, col_beta_bg);
        }
        table.add_column("alpha_null", {"number", "Constrained neutral rate (alpha = beta)", "substitutions/site"}, col_alpha_null);
        table.add_column("lrt", {"number", "Likelihood ratio test statistic for beta != alpha"}, col_lrt);
        table.add_column("p_value", {"number", "Asymptotic p-value from chi^2_1 distribution"}, col_pval);
        table.add_column("total_branch_length", {"number", "Total length of tested branches contributing to inference at this site", "substitutions/site"}, col_tbl);
        table.add_column("classification", {"string", "Selection classification at threshold (diversifying, purifying, neutral)"}, col_class);
        table.add_column("synonymous_substitutions", {"integer", "Count of synonymous substitutions on tested branches"}, col_syn);
        table.add_column("nonsynonymous_substitutions", {"integer", "Count of non-synonymous substitutions on tested branches"}, col_nonsyn);
        if (has_bg) {
            table.add_column("synonymous_substitutions_background", {"integer", "Count of synonymous substitutions on background branches"}, col_syn_bg);
            table.add_column("nonsynonymous_substitutions_background", {"integer", "Count of non-synonymous substitutions on background branches"}, col_nonsyn_bg);
        }
        table.add_column("substitutions_1nt", {"integer", "Count of single-nucleotide codon substitutions"}, col_1nt);
        table.add_column("substitutions_2nt", {"integer", "Count of double-nucleotide codon substitutions"}, col_2nt);
        table.add_column("substitutions_3nt", {"integer", "Count of triple-nucleotide codon substitutions"}, col_3nt);
        table.add_column("mean_substitution_surprise", {"number", "Mean Shannon self-information (-log2 P) of substitutions at this site", "bits"}, col_surprise);
        table.add_column("mean_substitution_branch_length", {"number", "Mean length of branches where substitutions occurred at this site", "substitutions/site"}, col_mut_bl);
        table.add_column("substitution_history", {"string", "Formatted substitution history summary (AA changes, hit counts, surprise)"}, col_history);
        j["site_results"] = table.to_json();

        // Compressed substitutions mapping (HyPhy 2.5 datamonkey / hyphy-vision parity)
        nlohmann::json subs_comp = nlohmann::json::object();
        for (size_t s = 0; s < site_results.size(); ++s) {
            const auto& r = site_results[s];
            std::string site_str = std::to_string(s);
            if (r.subs.events.empty()) {
                subs_comp[site_str] = nullptr;
            } else {
                nlohmann::json site_obj = nlohmann::json::object();
                site_obj["root"] = r.subs.root_codon;
                for (const auto& ev : r.subs.events) {
                    site_obj[ev.branch_name] = ev.child_codon;
                }
                subs_comp[site_str] = site_obj;
            }
        }
        j["substitutions"]["0"] = subs_comp;

        // Detailed substitution events per site
        nlohmann::json sub_events_json = nlohmann::json::array();
        for (size_t s = 0; s < site_results.size(); ++s) {
            const auto& r = site_results[s];
            if (!r.subs.events.empty()) {
                nlohmann::json site_ev_obj;
                site_ev_obj["site"] = static_cast<int>(s + 1);
                site_ev_obj["root_codon"] = r.subs.root_codon;
                site_ev_obj["root_aa"] = std::string(1, r.subs.root_aa);
                nlohmann::json ev_arr = nlohmann::json::array();
                for (const auto& ev : r.subs.events) {
                    nlohmann::json ev_json = {
                        {"branch", ev.branch_name},
                        {"tested", ev.is_tested},
                        {"branch_length", ev.branch_length},
                        {"parent_codon", ev.parent_codon},
                        {"child_codon", ev.child_codon},
                        {"parent_aa", std::string(1, ev.parent_aa)},
                        {"child_aa", std::string(1, ev.child_aa)},
                        {"nt_distance", ev.nt_distance},
                        {"is_synonymous", ev.is_synonymous},
                        {"transition_probability", ev.transition_prob},
                        {"surprise_bits", ev.surprise_bits}
                    };
                    ev_arr.push_back(ev_json);
                }
                site_ev_obj["events"] = ev_arr;
                sub_events_json.push_back(site_ev_obj);
            }
        }
        j["substitution_events"] = sub_events_json;

        return j;
    }

} // namespace hyphy::analyses
