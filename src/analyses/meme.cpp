// Out-of-line JSON emitter definitions for meme.
// Moved from include/hyphy/analyses/meme.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/meme.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json MEMEAnalyzer::to_legacy_json() const {
        nlohmann::json j;

        // Analysis metadata
        j["analysis"] = {
            {"authors", "Sergei L. Kosakovsky Pond, Steven Weaver"},
            {"citation", "Detecting Individual Sites Subject to Episodic Diversifying Selection. PLoS Genet 8(7): e1002764."},
            {"contact", "spond@temple.edu"},
            {"info", "MEME (Mixed Effects Model of Evolution) estimates a site-wise synonymous (alpha) and a two-category mixture of non-synonymous (beta1 with p1, and beta+ with p+) rates, and uses a likelihood ratio test to determine if beta+ > alpha at a site."},
            {"requirements", "in-frame codon alignment and a phylogenetic tree"},
            {"settings", {{"p-value", p_value_threshold}, {"code", aln.code ? aln.code->name : "Universal"}}},
            {"version", "4.1 (Modern C++ hyphy-next)"}
        };

        // Input
        j["input"]["file name"] = input_filepath;
        j["input"]["genetic code"] = aln.code ? aln.code->name : "Universal";
        j["input"]["number of sequences"] = aln.num_taxa;
        j["input"]["number of sites"] = aln.num_codons;
        j["input"]["partition count"] = 1;
        if (!tree_string.empty()) {
            j["input"]["trees"]["0"] = tree_string;
        }

        // Fits
        nlohmann::json fits_json;
        if (has_gtr_fit) {
            fits_json["Nucleotide GTR"] = {
                {"AIC-c", gtr_aicc},
                {"Log Likelihood", gtr_log_l},
                {"Equilibrium frequencies", {
                    aln.nuc_frequencies(0),
                    aln.nuc_frequencies(1),
                    aln.nuc_frequencies(2),
                    aln.nuc_frequencies(3)
                }},
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

        std::vector<Scalar> eq_freqs(aln.codon_frequencies_f3x4.data(), aln.codon_frequencies_f3x4.data() + aln.codon_frequencies_f3x4.size());
        fits_json["Global MG94xREV"] = {
            {"AIC-c", global_aicc},
            {"Log Likelihood", global_log_l},
            {"Equilibrium frequencies", eq_freqs},
            {"Rate Distributions", {
                {"non-synonymous/synonymous rate ratio", base_params.beta / std::max(base_params.alpha, 1e-12)}
            }}
        };
        j["fits"] = fits_json;

        // Branch Attributes
        nlohmann::json branch_attrs;
        nlohmann::json branch_dict;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                nlohmann::json b_info;
                b_info["Global MG94xREV"] = node.branch_length / conversion_factor;
                auto it = gtr_branch_lengths.find(node.name);
                if (it != gtr_branch_lengths.end()) {
                    b_info["Nucleotide GTR"] = it->second;
                }
                branch_dict[node.name] = b_info;
            }
        }
        branch_attrs["0"] = branch_dict;
        j["branch attributes"] = branch_attrs;

        // MLE headers and content
        j["MLE"]["headers"] = nlohmann::json::array({
            nlohmann::json::array({"&alpha;", "Synonymous substitution rate at a site"}),
            nlohmann::json::array({"&beta;<sup>1</sup>", "Non-synonymous substitution rate at a site for the negative/neutral evolution component 1"}),
            nlohmann::json::array({"p<sup>1</sup>", "Mixture distribution weight allocated to negative/neutral evolution component 1"}),
            nlohmann::json::array({"&beta;<sup>+</sup>", "Non-synonymous substitution rate at a site for the positive selection component"}),
            nlohmann::json::array({"p<sup>+</sup>", "Mixture distribution weight allocated to the positive selection component"}),
            nlohmann::json::array({"LRT", "Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> > 0 and &beta;<sup>+</sup> > &alpha;"}),
            nlohmann::json::array({"p-value", "Asymptotic p-value for episodic diversification, i.e., p<sup>+</sup> > 0 and &beta;<sup>+</sup> > &alpha;"}),
            nlohmann::json::array({"# branches under selection", "The count of branches with empirical Bayes factor >= 100 for the &beta;<sup>+</sup> rate"}),
            nlohmann::json::array({"Total branch length", "The total length of branches contributing to inference at this site"}),
            nlohmann::json::array({"MEME LogL", "Site Log-likelihood under the MEME model"}),
            nlohmann::json::array({"FEL LogL", "Site Log-likelihood under the FEL model"}),
            nlohmann::json::array({"LRT MEME vs FEL", "Likelihood ratio test statistic p-value for MEME vs FEL"}),
            nlohmann::json::array({"FEL &alpha;", "Synonymous substitution rate at a site under the FEL model"}),
            nlohmann::json::array({"FEL &beta;", "Non-synonymous substitution rate at a site under the FEL model"})
        });

        std::vector<std::vector<Scalar>> content(site_results.size());
        for (size_t s = 0; s < site_results.size(); ++s) {
            const auto& r = site_results[s];
            content[s] = {
                r.alpha,
                r.beta1,
                r.p1,
                r.beta_plus,
                r.p_plus,
                r.lrt,
                r.p_value,
                static_cast<Scalar>(r.branches_under_selection),
                r.total_branch_length,
                r.log_l_meme,
                r.log_l_fel,
                r.p_value_meme_vs_fel,
                r.fel_alpha,
                r.fel_beta
            };
        }
        j["MLE"]["content"]["0"] = content;

        // Tested branches
        nlohmann::json tested_dict;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                bool is_test = (static_cast<size_t>(node.id) < is_test_branch.size()) ? is_test_branch[node.id] : true;
                tested_dict[node.name] = is_test ? "test" : "background";
            }
        }
        j["tested"]["0"] = tested_dict;

        return j;
    }

nlohmann::json MEMEAnalyzer::to_modern_json(const Provenance& prov) const {
        AnalysisMetadata meta{
            "meme",
            "Mixed Effects Model of Evolution",
            "3.0.0",
            "episodic_site_selection",
            "MEME (Mixed Effects Model of Evolution) estimates site-wise synonymous (alpha) and two-category non-synonymous rates (beta1 with p1, and beta+ with p+), using a likelihood ratio test for episodic diversifying selection.",
            {citation("Murrell B, Wertheim JO, Moola S, Weighill T, Scheffler K, Kosakovsky Pond SL (2012). Detecting Individual Sites Subject to Episodic Diversifying Selection. PLoS Genet 8(7): e1002764.",
                      "10.1371/journal.pgen.1002764", "22844254")},
            {
                {"pvalue_threshold", p_value_threshold},
                {"code", aln.code ? aln.code->name : "Universal"},
                {"test_branches", test_branches}
            }
        };

        nlohmann::json j = new_analysis_json(meta, prov, input_filepath);
        j["dataset"] = dataset_block(aln);
        j["phylogeny"] = phylogeny_block(
            tree, tree_string,
            mg94_branch_lengths_fn(conversion_factor, has_gtr_fit, gtr_branch_lengths),
            [this](const TreeNode& node) {
                nlohmann::json extras;
                if (node.id != tree.root_id) {
                    extras["tested"] =
                        (static_cast<size_t>(node.id) < is_test_branch.size())
                            ? is_test_branch[node.id] : true;
                }
                if (!node.model_tag.empty()) {
                    extras["model_tags"] = nlohmann::json::array({node.model_tag});
                }
                return extras;
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
        int diversifying_sites = 0;
        for (const auto& r : site_results) {
            if (r.p_value <= p_value_threshold) diversifying_sites++;
        }
        j["statistical_tests"]["site_level_summary"] = {
            {"method", "Mixed Effects Model of Evolution (MEME)"},
            {"test_statistic", "Likelihood Ratio Test (LRT)"},
            {"distribution", "Asymptotic Chi-squared Mixture (0.333 chi^2_0 + 0.300 chi^2_1 + 0.367 chi^2_2)"},
            {"degrees_of_freedom", 2},
            {"null_hypothesis", "beta+ <= alpha or p+ = 0 (No episodic positive selection)"},
            {"alternative_hypothesis", "beta+ > alpha and p+ > 0 (Episodic positive selection on a subset of branches)"},
            {"threshold", p_value_threshold},
            {"total_sites", static_cast<int>(site_results.size())},
            {"diversifying_count", diversifying_sites}
        };

        // Site results (columnar; ColumnarResults enforces equal-length)
        ColumnarResults table;
        const size_t n = site_results.size();
        std::vector<int> col_site, col_branches, col_syn, col_nonsyn;
        std::vector<double> col_alpha, col_beta1, col_p1, col_beta_plus, col_p_plus,
                            col_lrt, col_pval, col_tbl, col_ll_meme, col_ll_fel,
                            col_surprise, col_perplexity;
        std::vector<std::string> col_class, col_history;
        nlohmann::json col_branch_ebf = nlohmann::json::array();

        col_site.reserve(n); col_alpha.reserve(n); col_beta1.reserve(n);
        col_p1.reserve(n); col_beta_plus.reserve(n); col_p_plus.reserve(n);
        col_lrt.reserve(n); col_pval.reserve(n); col_branches.reserve(n);
        col_tbl.reserve(n); col_ll_meme.reserve(n); col_ll_fel.reserve(n);
        col_class.reserve(n); col_history.reserve(n);
        col_surprise.reserve(n); col_perplexity.reserve(n);
        col_syn.reserve(n); col_nonsyn.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            const auto& r = site_results[i];
            col_site.push_back(static_cast<int>(i + 1));
            col_alpha.push_back(r.alpha);
            col_beta1.push_back(r.beta1);
            col_p1.push_back(r.p1);
            col_beta_plus.push_back(r.beta_plus);
            col_p_plus.push_back(r.p_plus);
            col_lrt.push_back(r.lrt);
            col_pval.push_back(r.p_value);
            col_branches.push_back(r.branches_under_selection);
            col_tbl.push_back(r.total_branch_length);
            col_ll_meme.push_back(r.log_l_meme);
            col_ll_fel.push_back(r.log_l_fel);
            col_class.push_back((r.p_value <= p_value_threshold) ? "diversifying" : "neutral");

            col_history.push_back(r.subs.summary_string);
            col_surprise.push_back(r.subs.mean_surprise_bits);
            col_perplexity.push_back(r.subs.site_perplexity);
            col_syn.push_back(static_cast<int>(r.subs.syn_test + r.subs.syn_bg));
            col_nonsyn.push_back(static_cast<int>(r.subs.nonsyn_test + r.subs.nonsyn_bg));

            nlohmann::json site_ebf = nlohmann::json::object();
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && node.id >= 0 &&
                    static_cast<size_t>(node.id) < r.branch_ebf.size()) {
                    Scalar ebf = r.branch_ebf[node.id];
                    if (ebf >= 1.0) {
                        site_ebf[node.name] = finite_or_null(ebf);
                    }
                }
            }
            col_branch_ebf.push_back(site_ebf);
        }

        table.add_column("site", {"integer", "Alignment codon site index (1-based)", "1-based codon position"}, col_site);
        table.add_column("alpha", {"number", "Synonymous substitution rate estimate (alpha)", "substitutions/site"}, col_alpha);
        table.add_column("beta1", {"number", "Non-synonymous rate for negative/neutral component", "substitutions/site"}, col_beta1);
        table.add_column("p1", {"number", "Weight allocated to negative/neutral component"}, col_p1);
        table.add_column("beta_plus", {"number", "Non-synonymous rate for positive selection component", "substitutions/site"}, col_beta_plus);
        table.add_column("p_plus", {"number", "Weight allocated to positive selection component"}, col_p_plus);
        table.add_column("lrt", {"number", "Likelihood ratio test statistic for episodic diversification"}, col_lrt);
        table.add_column("p_value", {"number", "Asymptotic p-value from mixture chi^2 distribution"}, col_pval);
        table.add_column("branches_under_selection", {"integer", "Number of branches with empirical Bayes factor >= 100"}, col_branches);
        table.add_column("total_branch_length", {"number", "Total length of branches contributing to inference at this site", "substitutions/site"}, col_tbl);
        table.add_column("log_l_meme", {"number", "Site log-likelihood under MEME model"}, col_ll_meme);
        table.add_column("log_l_fel", {"number", "Site log-likelihood under FEL model"}, col_ll_fel);
        table.add_column("classification", {"string", "Selection classification at threshold (diversifying, neutral)"}, col_class);
        table.add_column("substitution_history", {"string", "Formatted string of observed substitution events at this site"}, col_history);
        table.add_column("mean_substitution_surprise", {"number", "Mean Shannon perplexity surprise across substitution events at this site", "bits"}, col_surprise);
        table.add_column("site_perplexity", {"number", "Empirical site perplexity metric"}, col_perplexity);
        table.add_column("synonymous_substitutions", {"integer", "Count of synonymous substitutions at this site"}, col_syn);
        table.add_column("nonsynonymous_substitutions", {"integer", "Count of non-synonymous substitutions at this site"}, col_nonsyn);
        table.add_column("branch_ebf", {"object", "Map of branch names to Empirical Bayes Factors (EBF >= 1.0) under positive selection"}, col_branch_ebf);
        j["site_results"] = table.to_json();

        // Substitution events block (Root level)
        nlohmann::json sub_events_arr = nlohmann::json::array();
        for (size_t i = 0; i < site_results.size(); ++i) {
            const auto& r = site_results[i];
            if (r.subs.events.empty()) continue;

            nlohmann::json site_entry;
            site_entry["site"] = static_cast<int>(i + 1);
            site_entry["root_codon"] = r.subs.root_codon;
            site_entry["root_aa"] = std::string(1, r.subs.root_aa);

            nlohmann::json ev_list = nlohmann::json::array();
            for (const auto& ev : r.subs.events) {
                nlohmann::json ev_j;
                ev_j["branch"] = ev.branch_name;
                ev_j["parent_codon"] = ev.parent_codon;
                ev_j["child_codon"] = ev.child_codon;
                ev_j["parent_aa"] = std::string(1, ev.parent_aa);
                ev_j["child_aa"] = std::string(1, ev.child_aa);
                ev_j["is_synonymous"] = ev.is_synonymous;
                ev_j["nucleotide_differences"] = ev.nt_distance;
                ev_j["transition_probability"] = ev.transition_prob;
                ev_j["surprise_bits"] = ev.surprise_bits;
                ev_j["branch_length"] = ev.branch_length;
                ev_j["is_tested"] = ev.is_tested;
                ev_list.push_back(ev_j);
            }
            site_entry["events"] = ev_list;
            sub_events_arr.push_back(site_entry);
        }
        j["substitution_events"] = sub_events_arr;

        return j;
    }

} // namespace hyphy::analyses
