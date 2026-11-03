// Out-of-line JSON emitter definitions for busted.
// Moved from include/hyphy/analyses/busted.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/busted.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json BUSTEDAnalyzer::to_legacy_json(const BUSTEDResult& res) const {
        nlohmann::json j;

        if (res.is_failed) {
            j["status"] = "failed";
            j["error"] = res.error_message;
        } else {
            j["status"] = "completed";
        }

        j["analysis"] = {
            {"authors", "Sergei L Kosakovsky Pond"},
            {"citation", "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71"},
            {"contact", "spond@temple.edu"},
            {"version", "4.7-next"},
            {"settings", {
                {"error-sink", res.settings.error_sink ? 1 : 0},
                {"mss", res.settings.mss_type.empty() ? "No" : res.settings.mss_type},
                {"multiple-hit", res.settings.multiple_hits},
                {"srv", (!res.unconstrained.test_distribution.syn_rates.empty() && res.unconstrained.test_distribution.syn_rates.size() > 1) ? "Yes" : "No"},
                {"optimal_k", res.optimal_k}
            }}
        };
        if (!res.mss_rates.empty()) {
            nlohmann::json mss_block = nlohmann::json::object();
            for (const auto& r : res.mss_rates) {
                mss_block[r.class_name] = r.rate;
            }
            j["analysis"]["settings"]["mss_rates"] = mss_block;
            j["analysis"]["settings"]["mss_normalization"] = res.settings.mss_normalization;
        }

        j["test results"] = {
            {"LRT", res.lrt},
            {"p-value", res.p_value},
            {"optimal_k", res.optimal_k}
        };
        if (res.has_background_test) {
            j["test results background"] = {
                {"LRT", res.background_lrt},
                {"p-value", res.background_p_value},
                {"delta AIC-c", res.background_delta_aicc}
            };
            j["test results shared distributions"] = {
                {"LRT", res.shared_lrt},
                {"p-value", res.shared_p_value},
                {"df", res.shared_df}
            };
            j["BUSTED-PH"] = {
                {"Uncorrected P-value", {
                    {"FG", res.p_value},
                    {"BG", res.background_p_value},
                    {"Comparative", res.shared_p_value}
                }},
                {"union_intersection_p_value", res.union_intersection_p_value},
                {"background_delta_aicc", res.background_delta_aicc},
                {"background_neutral_supported", res.background_neutral_supported},
                {"ranking_score", res.ranking_score},
                {"Level", res.settings.p_value_threshold},
                {"Summary", res.phenotype_verdict}
            };
            j["phenotype association"] = {
                {"associated", res.phenotype_associated},
                {"union_intersection_p_value", res.union_intersection_p_value},
                {"background_delta_aicc", res.background_delta_aicc},
                {"background_neutral_supported", res.background_neutral_supported},
                {"ranking_score", res.ranking_score},
                {"verdict", res.phenotype_verdict}
            };
        }
        if (res.has_relax_test) {
            j["test results relaxation"] = {
                {"LRT", res.relax_lrt},
                {"p-value", res.relax_p_value},
                {"relaxation or intensification parameter", res.relax_k},
                {"direction", res.relax_direction}
            };
        }

        // Fits
        nlohmann::json fits;
        if (has_gtr_fit) {
            fits["Nucleotide GTR"] = {
                {"AIC-c", gtr_aicc},
                {"Log Likelihood", gtr_log_l},
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

        nlohmann::json mg94_rates = {
            {"non-synonymous/synonymous rate ratio for *test*", {{mg94_omega, 1.0}}}
        };
        if (mg94_background_omega >= 0.0) {
            mg94_rates["non-synonymous/synonymous rate ratio for *background*"] = {{mg94_background_omega, 1.0}};
        }
        fits["MG94xREV with separate rates for branch sets"] = {
            {"Log Likelihood", mg94_log_l},
            {"Rate Distributions", mg94_rates}
        };

        auto format_distro = [&](const BUSTEDRateDistribution& d) {
            nlohmann::json distro;
            Scalar bio_scale = res.settings.error_sink ? (1.0 - d.error_sink_weight) : 1.0;
            for (size_t k = 0; k < d.omegas.size(); ++k) {
                nlohmann::json entry = {
                    {"omega", d.omegas[k]},
                    {"proportion", bio_scale * d.weights[k]}
                };
                if (k < d.annotations.size() && !d.annotations[k].empty()) {
                    entry["annotation"] = d.annotations[k];
                }
                distro[std::to_string(k)] = entry;
            }
            if (res.settings.error_sink) {
                distro["error_sink"] = {
                    {"type", "saturation_replacement"},
                    {"proportion", d.error_sink_weight}
                };
            }
            return distro;
        };

        auto format_srv_distro = [](const BUSTEDRateDistribution& d) {
            nlohmann::json distro;
            for (size_t m = 0; m < d.syn_rates.size(); ++m) {
                distro[std::to_string(m)] = {
                    {"rate", d.syn_rates[m]},
                    {"proportion", d.syn_weights[m]}
                };
            }
            return distro;
        };

        auto format_mss_fit_rates = [&](const BUSTEDFit& fit) -> nlohmann::json {
            nlohmann::json mss_obj = nlohmann::json::object();
            if (!res.mss_rates.empty() && fit.mss_class_rates.size() == res.mss_rates.size()) {
                for (size_t c = 0; c < res.mss_rates.size(); ++c) {
                    mss_obj[res.mss_rates[c].class_name] = fit.mss_class_rates[c];
                }
            }
            return mss_obj;
        };

        nlohmann::json unc_rates = {{"Test", format_distro(res.unconstrained.test_distribution)}};
        if (res.unconstrained.has_background_distribution) {
            unc_rates["Background"] = format_distro(res.unconstrained.background_distribution);
        }
        if (res.unconstrained.has_unclassified_distribution) {
            unc_rates["Unclassified"] = format_distro(res.unconstrained.unclassified_distribution);
        }
        if (!res.unconstrained.test_distribution.syn_rates.empty() && res.unconstrained.test_distribution.syn_rates.size() > 1) {
            unc_rates["Synonymous site-to-site rates"] = format_srv_distro(res.unconstrained.test_distribution);
        }
        nlohmann::json unc_fit = {
            {"AIC-c", res.unconstrained.aicc},
            {"Log Likelihood", res.unconstrained.log_likelihood},
            {"Rate Distributions", unc_rates}
        };
        if (!res.mss_rates.empty() && res.unconstrained.mss_class_rates.size() == res.mss_rates.size()) {
            unc_fit["MSS Synonymous Rates"] = format_mss_fit_rates(res.unconstrained);
        }
        if (res.settings.multiple_hits != "None") {
            unc_fit["rate at which 2 nucleotides are changed instantly within a single codon"] = res.unconstrained.delta;
            unc_fit["Fraction of subs rate at which 2 nucleotides are changed instantly within a single codon"] = res.unconstrained.frac_delta;
            if (res.settings.multiple_hits == "Double+Triple") {
                unc_fit["rate at which 3 nucleotides are changed instantly within a single codon"] = res.unconstrained.psi;
                unc_fit["Fraction of subs rate at which 3 nucleotides are changed instantly within a single codon"] = res.unconstrained.frac_psi;
            }
        }
        fits["Unconstrained model"] = unc_fit;

        nlohmann::json con_rates = {{"Test", format_distro(res.constrained.test_distribution)}};
        if (res.constrained.has_background_distribution) {
            con_rates["Background"] = format_distro(res.constrained.background_distribution);
        }
        if (res.constrained.has_unclassified_distribution) {
            con_rates["Unclassified"] = format_distro(res.constrained.unclassified_distribution);
        }
        if (!res.constrained.test_distribution.syn_rates.empty() && res.constrained.test_distribution.syn_rates.size() > 1) {
            con_rates["Synonymous site-to-site rates"] = format_srv_distro(res.constrained.test_distribution);
        }
        nlohmann::json con_fit = {
            {"AIC-c", res.constrained.aicc},
            {"Log Likelihood", res.constrained.log_likelihood},
            {"Rate Distributions", con_rates}
        };
        if (!res.mss_rates.empty() && res.constrained.mss_class_rates.size() == res.mss_rates.size()) {
            con_fit["MSS Synonymous Rates"] = format_mss_fit_rates(res.constrained);
        }
        if (res.settings.multiple_hits != "None") {
            con_fit["rate at which 2 nucleotides are changed instantly within a single codon"] = res.constrained.delta;
            con_fit["Fraction of subs rate at which 2 nucleotides are changed instantly within a single codon"] = res.constrained.frac_delta;
            if (res.settings.multiple_hits == "Double+Triple") {
                con_fit["rate at which 3 nucleotides are changed instantly within a single codon"] = res.constrained.psi;
                con_fit["Fraction of subs rate at which 3 nucleotides are changed instantly within a single codon"] = res.constrained.frac_psi;
            }
        }
        fits["Constrained model"] = con_fit;

        if (res.has_background_test) {
            nlohmann::json con_bg_rates = {
                {"Test", format_distro(res.constrained_background.test_distribution)},
                {"Background", format_distro(res.constrained_background.background_distribution)}
            };
            if (res.constrained_background.has_unclassified_distribution) {
                con_bg_rates["Unclassified"] = format_distro(res.constrained_background.unclassified_distribution);
            }
            fits["Constrained model (background)"] = {
                {"AIC-c", res.constrained_background.aicc},
                {"Log Likelihood", res.constrained_background.log_likelihood},
                {"Rate Distributions", con_bg_rates}
            };
            if (!res.mss_rates.empty() && res.constrained_background.mss_class_rates.size() == res.mss_rates.size()) {
                fits["Constrained model (background)"]["MSS Synonymous Rates"] = format_mss_fit_rates(res.constrained_background);
            }
        }
        if (res.has_background_test || res.has_relax_test) {
            nlohmann::json shared_rates = {
                {"Test", format_distro(res.shared_distributions.test_distribution)},
                {"Background", format_distro(res.shared_distributions.background_distribution)}
            };
            if (res.shared_distributions.has_unclassified_distribution) {
                shared_rates["Unclassified"] = format_distro(res.shared_distributions.unclassified_distribution);
            }
            fits["Unconstrained model with shared distributions"] = {
                {"AIC-c", res.shared_distributions.aicc},
                {"Log Likelihood", res.shared_distributions.log_likelihood},
                {"Rate Distributions", shared_rates}
            };
            if (!res.mss_rates.empty() && res.shared_distributions.mss_class_rates.size() == res.mss_rates.size()) {
                fits["Unconstrained model with shared distributions"]["MSS Synonymous Rates"] = format_mss_fit_rates(res.shared_distributions);
            }
        }
        if (res.has_relax_test) {
            nlohmann::json relax_rates = {
                {"Test", format_distro(res.relax_alternative.test_distribution)},
                {"Background", format_distro(res.relax_alternative.background_distribution)}
            };
            if (res.relax_alternative.has_unclassified_distribution) {
                relax_rates["Unclassified"] = format_distro(res.relax_alternative.unclassified_distribution);
            }
            fits["RELAX alternative (coupled-k)"] = {
                {"AIC-c", res.relax_alternative.aicc},
                {"Log Likelihood", res.relax_alternative.log_likelihood},
                {"relaxation or intensification parameter", res.relax_k},
                {"Rate Distributions", relax_rates}
            };
            if (!res.mss_rates.empty() && res.relax_alternative.mss_class_rates.size() == res.mss_rates.size()) {
                fits["RELAX alternative (coupled-k)"]["MSS Synonymous Rates"] = format_mss_fit_rates(res.relax_alternative);
            }
        }

        j["fits"] = fits;

        // Branch attributes
        nlohmann::json branch_attr;
        nlohmann::json tested_map = nlohmann::json::object();
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                nlohmann::json b_entry;
                if (!res.unconstrained.branch_lengths.empty() && node.id >= 0 && static_cast<size_t>(node.id) < res.unconstrained.branch_lengths.size()) {
                    b_entry["unconstrained"] = res.unconstrained.branch_lengths[node.id];
                }
                if (!res.constrained.branch_lengths.empty() && node.id >= 0 && static_cast<size_t>(node.id) < res.constrained.branch_lengths.size()) {
                    b_entry["constrained"] = res.constrained.branch_lengths[node.id];
                }
                auto it_post = res.branch_posterior_probabilities.find(node.name);
                if (it_post != res.branch_posterior_probabilities.end()) {
                    b_entry["Posterior probability of positive selection"] = {it_post->second};
                }
                auto it_ebf = res.branch_empirical_bayes_factors.find(node.name);
                if (it_ebf != res.branch_empirical_bayes_factors.end()) {
                    b_entry["Empirical Bayes Factor for positive selection"] = {it_ebf->second};
                }
                branch_attr[node.name] = b_entry;
                int part = (static_cast<size_t>(node.id) < branch_partition.size()) ? branch_partition[node.id] : (is_test_branch[node.id] ? 0 : 1);
                tested_map[node.name] = (part == 0) ? "test" : ((part == 1) ? "background" : "unclassified");
            }
        }
        j["branch attributes"] = {{"0", branch_attr}};
        j["tested"] = {{"0", tested_map}};

        j["Evidence Ratios"] = {
            {"constrained", {res.constrained_evidence_ratios.empty() ? res.evidence_ratios : res.constrained_evidence_ratios}},
            {"optimized null", {res.evidence_ratios}}
        };
        if (!res.branch_er_2h.empty()) {
            nlohmann::json er_2h = nlohmann::json::object();
            for (const auto& [bname, vec] : res.branch_er_2h) {
                er_2h[bname] = vec;
            }
            j["Evidence Ratios"]["2H"] = er_2h;
        }
        if (!res.branch_er_3h.empty()) {
            nlohmann::json er_3h = nlohmann::json::object();
            for (const auto& [bname, vec] : res.branch_er_3h) {
                er_3h[bname] = vec;
            }
            j["Evidence Ratios"]["3H"] = er_3h;
        }
        if (!res.branch_er_23h.empty()) {
            nlohmann::json er_23h = nlohmann::json::object();
            for (const auto& [bname, vec] : res.branch_er_23h) {
                er_23h[bname] = vec;
            }
            j["Evidence Ratios"]["2H+3H"] = er_23h;
        }
        if (res.settings.error_sink) {
            j["Evidence Ratios"]["error_sink_max_bf"] = {res.site_error_bfs};
        }

        // Site Log Likelihood
        std::vector<Scalar> unc_site_ll(aln.num_codons), con_site_ll(aln.num_codons);
        for (size_t s = 0; s < aln.num_codons; ++s) {
            size_t p = aln.site_to_pattern[s];
            unc_site_ll[s] = res.unconstrained.site_log_likelihoods[p];
            con_site_ll[s] = res.constrained.site_log_likelihoods[p];
        }
        j["Site Log Likelihood"] = {
            {"unconstrained", {unc_site_ll}},
            {"constrained", {con_site_ll}}
        };
        if (res.has_background_test &&
            res.constrained_background.site_log_likelihoods.size() == aln.patterns.size() &&
            res.shared_distributions.site_log_likelihoods.size() == aln.patterns.size()) {
            std::vector<Scalar> con_bg_ll(aln.num_codons), shared_ll(aln.num_codons);
            for (size_t s = 0; s < aln.num_codons; ++s) {
                size_t p = aln.site_to_pattern[s];
                con_bg_ll[s] = res.constrained_background.site_log_likelihoods[p];
                shared_ll[s] = res.shared_distributions.site_log_likelihoods[p];
            }
            j["Site Log Likelihood"]["constrained_background"] = {con_bg_ll};
            j["Site Log Likelihood"]["shared_distributions"] = {shared_ll};
        }

        j["runtime"] = res.runtime_seconds;
        return j;
    }

nlohmann::json BUSTEDAnalyzer::to_modern_json(const BUSTEDResult& res, const Provenance& prov) const {
        auto json_settings = nlohmann::json{
            {"rates", res.settings.num_rate_classes},
            {"srv", res.settings.srv},
            {"syn_rates", res.settings.num_syn_rate_classes},
            {"auto_k", res.settings.auto_select_k},
            {"multiple_hits", res.settings.multiple_hits},
            {"error_sink", res.settings.error_sink},
            {"phenotype_mode", res.settings.phenotype_mode},
            {"relax_mode", res.has_relax_test},
            {"refine_branch_lengths", res.settings.refine_branch_lengths},
            {"pvalue_threshold", res.settings.p_value_threshold},
            {"optimal_k", res.optimal_k},
            {"test_branches", res.settings.test_branches}
        };
        if (!res.settings.comparison_branches.empty()) {
            json_settings["comparison_branches"] = res.settings.comparison_branches;
        }
        if (res.settings.error_sink) {
            json_settings["error_sink_max_weight"] = res.settings.error_sink_max_weight;
            json_settings["error_sink_bf_threshold"] = res.settings.error_sink_bf_threshold;
        }
        if (!res.settings.mss_type.empty()) {
            json_settings["mss_type"] = res.settings.mss_type;
            json_settings["mss_neutral"] = res.settings.mss_neutral;
            json_settings["mss_normalization"] = res.settings.mss_normalization;
        }

        AnalysisMetadata meta{
            "busted",
            "Branch-Site Unrestricted Statistical Test for Episodic Diversification",
            "3.0.0",
            "gene_selection",
            "BUSTED uses a random effects likelihood approach to test whether a gene has experienced positive selection at any site on specified (or all) branches of a phylogeny.",
            {citation("Murrell B, Weaver S, Smith MD, Wertheim JO, Murrell S, Aylward A, Eren K, Pollner T, Martin DP, Frost SDW, Kosakovsky Pond SL (2015). Gene-Wide Identification of Episodic Selection. Mol Biol Evol 32(5): 1365-1371.",
                      "10.1093/molbev/msv035", "25701167")},
            std::move(json_settings)
        };

        nlohmann::json j = new_analysis_json(meta, prov);
        j["status"] = res.is_failed ? "failed" : "completed";
        if (res.is_failed) {
            j["error"] = res.error_message;
        }
        // Preserve measured runtime when the caller supplied no real provenance:
        // the auto block reports the actual wall clock of the analysis.
        if (prov.invocation.cli_command.empty() && prov.inputs.empty() &&
            j["provenance"].contains("execution")) {
            j["provenance"]["execution"]["wall_time_seconds"] = res.runtime_seconds;
        }
        // BUSTED-PH shares the BUSTED schema, so meta.id stays "busted" and the
        // analysis identity is post-fixed here rather than derived from meta.
        if (res.settings.phenotype_mode) {
            j["analysis"]["id"] = "busted-ph";
            j["analysis"]["name"] =
                "BUSTED-Phenotype: Branch-Site Unrestricted Statistical Test for Episodic Diversification Associated with a Phenotype";
            j["analysis"]["description"] =
                "BUSTED-PH tests whether episodic positive selection is specifically associated with a foreground phenotype/trait using a Union-Intersection Test on foreground selection and regime difference, followed by an AIC-c lack-of-evidence check on background selection.";
        }
        j["dataset"] = dataset_block(aln);
        j["phylogeny"] = phylogeny_block(
            tree, "",
            [&res](const TreeNode& node) {
                std::map<std::string, Scalar> lengths;
                if (node.id >= 0 &&
                    static_cast<size_t>(node.id) < res.unconstrained.branch_lengths.size()) {
                    lengths["unconstrained"] = res.unconstrained.branch_lengths[node.id];
                }
                if (node.id >= 0 &&
                    static_cast<size_t>(node.id) < res.constrained.branch_lengths.size()) {
                    lengths["constrained"] = res.constrained.branch_lengths[node.id];
                }
                return lengths;
            },
            [this](const TreeNode& node) {
                nlohmann::json extras;
                if (node.id != tree.root_id) {
                    int part = (node.id >= 0 && static_cast<size_t>(node.id) < branch_partition.size())
                        ? branch_partition[node.id]
                        : (static_cast<size_t>(node.id) < is_test_branch.size() && is_test_branch[node.id] ? 0 : 1);
                    extras["tested"] = (part == 0);
                    extras["partition"] = (part == 0) ? "test" : ((part == 1) ? "background" : "unclassified");
                }
                return extras;
            });

        // Model fits
        auto format_fit_modern = [&](const BUSTEDFit& fit, const std::string& desc) {
            nlohmann::json f;
            f["description"] = desc;
            f["log_likelihood"] = finite_or_null(fit.log_likelihood);
            f["aicc"] = finite_or_null(fit.aicc);
            f["tree_scale"] = finite_or_null(fit.tree_scale);
            f["rate_classes"] = static_cast<int>(fit.num_rate_classes);

            const Scalar bio_scale =
                res.settings.error_sink ? (1.0 - fit.error_sink_weight) : 1.0;
            const auto& d = fit.test_distribution;
            if (d.omegas.size() != d.weights.size() ||
                (!d.annotations.empty() && d.annotations.size() != d.omegas.size())) {
                throw std::runtime_error(
                    "busted: rate distribution length mismatch (omegas=" +
                    std::to_string(d.omegas.size()) + ", weights=" +
                    std::to_string(d.weights.size()) + ", annotations=" +
                    std::to_string(d.annotations.size()) + ")");
            }
            nlohmann::json test_dist = nlohmann::json::array();
            for (size_t k = 0; k < d.omegas.size(); ++k) {
                nlohmann::json cls = {
                    {"class", static_cast<int>(k + 1)},
                    {"omega", finite_or_null(d.omegas[k])},
                    {"proportion", finite_or_null(bio_scale * d.weights[k])}
                };
                if (k < d.annotations.size() && !d.annotations[k].empty()) {
                    cls["annotation"] = d.annotations[k];
                }
                test_dist.push_back(cls);
            }
            f["rate_distributions"]["test"] = test_dist;

            if (fit.has_background_distribution && !fit.background_distribution.omegas.empty()) {
                Scalar bio_scale_bg = res.settings.error_sink ? (1.0 - fit.background_distribution.error_sink_weight) : 1.0;
                nlohmann::json bg_dist = nlohmann::json::array();
                for (size_t k = 0; k < fit.background_distribution.omegas.size(); ++k) {
                    bg_dist.push_back({
                        {"class", static_cast<int>(k + 1)},
                        {"omega", fit.background_distribution.omegas[k]},
                        {"proportion", bio_scale_bg * fit.background_distribution.weights[k]}
                    });
                }
                f["rate_distributions"]["background"] = bg_dist;
            }

            if (fit.has_unclassified_distribution && !fit.unclassified_distribution.omegas.empty()) {
                nlohmann::json unclass_dist = nlohmann::json::array();
                for (size_t k = 0; k < fit.unclassified_distribution.omegas.size(); ++k) {
                    unclass_dist.push_back({
                        {"class", static_cast<int>(k + 1)},
                        {"omega", finite_or_null(fit.unclassified_distribution.omegas[k])},
                        {"proportion", finite_or_null(fit.unclassified_distribution.weights[k])}
                    });
                }
                f["rate_distributions"]["unclassified"] = unclass_dist;
            }

            if (res.settings.error_sink) {
                f["error_sink"] = {
                    {"type", "saturation_replacement"},
                    {"kernel", "P_sink(i,j) = pi_j / (1 - pi_i) for i != j, 0 for i == j"},
                    {"proportion", finite_or_null(fit.error_sink_weight)},
                    {"max_weight_bound", res.settings.error_sink_max_weight}
                };
            }

            if (d.syn_rates.size() > 1 && d.syn_rates.size() == d.syn_weights.size()) {
                nlohmann::json syn_dist = nlohmann::json::array();
                for (size_t m = 0; m < d.syn_rates.size(); ++m) {
                    syn_dist.push_back({
                        {"class", static_cast<int>(m + 1)},
                        {"rate", finite_or_null(d.syn_rates[m])},
                        {"proportion", finite_or_null(d.syn_weights[m])}
                    });
                }
                f["rate_distributions"]["synonymous"] = syn_dist;
            }

            if (res.settings.multiple_hits != "None") {
                f["multiple_hits"]["double_hit_rate"] = finite_or_null(fit.delta);
                f["multiple_hits"]["double_hit_fraction"] = finite_or_null(fit.frac_delta);
                if (res.settings.multiple_hits == "Double+Triple") {
                    f["multiple_hits"]["triple_hit_rate"] = finite_or_null(fit.psi);
                    f["multiple_hits"]["triple_hit_fraction"] = finite_or_null(fit.frac_psi);
                }
            }
            if (!res.mss_rates.empty() && fit.mss_class_rates.size() == res.mss_rates.size()) {
                nlohmann::json mss_obj = nlohmann::json::object();
                for (size_t c = 0; c < res.mss_rates.size(); ++c) {
                    mss_obj[res.mss_rates[c].class_name] = finite_or_null(fit.mss_class_rates[c]);
                }
                f["mss_rates"] = mss_obj;
            }
            return f;
        };

        j["model_fits"]["unconstrained"] = format_fit_modern(res.unconstrained, "Unconstrained BUSTED mixture model allowing omega_k > 1");
        j["model_fits"]["constrained"] = format_fit_modern(res.constrained, "Constrained null BUSTED model fixing omega_k <= 1 on test branches");
        if (res.has_background_test) {
            j["model_fits"]["constrained_background"] = format_fit_modern(
                res.constrained_background, "Constrained background BUSTED-PH model fixing omega_k <= 1 on background/comparison branches"
            );
        }
        if (res.has_background_test || res.has_relax_test) {
            j["model_fits"]["shared_distributions"] = format_fit_modern(
                res.shared_distributions, "Shared-distribution BUSTED-PH / RELAX Null (k=1) model forcing identical omega distributions on test and background branches"
            );
        }
        if (res.has_relax_test) {
            auto relax_json = format_fit_modern(
                res.relax_alternative, "Coupled-k RELAX alternative model (omega_FG,c = (omega_BG,c)^k with shared category weights)"
            );
            relax_json["relaxation_parameter_k"] = res.relax_k;
            j["model_fits"]["relax_alternative"] = relax_json;
        }

        if (has_gtr_fit) {
            j["model_fits"]["nucleotide_gtr"] = nucleotide_gtr_fit(
                aln, gtr_rates, gtr_log_l, gtr_aicc,
                5 + static_cast<int>(tree.nodes.size()) - 1);
        }

        j["model_fits"]["codon_mg94"]["description"] = "Muse-Gaut 1994 x GTR codon model with shared global omega";
        j["model_fits"]["codon_mg94"]["log_likelihood"] = finite_or_null(mg94_log_l);
        j["model_fits"]["codon_mg94"]["global_parameters"] = {
            {"omega", finite_or_null(mg94_omega)}
        };

        // Statistical tests
        j["statistical_tests"]["gene_level"] = {
            {"method", "Likelihood Ratio Test (LRT)"},
            {"test_statistic", "LRT"},
            {"statistic_value", finite_or_null(res.lrt)},
            {"distribution", "0.5 * chi^2_0 + 0.5 * chi^2_2"},
            {"degrees_of_freedom", 2},
            {"null_hypothesis", "omega_k <= 1 for all k on test branches (Constrained model)"},
            {"alternative_hypothesis", "omega_k > 1 for at least one k on test branches (Unconstrained model)"},
            {"p_value", finite_or_null(res.p_value)},
            {"decision", res.is_failed ? ("Analysis failed: " + res.error_message) : ((res.p_value <= res.settings.p_value_threshold) ? "Episodic diversifying selection detected" : "No evidence of episodic diversifying selection")},
            {"optimal_k", static_cast<int>(res.optimal_k)}
        };
        if (res.has_background_test) {
            j["statistical_tests"]["background_selection"] = {
                {"method", "AIC-c Lack-of-Evidence Screen & LRT"},
                {"test_statistic", "LRT"},
                {"statistic_value", res.background_lrt},
                {"delta_aicc", res.background_delta_aicc},
                {"neutral_supported_by_aicc", res.background_neutral_supported},
                {"distribution", "0.5 * chi^2_0 + 0.5 * chi^2_2"},
                {"degrees_of_freedom", 2},
                {"null_hypothesis", "omega_k <= 1 for all k on background/comparison branches"},
                {"alternative_hypothesis", "omega_k > 1 for at least one k on background/comparison branches"},
                {"p_value", res.background_p_value},
                {"decision", res.background_neutral_supported
                    ? "Neutral/constrained model for background branches is sufficiently supported by AIC-c"
                    : "Neutral/constrained model for background branches is not supported by AIC-c (evidence of background selection)"}
            };
            j["statistical_tests"]["shared_distributions"] = {
                {"method", "Likelihood Ratio Test (LRT)"},
                {"test_statistic", "LRT"},
                {"statistic_value", res.shared_lrt},
                {"distribution", "chi^2_" + std::to_string(res.shared_df)},
                {"degrees_of_freedom", res.shared_df},
                {"null_hypothesis", "Identical omega rate distributions on test and background branches (D_T == D_B)"},
                {"alternative_hypothesis", "Distinct omega rate distributions on test and background branches (D_T != D_B)"},
                {"p_value", res.shared_p_value},
                {"decision", (res.shared_p_value <= res.settings.p_value_threshold) ? "Selective regimes differ significantly between test and background" : "No significant difference between test and background selective regimes"}
            };
            j["statistical_tests"]["union_intersection"] = {
                {"method", "Union-Intersection Test (UIT)"},
                {"null_hypothesis", "(omega_K^(FG) <= 1) U (D_FG == D_BG)"},
                {"alternative_hypothesis", "(omega_K^(FG) > 1) n (D_FG != D_BG)"},
                {"p_value_fg", res.p_value},
                {"p_value_diff", res.shared_p_value},
                {"p_value", res.union_intersection_p_value},
                {"rejected", (res.union_intersection_p_value <= res.settings.p_value_threshold)}
            };
            j["statistical_tests"]["phenotype_association"] = {
                {"associated", res.phenotype_associated},
                {"verdict", res.phenotype_verdict},
                {"union_intersection_p_value", res.union_intersection_p_value},
                {"background_delta_aicc", res.background_delta_aicc},
                {"background_neutral_supported", res.background_neutral_supported},
                {"ranking_score", res.ranking_score},
                {"p_value_test", res.p_value},
                {"p_value_background", res.background_p_value},
                {"p_value_shared", res.shared_p_value},
                {"threshold", res.settings.p_value_threshold}
            };
        }
        if (res.has_relax_test) {
            j["statistical_tests"]["relaxation"] = {
                {"method", "RELAX Coupled-k Likelihood Ratio Test"},
                {"test_statistic", "LRT"},
                {"statistic_value", res.relax_lrt},
                {"distribution", "chi^2_1"},
                {"degrees_of_freedom", 1},
                {"relaxation_parameter_k", res.relax_k},
                {"direction", res.relax_direction},
                {"null_hypothesis", "k = 1 (identical omega distributions on test and background branches)"},
                {"alternative_hypothesis", "k != 1 (omega_FG,c = (omega_BG,c)^k with shared category weights)"},
                {"p_value", res.relax_p_value},
                {"decision", (res.relax_p_value <= res.settings.p_value_threshold)
                    ? ((res.relax_k > 1.0) ? "Significant intensification of selection (k > 1)" : "Significant relaxation of selection (k < 1)")
                    : "No significant relaxation or intensification of selection"}
            };
        }

        if (res.auto_diagnostics.enabled) {
            j["auto_model_selection"] = {
                {"error_sink_gate", {
                    {"max_surprisal_bits", res.auto_diagnostics.max_sink_surprisal_bits},
                    {"max_baseline_bf", res.auto_diagnostics.max_sink_bf},
                    {"candidate_sites_gt_8bits", static_cast<int>(res.auto_diagnostics.sink_candidate_sites)},
                    {"adjacent_clusters", static_cast<int>(res.auto_diagnostics.sink_adjacent_clusters)},
                    {"delta_log_l", res.auto_diagnostics.sink_delta_ll},
                    {"activated", res.auto_diagnostics.activated_error_sink}
                }},
                {"multiple_hits_gate", {
                    {"quarantined_double_hit_score", res.auto_diagnostics.double_hit_score},
                    {"quarantined_triple_hit_score", res.auto_diagnostics.triple_hit_score},
                    {"double_hit_delta_log_l", res.auto_diagnostics.double_hit_delta_ll},
                    {"triple_hit_delta_log_l", res.auto_diagnostics.triple_hit_delta_ll},
                    {"selected", res.auto_diagnostics.selected_multiple_hits}
                }},
                {"srv_gate", {
                    {"overdispersion_z", res.auto_diagnostics.srv_overdispersion_z},
                    {"delta_log_l", res.auto_diagnostics.srv_delta_ll},
                    {"activated", res.auto_diagnostics.activated_srv},
                    {"selected_syn_rates", static_cast<int>(res.auto_diagnostics.selected_syn_rates)}
                }},
                {"omega_classes_gate", {
                    {"selected_k", static_cast<int>(res.auto_diagnostics.selected_k)},
                    {"rationale", res.auto_diagnostics.k_rationale}
                }}
            };
        }

        if (res.settings.error_sink) {
            nlohmann::json flagged = nlohmann::json::array();
            for (const auto& ev : res.flagged_error_events) {
                flagged.push_back({
                    {"site", static_cast<int>(ev.site)},
                    {"branch", ev.branch_name},
                    {"posterior", finite_or_null(ev.posterior)},
                    {"bayes_factor", finite_or_null(ev.bayes_factor)}
                });
            }
            j["error_sink_diagnostics"] = {
                {"unconstrained_epsilon", finite_or_null(res.unconstrained.error_sink_weight)},
                {"constrained_epsilon", finite_or_null(res.constrained.error_sink_weight)},
                {"bf_threshold", res.settings.error_sink_bf_threshold},
                {"flagged_branch_sites_count", static_cast<int>(res.flagged_error_events.size())},
                {"flagged_branch_sites", flagged}
            };
        }

        if (!res.mss_rates.empty()) {
            nlohmann::json mss_arr = nlohmann::json::array();
            for (const auto& r : res.mss_rates) {
                mss_arr.push_back({
                    {"name", r.class_name},
                    {"rate", finite_or_null(r.rate)},
                    {"pair_count", static_cast<int>(r.num_pairs)},
                    {"is_reference", r.is_reference}
                });
            }
            j["mss_rates"] = mss_arr;
        }

        // Site results (columnar; ColumnarResults enforces equal-length)
        ColumnarResults table;
        std::vector<int> col_site;
        std::vector<double> col_er, col_con_er, col_unc_ll, col_con_ll, col_err_bf;
        col_site.reserve(aln.num_codons);
        col_er.reserve(aln.num_codons);
        col_con_er.reserve(aln.num_codons);
        col_unc_ll.reserve(aln.num_codons);
        col_con_ll.reserve(aln.num_codons);
        if (res.settings.error_sink) col_err_bf.reserve(aln.num_codons);

        for (size_t s = 0; s < aln.num_codons; ++s) {
            col_site.push_back(static_cast<int>(s + 1));
            col_er.push_back(res.evidence_ratios[s]);
            col_con_er.push_back(s < res.constrained_evidence_ratios.size()
                ? res.constrained_evidence_ratios[s]
                : res.evidence_ratios[s]);
            const size_t p = aln.site_to_pattern[s];
            col_unc_ll.push_back(res.unconstrained.site_log_likelihoods[p]);
            col_con_ll.push_back(res.constrained.site_log_likelihoods[p]);
            if (res.settings.error_sink) {
                col_err_bf.push_back(s < res.site_error_bfs.size()
                    ? res.site_error_bfs[s]
                    : std::numeric_limits<double>::quiet_NaN());
            }
        }

        table.add_column("site", {"integer", "Alignment codon site index (1-based)", "1-based codon position"}, col_site);
        table.add_column("evidence_ratio", {"number", "Per-site Evidence Ratio favoring positive selection against optimized null model"}, col_er);
        table.add_column("constrained_evidence_ratio", {"number", "Per-site Evidence Ratio favoring positive selection against constrained (omega_K=1) alternative model"}, col_con_er);
        table.add_column("unconstrained_log_l", {"number", "Site log-likelihood under unconstrained model"}, col_unc_ll);
        table.add_column("constrained_log_l", {"number", "Site log-likelihood under constrained null model"}, col_con_ll);
        if (res.settings.error_sink) {
            table.add_column("error_sink_max_bf", {"number", "Maximum branch-site Empirical Bayes Factor for the saturation error sink across branches at this site"}, col_err_bf);
        }
        j["site_results"] = table.to_json();

        if (!res.branch_posterior_probabilities.empty()) {
            nlohmann::json b_post = nlohmann::json::object();
            for (const auto& [bname, vec] : res.branch_posterior_probabilities) {
                b_post[bname]["posterior_probability"] = vec;
                auto it_ebf = res.branch_empirical_bayes_factors.find(bname);
                if (it_ebf != res.branch_empirical_bayes_factors.end()) {
                    b_post[bname]["empirical_bayes_factor"] = it_ebf->second;
                }
            }
            j["branch_site_posteriors"] = std::move(b_post);
        }
        if (!res.branch_er_2h.empty() || !res.branch_er_3h.empty() || !res.branch_er_23h.empty()) {
            nlohmann::json mh_er = nlohmann::json::object();
            if (!res.branch_er_2h.empty()) mh_er["2H"] = res.branch_er_2h;
            if (!res.branch_er_3h.empty()) mh_er["3H"] = res.branch_er_3h;
            if (!res.branch_er_23h.empty()) mh_er["2H+3H"] = res.branch_er_23h;
            j["multiple_hit_branch_evidence_ratios"] = std::move(mh_er);
        }

        return j;
    }

} // namespace hyphy::analyses
