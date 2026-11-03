// Out-of-line JSON emitter definitions for relax.
// Moved from include/hyphy/analyses/relax.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/relax.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json RELAXResult::to_legacy_json(const Tree& tree, const Alignment& aln) const {
        nlohmann::json j;

        // 1. Analysis metadata
        j["analysis"] = {
            {"info", "RELAX (a random effects test of selection relaxation) uses a random effects branch-site model framework to test whether a set of 'Test' branches evolves under relaxed selection relative to a set of 'Reference' branches (R), as measured by the relaxation parameter (K)."},
            {"citation", "RELAX: Detecting Relaxed Selection in a Phylogenetic Framework (2015). Mol Biol Evol 32 (3): 820-832"},
            {"authors", "Sergei L Kosakovsky Pond, Ben Murrell, Steven Weaver and Temple iGEM / UCSD viral evolution group"},
            {"contact", "spond@temple.edu"},
            {"version", "3.0.0"}
        };

        // 2. Input
        j["input"] = {
            {"number of sequences", aln.num_taxa},
            {"number of sites", aln.num_codons},
            {"partition count", 1},
            {"trees", {{"0", aln.embedded_tree_newick}}}
        };

        // 3. Fits
        nlohmann::json fits;

        // GTR fit
        nlohmann::json gtr_freqs = nlohmann::json::array();
        for (int i = 0; i < 4; ++i) {
            gtr_freqs.push_back({aln.nuc_frequencies(i)});
        }
        fits["Nucleotide GTR"] = {
            {"AIC-c", gtr_aicc},
            {"Log Likelihood", gtr_log_likelihood},
            {"display order", 0},
            {"estimated parameters", gtr_parameters},
            {"Equilibrium frequencies", gtr_freqs},
            {"Rate Distributions", {
                {"Substitution rate from nucleotide A to nucleotide C", gtr_rates.theta_AC},
                {"Substitution rate from nucleotide A to nucleotide G", 1.0},
                {"Substitution rate from nucleotide A to nucleotide T", gtr_rates.theta_AT},
                {"Substitution rate from nucleotide C to nucleotide G", gtr_rates.theta_CG},
                {"Substitution rate from nucleotide C to nucleotide T", gtr_rates.theta_CT},
                {"Substitution rate from nucleotide G to nucleotide T", gtr_rates.theta_GT}
            }}
        };

        // MG94 separate rates fit
        nlohmann::json mg_freqs = nlohmann::json::array();
        for (int i = 0; i < aln.codon_frequencies_f3x4.size(); ++i) {
            mg_freqs.push_back({aln.codon_frequencies_f3x4(i)});
        }
        fits["MG94xREV with separate rates for branch sets"] = {
            {"AIC-c", mg94_aicc},
            {"Log Likelihood", mg94_log_likelihood},
            {"display order", 1},
            {"estimated parameters", mg94_parameters},
            {"Equilibrium frequencies", mg_freqs},
            {"Rate Distributions", {
                {"non-synonymous/synonymous rate ratio for *Reference*", {{mg94_omega_R, 1.0}}},
                {"non-synonymous/synonymous rate ratio for *Test*", {{mg94_omega_T, 1.0}}}
            }}
        };

        auto format_distro = [](const RELAXRateDistribution& ref_d, const RELAXRateDistribution& test_d, const RELAXRateDistribution& uncl_d) {
            nlohmann::json res;
            nlohmann::json r_obj, t_obj;
            for (size_t k = 0; k < ref_d.omegas.size(); ++k) {
                r_obj[std::to_string(k)] = {
                    {"omega", ref_d.omegas[k]},
                    {"proportion", ref_d.weights[k]}
                };
                t_obj[std::to_string(k)] = {
                    {"omega", test_d.omegas[k]},
                    {"proportion", test_d.weights[k]}
                };
            }
            res["Reference"] = r_obj;
            res["Test"] = t_obj;
            if (!uncl_d.omegas.empty() && uncl_d.omegas.size() == uncl_d.weights.size()) {
                nlohmann::json u_obj;
                for (size_t k = 0; k < uncl_d.omegas.size(); ++k) {
                    u_obj[std::to_string(k)] = {
                        {"omega", uncl_d.omegas[k]},
                        {"proportion", uncl_d.weights[k]}
                    };
                }
                res["Unclassified"] = u_obj;
            }
            return res;
        };

        // Alternative fit
        fits["RELAX alternative"] = {
            {"AIC-c", alternative_fit.aicc},
            {"Log Likelihood", alternative_fit.log_likelihood},
            {"display order", 2},
            {"estimated parameters", alternative_fit.parameters},
            {"realized_cycles", alternative_fit.realized_cycles},
            {"converged", alternative_fit.converged},
            {"Rate Distributions", format_distro(alternative_fit.reference_distribution, alternative_fit.test_distribution, alternative_fit.unclassified_distribution)}
        };

        // Null fit
        fits["RELAX null"] = {
            {"AIC-c", null_fit.aicc},
            {"Log Likelihood", null_fit.log_likelihood},
            {"display order", 3},
            {"estimated parameters", null_fit.parameters},
            {"realized_cycles", null_fit.realized_cycles},
            {"converged", null_fit.converged},
            {"Rate Distributions", format_distro(null_fit.reference_distribution, null_fit.test_distribution, null_fit.unclassified_distribution)}
        };

        j["fits"] = fits;

        // 4. Test Results
        j["test results"] = {
            {"LRT", lrt},
            {"p-value", p_value},
            {"relaxation or intensification parameter", k}
        };
        if (is_failed) {
            j["error"] = error_message;
            j["status"] = "failed";
        }

        // 5. Tested branch classification
        nlohmann::json tested_map;
        for (const auto& [name, cls] : branch_class) {
            tested_map[name] = cls;
        }
        j["tested"] = {{"0", tested_map}};

        // 6. Branch attributes
        nlohmann::json branch_attr;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                nlohmann::json b_entry;
                b_entry["original name"] = node.name;
                if (gtr_branch_lengths.count(node.name)) {
                    b_entry["Nucleotide GTR"] = gtr_branch_lengths.at(node.name);
                }
                if (mg94_branch_lengths.count(node.name)) {
                    b_entry["MG94xREV with separate rates for branch sets"] = mg94_branch_lengths.at(node.name);
                }
                if (!alternative_fit.branch_lengths.empty() && node.id >= 0 && static_cast<size_t>(node.id) < alternative_fit.branch_lengths.size()) {
                    b_entry["RELAX alternative"] = alternative_fit.branch_lengths[node.id];
                }
                if (!null_fit.branch_lengths.empty() && node.id >= 0 && static_cast<size_t>(node.id) < null_fit.branch_lengths.size()) {
                    b_entry["RELAX null"] = null_fit.branch_lengths[node.id];
                }
                branch_attr[node.name] = b_entry;
            }
        }
        j["branch attributes"] = {{"0", branch_attr}};

        // 7. Site Log Likelihoods
        if (!alternative_fit.site_log_likelihoods.empty()) {
            const auto& constrained_vals = !branch_switched_evidence_ratios.empty()
                ? branch_switched_evidence_ratios
                : null_fit.site_log_likelihoods;
            if (!constrained_vals.empty()) {
                j["Site Log Likelihood"] = {
                    {"0", {
                        {"unconstrained", {alternative_fit.site_log_likelihoods}},
                        {"constrained", {constrained_vals}}
                    }}
                };
            }
        }

        j["runtime"] = "3.0.0";
        return j;
    }

nlohmann::json RELAXResult::to_modern_json(const Tree& tree, const Alignment& aln, const Provenance& prov) const {
        AnalysisMetadata meta{
            "relax",
            "Relaxed Selection Analysis (RELAX)",
            "3.0.0",
            "comparative_selection",
            "RELAX uses a random effects branch-site model framework to test whether a set of 'Test' branches evolves under relaxed selection relative to a set of 'Reference' branches, as measured by the selection intensity parameter K.",
            {citation("Wertheim JO, Murrell B, Smith MD, Kosakovsky Pond SL, Scheffler K (2015). RELAX: Detecting Relaxed Selection in a Phylogenetic Framework. Mol Biol Evol 32(3): 820-832.",
                      "10.1093/molbev/msu400", "25540451")},
            {
                {"test_branches_count", test_branches.size()},
                {"reference_branches_count", reference_branches.size()},
                {"normalization",
                 (settings.normalization == RELAXNormalization::SharedMixture) ? "shared_mixture" : "per_class"},
                {"mixture_level",
                 (settings.mixture_level == RELAXMixtureLevel::BranchSite) ? "branch_site" : "site"},
                {"refine_branch_lengths", settings.refine_branch_lengths},
                {"max_joint_cycles", settings.max_joint_cycles}
            }
        };

        nlohmann::json j = new_analysis_json(meta, prov);
        if (prov.invocation.cli_command.empty() && prov.inputs.empty() &&
            j["provenance"].contains("execution")) {
            j["provenance"]["execution"]["wall_time_seconds"] = runtime_seconds;
        }
        j["dataset"] = dataset_block(aln);
        j["phylogeny"] = phylogeny_block(
            tree, "",
            [this](const TreeNode& node) {
                std::map<std::string, Scalar> lengths;
                if (node.id >= 0 &&
                    static_cast<size_t>(node.id) < alternative_fit.branch_lengths.size()) {
                    lengths["alternative"] = alternative_fit.branch_lengths[node.id];
                }
                if (node.id >= 0 &&
                    static_cast<size_t>(node.id) < null_fit.branch_lengths.size()) {
                    lengths["null"] = null_fit.branch_lengths[node.id];
                }
                return lengths;
            },
            [this](const TreeNode& node) -> nlohmann::json {
                auto it = branch_class.find(node.name);
                if (it == branch_class.end()) return nlohmann::json::object();
                return {{"branch_class", it->second}};
            });

        // Model fits
        auto format_distro_modern = [](const RELAXRateDistribution& ref_d,
                                       const RELAXRateDistribution& test_d,
                                       const RELAXRateDistribution& uncl_d) {
            if (ref_d.omegas.size() != ref_d.weights.size() ||
                test_d.omegas.size() != test_d.weights.size() ||
                ref_d.omegas.size() != test_d.omegas.size()) {
                throw std::runtime_error(
                    "relax: rate distribution length mismatch (ref omegas=" +
                    std::to_string(ref_d.omegas.size()) + ", ref weights=" +
                    std::to_string(ref_d.weights.size()) + ", test omegas=" +
                    std::to_string(test_d.omegas.size()) + ", test weights=" +
                    std::to_string(test_d.weights.size()) + ")");
            }
            nlohmann::json d;
            nlohmann::json ref_arr = nlohmann::json::array();
            nlohmann::json test_arr = nlohmann::json::array();
            for (size_t k = 0; k < ref_d.omegas.size(); ++k) {
                ref_arr.push_back({
                    {"class", static_cast<int>(k + 1)},
                    {"omega", finite_or_null(ref_d.omegas[k])},
                    {"proportion", finite_or_null(ref_d.weights[k])}
                });
                test_arr.push_back({
                    {"class", static_cast<int>(k + 1)},
                    {"omega", finite_or_null(test_d.omegas[k])},
                    {"proportion", finite_or_null(test_d.weights[k])}
                });
            }
            d["reference"] = ref_arr;
            d["test"] = test_arr;
            if (!uncl_d.omegas.empty() && uncl_d.omegas.size() == uncl_d.weights.size()) {
                nlohmann::json uncl_arr = nlohmann::json::array();
                for (size_t k = 0; k < uncl_d.omegas.size(); ++k) {
                    uncl_arr.push_back({
                        {"class", static_cast<int>(k + 1)},
                        {"omega", finite_or_null(uncl_d.omegas[k])},
                        {"proportion", finite_or_null(uncl_d.weights[k])}
                    });
                }
                d["unclassified"] = uncl_arr;
            }
            return d;
        };

        if (gtr_parameters > 0 || gtr_log_likelihood != 0.0) {
            j["model_fits"]["nucleotide_gtr"] = nucleotide_gtr_fit(
                aln, gtr_rates, gtr_log_likelihood, gtr_aicc,
                static_cast<int>(gtr_parameters));
        }

        nlohmann::json mg94_fit = fit_entry(
            "Muse-Gaut 1994 x GTR with separate rates for Reference and Test branches",
            mg94_log_likelihood, mg94_aicc, static_cast<int>(mg94_parameters));
        mg94_fit["rate_distributions"] = {
            {"omega_reference", finite_or_null(mg94_omega_R)},
            {"omega_test", finite_or_null(mg94_omega_T)}
        };
        j["model_fits"]["codon_mg94"] = mg94_fit;

        nlohmann::json alt = fit_entry(
            "RELAX alternative model with unconstrained relaxation parameter K",
            alternative_fit.log_likelihood, alternative_fit.aicc,
            static_cast<int>(alternative_fit.parameters));
        alt["k_parameter"] = finite_or_null(k);
        alt["realized_cycles"] = alternative_fit.realized_cycles;
        alt["converged"] = alternative_fit.converged;
        alt["rate_distributions"] = format_distro_modern(
            alternative_fit.reference_distribution, alternative_fit.test_distribution,
            alternative_fit.unclassified_distribution);
        j["model_fits"]["alternative"] = alt;

        nlohmann::json nul = fit_entry(
            "RELAX null model fixing K = 1",
            null_fit.log_likelihood, null_fit.aicc,
            static_cast<int>(null_fit.parameters));
        nul["k_parameter"] = 1.0;
        nul["realized_cycles"] = null_fit.realized_cycles;
        nul["converged"] = null_fit.converged;
        nul["rate_distributions"] = format_distro_modern(
            null_fit.reference_distribution, null_fit.test_distribution,
            null_fit.unclassified_distribution);
        j["model_fits"]["null"] = nul;

        // Statistical tests
        std::string decision_str;
        if (is_failed) {
            decision_str = "Optimization failed: " + error_message;
            j["status"] = "failed";
            j["error"] = error_message;
        } else if (p_value <= 0.05) {
            decision_str = (k < 1.0) ? "Relaxed selection detected on test branches" : "Intensified selection detected on test branches";
        } else {
            decision_str = "No evidence of relaxed or intensified selection";
        }
        j["statistical_tests"]["hypothesis_test"] = {
            {"method", "Likelihood Ratio Test (LRT)"},
            {"test_statistic", "LRT"},
            {"statistic_value", finite_or_null(lrt)},
            {"distribution", "Asymptotic Chi-squared"},
            {"degrees_of_freedom", 1},
            {"null_hypothesis", "K = 1 (Selection intensity is equal between Test and Reference)"},
            {"alternative_hypothesis", "K != 1 (Selection is relaxed K < 1 or intensified K > 1)"},
            {"p_value", finite_or_null(p_value)},
            {"k_parameter", finite_or_null(k)},
            {"decision", decision_str}
        };

        return j;
    }

} // namespace hyphy::analyses
