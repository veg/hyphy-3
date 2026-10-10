#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "nlohmann/json.hpp"
#include "hyphy/core/provenance.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <regex>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

struct RELAXRateDistribution {
    std::vector<Scalar> omegas;
    std::vector<Scalar> weights;
};

struct RELAXModelFit {
    Scalar log_likelihood = 0.0;
    Scalar aicc = 0.0;
    size_t parameters = 0;
    Scalar tree_scale = 1.0;
    Scalar k = 1.0;
    RELAXRateDistribution reference_distribution;
    RELAXRateDistribution test_distribution;
    std::vector<Scalar> branch_lengths;
    std::vector<Scalar> site_log_likelihoods;
};

struct RELAXSettings {
    Scalar p_value_threshold = 0.05;
    std::string test_branch_regex = "";
    std::vector<std::string> test_branch_names = {};
    std::string reference_branch_regex = "";
    std::vector<std::string> reference_branch_names = {};
    bool refine_branch_lengths = true;
    bool verbose = false;
};

struct RELAXResult {
    RELAXSettings settings;

    // Phase 1: Nucleotide GTR
    Scalar gtr_log_likelihood = 0.0;
    Scalar gtr_aicc = 0.0;
    size_t gtr_parameters = 0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    // Phase 2: MG94 with separate rates for branch sets
    Scalar mg94_log_likelihood = 0.0;
    Scalar mg94_aicc = 0.0;
    size_t mg94_parameters = 0;
    Scalar mg94_omega_R = 1.0;
    Scalar mg94_omega_T = 1.0;
    std::unordered_map<std::string, Scalar> mg94_branch_lengths;

    // Phase 3 & 4: RELAX fits
    RELAXModelFit alternative_fit;
    RELAXModelFit null_fit;

    // Phase 5: Hypothesis testing
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
    Scalar k = 1.0;
    bool is_relaxed = false;
    bool is_intensified = false;
    bool is_significant = false;

    // Branch classification
    std::vector<std::string> test_branches;
    std::vector<std::string> reference_branches;
    std::unordered_map<std::string, std::string> branch_class; // name -> "Test" / "Reference"

    double runtime_seconds = 0.0;

    nlohmann::json to_legacy_json(const Tree& tree, const Alignment& aln) const {
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

        auto format_distro = [](const RELAXRateDistribution& ref_d, const RELAXRateDistribution& test_d) {
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
            return res;
        };

        // Alternative fit
        fits["RELAX alternative"] = {
            {"AIC-c", alternative_fit.aicc},
            {"Log Likelihood", alternative_fit.log_likelihood},
            {"display order", 2},
            {"estimated parameters", alternative_fit.parameters},
            {"Rate Distributions", format_distro(alternative_fit.reference_distribution, alternative_fit.test_distribution)}
        };

        // Null fit
        fits["RELAX null"] = {
            {"AIC-c", null_fit.aicc},
            {"Log Likelihood", null_fit.log_likelihood},
            {"display order", 3},
            {"estimated parameters", null_fit.parameters},
            {"Rate Distributions", format_distro(null_fit.reference_distribution, null_fit.test_distribution)}
        };

        j["fits"] = fits;

        // 4. Test Results
        j["test results"] = {
            {"LRT", lrt},
            {"p-value", p_value},
            {"relaxation or intensification parameter", k}
        };

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
        if (!alternative_fit.site_log_likelihoods.empty() && !null_fit.site_log_likelihoods.empty()) {
            j["Site Log Likelihood"] = {
                {"0", {
                    {"unconstrained", {alternative_fit.site_log_likelihoods}},
                    {"constrained", {null_fit.site_log_likelihoods}}
                }}
            };
        }

        j["runtime"] = "3.0.0";
        return j;
    }

    nlohmann::json to_modern_json(const Tree& tree, const Alignment& aln, const Provenance& prov = {}) const {
        nlohmann::json j;

        j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/relax.v3.schema.json";
        j["schema_version"] = "3.0.0";

        // Analysis block
        j["analysis"]["id"] = "relax";
        j["analysis"]["name"] = "Relaxed Selection Analysis (RELAX)";
        j["analysis"]["version"] = "3.0.0";
        j["analysis"]["category"] = "comparative_selection";
        j["analysis"]["description"] = "RELAX uses a random effects branch-site model framework to test whether a set of 'Test' branches evolves under relaxed selection relative to a set of 'Reference' branches, as measured by the selection intensity parameter K.";
        j["analysis"]["citations"] = nlohmann::json::array({
            {
                {"citation", "Wertheim JO, Murrell B, Smith MD, Kosakovsky Pond SL, Scheffler K (2015). RELAX: Detecting Relaxed Selection in a Phylogenetic Framework. Mol Biol Evol 32(3): 820-832."},
                {"doi", "10.1093/molbev/msu400"},
                {"pmid", "25540451"}
            }
        });
        j["analysis"]["settings"]["test_branches_count"] = test_branches.size();
        j["analysis"]["settings"]["reference_branches_count"] = reference_branches.size();

        // Software block
        j["software"]["name"] = "hyphy";
        j["software"]["version"] = "3.0.0";
        j["software"]["git_commit"] = "6e22db3";
        j["software"]["build_type"] = "Release";
        j["software"]["compiler"] = Provenance::detect_compiler();

        // Provenance block
        if (!prov.invocation.cli_command.empty() || !prov.inputs.empty()) {
            j["provenance"] = prov.to_json();
        } else {
            Provenance auto_prov;
            auto_prov.invocation.cli_command = "hyphy3 relax";
            auto_prov.invocation.working_directory = Provenance::get_cwd();
            auto_prov.execution.start_time = Provenance::current_iso8601();
            auto_prov.execution.end_time = Provenance::current_iso8601();
            auto_prov.execution.wall_time_seconds = runtime_seconds;
            auto_prov.execution.cpu_threads = 1;
            auto_prov.execution.hostname = Provenance::get_hostname();
            auto_prov.execution.os = Provenance::detect_os();
            auto_prov.execution.compiler = Provenance::detect_compiler();
            j["provenance"] = auto_prov.to_json();
        }

        // Dataset block
        j["dataset"]["taxa_count"] = aln.num_taxa;
        j["dataset"]["codon_sites"] = aln.num_codons;
        j["dataset"]["nucleotide_sites"] = aln.num_codons * 3;
        j["dataset"]["unique_patterns"] = aln.patterns.size();
        j["dataset"]["taxa"] = aln.taxon_names;

        nlohmann::json gcode;
        gcode["id"] = aln.code ? aln.code->name : "Universal";
        gcode["name"] = aln.code ? aln.code->name : "Universal";
        gcode["sense_codons"] = aln.code ? static_cast<int>(aln.code->sense_codons.size()) : 61;
        if (aln.code) {
            gcode["stop_codons"] = aln.code->stop_codons;
        } else {
            gcode["stop_codons"] = {"TAA", "TAG", "TGA"};
        }
        j["dataset"]["genetic_code"] = gcode;

        j["dataset"]["partitions"] = nlohmann::json::array({
            {
                {"id", "default"},
                {"name", "Full Alignment"},
                {"span", nlohmann::json::array({nlohmann::json::array({1, aln.num_codons})})},
                {"sites_count", aln.num_codons},
                {"patterns_count", aln.patterns.size()}
            }
        });

        // Phylogeny block
        j["phylogeny"]["newick"] = tree.to_newick();
        nlohmann::json nodes_json;
        double sum_bl_alt = 0.0;
        double sum_bl_null = 0.0;

        for (const auto& node : tree.nodes) {
            nlohmann::json node_entry;
            node_entry["type"] = node.children.empty() ? "leaf" : "internal";

            auto it = branch_class.find(node.name);
            if (it != branch_class.end()) {
                node_entry["branch_class"] = it->second;
            }

            if (!alternative_fit.branch_lengths.empty() && node.id >= 0 && static_cast<size_t>(node.id) < alternative_fit.branch_lengths.size()) {
                node_entry["branch_lengths"]["alternative"] = alternative_fit.branch_lengths[node.id];
                if (node.id != tree.root_id) sum_bl_alt += alternative_fit.branch_lengths[node.id];
            }
            if (!null_fit.branch_lengths.empty() && node.id >= 0 && static_cast<size_t>(node.id) < null_fit.branch_lengths.size()) {
                node_entry["branch_lengths"]["null"] = null_fit.branch_lengths[node.id];
                if (node.id != tree.root_id) sum_bl_null += null_fit.branch_lengths[node.id];
            }

            nodes_json[node.name] = node_entry;
        }

        j["phylogeny"]["nodes"] = nodes_json;
        j["phylogeny"]["tree_lengths"]["alternative"] = sum_bl_alt;
        j["phylogeny"]["tree_lengths"]["null"] = sum_bl_null;

        // Model fits block
        auto format_distro_modern = [](const RELAXRateDistribution& ref_d, const RELAXRateDistribution& test_d) {
            nlohmann::json d;
            nlohmann::json ref_arr = nlohmann::json::array();
            nlohmann::json test_arr = nlohmann::json::array();
            for (size_t k = 0; k < ref_d.omegas.size(); ++k) {
                ref_arr.push_back({
                    {"class", static_cast<int>(k + 1)},
                    {"omega", ref_d.omegas[k]},
                    {"proportion", ref_d.weights[k]}
                });
                test_arr.push_back({
                    {"class", static_cast<int>(k + 1)},
                    {"omega", test_d.omegas[k]},
                    {"proportion", test_d.weights[k]}
                });
            }
            d["reference"] = ref_arr;
            d["test"] = test_arr;
            return d;
        };

        if (gtr_parameters > 0 || gtr_log_likelihood != 0.0) {
            j["model_fits"]["nucleotide_gtr"] = {
                {"description", "General Time Reversible nucleotide model with empirical frequencies"},
                {"log_likelihood", gtr_log_likelihood},
                {"parameters_count", gtr_parameters},
                {"aicc", gtr_aicc}
            };
        }

        j["model_fits"]["codon_mg94"] = {
            {"description", "Muse-Gaut 1994 x GTR with separate rates for Reference and Test branches"},
            {"log_likelihood", mg94_log_likelihood},
            {"parameters_count", mg94_parameters},
            {"aicc", mg94_aicc},
            {"rate_distributions", {
                {"omega_reference", mg94_omega_R},
                {"omega_test", mg94_omega_T}
            }}
        };

        j["model_fits"]["alternative"] = {
            {"description", "RELAX alternative model with unconstrained relaxation parameter K"},
            {"log_likelihood", alternative_fit.log_likelihood},
            {"parameters_count", alternative_fit.parameters},
            {"aicc", alternative_fit.aicc},
            {"k_parameter", k},
            {"rate_distributions", format_distro_modern(alternative_fit.reference_distribution, alternative_fit.test_distribution)}
        };

        j["model_fits"]["null"] = {
            {"description", "RELAX null model fixing K = 1"},
            {"log_likelihood", null_fit.log_likelihood},
            {"parameters_count", null_fit.parameters},
            {"aicc", null_fit.aicc},
            {"k_parameter", 1.0},
            {"rate_distributions", format_distro_modern(null_fit.reference_distribution, null_fit.test_distribution)}
        };

        // Statistical tests block
        std::string decision_str;
        if (p_value <= 0.05) {
            decision_str = (k < 1.0) ? "Relaxed selection detected on test branches" : "Intensified selection detected on test branches";
        } else {
            decision_str = "No evidence of relaxed or intensified selection";
        }

        j["statistical_tests"]["hypothesis_test"] = {
            {"method", "Likelihood Ratio Test (LRT)"},
            {"test_statistic", "LRT"},
            {"statistic_value", lrt},
            {"distribution", "Asymptotic Chi-squared"},
            {"degrees_of_freedom", 1},
            {"null_hypothesis", "K = 1 (Selection intensity is equal between Test and Reference)"},
            {"alternative_hypothesis", "K != 1 (Selection is relaxed K < 1 or intensified K > 1)"},
            {"p_value", p_value},
            {"k_parameter", k},
            {"decision", decision_str}
        };

        return j;
    }

    nlohmann::json to_json(
        const Tree& tree,
        const Alignment& aln,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        if (format == JSONFormat::Legacy) {
            return to_legacy_json(tree, aln);
        }
        return to_modern_json(tree, aln, prov);
    }

    void save_json(
        const std::string& filepath,
        const Tree& tree,
        const Alignment& aln,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open file for writing: " + filepath);
        }
        nlohmann::json j = to_json(tree, aln, format, prov);
        out << j.dump(2) << "\n";
    }

};

class RELAXAnalyzer {
public:
    Tree tree;
    Alignment aln;
    RELAXSettings settings;
    std::shared_ptr<const GeneticCode> code;
    MG94Parameters base_params;

    std::vector<bool> is_test_branch;
    std::vector<size_t> leaf_to_taxon;
    std::vector<std::string> test_branch_names;
    std::vector<std::string> ref_branch_names;
    std::unordered_map<std::string, std::string> branch_class;

    RELAXAnalyzer(
        Tree t,
        Alignment a,
        RELAXSettings s = RELAXSettings{}
    ) : tree(std::move(t)), aln(std::move(a)), settings(std::move(s)) {
        code = aln.code ? aln.code : GeneticCode::universal();
        size_t num_nodes = tree.num_nodes();
        is_test_branch.assign(num_nodes, false);
        leaf_to_taxon.assign(num_nodes, static_cast<size_t>(-1));

        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        partition_branches();
    }

    static RELAXAnalyzer create(
        Tree tree,
        Alignment aln,
        RELAXSettings settings = RELAXSettings{}
    ) {
        return RELAXAnalyzer(std::move(tree), std::move(aln), std::move(settings));
    }

    void partition_branches() {
        bool has_tagged_test = false;

        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                if (node.model_tag == "T" || node.model_tag == "Test" || node.model_tag == "test") {
                    has_tagged_test = true;
                }
            }
        }

        std::regex test_regex;
        bool use_test_regex = false;
        if (!settings.test_branch_regex.empty()) {
            try {
                test_regex = std::regex(settings.test_branch_regex, std::regex::icase);
                use_test_regex = true;
            } catch (const std::exception& e) {
                std::cerr << "Warning: Invalid test branch regex: " << e.what() << "\n";
            }
        }

        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;

            bool is_test = false;

            if (use_test_regex) {
                is_test = std::regex_search(node.name, test_regex);
            } else if (!settings.test_branch_names.empty()) {
                is_test = (std::find(settings.test_branch_names.begin(), settings.test_branch_names.end(), node.name) != settings.test_branch_names.end());
            } else if (has_tagged_test) {
                is_test = (node.model_tag == "T" || node.model_tag == "Test" || node.model_tag == "test");
            }

            is_test_branch[node.id] = is_test;
            if (is_test) {
                test_branch_names.push_back(node.name);
                branch_class[node.name] = "Test";
            } else {
                ref_branch_names.push_back(node.name);
                branch_class[node.name] = "Reference";
            }
        }
    }

    // Evaluates 3-category RELAX mixture log-likelihood
    Scalar evaluate_relax_mixture(
        const Tree& t,
        const std::vector<Scalar>& omegas_R,
        const std::vector<Scalar>& weights,
        Scalar k,
        Scalar tree_scale,
        std::vector<Scalar>* out_site_ll = nullptr,
        const std::vector<Scalar>* custom_branch_lengths = nullptr
    ) const {
        size_t num_patterns = aln.patterns.size();
        size_t num_nodes = t.num_nodes();
        size_t S = code->num_sense_codons;
        const Vector& eq_freqs = aln.codon_frequencies_f3x4;

        // Derive Test branch omegas: omega_T,i = (omega_R,i)^k
        std::vector<Scalar> omegas_T(3);
        for (int c = 0; c < 3; ++c) {
            if (omegas_R[c] <= 1e-12) {
                omegas_T[c] = (k > 0) ? 0.0 : 1.0;
            } else {
                omegas_T[c] = std::clamp(std::pow(omegas_R[c], k), 1e-6, 10000.0);
            }
        }

        // Precompute transition matrices for all 3 classes across all tree branches
        std::vector<std::vector<Matrix>> P_branches(3, std::vector<Matrix>(num_nodes));

        for (int c = 0; c < 3; ++c) {
            MG94Parameters p_R = base_params;
            p_R.alpha = 1.0;
            p_R.beta = omegas_R[c];
            MG94Matrix mat_R;
            mat_R.update(p_R, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            Scalar conv_R = (mat_R.scale_factor > 1e-12) ? (3.0 / mat_R.scale_factor) : 1.0;

            MG94Parameters p_T = base_params;
            p_T.alpha = 1.0;
            p_T.beta = omegas_T[c];
            MG94Matrix mat_T;
            mat_T.update(p_T, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            Scalar conv_T = (mat_T.scale_factor > 1e-12) ? (3.0 / mat_T.scale_factor) : 1.0;

            for (const auto& node : t.nodes) {
                if (node.id != t.root_id) {
                    int nid = node.id;
                    Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[nid] : (node.branch_length * tree_scale);
                    if (is_test_branch[nid]) {
                        P_branches[c][nid] = mat_T.transition_matrix(bl * conv_T);
                    } else {
                        P_branches[c][nid] = mat_R.transition_matrix(bl * conv_R);
                    }
                }
            }
        }

        std::vector<Scalar> pattern_ll(num_patterns, 0.0);
        Scalar total_lnl = 0.0;

        #pragma omp parallel
        {
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));

            #pragma omp for reduction(+:total_lnl) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pat_L = 0.0;

                for (int c = 0; c < 3; ++c) {
                    if (weights[c] <= 1e-12) continue;

                    for (int32_t nid : t.post_order) {
                        const auto& node = t.nodes[nid];
                        if (node.is_leaf) {
                            size_t t_idx = leaf_to_taxon[nid];
                            if (t_idx != static_cast<size_t>(-1)) {
                                int code_val = pattern.states[t_idx];
                                if (code_val >= 0 && code_val < static_cast<int>(S)) {
                                    node_L[nid].setZero();
                                    node_L[nid](code_val) = 1.0;
                                } else {
                                    node_L[nid].setOnes();
                                }
                            } else {
                                node_L[nid].setOnes();
                            }
                        } else {
                            node_L[nid].setOnes();
                            for (int32_t child_id : node.children) {
                                node_L[nid].array() *= (P_branches[c][child_id] * node_L[child_id]).array();
                            }
                        }
                    }

                    Scalar class_L = eq_freqs.dot(node_L[t.root_id]);
                    if (class_L > 0.0) {
                        pat_L += weights[c] * class_L;
                    }
                }

                Scalar pll = (pat_L > 0.0) ? std::log(pat_L) : -1e20;
                total_lnl += pattern.weight * pll;
                if (out_site_ll) {
                    pattern_ll[p] = pll;
                }
            }
        }

        if (out_site_ll) {
            out_site_ll->resize(aln.num_codons);
            for (size_t s = 0; s < aln.num_codons; ++s) {
                (*out_site_ll)[s] = pattern_ll[aln.site_to_pattern[s]];
            }
        }

        return total_lnl;
    }

    // Evaluates log-likelihood and exact analytical branch gradients d ln L / d t_b
    // under the 3-category RELAX mixture model using Inside-Outside adjoints
    std::pair<Scalar, std::vector<Scalar>> compute_relax_branch_gradients(
        const std::vector<Scalar>& branch_lengths,
        const std::vector<Scalar>& omegas_R,
        const std::vector<Scalar>& weights,
        Scalar k
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        const Vector& eq_freqs = aln.codon_frequencies_f3x4;

        // Derive Test branch omegas: omega_T,i = (omega_R,i)^k
        std::vector<Scalar> omegas_T(3);
        for (int c = 0; c < 3; ++c) {
            if (omegas_R[c] <= 1e-12) {
                omegas_T[c] = (k > 0) ? 0.0 : 1.0;
            } else {
                omegas_T[c] = std::clamp(std::pow(omegas_R[c], k), 1e-6, 10000.0);
            }
        }

        // Precompute transition matrices P_c and derivatives dP_c for all 3 classes
        std::vector<std::vector<Matrix>> P_c(3, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<Matrix>> dP_c(3, std::vector<Matrix>(num_nodes));

        for (int c = 0; c < 3; ++c) {
            if (weights[c] <= 1e-12) continue;

            MG94Parameters p_R = base_params;
            p_R.alpha = 1.0;
            p_R.beta = omegas_R[c];
            MG94Matrix mat_R;
            mat_R.update(p_R, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            Scalar conv_R = (mat_R.scale_factor > 1e-12) ? (3.0 / mat_R.scale_factor) : 1.0;

            MG94Parameters p_T = base_params;
            p_T.alpha = 1.0;
            p_T.beta = omegas_T[c];
            MG94Matrix mat_T;
            mat_T.update(p_T, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            Scalar conv_T = (mat_T.scale_factor > 1e-12) ? (3.0 / mat_T.scale_factor) : 1.0;

            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = branch_lengths[nid];
                    if (is_test_branch[nid]) {
                        Matrix P = mat_T.transition_matrix(bl * conv_T);
                        P_c[c][nid] = P;
                        dP_c[c][nid] = conv_T * (mat_T.Q * P);
                    } else {
                        Matrix P = mat_R.transition_matrix(bl * conv_R);
                        P_c[c][nid] = P;
                        dP_c[c][nid] = conv_R * (mat_R.Q * P);
                    }
                }
            }
        }

        Scalar total_log_l = 0.0;
        std::vector<Scalar> grad_b(num_nodes, 0.0);

        #pragma omp parallel
        {
            std::vector<Scalar> local_grad(num_nodes, 0.0);
            Scalar local_ll = 0.0;

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pat_L = 0.0;
                std::array<LikelihoodEngine::InsideOutsideResult, 3> io_c;

                for (int c = 0; c < 3; ++c) {
                    if (weights[c] <= 1e-12) continue;
                    io_c[c] = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_c[c], eq_freqs, S
                    );
                    pat_L += weights[c] * io_c[c].likelihood;
                }

                if (pat_L > 0.0) {
                    local_ll += pattern.weight * std::log(pat_L);
                    Scalar inv_L = pattern.weight / pat_L;
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t vid = node.id;
                            Scalar dL_dt = 0.0;
                            for (int c = 0; c < 3; ++c) {
                                if (weights[c] > 1e-12) {
                                    dL_dt += weights[c] * io_c[c].V[vid].dot(dP_c[c][vid] * io_c[c].D[vid]);
                                }
                            }
                            local_grad[vid] += inv_L * dL_dt;
                        }
                    }
                }
            }

            #pragma omp critical
            {
                total_log_l += local_ll;
                for (size_t i = 0; i < num_nodes; ++i) {
                    grad_b[i] += local_grad[i];
                }
            }
        }

        return {total_log_l, grad_b};
    }

    // High-performance joint branch length refinement via exact Inside-Outside adjoint gradients and L-BFGS
    std::vector<Scalar> refine_branch_lengths(
        const std::vector<Scalar>& init_bl,
        const std::vector<Scalar>& omegas_R,
        const std::vector<Scalar>& weights,
        Scalar k,
        int max_iters = 10
    ) const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Scalar> bl = init_bl;
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id)) {
                bl[i] = std::clamp(bl[i], 1e-6, 10.0);
            }
        }

        const size_t m_history = 5;
        std::vector<std::vector<Scalar>> s_hist;
        std::vector<std::vector<Scalar>> y_hist;
        std::vector<Scalar> rho_hist;

        auto [curr_ll, curr_grad] = compute_relax_branch_gradients(bl, omegas_R, weights, k);

        Scalar best_ll = curr_ll;
        std::vector<Scalar> best_bl = bl;

        for (int iter = 0; iter < max_iters; ++iter) {
            Scalar max_g = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_g = std::max(max_g, std::abs(curr_grad[i]));
                }
            }
            if (max_g < 1e-3) break;

            // Two-loop L-BFGS recursion
            std::vector<Scalar> q = curr_grad;
            std::vector<Scalar> alpha(s_hist.size());

            for (int j = static_cast<int>(s_hist.size()) - 1; j >= 0; --j) {
                Scalar s_dot_q = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    s_dot_q += s_hist[j][i] * q[i];
                }
                alpha[j] = rho_hist[j] * s_dot_q;
                for (size_t i = 0; i < num_nodes; ++i) {
                    q[i] -= alpha[j] * y_hist[j][i];
                }
            }

            Scalar gamma = 1.0;
            if (!s_hist.empty()) {
                Scalar s_dot_y = 0.0;
                Scalar y_dot_y = 0.0;
                const auto& s_last = s_hist.back();
                const auto& y_last = y_hist.back();
                for (size_t i = 0; i < num_nodes; ++i) {
                    s_dot_y += s_last[i] * y_last[i];
                    y_dot_y += y_last[i] * y_last[i];
                }
                if (y_dot_y > 1e-12) {
                    gamma = std::clamp(s_dot_y / y_dot_y, 1e-4, 10.0);
                }
            }

            std::vector<Scalar> r(num_nodes);
            for (size_t i = 0; i < num_nodes; ++i) {
                r[i] = gamma * q[i];
            }

            for (size_t j = 0; j < s_hist.size(); ++j) {
                Scalar y_dot_r = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    y_dot_r += y_hist[j][i] * r[i];
                }
                Scalar beta = rho_hist[j] * y_dot_r;
                for (size_t i = 0; i < num_nodes; ++i) {
                    r[i] += s_hist[j][i] * (alpha[j] - beta);
                }
            }

            std::vector<Scalar> d = r;
            d[tree.root_id] = 0.0;

            Scalar dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    dir_dot_grad += d[i] * curr_grad[i];
                }
            }
            if (dir_dot_grad <= 0.0) {
                d = curr_grad;
                d[tree.root_id] = 0.0;
                dir_dot_grad = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) {
                        dir_dot_grad += d[i] * curr_grad[i];
                    }
                }
                s_hist.clear();
                y_hist.clear();
                rho_hist.clear();
            }

            // Backtracking Armijo line search
            Scalar max_step = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_step = std::max(max_step, std::abs(d[i]));
                }
            }
            Scalar step_size = s_hist.empty() ? std::min(1.0, 0.05 / std::max(1.0, max_step)) : 1.0;
            std::vector<Scalar> new_bl(num_nodes);
            Scalar new_ll = curr_ll;
            std::vector<Scalar> new_grad;

            bool line_search_ok = false;
            for (int ls = 0; ls < 8; ++ls) {
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) {
                        new_bl[i] = std::clamp(bl[i] + step_size * d[i], 1e-6, 10.0);
                    } else {
                        new_bl[i] = 0.0;
                    }
                }
                auto [test_ll, test_grad] = compute_relax_branch_gradients(new_bl, omegas_R, weights, k);
                if (test_ll > curr_ll + 1e-4 * step_size * dir_dot_grad) {
                    new_ll = test_ll;
                    new_grad = test_grad;
                    line_search_ok = true;
                    break;
                }
                step_size *= 0.5;
            }

            if (!line_search_ok) {
                break;
            }

            // Update L-BFGS history
            std::vector<Scalar> s_k(num_nodes);
            std::vector<Scalar> y_k(num_nodes);
            Scalar s_dot_y = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                s_k[i] = new_bl[i] - bl[i];
                y_k[i] = curr_grad[i] - new_grad[i];
                s_dot_y += s_k[i] * y_k[i];
            }

            if (s_dot_y > 1e-8) {
                if (s_hist.size() >= m_history) {
                    s_hist.erase(s_hist.begin());
                    y_hist.erase(y_hist.begin());
                    rho_hist.erase(rho_hist.begin());
                }
                s_hist.push_back(s_k);
                y_hist.push_back(y_k);
                rho_hist.push_back(1.0 / s_dot_y);
            }

            bl = new_bl;
            curr_ll = new_ll;
            curr_grad = new_grad;

            if (curr_ll > best_ll) {
                best_ll = curr_ll;
                best_bl = bl;
            }
        }

        return best_bl;
    }

    RELAXResult run(
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        bool show_progress = false,
        bool force_progress = false
    ) {
        auto start_time = std::chrono::high_resolution_clock::now();
        RELAXResult result;
        result.settings = settings;
        result.test_branches = test_branch_names;
        result.reference_branches = ref_branch_names;
        result.branch_class = branch_class;

        std::unique_ptr<ProgressBar> pbar;
        if (show_progress || force_progress) {
            pbar = std::make_unique<ProgressBar>(100, "RELAX Analysis");
        }

        auto update_progress = [&](const std::string& stage, double frac) {
            if (pbar) {
                pbar->set_status(stage);
                pbar->update(static_cast<int>(frac * 100.0));
            } else if (progress_cb) {
                progress_cb(stage, frac);
            }
        };

        // ---------------------------------------------------------
        // Phase 1: Nucleotide GTR Model Fit
        // ---------------------------------------------------------
        update_progress("Phase 1: Fitting Nucleotide GTR model", 0.05);
        GTRFitter gtr_fitter(tree, aln);
        auto gtr_res = gtr_fitter.fit();

        result.gtr_log_likelihood = gtr_res.log_likelihood;
        result.gtr_aicc = gtr_res.aicc;
        result.gtr_parameters = 5 + (tree.num_nodes() - 1) + 2;
        result.gtr_rates = gtr_res.params;

        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                result.gtr_branch_lengths[node.name] = node.branch_length;
            }
        }

        base_params.theta_AC = gtr_res.params.theta_AC;
        base_params.theta_AT = gtr_res.params.theta_AT;
        base_params.theta_CG = gtr_res.params.theta_CG;
        base_params.theta_CT = gtr_res.params.theta_CT;
        base_params.theta_GT = gtr_res.params.theta_GT;

        update_progress("Phase 1: GTR fit complete", 0.20);

        // ---------------------------------------------------------
        // Phase 2: MG94 with Separate Rates for Branch Sets
        // ---------------------------------------------------------
        update_progress("Phase 2: Fitting MG94 with separate branch rates", 0.25);

        auto eval_sep_mg94 = [&](const std::array<Scalar, 3>& x) -> Scalar {
            Scalar w_R = x[0];
            Scalar w_T = x[1];
            Scalar s = x[2];
            Scalar dummy_k = (std::abs(std::log(std::max(1e-5, w_R))) > 1e-6)
                ? (std::log(std::max(1e-5, w_T)) / std::log(std::max(1e-5, w_R)))
                : 1.0;
            return -evaluate_relax_mixture(gtr_res.tree, {w_R, w_R, w_R}, {1.0, 0.0, 0.0}, dummy_k, s);
        };

        auto opt_sep = NelderMeadND<3>::minimize(
            eval_sep_mg94,
            {0.1, 0.5, 1.0},
            {1e-4, 1e-4, 0.1},
            {10.0, 10.0, 10.0},
            1e-4, 100
        );

        result.mg94_log_likelihood = -opt_sep.f;
        result.mg94_omega_R = opt_sep.x[0];
        result.mg94_omega_T = opt_sep.x[1];
        Scalar sep_scale = opt_sep.x[2];

        size_t B = tree.num_nodes() - 1;
        result.mg94_parameters = 5 + B + 2 + 6; // GTR rates + branches + 2 omegas + CF3x4
        Scalar n_codons = static_cast<Scalar>(aln.num_codons);
        Scalar k_mg = static_cast<Scalar>(result.mg94_parameters);
        result.mg94_aicc = 2.0 * k_mg - 2.0 * result.mg94_log_likelihood +
            (2.0 * k_mg * (k_mg + 1.0)) / std::max(1.0, n_codons - k_mg - 1.0);

        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                result.mg94_branch_lengths[node.name] = node.branch_length * sep_scale;
            }
        }

        update_progress("Phase 2: MG94 separate rates fit complete", 0.45);

        // ---------------------------------------------------------
        // Phase 3: RELAX Alternative Model Fit (k free)
        // ---------------------------------------------------------
        update_progress("Phase 3: Fitting RELAX alternative model (K free)", 0.50);

        // Initial guess for k from MG94 separate rates
        Scalar k_init = 1.0;
        if (std::abs(std::log(std::max(1e-5, result.mg94_omega_R))) > 1e-5) {
            k_init = std::clamp(
                std::log(std::max(1e-5, result.mg94_omega_T)) / std::log(std::max(1e-5, result.mg94_omega_R)),
                0.01, 10.0
            );
        }

        // Optimization over 7 parameters: [s, w0, w1, w2, p0, p1, k]
        // Reference omegas: w0 <= w1 <= 1 <= w2
        // Weights: p0 in [0, 1], p1 in [0, 1-p0], p2 = 1 - p0 - p1
        auto eval_alt = [&](const std::array<Scalar, 7>& x) -> Scalar {
            Scalar s = x[0];
            Scalar w0 = x[1];
            Scalar w1 = x[2];
            Scalar w2 = x[3];
            Scalar p0 = x[4];
            Scalar p1 = x[5];
            if (w0 > w1 || w1 > 1.0 || w2 < 1.0) return 1e20;
            if (p0 + p1 >= 1.0) return 1e20;
            Scalar p2 = 1.0 - p0 - p1;
            Scalar k = x[6];

            std::vector<Scalar> omegas_R = {w0, w1, w2};
            std::vector<Scalar> weights = {p0, p1, p2};
            return -evaluate_relax_mixture(gtr_res.tree, omegas_R, weights, k, s);
        };

        std::vector<std::array<Scalar, 7>> grid_alt = {
            {sep_scale, 0.0, result.mg94_omega_R, 5.0, 0.35, 0.60, k_init},
            {sep_scale, 0.0, 0.16, 30.0, 0.35, 0.65, 0.0001},
            {sep_scale, 0.05, 0.3, 2.0, 0.40, 0.50, 0.5},
            {sep_scale, 0.0, 0.5, 5.0, 0.3, 0.4, 2.0}
        };

        auto opt_alt = NelderMeadND<7>::minimize(
            eval_alt,
            grid_alt[0],
            {0.1, 0.0, 1e-4, 1.0, 1e-4, 1e-4, 0.0},
            {5.0, 1.0, 1.0, 50.0, 0.99, 0.99, 50.0},
            1e-4, 300,
            grid_alt
        );

        result.alternative_fit.log_likelihood = -opt_alt.f;
        result.alternative_fit.tree_scale = opt_alt.x[0];
        Scalar alt_s = opt_alt.x[0];
        Scalar alt_w0 = opt_alt.x[1];
        Scalar alt_w1 = opt_alt.x[2];
        Scalar alt_w2 = opt_alt.x[3];
        Scalar alt_p0 = opt_alt.x[4];
        Scalar alt_p1 = opt_alt.x[5];
        Scalar alt_p2 = std::max(0.0, 1.0 - alt_p0 - alt_p1);
        Scalar alt_k = opt_alt.x[6];

        result.alternative_fit.k = alt_k;
        result.alternative_fit.reference_distribution.omegas = {alt_w0, alt_w1, alt_w2};
        result.alternative_fit.reference_distribution.weights = {alt_p0, alt_p1, alt_p2};

        std::vector<Scalar> alt_omegas_T = {
            (alt_w0 <= 1e-12) ? ((alt_k > 0) ? 0.0 : 1.0) : std::pow(alt_w0, alt_k),
            (alt_w1 <= 1e-12) ? ((alt_k > 0) ? 0.0 : 1.0) : std::pow(alt_w1, alt_k),
            std::pow(alt_w2, alt_k)
        };
        result.alternative_fit.test_distribution.omegas = alt_omegas_T;
        result.alternative_fit.test_distribution.weights = {alt_p0, alt_p1, alt_p2};

        result.alternative_fit.parameters = 5 + B + 5 + 1 + 6; // GTR rates + branches + 3 omegas + 2 weights + k + CF3x4
        Scalar k_alt_params = static_cast<Scalar>(result.alternative_fit.parameters);
        result.alternative_fit.aicc = 2.0 * k_alt_params - 2.0 * result.alternative_fit.log_likelihood +
            (2.0 * k_alt_params * (k_alt_params + 1.0)) / std::max(1.0, n_codons - k_alt_params - 1.0);

        result.alternative_fit.branch_lengths.assign(tree.num_nodes(), 0.0);
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                result.alternative_fit.branch_lengths[node.id] = node.branch_length * alt_s;
            }
        }

        if (settings.refine_branch_lengths) {
            update_progress("Phase 3: Refining alternative branch lengths", 0.65);
            result.alternative_fit.branch_lengths = refine_branch_lengths(
                result.alternative_fit.branch_lengths,
                result.alternative_fit.reference_distribution.omegas,
                result.alternative_fit.reference_distribution.weights,
                alt_k,
                15
            );
            result.alternative_fit.log_likelihood = evaluate_relax_mixture(
                gtr_res.tree, result.alternative_fit.reference_distribution.omegas,
                result.alternative_fit.reference_distribution.weights, alt_k, 1.0,
                &result.alternative_fit.site_log_likelihoods,
                &result.alternative_fit.branch_lengths
            );
            result.alternative_fit.aicc = 2.0 * k_alt_params - 2.0 * result.alternative_fit.log_likelihood +
                (2.0 * k_alt_params * (k_alt_params + 1.0)) / std::max(1.0, n_codons - k_alt_params - 1.0);
        } else {
            // Compute site log-likelihoods for alternative model
            evaluate_relax_mixture(gtr_res.tree, result.alternative_fit.reference_distribution.omegas,
                                   result.alternative_fit.reference_distribution.weights, alt_k, alt_s,
                                   &result.alternative_fit.site_log_likelihoods);
        }

        update_progress("Phase 3: RELAX alternative fit complete", 0.70);

        // ---------------------------------------------------------
        // Phase 4: RELAX Null Model Fit (k = 1.0 fixed)
        // ---------------------------------------------------------
        update_progress("Phase 4: Fitting RELAX null model (K = 1 fixed)", 0.75);

        auto eval_null = [&](const std::array<Scalar, 6>& x) -> Scalar {
            Scalar s = x[0];
            Scalar w0 = x[1];
            Scalar w1 = x[2];
            Scalar w2 = x[3];
            Scalar p0 = x[4];
            Scalar p1 = x[5];
            if (w0 > w1 || w1 > 1.0 || w2 < 1.0) return 1e20;
            if (p0 + p1 >= 1.0) return 1e20;
            Scalar p2 = 1.0 - p0 - p1;

            std::vector<Scalar> omegas_R = {w0, w1, w2};
            std::vector<Scalar> weights = {p0, p1, p2};
            return -evaluate_relax_mixture(gtr_res.tree, omegas_R, weights, 1.0, s);
        };

        std::vector<std::array<Scalar, 6>> grid_null = {
            {sep_scale, 0.0, 0.5, 1.5, 0.6, 0.3},
            {sep_scale, 0.0, 0.6, 1.2, 0.67, 0.27},
            {alt_s, alt_w0, alt_w1, alt_w2, alt_p0, alt_p1}
        };

        auto opt_null = NelderMeadND<6>::minimize(
            eval_null,
            grid_null[0],
            {0.1, 0.0, 1e-4, 1.0, 1e-4, 1e-4},
            {5.0, 1.0, 1.0, 50.0, 0.99, 0.99},
            1e-4, 300,
            grid_null
        );

        result.null_fit.log_likelihood = -opt_null.f;
        result.null_fit.tree_scale = opt_null.x[0];
        Scalar null_s = opt_null.x[0];
        Scalar null_w0 = opt_null.x[1];
        Scalar null_w1 = opt_null.x[2];
        Scalar null_w2 = opt_null.x[3];
        Scalar null_p0 = opt_null.x[4];
        Scalar null_p1 = opt_null.x[5];
        Scalar null_p2 = std::max(0.0, 1.0 - null_p0 - null_p1);

        result.null_fit.k = 1.0;
        result.null_fit.reference_distribution.omegas = {null_w0, null_w1, null_w2};
        result.null_fit.reference_distribution.weights = {null_p0, null_p1, null_p2};
        result.null_fit.test_distribution = result.null_fit.reference_distribution;

        result.null_fit.parameters = 5 + B + 5 + 6; // 1 fewer parameter than alternative (k fixed)
        Scalar k_null_params = static_cast<Scalar>(result.null_fit.parameters);
        result.null_fit.aicc = 2.0 * k_null_params - 2.0 * result.null_fit.log_likelihood +
            (2.0 * k_null_params * (k_null_params + 1.0)) / std::max(1.0, n_codons - k_null_params - 1.0);

        result.null_fit.branch_lengths.assign(tree.num_nodes(), 0.0);
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                result.null_fit.branch_lengths[node.id] = node.branch_length * null_s;
            }
        }

        if (settings.refine_branch_lengths) {
            update_progress("Phase 4: Refining null branch lengths", 0.85);
            result.null_fit.branch_lengths = refine_branch_lengths(
                result.null_fit.branch_lengths,
                result.null_fit.reference_distribution.omegas,
                result.null_fit.reference_distribution.weights,
                1.0,
                15
            );
            result.null_fit.log_likelihood = evaluate_relax_mixture(
                gtr_res.tree, result.null_fit.reference_distribution.omegas,
                result.null_fit.reference_distribution.weights, 1.0, 1.0,
                &result.null_fit.site_log_likelihoods,
                &result.null_fit.branch_lengths
            );
            result.null_fit.aicc = 2.0 * k_null_params - 2.0 * result.null_fit.log_likelihood +
                (2.0 * k_null_params * (k_null_params + 1.0)) / std::max(1.0, n_codons - k_null_params - 1.0);
        } else {
            evaluate_relax_mixture(gtr_res.tree, result.null_fit.reference_distribution.omegas,
                                   result.null_fit.reference_distribution.weights, 1.0, null_s,
                                   &result.null_fit.site_log_likelihoods);
        }

        update_progress("Phase 4: RELAX null fit complete", 0.90);

        // ---------------------------------------------------------
        // Phase 5: Hypothesis Testing & Interpretation
        // ---------------------------------------------------------
        update_progress("Phase 5: Calculating LRT & asymptotic p-value", 0.95);

        result.lrt = 2.0 * std::max(0.0, result.alternative_fit.log_likelihood - result.null_fit.log_likelihood);
        result.p_value = std::erfc(std::sqrt(result.lrt / 2.0));
        result.k = alt_k;

        result.is_significant = (result.p_value <= settings.p_value_threshold);
        if (result.is_significant) {
            if (result.k > 1.0) {
                result.is_intensified = true;
            } else {
                result.is_relaxed = true;
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        result.runtime_seconds = std::chrono::duration<double>(end_time - start_time).count();

        update_progress("RELAX analysis complete", 1.0);

        return result;
    }
};

} // namespace hyphy::analyses
