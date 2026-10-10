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
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

struct ABSRELRateDistribution {
    std::vector<Scalar> rates;   // omega values
    std::vector<Scalar> weights; // mixture proportions summing to 1.0
};

struct ABSRELBranchResult {
    std::string branch_name;
    int32_t node_id = -1;
    Scalar baseline_omega = 1.0;
    Scalar baseline_branch_length = 0.0;
    int rate_classes = 1;
    ABSRELRateDistribution rate_distribution;
    Scalar full_branch_length = 0.0;
    Scalar full_es = 0.0; // synonymous subs/site
    Scalar full_en = 0.0; // nonsynonymous subs/site
    Scalar lrt = 0.0;
    Scalar uncorrected_p_value = 1.0;
    Scalar corrected_p_value = 1.0;
    bool is_tested = true;
    bool is_positive = false;
    int sites_ebf_100 = 0;
    std::vector<Scalar> site_ebf;
};

struct ABSRELFitSummary {
    Scalar log_likelihood = 0.0;
    size_t parameters = 0;
    Scalar aicc = 0.0;
};

struct ABSRELSettings {
    int max_rate_classes = 3;
    Scalar p_threshold = 0.05;
    std::string test_branches = "All"; // "All", "Internal", "Leaves"
    bool do_srv = false;
    int syn_rate_classes = 3;
    bool verbose = false;
};

struct ABSRELResult {
    ABSRELFitSummary gtr_fit;
    ABSRELFitSummary baseline_fit;
    ABSRELFitSummary full_adaptive_fit;
    std::unordered_map<std::string, ABSRELBranchResult> branches;
    std::vector<std::string> tested_branches;
    std::vector<std::string> positive_branches;
    Scalar p_threshold = 0.05;
    double runtime_seconds = 0.0;

    nlohmann::json to_legacy_json(const Tree& tree, const Alignment& aln) const {
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
            if (bres.is_tested) {
                tested_map[name] = "test";
            }
        }
        j["tested"] = {{"0", tested_map}};

        // Timers
        j["timers"] = {
            {"Overall", {{"timer", static_cast<int>(runtime_seconds)}, {"order", 0}}}
        };

        return j;
    }

    nlohmann::json to_modern_json(const Tree& tree, const Alignment& aln, const Provenance& prov = {}) const {
        nlohmann::json j;

        j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/absrel.v3.schema.json";
        j["schema_version"] = "3.0.0";

        // Analysis block
        j["analysis"]["id"] = "absrel";
        j["analysis"]["name"] = "Adaptive Branch-Site Random Effects Likelihood";
        j["analysis"]["version"] = "3.0.0";
        j["analysis"]["category"] = "branch_site_selection";
        j["analysis"]["description"] = "aBSREL (Adaptive branch-site random effects likelihood) uses an adaptive random effects branch-site model framework to test whether each branch has evolved under positive selection, inferring an optimal number of rate categories per branch.";
        j["analysis"]["citations"] = nlohmann::json::array({
            {
                {"citation", "Smith MD, Wertheim JO, Weaver S, Murrell B, Scheffler K, Kosakovsky Pond SL (2015). Less Is More: An Adaptive Branch-Site Random Effects Model for Efficient Detection of Episodic Diversifying Selection. Mol Biol Evol 32(5): 1342-1353."},
                {"doi", "10.1093/molbev/msv022"},
                {"pmid", "25697341"}
            }
        });

        j["analysis"]["settings"]["p_threshold"] = p_threshold;
        j["analysis"]["settings"]["tested_branches"] = tested_branches.size();

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
            auto_prov.invocation.cli_command = "hyphy3 absrel";
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
        double sum_bl_base = 0.0;
        double sum_bl_full = 0.0;

        for (const auto& node : tree.nodes) {
            nlohmann::json node_entry;
            node_entry["type"] = node.children.empty() ? "leaf" : "internal";

            auto it = branches.find(node.name);
            if (it != branches.end()) {
                const auto& bres = it->second;
                node_entry["branch_lengths"]["baseline_mg94"] = bres.baseline_branch_length;
                node_entry["branch_lengths"]["full_adaptive"] = bres.full_branch_length;
                node_entry["rate_classes"] = bres.rate_classes;
                node_entry["baseline_omega"] = bres.baseline_omega;
                node_entry["is_tested"] = bres.is_tested;
                node_entry["is_positive"] = bres.is_positive;

                if (node.id != tree.root_id) {
                    sum_bl_base += bres.baseline_branch_length;
                    sum_bl_full += bres.full_branch_length;
                }

                nlohmann::json rdist = nlohmann::json::array();
                for (size_t k = 0; k < bres.rate_distribution.rates.size(); ++k) {
                    rdist.push_back({
                        {"omega", bres.rate_distribution.rates[k]},
                        {"proportion", bres.rate_distribution.weights[k]}
                    });
                }
                node_entry["rate_distribution"] = rdist;

                if (bres.is_tested) {
                    node_entry["lrt"] = bres.lrt;
                    node_entry["uncorrected_p_value"] = bres.uncorrected_p_value;
                    node_entry["corrected_p_value"] = bres.corrected_p_value;
                    node_entry["sites_ebf_100"] = bres.sites_ebf_100;
                }
            }

            nodes_json[node.name] = node_entry;
        }

        j["phylogeny"]["nodes"] = nodes_json;
        j["phylogeny"]["tree_lengths"]["baseline_mg94"] = sum_bl_base;
        j["phylogeny"]["tree_lengths"]["full_adaptive"] = sum_bl_full;

        // Model fits block
        if (gtr_fit.parameters > 0 || gtr_fit.log_likelihood != 0.0) {
            j["model_fits"]["nucleotide_gtr"] = {
                {"description", "General Time Reversible nucleotide model with empirical frequencies"},
                {"log_likelihood", gtr_fit.log_likelihood},
                {"parameters_count", gtr_fit.parameters},
                {"aicc", gtr_fit.aicc}
            };
        }

        j["model_fits"]["baseline_mg94"] = {
            {"description", "Baseline MG94xREV model with individual branch omegas and branch lengths"},
            {"log_likelihood", baseline_fit.log_likelihood},
            {"parameters_count", baseline_fit.parameters},
            {"aicc", baseline_fit.aicc}
        };

        j["model_fits"]["full_adaptive"] = {
            {"description", "Full adaptive aBSREL mixture model with branch-specific rate categories"},
            {"log_likelihood", full_adaptive_fit.log_likelihood},
            {"parameters_count", full_adaptive_fit.parameters},
            {"aicc", full_adaptive_fit.aicc}
        };

        // Statistical tests block
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

        // Branch results (Columnar format)
        j["branch_results"]["columns"] = {
            {"branch", {{"type", "string"}, {"description", "Branch name in phylogeny"}}},
            {"rate_classes", {{"type", "integer"}, {"description", "Number of inferred omega rate classes"}}},
            {"baseline_omega", {{"type", "float"}, {"description", "Baseline single omega estimate"}}},
            {"full_branch_length", {{"type", "float"}, {"unit", "substitutions/site"}, {"description", "Estimated branch length under full adaptive model"}}},
            {"lrt", {{"type", "float"}, {"description", "Likelihood ratio test statistic for episodic positive selection"}}},
            {"uncorrected_p_value", {{"type", "float"}, {"description", "Asymptotic uncorrected p-value"}}},
            {"corrected_p_value", {{"type", "float"}, {"description", "Holm-Bonferroni corrected p-value"}}},
            {"sites_ebf_100", {{"type", "integer"}, {"description", "Number of sites with Empirical Bayes Factor (EBF) >= 100 favoring selection"}}},
            {"positive", {{"type", "boolean"}, {"description", "Flag indicating episodic diversifying selection at threshold"}}}
        };

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

        j["branch_results"]["data"]["branch"] = col_branch;
        j["branch_results"]["data"]["rate_classes"] = col_k;
        j["branch_results"]["data"]["baseline_omega"] = col_base_omega;
        j["branch_results"]["data"]["full_branch_length"] = col_len;
        j["branch_results"]["data"]["lrt"] = col_lrt;
        j["branch_results"]["data"]["uncorrected_p_value"] = col_raw_p;
        j["branch_results"]["data"]["corrected_p_value"] = col_adj_p;
        j["branch_results"]["data"]["sites_ebf_100"] = col_ebf;
        j["branch_results"]["data"]["positive"] = col_pos;

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

class ABSRELAnalyzer {
public:
    Tree tree;
    Alignment aln;
    ABSRELSettings settings;

    const GeneticCode* gcode = nullptr;
    Eigen::Matrix<Scalar, 3, 4> pos_nuc_freqs;
    Vector codon_freqs;
    GTRParameters gtr_params;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;

    std::vector<size_t> leaf_to_taxon;

    // Per-branch parameters:
    // branch_syn_alpha[id]: synonymous rate multiplier alpha_b
    std::vector<Scalar> branch_alpha;
    // branch_omegas[id]: omega rate categories for branch id
    std::vector<std::vector<Scalar>> branch_omegas;
    // branch_weights[id]: mixture proportions for branch id
    std::vector<std::vector<Scalar>> branch_weights;

    // Scale factors for Q matrices
    Scalar syn_scale_factor = 1.0;

    ABSRELAnalyzer(Tree t, Alignment a, ABSRELSettings s = ABSRELSettings{})
        : tree(std::move(t)), aln(std::move(a)), settings(s) {
        gcode = aln.code ? aln.code.get() : GeneticCode::universal().get();
        pos_nuc_freqs = aln.pos_nuc_frequencies;
        codon_freqs = aln.codon_frequencies_f3x4;

        size_t num_nodes = tree.num_nodes();
        leaf_to_taxon.assign(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        branch_alpha.assign(num_nodes, 0.01);
        branch_omegas.resize(num_nodes);
        branch_weights.resize(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i) {
            branch_omegas[i] = {0.5};
            branch_weights[i] = {1.0};
        }

        compute_syn_scale();
    }

    static ABSRELAnalyzer create(Tree t, Alignment a, ABSRELSettings s = ABSRELSettings{}) {
        return ABSRELAnalyzer(std::move(t), std::move(a), s);
    }

    void compute_syn_scale() {
        MG94Matrix mat;
        MG94Parameters p;
        p.alpha = 1.0;
        p.beta = 0.0; // synonymous only
        p.theta_AC = gtr_params.theta_AC;
        p.theta_AT = gtr_params.theta_AT;
        p.theta_CG = gtr_params.theta_CG;
        p.theta_CT = gtr_params.theta_CT;
        p.theta_GT = gtr_params.theta_GT;
        mat.update(p, pos_nuc_freqs, codon_freqs, *gcode);
        syn_scale_factor = mat.scale_factor;
    }

    // Build rate matrix Q(1, omega)
    MG94Matrix build_q_matrix(Scalar omega) const {
        MG94Matrix mat;
        MG94Parameters p;
        p.alpha = 1.0;
        p.beta = omega;
        p.theta_AC = gtr_params.theta_AC;
        p.theta_AT = gtr_params.theta_AT;
        p.theta_CG = gtr_params.theta_CG;
        p.theta_CT = gtr_params.theta_CT;
        p.theta_GT = gtr_params.theta_GT;
        mat.update(p, pos_nuc_freqs, codon_freqs, *gcode);
        return mat;
    }

    // Evaluate transition matrix for branch node_id given its current parameters
    Matrix compute_branch_transition_matrix(int32_t node_id) const {
        int S = gcode->num_sense_codons;
        Scalar alpha = branch_alpha[node_id];
        const auto& omegas = branch_omegas[node_id];
        const auto& weights = branch_weights[node_id];
        size_t K = omegas.size();

        Matrix P_mix = Matrix::Zero(S, S);
        for (size_t k = 0; k < K; ++k) {
            MG94Matrix mat = build_q_matrix(omegas[k]);
            Matrix P_k = mat.transition_matrix(alpha);
            P_mix += weights[k] * P_k;
        }
        return P_mix;
    }

    // Precompute all branch transition matrices across the tree
    std::vector<Matrix> compute_all_branch_transition_matrices() const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Matrix> P_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P_branches[node.id] = compute_branch_transition_matrix(node.id);
            }
        }
        return P_branches;
    }

    // Compute full tree log-likelihood
    Scalar compute_tree_log_likelihood(const std::vector<Matrix>& P_branches) const {
        int S = gcode->num_sense_codons;
        return LikelihoodEngine::compute_tree_log_likelihood_from_P(
            tree, aln, P_branches, codon_freqs, S, &leaf_to_taxon
        );
    }

    // Compute expected branch length (substitutions per nucleotide site) for branch node_id
    Scalar compute_expected_branch_length(int32_t node_id) const {
        Scalar alpha = branch_alpha[node_id];
        const auto& omegas = branch_omegas[node_id];
        const auto& weights = branch_weights[node_id];
        Scalar total_sub_rate = 0.0;
        for (size_t k = 0; k < omegas.size(); ++k) {
            MG94Matrix mat = build_q_matrix(omegas[k]);
            total_sub_rate += weights[k] * (mat.scale_factor / 3.0);
        }
        return alpha * total_sub_rate;
    }

    // Compute synonymous and nonsynonymous component branch lengths
    std::pair<Scalar, Scalar> compute_es_en(int32_t node_id) const {
        Scalar alpha = branch_alpha[node_id];
        Scalar syn_rate = syn_scale_factor / 3.0;
        Scalar es = alpha * syn_rate;

        Scalar total_len = compute_expected_branch_length(node_id);
        Scalar en = std::max(0.0, total_len - es);
        return {es, en};
    }

    // AICc helper: N = num_sequences * num_codons
    Scalar compute_aicc(Scalar log_l, size_t num_params) const {
        Scalar N = static_cast<Scalar>(aln.num_taxa * aln.num_codons);
        Scalar P = static_cast<Scalar>(num_params);
        if (N - P - 1.0 > 0.0) {
            return -2.0 * log_l + 2.0 * P + (2.0 * P * (P + 1.0)) / (N - P - 1.0);
        }
        return -2.0 * log_l + 2.0 * P;
    }

    // Asymptotic p-value for mixture distribution: 0.5 * [0.4 * chi2_1 + 0.6 * chi2_2]
    static Scalar compute_mixture_p_value(Scalar lrt) {
        if (lrt <= 0.0) return 1.0;
        Scalar p_chi1 = std::erfc(std::sqrt(lrt / 2.0));
        Scalar p_chi2 = std::exp(-lrt / 2.0);
        Scalar p = 0.5 * (0.4 * p_chi1 + 0.6 * p_chi2);
        return std::clamp(p, 0.0, 1.0);
    }

    // Phase 1: Fit GTR
    void fit_gtr() {
        GTRFitter gtr_fitter(tree, aln);
        auto gtr_res = gtr_fitter.fit();
        gtr_params = gtr_res.params;
        gtr_log_l = gtr_res.log_likelihood;
        gtr_aicc = gtr_res.aicc;

        compute_syn_scale();

        // Initialize branch alpha using nucleotide branch lengths
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != tree.root_id) {
                // Initial conversion from nucleotide substitutions to synonymous time
                Scalar bl = std::max(node.branch_length, 1e-4);
                branch_alpha[node.id] = bl * (3.0 / std::max(syn_scale_factor, 0.1));
                branch_omegas[node.id] = {0.3};
                branch_weights[node.id] = {1.0};
            }
        }
    }

    // Phase 2: Fit Baseline MG94xREV Model (each branch has 1 omega)
    ABSRELFitSummary fit_baseline(int max_cycles = 4) {
        int S = gcode->num_sense_codons;
        size_t num_patterns = aln.patterns.size();

        std::vector<Matrix> P_branches = compute_all_branch_transition_matrices();
        Scalar prev_ll = compute_tree_log_likelihood(P_branches);

        for (int cycle = 0; cycle < max_cycles; ++cycle) {
            // Precompute V and D for each branch across all patterns
            // To be memory efficient, we iterate over branches
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                int32_t bid = node.id;

                // Collect V and D for branch bid across all patterns
                std::vector<Vector> V_patterns(num_patterns);
                std::vector<Vector> D_patterns(num_patterns);

                #pragma omp parallel for schedule(dynamic)
                for (size_t p = 0; p < num_patterns; ++p) {
                    const auto& pattern = aln.patterns[p];
                    auto io = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_branches, codon_freqs, S
                    );
                    V_patterns[p] = io.V[bid];
                    D_patterns[p] = io.D[bid];
                }

                // Objective function for branch bid: optimize log(alpha) and log(omega)
                auto branch_obj = [&](Scalar log_alpha, Scalar log_omega) -> Scalar {
                    Scalar cur_alpha = std::exp(log_alpha);
                    Scalar cur_omega = std::exp(log_omega);
                    if (cur_alpha > 50.0 || cur_omega > 5000.0) return 1e20;

                    MG94Matrix mat = build_q_matrix(cur_omega);
                    Matrix P_b = mat.transition_matrix(cur_alpha);

                    Scalar ll = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Scalar Lp = V_patterns[p].dot(P_b * D_patterns[p]);
                        if (Lp > 0.0) {
                            ll += aln.patterns[p].weight * std::log(Lp);
                        } else {
                            ll += aln.patterns[p].weight * (-1e20);
                        }
                    }
                    return -ll;
                };

                Scalar init_log_alpha = std::log(std::clamp(branch_alpha[bid], 1e-4, 20.0));
                Scalar init_log_omega = std::log(std::clamp(branch_omegas[bid][0], 1e-3, 50.0));

                auto opt_res = NelderMead2D::minimize(branch_obj, init_log_alpha, init_log_omega, -10.0, 10.0, 1e-4, 150);
                Scalar opt_alpha = std::exp(opt_res.a);
                Scalar opt_omega = std::exp(opt_res.b);

                branch_alpha[bid] = opt_alpha;
                branch_omegas[bid][0] = opt_omega;
                branch_weights[bid][0] = 1.0;

                // Update this branch's transition matrix immediately
                MG94Matrix updated_mat = build_q_matrix(opt_omega);
                P_branches[bid] = updated_mat.transition_matrix(opt_alpha);
            }

            Scalar cur_ll = compute_tree_log_likelihood(P_branches);
            if (std::abs(cur_ll - prev_ll) < 1e-3) {
                prev_ll = cur_ll;
                break;
            }
            prev_ll = cur_ll;
        }

        size_t num_branches = tree.num_nodes() - 1;
        size_t baseline_params = 5 + 9 + 2 * num_branches;
        Scalar baseline_aicc = compute_aicc(prev_ll, baseline_params);

        return {prev_ll, baseline_params, baseline_aicc};
    }

    // Optimize a single branch for M rate classes using local Inside-Outside cache
    // Returns: {best_log_l, best_omegas, best_weights, best_alpha}
    struct LocalBranchFit {
        Scalar log_l = -1e20;
        std::vector<Scalar> omegas;
        std::vector<Scalar> weights;
        Scalar alpha = 0.01;
    };

    LocalBranchFit optimize_branch_mixture(
        int32_t bid,
        int M,
        const std::vector<Vector>& V_patterns,
        const std::vector<Vector>& D_patterns,
        bool constrain_null = false
    ) const {
        size_t num_patterns = aln.patterns.size();
        Scalar cur_alpha = std::max(1e-5, branch_alpha[bid]);

        if (M == 1) {
            auto obj_1d = [&](Scalar log_a, Scalar log_w) -> Scalar {
                Scalar a = std::exp(log_a);
                Scalar w = constrain_null ? 1.0 : std::exp(log_w);
                MG94Matrix mat = build_q_matrix(w);
                Matrix P_b = mat.transition_matrix(a);

                Scalar ll = 0.0;
                for (size_t p = 0; p < num_patterns; ++p) {
                    Scalar Lp = V_patterns[p].dot(P_b * D_patterns[p]);
                    if (Lp > 0.0) {
                        ll += aln.patterns[p].weight * std::log(Lp);
                    } else {
                        ll += aln.patterns[p].weight * (-1e20);
                    }
                }
                return -ll;
            };

            LocalBranchFit res;
            if (constrain_null) {
                auto obj_a = [&](Scalar log_a) -> Scalar {
                    return obj_1d(log_a, 0.0);
                };
                Scalar init_log_a = std::log(cur_alpha);
                auto [best_log_a, best_nll] = Brent1D::minimize(obj_a, init_log_a - 1.5, init_log_a, init_log_a + 1.5, 1e-4, 50);
                res.log_l = -best_nll;
                res.omegas = {1.0};
                res.weights = {1.0};
                res.alpha = std::exp(best_log_a);
            } else {
                Scalar init_log_a = std::log(cur_alpha);
                Scalar init_log_w = std::log(std::max(1e-4, branch_omegas[bid].empty() ? 0.5 : branch_omegas[bid][0]));
                auto opt = NelderMead2D::minimize(obj_1d, init_log_a, init_log_w, -12.0, 10.0, 1e-4, 80);
                res.log_l = -opt.f;
                res.omegas = {std::exp(opt.b)};
                res.weights = {1.0};
                res.alpha = std::exp(opt.a);
            }
            return res;
        }

        // M == 2 classes:
        // Parameters:
        // x[0] = log(alpha)
        // x[1] = logit(p1): p1 = 1 / (1 + exp(-x[1]))
        // x[2] = logit(omega1): omega1 = 1 / (1 + exp(-x[2])) in [0, 1]
        // x[3] = log(omega2): omega2 = exp(x[3]) in [0, 10000]
        if (M == 2) {
            if (constrain_null) {
                // 3D optimization: log(alpha), logit(p1), logit(omega1) with omega2 = 1.0 fixed
                auto obj_3d = [&](const std::array<Scalar, 3>& x) -> Scalar {
                    Scalar a = std::exp(x[0]);
                    Scalar p1 = 1.0 / (1.0 + std::exp(-x[1]));
                    Scalar p2 = 1.0 - p1;
                    Scalar w1 = 1.0 / (1.0 + std::exp(-x[2]));
                    Scalar w2 = 1.0;

                    MG94Matrix mat1 = build_q_matrix(w1);
                    MG94Matrix mat2 = build_q_matrix(w2);
                    Matrix P_mix = p1 * mat1.transition_matrix(a) + p2 * mat2.transition_matrix(a);

                    Scalar ll = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Scalar Lp = V_patterns[p].dot(P_mix * D_patterns[p]);
                        if (Lp > 0.0) {
                            ll += aln.patterns[p].weight * std::log(Lp);
                        } else {
                            ll += aln.patterns[p].weight * (-1e20);
                        }
                    }
                    return -ll;
                };

                Scalar init_log_a = std::log(cur_alpha);
                std::array<Scalar, 3> lb = {init_log_a - 5.0, -6.0, -8.0};
                std::array<Scalar, 3> ub = {init_log_a + 6.0,  6.0,  5.0};

                // Warm-start from current branch state if available
                Scalar cur_p1 = branch_weights[bid].empty() ? 0.90 : branch_weights[bid][0];
                Scalar cur_w1 = branch_omegas[bid].empty() ? 0.20 : branch_omegas[bid][0];
                cur_p1 = std::clamp(cur_p1, 0.01, 0.99);
                cur_w1 = std::clamp(cur_w1, 1e-4, 0.99);
                std::array<Scalar, 3> best_pt = {
                    init_log_a,
                    std::log(cur_p1 / (1.0 - cur_p1)),
                    std::log(cur_w1 / (1.0 - cur_w1))
                };
                Scalar best_sc = obj_3d(best_pt);

                for (Scalar a_mult : {0.5, 1.0, 2.0}) {
                    Scalar cand_log_a = std::log(cur_alpha * a_mult);
                    for (Scalar p1 : {0.98, 0.90, 0.75, 0.50}) {
                        Scalar cand_x1 = std::log(p1 / (1.0 - p1));
                        for (Scalar w1 : {0.005, 0.05, 0.20, 0.50}) {
                            Scalar cand_x2 = std::log(w1 / (1.0 - w1));
                            std::array<Scalar, 3> pt = {cand_log_a, cand_x1, cand_x2};
                            Scalar sc = obj_3d(pt);
                            if (sc < best_sc) {
                                best_sc = sc;
                                best_pt = pt;
                            }
                        }
                    }
                }

                auto opt = NelderMeadND<3>::minimize(obj_3d, best_pt, lb, ub, 1e-4, 120);
                Scalar opt_a = std::exp(opt.x[0]);
                Scalar opt_p1 = 1.0 / (1.0 + std::exp(-opt.x[1]));
                Scalar opt_w1 = 1.0 / (1.0 + std::exp(-opt.x[2]));

                LocalBranchFit res;
                res.log_l = -opt.f;
                res.alpha = opt_a;
                res.omegas = {opt_w1, 1.0};
                res.weights = {opt_p1, 1.0 - opt_p1};
                return res;
            } else {
                // 4D optimization: log(alpha), logit(p1), logit(omega1), log(omega2)
                auto obj_4d = [&](const std::array<Scalar, 4>& x) -> Scalar {
                    Scalar a = std::exp(x[0]);
                    Scalar p1 = 1.0 / (1.0 + std::exp(-x[1]));
                    Scalar p2 = 1.0 - p1;
                    Scalar w1 = 1.0 / (1.0 + std::exp(-x[2]));
                    Scalar w2 = std::exp(x[3]);

                    MG94Matrix mat1 = build_q_matrix(w1);
                    MG94Matrix mat2 = build_q_matrix(w2);
                    Matrix P_mix = p1 * mat1.transition_matrix(a) + p2 * mat2.transition_matrix(a);

                    Scalar ll = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Scalar Lp = V_patterns[p].dot(P_mix * D_patterns[p]);
                        if (Lp > 0.0) {
                            ll += aln.patterns[p].weight * std::log(Lp);
                        } else {
                            ll += aln.patterns[p].weight * (-1e20);
                        }
                    }
                    return -ll;
                };

                Scalar init_log_a = std::log(cur_alpha);
                std::array<Scalar, 4> lb = {init_log_a - 5.0, -6.0, -8.0, -2.3};
                std::array<Scalar, 4> ub = {init_log_a + 8.0,  6.0,  6.0, 10.0};

                // Warm-start from current branch state if available
                Scalar cur_p1 = branch_weights[bid].empty() ? 0.90 : branch_weights[bid][0];
                Scalar cur_w1 = branch_omegas[bid].empty() ? 0.20 : branch_omegas[bid][0];
                Scalar cur_w2 = branch_omegas[bid].size() > 1 ? branch_omegas[bid][1] : 5.0;
                cur_p1 = std::clamp(cur_p1, 0.005, 0.995);
                cur_w1 = std::clamp(cur_w1, 1e-4, 0.995);
                cur_w2 = std::clamp(cur_w2, 0.1, 20000.0);

                std::array<Scalar, 4> best_pt = {
                    init_log_a,
                    std::log(cur_p1 / (1.0 - cur_p1)),
                    std::log(cur_w1 / (1.0 - cur_w1)),
                    std::log(cur_w2)
                };
                Scalar best_sc = obj_4d(best_pt);

                for (Scalar a_mult : {0.5, 1.0, 2.0, 6.0}) {
                    Scalar cand_log_a = std::log(cur_alpha * a_mult);
                    for (Scalar p1 : {0.98, 0.92, 0.80, 0.50}) {
                        Scalar cand_x1 = std::log(p1 / (1.0 - p1));
                        for (Scalar w1 : {0.005, 0.05, 0.20}) {
                            Scalar cand_x2 = std::log(w1 / (1.0 - w1));
                            for (Scalar w2 : {2.0, 15.0, 80.0, 500.0, 2500.0}) {
                                Scalar cand_x3 = std::log(w2);
                                std::array<Scalar, 4> pt = {cand_log_a, cand_x1, cand_x2, cand_x3};
                                Scalar sc = obj_4d(pt);
                                if (sc < best_sc) {
                                    best_sc = sc;
                                    best_pt = pt;
                                }
                            }
                        }
                    }
                }

                auto opt = NelderMeadND<4>::minimize(obj_4d, best_pt, lb, ub, 1e-4, 150);
                Scalar opt_a = std::exp(opt.x[0]);
                Scalar opt_p1 = 1.0 / (1.0 + std::exp(-opt.x[1]));
                Scalar opt_w1 = 1.0 / (1.0 + std::exp(-opt.x[2]));
                Scalar opt_w2 = std::exp(opt.x[3]);

                LocalBranchFit res;
                res.log_l = -opt.f;
                res.alpha = opt_a;
                res.omegas = {opt_w1, opt_w2};
                res.weights = {opt_p1, 1.0 - opt_p1};
                return res;
            }
        }

        // Fallback for M >= 3
        LocalBranchFit res;
        res.log_l = -1e20;
        return res;
    }

    // Run complete aBSREL analysis
    ABSRELResult run(
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        bool show_progress = false,
        bool force_progress = false
    ) {
        auto t_start = std::chrono::high_resolution_clock::now();
        ABSRELResult result;
        result.p_threshold = settings.p_threshold;

        // Phase 1: GTR
        if (progress_cb) progress_cb("Phase 1: Fitting Nucleotide GTR Model", 0.05);
        fit_gtr();
        result.gtr_fit = {gtr_log_l, 39, gtr_aicc};

        // Phase 2: Baseline MG94xREV
        if (progress_cb) progress_cb("Phase 2: Fitting Baseline MG94xREV Model", 0.20);
        auto base_summary = fit_baseline(4);
        result.baseline_fit = base_summary;

        // Record baseline results per branch
        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            ABSRELBranchResult bres;
            bres.branch_name = node.name;
            bres.node_id = node.id;
            bres.baseline_omega = branch_omegas[node.id][0];
            bres.baseline_branch_length = compute_expected_branch_length(node.id);
            bres.rate_classes = 1;
            bres.rate_distribution.rates = branch_omegas[node.id];
            bres.rate_distribution.weights = branch_weights[node.id];
            result.branches[node.name] = bres;
        }

        // Phase 3: Step-up complexity analysis
        if (progress_cb) progress_cb("Phase 3: Branch Complexity Model Selection", 0.40);
        size_t num_patterns = aln.patterns.size();
        int S = gcode->num_sense_codons;

        // Sort branches by baseline length descending
        std::vector<int32_t> sorted_branches;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                sorted_branches.push_back(node.id);
            }
        }
        std::sort(sorted_branches.begin(), sorted_branches.end(), [&](int32_t a, int32_t b) {
            return result.branches[tree.nodes[a].name].baseline_branch_length >
                   result.branches[tree.nodes[b].name].baseline_branch_length;
        });

        Scalar current_best_aicc = base_summary.aicc;
        size_t current_params = base_summary.parameters;

        std::vector<Matrix> P_branches = compute_all_branch_transition_matrices();

        std::unique_ptr<ProgressBar> pbar_p3;
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar_p3 = std::make_unique<ProgressBar>(
                sorted_branches.size(),
                "[aBSREL] Complexity Selection",
                "branches",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        for (size_t b_idx = 0; b_idx < sorted_branches.size(); ++b_idx) {
            int32_t bid = sorted_branches[b_idx];
            const std::string& bname = tree.nodes[bid].name;

            // Compute local (V, D) context for branch bid
            std::vector<Vector> V_patterns(num_patterns);
            std::vector<Vector> D_patterns(num_patterns);
            #pragma omp parallel for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                auto io = LikelihoodEngine::compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_branches, codon_freqs, S
                );
                V_patterns[p] = io.V[bid];
                D_patterns[p] = io.D[bid];
            }

            int current_classes = 1;
            LocalBranchFit best_fit;
            best_fit.omegas = branch_omegas[bid];
            best_fit.weights = branch_weights[bid];
            best_fit.alpha = branch_alpha[bid];

            for (int M = 2; M <= settings.max_rate_classes; ++M) {
                LocalBranchFit cand_fit = optimize_branch_mixture(bid, M, V_patterns, D_patterns, false);
                size_t cand_params = current_params + 2;
                Scalar cand_aicc = compute_aicc(cand_fit.log_l, cand_params);

                if (cand_aicc < current_best_aicc) {
                    current_best_aicc = cand_aicc;
                    current_params = cand_params;
                    current_classes = M;
                    best_fit = cand_fit;
                } else {
                    break;
                }
            }

            if (current_classes > 1) {
                branch_alpha[bid] = best_fit.alpha;
                branch_omegas[bid] = best_fit.omegas;
                branch_weights[bid] = best_fit.weights;
                P_branches[bid] = compute_branch_transition_matrix(bid);
                result.branches[bname].rate_classes = current_classes;
                result.branches[bname].rate_distribution.rates = best_fit.omegas;
                result.branches[bname].rate_distribution.weights = best_fit.weights;
            }

            if (pbar_p3) {
                std::string stat = bname + ": " + std::to_string(current_classes) + " rate class(es)";
                pbar_p3->set_status(stat);
                pbar_p3->tick();
            }
        }

        if (pbar_p3) {
            size_t multi_rate_count = 0;
            for (const auto& pair : result.branches) {
                if (pair.second.rate_classes > 1) multi_rate_count++;
            }
            std::ostringstream summary;
            summary << "Model selection complete (" << multi_rate_count << " branch(es) with >1 rate class)";
            pbar_p3->finish(summary.str());
        }

        // Phase 4: Full Adaptive Model Fitting
        if (progress_cb) progress_cb("Phase 4: Full Adaptive Model Fitting", 0.65);

        // Joint full model optimization:
        // 1. Joint L-BFGS on all branch lengths using analytical Inside-Outside gradients
        // 2. GTR nucleotide exchangeability rate optimization (Brent 1D on 5 rates)
        // 3. Branch mixture parameter refinement (Nelder-Mead on local mixtures)
        std::vector<int32_t> all_branches;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                all_branches.push_back(node.id);
            }
        }
        size_t n_b = all_branches.size();

        auto run_joint_branch_length_lbfgs = [&]() {
            struct BranchObjective {
                ABSRELAnalyzer& abs;
                const std::vector<int32_t>& nodes;
                int S_val;
                size_t n_pat;

                Scalar operator()(const Vector& y, Vector& grad) {
                    size_t n = nodes.size();
                    for (size_t i = 0; i < n; ++i) {
                        abs.branch_alpha[nodes[i]] = std::exp(y(i));
                    }

                    auto P_b = abs.compute_all_branch_transition_matrices();
                    Scalar total_ll = 0.0;

                    size_t num_nodes = abs.tree.num_nodes();
                    std::vector<Matrix> dP_branches(num_nodes, Matrix::Zero(S_val, S_val));
                    for (int32_t bid : nodes) {
                        Scalar alpha = abs.branch_alpha[bid];
                        for (size_t k = 0; k < abs.branch_omegas[bid].size(); ++k) {
                            Scalar w = abs.branch_omegas[bid][k];
                            Scalar p = abs.branch_weights[bid][k];
                            MG94Matrix mat = abs.build_q_matrix(w);
                            Matrix P_k = mat.transition_matrix(alpha);
                            dP_branches[bid] += p * (mat.Q * P_k);
                        }
                    }

                    std::vector<Scalar> dL_da(num_nodes, 0.0);

                    #pragma omp parallel
                    {
                        Scalar local_ll = 0.0;
                        std::vector<Scalar> local_dL(num_nodes, 0.0);

                        #pragma omp for schedule(dynamic)
                        for (size_t p = 0; p < n_pat; ++p) {
                            const auto& pattern = abs.aln.patterns[p];
                            auto io = LikelihoodEngine::compute_inside_outside(
                                abs.tree, pattern, abs.leaf_to_taxon, P_b, abs.codon_freqs, S_val
                            );
                            Scalar L = io.likelihood;
                            if (L > 0.0) {
                                local_ll += pattern.weight * std::log(L);
                                Scalar inv_L = pattern.weight / L;
                                for (int32_t bid : nodes) {
                                    Scalar term = io.V[bid].dot(dP_branches[bid] * io.D[bid]);
                                    local_dL[bid] += inv_L * term;
                                }
                            }
                        }

                        #pragma omp critical
                        {
                            total_ll += local_ll;
                            for (int32_t bid : nodes) {
                                dL_da[bid] += local_dL[bid];
                            }
                        }
                    }

                    grad.resize(n);
                    for (size_t i = 0; i < n; ++i) {
                        int32_t bid = nodes[i];
                        Scalar alpha = abs.branch_alpha[bid];
                        grad(i) = -alpha * dL_da[bid];
                    }

                    return -total_ll;
                }
            };

            BranchObjective b_obj{*this, all_branches, S, num_patterns};
            Vector y(n_b), lb(n_b), ub(n_b);
            for (size_t i = 0; i < n_b; ++i) {
                y(i) = std::log(std::clamp(branch_alpha[all_branches[i]], 1e-6, 50.0));
                lb(i) = -14.0;
                ub(i) =   5.0;
            }

            LBFGSpp::LBFGSBParam<Scalar> opt_param;
            opt_param.m = 6;
            opt_param.epsilon = 1e-4;
            opt_param.max_iterations = 35;

            LBFGSpp::LBFGSBSolver<Scalar> solver(opt_param);
            Scalar fx = 0.0;
            try {
                solver.minimize(b_obj, y, fx, lb, ub);
            } catch (...) {}
            return -fx;
        };

        auto run_gtr_optimization = [&]() {
            auto opt_rate = [&](Scalar& r_ref) {
                auto rate_obj = [&](Scalar val) -> Scalar {
                    Scalar old_val = r_ref;
                    r_ref = val;
                    compute_syn_scale();
                    auto P_b = compute_all_branch_transition_matrices();
                    Scalar lnl = compute_tree_log_likelihood(P_b);
                    r_ref = old_val;
                    return -lnl;
                };
                auto [best_r, _] = Brent1D::minimize(rate_obj, 1e-4, r_ref, 20.0, 1e-3, 15);
                r_ref = best_r;
                compute_syn_scale();
            };

            opt_rate(gtr_params.theta_AC);
            opt_rate(gtr_params.theta_AT);
            opt_rate(gtr_params.theta_CG);
            opt_rate(gtr_params.theta_CT);
            opt_rate(gtr_params.theta_GT);

            auto P_b = compute_all_branch_transition_matrices();
            return compute_tree_log_likelihood(P_b);
        };

        auto run_mixture_refinement = [&]() {
            auto P_b = compute_all_branch_transition_matrices();
            for (int32_t bid : all_branches) {
                int M = branch_omegas[bid].size();
                std::vector<Vector> V_patterns(num_patterns);
                std::vector<Vector> D_patterns(num_patterns);
                #pragma omp parallel for schedule(dynamic)
                for (size_t p = 0; p < num_patterns; ++p) {
                    const auto& pattern = aln.patterns[p];
                    auto io = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_b, codon_freqs, S
                    );
                    V_patterns[p] = io.V[bid];
                    D_patterns[p] = io.D[bid];
                }

                auto refined_fit = optimize_branch_mixture(bid, M, V_patterns, D_patterns, false);
                branch_alpha[bid] = refined_fit.alpha;
                branch_omegas[bid] = refined_fit.omegas;
                branch_weights[bid] = refined_fit.weights;
                P_b[bid] = compute_branch_transition_matrix(bid);

                result.branches[tree.nodes[bid].name].rate_distribution.rates = refined_fit.omegas;
                result.branches[tree.nodes[bid].name].rate_distribution.weights = refined_fit.weights;
            }
            return compute_tree_log_likelihood(P_b);
        };

        Scalar prev_full_ll = compute_tree_log_likelihood(P_branches);
        const int max_joint_passes = 6;
        std::unique_ptr<ProgressBar> pbar_p4;
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar_p4 = std::make_unique<ProgressBar>(
                max_joint_passes,
                "[aBSREL] Full Model Refinement",
                "passes",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        for (int iter = 0; iter < max_joint_passes; ++iter) {
            run_joint_branch_length_lbfgs();
            run_gtr_optimization();
            Scalar cur_full_ll = run_mixture_refinement();

            if (settings.verbose) {
                std::cout << "  [Phase 4 Joint Pass " << iter << "] Log(L) = " << cur_full_ll << " (delta = " << cur_full_ll - prev_full_ll << ")\n";
            }
            if (pbar_p4) {
                std::ostringstream stat;
                stat << "Pass " << (iter + 1) << "/" << max_joint_passes << ": lnL = "
                     << std::fixed << std::setprecision(2) << cur_full_ll;
                pbar_p4->set_status(stat.str());
                pbar_p4->tick();
            }
            if (std::abs(cur_full_ll - prev_full_ll) < 0.05) {
                break;
            }
            prev_full_ll = cur_full_ll;
        }

        P_branches = compute_all_branch_transition_matrices();
        Scalar full_ll = compute_tree_log_likelihood(P_branches);
        Scalar full_aicc = compute_aicc(full_ll, current_params);
        result.full_adaptive_fit = {full_ll, current_params, full_aicc};

        if (pbar_p4) {
            std::ostringstream summary;
            summary << "Refinement complete: lnL = " << std::fixed << std::setprecision(2) << full_ll
                    << " (AICc = " << full_aicc << ")";
            pbar_p4->finish(summary.str());
        }

        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            int32_t bid = node.id;
            const std::string& bname = node.name;
            result.branches[bname].full_branch_length = compute_expected_branch_length(bid);
            auto [es, en] = compute_es_en(bid);
            result.branches[bname].full_es = es;
            result.branches[bname].full_en = en;
        }

        // Phase 5: Testing Selected Branches for Positive Selection
        if (progress_cb) progress_cb("Phase 5: Hypothesis Testing for Selection", 0.80);
        std::vector<std::string> to_test;
        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            bool test_this = true;
            if (settings.test_branches == "Internal" && node.is_leaf) test_this = false;
            if (settings.test_branches == "Leaves" && !node.is_leaf) test_this = false;
            if (test_this) {
                to_test.push_back(node.name);
            } else {
                result.branches[node.name].is_tested = false;
            }
        }
        result.tested_branches = to_test;

        std::unique_ptr<ProgressBar> pbar_p5;
        std::atomic<size_t> cand_pos_count{0};
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar_p5 = std::make_unique<ProgressBar>(
                to_test.size(),
                "[aBSREL] Selection Testing",
                "branches",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        for (const std::string& bname : to_test) {
            auto& bres = result.branches[bname];
            int32_t bid = bres.node_id;

            Scalar max_omega = *std::max_element(branch_omegas[bid].begin(), branch_omegas[bid].end());
            if (max_omega <= 1.0) {
                bres.lrt = 0.0;
                bres.uncorrected_p_value = 1.0;
                bres.sites_ebf_100 = 0;
                if (pbar_p5) {
                    pbar_p5->tick();
                }
                continue;
            }

            // Compute local (V, D) context under Full model
            std::vector<Vector> V_patterns(num_patterns);
            std::vector<Vector> D_patterns(num_patterns);
            #pragma omp parallel for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                auto io = LikelihoodEngine::compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_branches, codon_freqs, S
                );
                V_patterns[p] = io.V[bid];
                D_patterns[p] = io.D[bid];
            }

            // Optimize under null constraint omega_max = 1.0
            LocalBranchFit null_fit = optimize_branch_mixture(
                bid, bres.rate_classes, V_patterns, D_patterns, true
            );

            Scalar lrt = 2.0 * (full_ll - null_fit.log_l);
            if (lrt < 0.0) lrt = 0.0;
            bres.lrt = lrt;
            bres.uncorrected_p_value = compute_mixture_p_value(lrt);

            // Compute Empirical Bayes Factors per site
            size_t K = bres.rate_classes;
            Scalar p_pos = bres.rate_distribution.weights.back();
            Scalar prior_odds = (p_pos > 0.0 && p_pos < 1.0) ? (p_pos / (1.0 - p_pos)) : 1.0;

            int sites_ebf_count = 0;
            bres.site_ebf.assign(aln.num_codons, 1.0);

            // Compute pattern EBFs
            std::vector<Scalar> pattern_ebf(num_patterns, 1.0);
            Scalar alpha = branch_alpha[bid];
            std::vector<Matrix> P_comps(K);
            for (size_t k = 0; k < K; ++k) {
                MG94Matrix mat = build_q_matrix(bres.rate_distribution.rates[k]);
                P_comps[k] = mat.transition_matrix(alpha);
            }

            for (size_t p = 0; p < num_patterns; ++p) {
                std::vector<Scalar> comp_L(K);
                Scalar denom = 0.0;
                for (size_t k = 0; k < K; ++k) {
                    Scalar L_k = V_patterns[p].dot(P_comps[k] * D_patterns[p]);
                    comp_L[k] = bres.rate_distribution.weights[k] * L_k;
                    denom += comp_L[k];
                }
                if (denom > 1e-300 && p_pos < 1.0) {
                    Scalar post_pos = comp_L.back() / denom;
                    if (post_pos < 1.0) {
                        Scalar post_odds = post_pos / (1.0 - post_pos);
                        pattern_ebf[p] = post_odds / prior_odds;
                    } else {
                        pattern_ebf[p] = 1e6;
                    }
                }
            }

            // Map pattern EBF to site EBF
            for (size_t s = 0; s < aln.num_codons; ++s) {
                size_t p = aln.site_to_pattern[s];
                Scalar ebf = pattern_ebf[p];
                bres.site_ebf[s] = ebf;
                if (ebf >= 100.0) {
                    sites_ebf_count++;
                }
            }
            bres.sites_ebf_100 = sites_ebf_count;

            if (bres.uncorrected_p_value <= settings.p_threshold) {
                cand_pos_count.fetch_add(1, std::memory_order_relaxed);
            }
            if (pbar_p5) {
                size_t cpos = cand_pos_count.load(std::memory_order_relaxed);
                if (cpos > 0) {
                    pbar_p5->set_status("\033[1;32m+" + std::to_string(cpos) + "\033[0m candidate positive");
                }
                pbar_p5->tick();
            }
        }

        // Phase 6: Holm-Bonferroni correction
        if (progress_cb) progress_cb("Phase 6: Multiple Testing Correction", 0.95);
        size_t m = to_test.size();
        std::vector<size_t> p_indices(m);
        std::iota(p_indices.begin(), p_indices.end(), 0);
        std::sort(p_indices.begin(), p_indices.end(), [&](size_t i, size_t j) {
            return result.branches[to_test[i]].uncorrected_p_value <
                   result.branches[to_test[j]].uncorrected_p_value;
        });

        Scalar running_max = 0.0;
        for (size_t rank = 0; rank < m; ++rank) {
            size_t idx = p_indices[rank];
            const std::string& bname = to_test[idx];
            Scalar raw_p = result.branches[bname].uncorrected_p_value;
            Scalar adj_p = (m - rank) * raw_p;
            running_max = std::max(running_max, adj_p);
            Scalar final_p = std::min(1.0, running_max);
            result.branches[bname].corrected_p_value = final_p;
            if (final_p <= settings.p_threshold) {
                result.branches[bname].is_positive = true;
                result.positive_branches.push_back(bname);
            }
        }

        if (pbar_p5) {
            std::ostringstream summary;
            summary << "\033[1;32m" << result.positive_branches.size()
                    << " branch(es) under episodic diversifying selection\033[0m (p ≤ "
                    << std::fixed << std::setprecision(2) << settings.p_threshold << ")";
            pbar_p5->finish(summary.str());
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        result.runtime_seconds = std::chrono::duration<double>(t_end - t_start).count();
        if (progress_cb) progress_cb("Analysis Complete", 1.0);

        return result;
    }
};

} // namespace hyphy::analyses
