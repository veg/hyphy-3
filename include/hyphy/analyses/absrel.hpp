#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "nlohmann/json.hpp"

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

    nlohmann::json to_json(const Tree& tree, const Alignment& aln) const {
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
        size_t num_patterns = aln.patterns.size();
        Scalar total_log_l = 0.0;

        #pragma omp parallel for reduction(+:total_log_l) schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            const auto& pattern = aln.patterns[p];
            auto io = LikelihoodEngine::compute_inside_outside(
                tree, pattern, leaf_to_taxon, P_branches, codon_freqs, S
            );
            if (io.likelihood > 0.0) {
                total_log_l += pattern.weight * std::log(io.likelihood);
            } else {
                total_log_l += pattern.weight * (-1e20);
            }
        }
        return total_log_l;
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
        const auto& omegas = branch_omegas[node_id];
        const auto& weights = branch_weights[node_id];
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
        Scalar cur_alpha = branch_alpha[bid];

        if (M == 1) {
            auto obj_1d = [&](Scalar log_w) -> Scalar {
                Scalar w = constrain_null ? 1.0 : std::exp(log_w);
                MG94Matrix mat = build_q_matrix(w);
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

            LocalBranchFit res;
            if (constrain_null) {
                res.log_l = -obj_1d(0.0);
                res.omegas = {1.0};
                res.weights = {1.0};
                res.alpha = cur_alpha;
            } else {
                Scalar best_obj = 1e20;
                Scalar best_log_w = 0.0;
                for (Scalar cand_w : {0.05, 0.2, 0.5, 1.0, 2.5, 5.0}) {
                    Scalar val = obj_1d(std::log(cand_w));
                    if (val < best_obj) {
                        best_obj = val;
                        best_log_w = std::log(cand_w);
                    }
                }
                auto [best_opt_w, best_opt_nll] = Brent1D::minimize(obj_1d, best_log_w - 1.5, best_log_w, best_log_w + 1.5, 1e-4, 50);
                res.log_l = -best_opt_nll;
                res.omegas = {std::exp(best_opt_w)};
                res.weights = {1.0};
                res.alpha = cur_alpha;
            }
            return res;
        }

        // M == 2 classes: omega1 in [0, 1], omega2 in [0, 5000] (or 1.0 if constrained)
        // weight p1 in (0, 1), p2 = 1 - p1
        if (M == 2) {
            auto obj_2class = [&](Scalar x0, Scalar x1, Scalar x2) -> Scalar {
                // x0: transformed p1 via logit: p1 = sigmoid(x0)
                // x1: transformed omega1: omega1 = sigmoid(x1) in [0, 1]
                // x2: log(omega2)
                Scalar p1 = 1.0 / (1.0 + std::exp(-x0));
                Scalar p2 = 1.0 - p1;
                Scalar w1 = 1.0 / (1.0 + std::exp(-x1));
                Scalar w2 = constrain_null ? 1.0 : std::exp(x2);

                MG94Matrix mat1 = build_q_matrix(w1);
                MG94Matrix mat2 = build_q_matrix(w2);
                Matrix P_mix = p1 * mat1.transition_matrix(cur_alpha) + p2 * mat2.transition_matrix(cur_alpha);

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

            // Grid search for initialization
            Scalar best_grid_score = 1e20;
            Vector best_x = Vector::Zero(3);

            std::vector<Scalar> grid_p1 = {0.95, 0.85, 0.65};
            std::vector<Scalar> grid_w1 = {0.05, 0.2, 0.5};
            std::vector<Scalar> grid_w2 = constrain_null ? std::vector<Scalar>{1.0} : std::vector<Scalar>{1.5, 5.0, 20.0, 100.0};

            for (Scalar p1 : grid_p1) {
                for (Scalar w1 : grid_w1) {
                    for (Scalar w2 : grid_w2) {
                        Scalar cand_x0 = std::log(p1 / (1.0 - p1));
                        Scalar cand_x1 = std::log(w1 / (1.0 - w1));
                        Scalar cand_x2 = std::log(w2);

                        Scalar sc = obj_2class(cand_x0, cand_x1, cand_x2);
                        if (sc < best_grid_score) {
                            best_grid_score = sc;
                            best_x(0) = cand_x0;
                            best_x(1) = cand_x1;
                            best_x(2) = cand_x2;
                        }
                    }
                }
            }

            // Optimize with Nelder-Mead
            if (constrain_null) {
                // 2D Nelder-Mead on x(0) and x(1)
                auto obj_null_2d = [&](Scalar x0, Scalar x1) -> Scalar {
                    return obj_2class(x0, x1, 0.0);
                };
                auto opt = NelderMead2D::minimize(obj_null_2d, best_x(0), best_x(1), -8.0, 8.0, 1e-4, 100);
                Scalar opt_p1 = 1.0 / (1.0 + std::exp(-opt.a));
                Scalar opt_w1 = 1.0 / (1.0 + std::exp(-opt.b));

                LocalBranchFit res;
                res.log_l = -opt.f;
                res.omegas = {opt_w1, 1.0};
                res.weights = {opt_p1, 1.0 - opt_p1};
                res.alpha = cur_alpha;
                return res;
            } else {
                // 2D Nelder-Mead on x(0) and x(2)
                auto obj_free_2d = [&](Scalar x0, Scalar x2) -> Scalar {
                    return obj_2class(x0, best_x(1), x2);
                };
                auto opt = NelderMead2D::minimize(obj_free_2d, best_x(0), best_x(2), -8.0, 10.0, 1e-4, 120);
                Scalar opt_p1 = 1.0 / (1.0 + std::exp(-opt.a));
                Scalar opt_w1 = 1.0 / (1.0 + std::exp(-best_x(1)));
                Scalar opt_w2 = std::exp(opt.b);

                LocalBranchFit res;
                res.log_l = -opt.f;
                res.omegas = {opt_w1, opt_w2};
                res.weights = {opt_p1, 1.0 - opt_p1};
                res.alpha = cur_alpha;
                return res;
            }
        }

        // Fallback for M >= 3
        LocalBranchFit res;
        res.log_l = -1e20;
        return res;
    }

    // Run complete aBSREL analysis
    ABSRELResult run(std::function<void(const std::string&, double)> progress_cb = nullptr) {
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
                branch_omegas[bid] = best_fit.omegas;
                branch_weights[bid] = best_fit.weights;
                P_branches[bid] = compute_branch_transition_matrix(bid);
                result.branches[bname].rate_classes = current_classes;
                result.branches[bname].rate_distribution.rates = best_fit.omegas;
                result.branches[bname].rate_distribution.weights = best_fit.weights;
            }
        }

        // Phase 4: Full Adaptive Model Fitting
        if (progress_cb) progress_cb("Phase 4: Full Adaptive Model Fitting", 0.65);
        Scalar full_ll = compute_tree_log_likelihood(P_branches);
        result.full_adaptive_fit = {full_ll, current_params, current_best_aicc};

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

        for (const std::string& bname : to_test) {
            auto& bres = result.branches[bname];
            int32_t bid = bres.node_id;

            Scalar max_omega = *std::max_element(branch_omegas[bid].begin(), branch_omegas[bid].end());
            if (max_omega <= 1.0) {
                bres.lrt = 0.0;
                bres.uncorrected_p_value = 1.0;
                bres.sites_ebf_100 = 0;
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

        auto t_end = std::chrono::high_resolution_clock::now();
        result.runtime_seconds = std::chrono::duration<double>(t_end - t_start).count();
        if (progress_cb) progress_cb("Analysis Complete", 1.0);

        return result;
    }
};

} // namespace hyphy::analyses
