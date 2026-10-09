#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/squarem.hpp"
#include "nlohmann/json.hpp"

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

struct BUSTEDRateDistribution {
    std::vector<Scalar> omegas;
    std::vector<Scalar> weights;
    std::vector<Scalar> syn_rates;
    std::vector<Scalar> syn_weights;
    std::vector<std::string> annotations; // e.g. "Collapsed rate class"
};

struct BUSTEDFit {
    Scalar log_likelihood = 0.0;
    Scalar aicc = 0.0;
    Scalar tree_scale = 1.0;
    size_t num_rate_classes = 3;
    BUSTEDRateDistribution test_distribution;
    std::vector<Scalar> branch_lengths; // Branch lengths per node id
    std::vector<Scalar> site_log_likelihoods;
};

struct BUSTEDSettings {
    size_t num_rate_classes = 3;       // Default K = 3
    bool srv = false;                  // Enable synonymous rate variation (BUSTED-S)
    size_t num_syn_rate_classes = 3;   // Default M = 3 for alpha
    bool auto_select_k = false;        // Automatically select K in {1, ..., max_k} via AICc
    size_t max_k = 3;                  // Max K to test in auto_select_k
    bool refine_branch_lengths = true; // Re-optimize individual branch lengths
    Scalar p_value_threshold = 0.05;
};

struct BUSTEDResult {
    size_t optimal_k = 3;
    std::vector<BUSTEDFit> k_null_fits; // Null model fits for K=1, 2, ...
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
    BUSTEDFit unconstrained;
    BUSTEDFit constrained;
    std::vector<Scalar> evidence_ratios; // Bayes factors per site for omega_K > 1
    double runtime_seconds = 0.0;
};

class BUSTEDAnalyzer {
public:
    Tree tree;
    Alignment aln;
    std::shared_ptr<const GeneticCode> code;
    MG94Parameters base_params;

    // Relative branch lengths tau_b (normalized from MG94 baseline)
    std::vector<Scalar> branch_tau;
    std::vector<bool> is_test_branch;
    size_t num_test_branches = 0;

    // GTR baseline results
    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    // Global MG94 baseline results
    Scalar mg94_log_l = 0.0;
    Scalar mg94_omega = 1.0;
    Scalar mg94_tree_length = 0.0;

    BUSTEDAnalyzer(
        Tree t,
        Alignment a,
        MG94Parameters base_p = MG94Parameters{}
    ) : tree(std::move(t)), aln(std::move(a)), base_params(base_p) {
        code = aln.code ? aln.code : GeneticCode::universal();
        size_t num_nodes = tree.num_nodes();
        branch_tau.assign(num_nodes, 0.0);
        is_test_branch.assign(num_nodes, true); // Default: all branches are Test
        num_test_branches = num_nodes - 1; // Non-root branches
        prepare_baseline_lengths();
    }

    static BUSTEDAnalyzer create_and_fit(
        Tree raw_tree,
        Alignment aln
    ) {
        GTRFitter gtr_fitter(raw_tree, aln);
        auto gtr_res = gtr_fitter.fit();

        MG94Parameters base_p;
        base_p.theta_AC = gtr_res.params.theta_AC;
        base_p.theta_AT = gtr_res.params.theta_AT;
        base_p.theta_CG = gtr_res.params.theta_CG;
        base_p.theta_CT = gtr_res.params.theta_CT;
        base_p.theta_GT = gtr_res.params.theta_GT;

        MG94Fitter mg_fitter(gtr_res.tree, aln);
        auto mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
        base_p.alpha = 1.0;
        base_p.beta = mg_res.x_opt(0);

        BUSTEDAnalyzer analyzer(mg_res.tree, std::move(aln), base_p);
        analyzer.has_gtr_fit = true;
        analyzer.gtr_log_l = gtr_res.log_likelihood;
        analyzer.gtr_aicc = gtr_res.aicc;
        analyzer.gtr_rates = gtr_res.params;
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                analyzer.gtr_branch_lengths[node.name] = node.branch_length;
            }
        }

        analyzer.mg94_log_l = mg_res.log_likelihood;
        analyzer.mg94_omega = mg_res.x_opt(0);
        return analyzer;
    }

    void prepare_baseline_lengths() {
        MG94Matrix mat;
        mat.update(base_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        Scalar Q_scale = mat.scale_factor;
        Scalar conversion_factor = (Q_scale > 1e-12) ? (3.0 / Q_scale) : 1.0;

        mg94_tree_length = 0.0;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar scaled_len = node.branch_length * conversion_factor;
                branch_tau[node.id] = scaled_len;
                mg94_tree_length += scaled_len;
            }
        }
    }

    // Evaluates tree log-likelihood for arbitrary K rate categories, optional SRV (M classes), and custom branch lengths
    Scalar evaluate_log_l(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        std::vector<Scalar>* out_pattern_ll = nullptr
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();
        size_t M = syn_rates.size();

        // Build K rate matrices
        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        // Precompute branch transition matrices P_mix[m][node_id]
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[node.id] : (s * branch_tau[node.id]);
                    P_mix[m][node.id] = Matrix::Zero(S, S);
                    for (size_t k = 0; k < K; ++k) {
                        P_mix[m][node.id] += weights[k] * mg[k].transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        std::vector<size_t> leaf_to_taxon(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        if (out_pattern_ll) {
            out_pattern_ll->assign(num_patterns, 0.0);
        }

        Scalar total_log_l = 0.0;
        #pragma omp parallel
        {
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));
            #pragma omp for reduction(+:total_log_l) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pat_L = 0.0;

                for (size_t m = 0; m < M; ++m) {
                    for (int32_t node_id : tree.post_order) {
                        const auto& node = tree.nodes[node_id];
                        if (node.is_leaf) {
                            size_t t_idx = leaf_to_taxon[node_id];
                            if (t_idx != static_cast<size_t>(-1)) {
                                int8_t state = pattern.states[t_idx];
                                if (state >= 0 && state < S) {
                                    node_L[node_id].setZero();
                                    node_L[node_id](state) = 1.0;
                                } else {
                                    node_L[node_id].setOnes();
                                }
                            } else {
                                node_L[node_id].setOnes();
                            }
                        } else {
                            node_L[node_id].setOnes();
                            for (int32_t child_id : node.children) {
                                Vector child_msg = P_mix[m][child_id] * node_L[child_id];
                                node_L[node_id] = node_L[node_id].cwiseProduct(child_msg);
                            }
                        }
                    }
                    Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                    pat_L += syn_weights[m] * cond_L;
                }

                Scalar pll = (pat_L > 0.0) ? std::log(pat_L) : -1e20;
                total_log_l += pattern.weight * pll;
                if (out_pattern_ll) {
                    (*out_pattern_ll)[p] = pll;
                }
            }
        }
        return total_log_l;
    }

    // Expectation Step: Computes new mixture weights via inside-outside traversal for arbitrary K >= 1 and optional SRV
    std::vector<Scalar> compute_expected_weights(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr
    ) const {
        size_t K = omegas.size();
        if (K <= 1) return {1.0};

        size_t M = syn_rates.size();
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;

        std::vector<std::vector<std::vector<Matrix>>> P_comp(M, std::vector<std::vector<Matrix>>(K, std::vector<Matrix>(num_nodes)));
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            MG94Matrix mg;
            mg.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            for (size_t m = 0; m < M; ++m) {
                Scalar alpha_m = syn_rates[m];
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) {
                        Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[node.id] : (s * branch_tau[node.id]);
                        P_comp[m][k][node.id] = mg.transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        for (size_t m = 0; m < M; ++m) {
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    P_mix[m][node.id] = Matrix::Zero(S, S);
                    for (size_t k = 0; k < K; ++k) {
                        P_mix[m][node.id] += weights[k] * P_comp[m][k][node.id];
                    }
                }
            }
        }

        std::vector<size_t> leaf_to_taxon(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        Vector weight_accum = Vector::Zero(K);
        Scalar total_pattern_weight = 0.0;

        #pragma omp parallel
        {
            Vector local_accum = Vector::Zero(K);
            Scalar local_pat_wt = 0.0;
            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pat_L = 0.0;
                std::vector<LikelihoodEngine::InsideOutsideResult> io_m(M);

                for (size_t m = 0; m < M; ++m) {
                    io_m[m] = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S
                    );
                    pat_L += syn_weights[m] * io_m[m].likelihood;
                }

                if (pat_L > 0.0) {
                    Scalar inv_L = pattern.weight / pat_L;
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id && is_test_branch[node.id]) {
                            int nid = node.id;
                            for (size_t k = 0; k < K; ++k) {
                                Scalar num = 0.0;
                                for (size_t m = 0; m < M; ++m) {
                                    num += syn_weights[m] * io_m[m].V[nid].dot(P_comp[m][k][nid] * io_m[m].D[nid]);
                                }
                                local_accum(k) += weights[k] * inv_L * num;
                            }
                            local_pat_wt += pattern.weight;
                        }
                    }
                }
            }
            #pragma omp critical
            {
                weight_accum += local_accum;
                total_pattern_weight += local_pat_wt;
            }
        }

        std::vector<Scalar> new_weights(K);
        if (total_pattern_weight > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t k = 0; k < K; ++k) {
                new_weights[k] = std::clamp(weight_accum(k) / total_pattern_weight, 1e-6, 1.0 - 1e-6);
                sum_w += new_weights[k];
            }
            for (size_t k = 0; k < K; ++k) {
                new_weights[k] /= sum_w;
            }
        } else {
            new_weights = weights;
        }
        return new_weights;
    }

    // Expectation Step: Computes new synonymous rate weights via site posteriors
    std::vector<Scalar> compute_expected_syn_weights(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates,
        const std::vector<Scalar>& syn_weights,
        const std::vector<Scalar>* custom_branch_lengths = nullptr
    ) const {
        size_t M = syn_rates.size();
        if (M <= 1) return {1.0};

        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();

        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[node.id] : (s * branch_tau[node.id]);
                    P_mix[m][node.id] = Matrix::Zero(S, S);
                    for (size_t k = 0; k < K; ++k) {
                        P_mix[m][node.id] += weights[k] * mg[k].transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        std::vector<size_t> leaf_to_taxon(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        Vector syn_accum = Vector::Zero(M);
        Scalar total_pattern_weight = 0.0;

        #pragma omp parallel
        {
            Vector local_accum = Vector::Zero(M);
            Scalar local_pat_wt = 0.0;
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pat_L = 0.0;
                std::vector<Scalar> cond_L(M, 0.0);

                for (size_t m = 0; m < M; ++m) {
                    for (int32_t node_id : tree.post_order) {
                        const auto& node = tree.nodes[node_id];
                        if (node.is_leaf) {
                            size_t t_idx = leaf_to_taxon[node_id];
                            if (t_idx != static_cast<size_t>(-1)) {
                                int8_t state = pattern.states[t_idx];
                                if (state >= 0 && state < S) {
                                    node_L[node_id].setZero();
                                    node_L[node_id](state) = 1.0;
                                } else {
                                    node_L[node_id].setOnes();
                                }
                            } else {
                                node_L[node_id].setOnes();
                            }
                        } else {
                            node_L[node_id].setOnes();
                            for (int32_t child_id : node.children) {
                                Vector child_msg = P_mix[m][child_id] * node_L[child_id];
                                node_L[node_id] = node_L[node_id].cwiseProduct(child_msg);
                            }
                        }
                    }
                    cond_L[m] = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                    pat_L += syn_weights[m] * cond_L[m];
                }

                if (pat_L > 0.0) {
                    for (size_t m = 0; m < M; ++m) {
                        Scalar rho_m = (syn_weights[m] * cond_L[m]) / pat_L;
                        local_accum(m) += pattern.weight * rho_m;
                    }
                    local_pat_wt += pattern.weight;
                }
            }

            #pragma omp critical
            {
                syn_accum += local_accum;
                total_pattern_weight += local_pat_wt;
            }
        }

        std::vector<Scalar> new_syn_weights(M);
        if (total_pattern_weight > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t m = 0; m < M; ++m) {
                new_syn_weights[m] = std::clamp(syn_accum(m) / total_pattern_weight, 1e-6, 1.0 - 1e-6);
                sum_w += new_syn_weights[m];
            }
            for (size_t m = 0; m < M; ++m) {
                new_syn_weights[m] /= sum_w;
            }
        } else {
            new_syn_weights = syn_weights;
        }
        return new_syn_weights;
    }

    // Evaluates log-likelihood and exact analytical branch gradients d ln L / d t_b
    // under the BUSTED mixture model using Inside-Outside adjoints (with optional SRV)
    std::pair<Scalar, std::vector<Scalar>> compute_busted_branch_gradients(
        const std::vector<Scalar>& branch_lengths,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0}
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();
        size_t M = syn_rates.size();

        // 1. Build K rate matrices
        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        // 2. Precompute P_mix[m] and dP_mix[m] for each synonymous rate category and branch
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<Matrix>> dP_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = branch_lengths[nid];
                    P_mix[m][nid] = Matrix::Zero(S, S);
                    dP_mix[m][nid] = Matrix::Zero(S, S);
                    for (size_t k = 0; k < K; ++k) {
                        Matrix P_k = mg[k].transition_matrix(alpha_m * bl);
                        P_mix[m][nid] += weights[k] * P_k;
                        dP_mix[m][nid] += weights[k] * alpha_m * (mg[k].Q * P_k);
                    }
                }
            }
        }

        // 3. Leaf mapping
        std::vector<size_t> leaf_to_taxon(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        // 4. Parallel inside-outside traversal
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
                std::vector<LikelihoodEngine::InsideOutsideResult> io_m(M);

                for (size_t m = 0; m < M; ++m) {
                    io_m[m] = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S
                    );
                    pat_L += syn_weights[m] * io_m[m].likelihood;
                }

                if (pat_L > 0.0) {
                    local_ll += pattern.weight * std::log(pat_L);
                    Scalar inv_L = pattern.weight / pat_L;
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t vid = node.id;
                            Scalar dL_dt = 0.0;
                            for (size_t m = 0; m < M; ++m) {
                                dL_dt += syn_weights[m] * io_m[m].V[vid].dot(dP_mix[m][vid] * io_m[m].D[vid]);
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
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        int max_iters = 10
    ) const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Scalar> bl = init_bl;
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id)) {
                bl[i] = std::clamp(bl[i], 1e-4, 10.0);
            }
        }

        // L-BFGS state
        const size_t m_history = 5;
        std::vector<std::vector<Scalar>> s_hist;
        std::vector<std::vector<Scalar>> y_hist;
        std::vector<Scalar> rho_hist;

        auto [curr_ll, curr_grad] = compute_busted_branch_gradients(bl, omegas, weights, syn_rates, syn_weights);

        Scalar best_ll = curr_ll;
        std::vector<Scalar> best_bl = bl;

        for (int iter = 0; iter < max_iters; ++iter) {
            // Check convergence
            Scalar max_g = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_g = std::max(max_g, std::abs(curr_grad[i]));
                }
            }
            if (max_g < 1e-3) break;

            // Two-loop L-BFGS recursion to compute search direction d = H * grad
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

            // Initial Hessian scaling H0 = gamma * I
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

            // Direction d (for maximization, d points along r)
            std::vector<Scalar> d = r;
            d[tree.root_id] = 0.0;

            // Ensure ascent direction (d . grad > 0)
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
            Scalar step_size = 1.0;
            std::vector<Scalar> new_bl(num_nodes);
            Scalar new_ll = curr_ll;
            std::vector<Scalar> new_grad;

            bool line_search_ok = false;
            for (int ls = 0; ls < 8; ++ls) {
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) {
                        new_bl[i] = std::clamp(bl[i] + step_size * d[i], 1e-4, 10.0);
                    } else {
                        new_bl[i] = 0.0;
                    }
                }
                auto [test_ll, test_grad] = compute_busted_branch_gradients(new_bl, omegas, weights, syn_rates, syn_weights);
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
                y_k[i] = curr_grad[i] - new_grad[i]; // for maximization
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

    // Fits BUSTED model for specified number of rate classes K in {1, 2, 3} and optional SRV
    BUSTEDFit fit_model(
        bool is_constrained,
        size_t K = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3
    ) const {
        BUSTEDFit fit;
        fit.num_rate_classes = K;
        Scalar s = 1.0;

        // Initialize rates and weights
        std::vector<Scalar> omegas(K);
        std::vector<Scalar> weights(K);

        if (K == 1) {
            weights[0] = 1.0;
            omegas[0] = is_constrained ? std::min(0.5, mg94_omega) : mg94_omega;
        } else if (K == 2) {
            weights = {0.75, 0.25};
            if (is_constrained) {
                omegas = {0.10, 1.0};
            } else {
                omegas = {0.10, 3.0};
            }
        } else {
            // Default K = 3
            if (is_constrained) {
                omegas = {0.05, 0.50, 1.0};
                weights = {0.45, 0.45, 0.10};
            } else {
                omegas = {0.05, 0.50, 4.0};
                weights = {0.50, 0.40, 0.10};
            }
        }

        // Initialize SRV synonymous rate distribution
        size_t M = srv ? num_syn_rates : 1;
        std::vector<Scalar> syn_rates(M, 1.0);
        std::vector<Scalar> syn_weights(M, 1.0 / M);

        if (srv && M == 3) {
            syn_rates = {0.20, 0.60, 2.20};
            syn_weights = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        }

        // ECM Optimization Loop
        int max_ecm = (K == 1 && !srv) ? 1 : 6;
        for (int ecm_iter = 0; ecm_iter < max_ecm; ++ecm_iter) {
            // Block 1: SQUAREM acceleration for omega weights
            if (K == 2) {
                for (int sq = 0; sq < 2; ++sq) {
                    auto w0 = weights;
                    auto w1 = compute_expected_weights(s, omegas, w0, syn_rates, syn_weights);
                    auto w2 = compute_expected_weights(s, omegas, w1, syn_rates, syn_weights);
                    Scalar r = w1[0] - w0[0];
                    Scalar v = (w2[0] - w1[0]) - r;
                    if (std::abs(v) > 1e-12) {
                        Scalar step = -std::abs(r) / std::abs(v);
                        Scalar p_acc = std::clamp(w0[0] - 2.0 * step * r + step * step * v, 1e-5, 1.0 - 1e-5);
                        std::vector<Scalar> w_acc = {p_acc, 1.0 - p_acc};
                        auto w_stab = compute_expected_weights(s, omegas, w_acc, syn_rates, syn_weights);
                        if (evaluate_log_l(s, omegas, w_stab, syn_rates, syn_weights) >= evaluate_log_l(s, omegas, w2, syn_rates, syn_weights)) {
                            weights = w_stab;
                        } else {
                            weights = w2;
                        }
                    } else {
                        weights = w2;
                    }
                }
            } else if (K == 3) {
                for (int sq = 0; sq < 2; ++sq) {
                    auto w0 = weights;
                    auto w1 = compute_expected_weights(s, omegas, w0, syn_rates, syn_weights);
                    auto w2 = compute_expected_weights(s, omegas, w1, syn_rates, syn_weights);
                    Scalar r0 = w1[0] - w0[0], r1 = w1[1] - w0[1];
                    Scalar v0 = (w2[0] - w1[0]) - r0, v1 = (w2[1] - w1[1]) - r1;
                    Scalar r2 = r0 * r0 + r1 * r1;
                    Scalar v2 = v0 * v0 + v1 * v1;
                    if (v2 > 1e-12) {
                        Scalar step = -std::sqrt(r2 / v2);
                        Scalar p1_acc = std::clamp(w0[0] - 2.0 * step * r0 + step * step * v0, 1e-5, 1.0 - 1e-5);
                        Scalar p2_acc = std::clamp(w0[1] - 2.0 * step * r1 + step * step * v1, 1e-5, 1.0 - 1e-5);
                        if (p1_acc + p2_acc > 1.0 - 1e-5) {
                            Scalar sum_p = p1_acc + p2_acc;
                            p1_acc = (p1_acc / sum_p) * 0.99;
                            p2_acc = (p2_acc / sum_p) * 0.99;
                        }
                        Scalar p3_acc = 1.0 - p1_acc - p2_acc;
                        std::vector<Scalar> w_acc = {p1_acc, p2_acc, p3_acc};
                        auto w_stab = compute_expected_weights(s, omegas, w_acc, syn_rates, syn_weights);
                        if (evaluate_log_l(s, omegas, w_stab, syn_rates, syn_weights) >= evaluate_log_l(s, omegas, w2, syn_rates, syn_weights)) {
                            weights = w_stab;
                        } else {
                            weights = w2;
                        }
                    } else {
                        weights = w2;
                    }
                }
            }

            // Block 1b: SQUAREM acceleration for synonymous weights (if SRV)
            if (srv && M >= 2) {
                for (int sq = 0; sq < 2; ++sq) {
                    auto q0 = syn_weights;
                    auto q1 = compute_expected_syn_weights(s, omegas, weights, syn_rates, q0);
                    auto q2 = compute_expected_syn_weights(s, omegas, weights, syn_rates, q1);
                    if (M == 2) {
                        Scalar r = q1[0] - q0[0];
                        Scalar v = (q2[0] - q1[0]) - r;
                        if (std::abs(v) > 1e-12) {
                            Scalar step = -std::abs(r) / std::abs(v);
                            Scalar q_acc = std::clamp(q0[0] - 2.0 * step * r + step * step * v, 1e-5, 1.0 - 1e-5);
                            std::vector<Scalar> q_cand = {q_acc, 1.0 - q_acc};
                            auto q_stab = compute_expected_syn_weights(s, omegas, weights, syn_rates, q_cand);
                            if (evaluate_log_l(s, omegas, weights, syn_rates, q_stab) >= evaluate_log_l(s, omegas, weights, syn_rates, q2)) {
                                syn_weights = q_stab;
                            } else {
                                syn_weights = q2;
                            }
                        } else {
                            syn_weights = q2;
                        }
                    } else if (M == 3) {
                        Scalar r0 = q1[0] - q0[0], r1 = q1[1] - q0[1];
                        Scalar v0 = (q2[0] - q1[0]) - r0, v1 = (q2[1] - q1[1]) - r1;
                        Scalar r2 = r0 * r0 + r1 * r1;
                        Scalar v2 = v0 * v0 + v1 * v1;
                        if (v2 > 1e-12) {
                            Scalar step = -std::sqrt(r2 / v2);
                            Scalar q1_acc = std::clamp(q0[0] - 2.0 * step * r0 + step * step * v0, 1e-5, 1.0 - 1e-5);
                            Scalar q2_acc = std::clamp(q0[1] - 2.0 * step * r1 + step * step * v1, 1e-5, 1.0 - 1e-5);
                            if (q1_acc + q2_acc > 1.0 - 1e-5) {
                                Scalar sum_q = q1_acc + q2_acc;
                                q1_acc = (q1_acc / sum_q) * 0.99;
                                q2_acc = (q2_acc / sum_q) * 0.99;
                            }
                            Scalar q3_acc = 1.0 - q1_acc - q2_acc;
                            std::vector<Scalar> q_cand = {q1_acc, q2_acc, q3_acc};
                            auto q_stab = compute_expected_syn_weights(s, omegas, weights, syn_rates, q_cand);
                            if (evaluate_log_l(s, omegas, weights, syn_rates, q_stab) >= evaluate_log_l(s, omegas, weights, syn_rates, q2)) {
                                syn_weights = q_stab;
                            } else {
                                syn_weights = q2;
                            }
                        } else {
                            syn_weights = q2;
                        }
                    } else {
                        syn_weights = q2;
                    }
                }
            }

            // Block 2: Continuous parameter search on rates + scale s
            if (is_constrained) {
                if (K == 1) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_log_l(s_c, {w1}, weights, syn_rates, syn_weights);
                    }, omegas[0], s, 1e-4, 1.0, 1e-4, 40);
                    omegas[0] = opt.a;
                    s = opt.b;
                } else if (K == 2) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_log_l(s_c, {w1, 1.0}, weights, syn_rates, syn_weights);
                    }, omegas[0], s, 1e-4, 1.0, 1e-4, 40);
                    omegas[0] = opt.a;
                    omegas[1] = 1.0;
                    s = opt.b;
                } else if (K == 3) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        if (x[1] > x[2]) return 1e20;
                        return -evaluate_log_l(x[0], {x[1], x[2], 1.0}, weights, syn_rates, syn_weights);
                    }, {s, omegas[0], omegas[1]}, {0.01, 1e-4, 1e-4}, {50.0, 1.0, 1.0}, 1e-4, 40);
                    s = opt.x[0];
                    omegas[0] = opt.x[1];
                    omegas[1] = opt.x[2];
                    omegas[2] = 1.0;
                }
            } else {
                if (K == 1) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_log_l(s_c, {w1}, weights, syn_rates, syn_weights);
                    }, omegas[0], s, 1e-4, 2000.0, 1e-4, 40);
                    omegas[0] = opt.a;
                    s = opt.b;
                } else if (K == 2) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        return -evaluate_log_l(x[0], {x[1], x[2]}, weights, syn_rates, syn_weights);
                    }, {s, omegas[0], omegas[1]}, {0.01, 1e-4, 1.0}, {50.0, 1.0, 2000.0}, 1e-4, 50);
                    s = opt.x[0];
                    omegas[0] = opt.x[1];
                    omegas[1] = opt.x[2];
                } else if (K == 3) {
                    auto opt = NelderMeadND<4>::minimize([&](const std::array<Scalar, 4>& x) {
                        if (x[1] > x[2]) return 1e20;
                        return -evaluate_log_l(x[0], {x[1], x[2], x[3]}, weights, syn_rates, syn_weights);
                    }, {s, omegas[0], omegas[1], omegas[2]}, {0.01, 1e-4, 1e-4, 1.0}, {50.0, 1.0, 1.0, 2000.0}, 1e-4, 50);
                    s = opt.x[0];
                    omegas[0] = opt.x[1];
                    omegas[1] = opt.x[2];
                    omegas[2] = opt.x[3];
                }
            }

            // Block 2b: Continuous search on synonymous rates (if SRV && M == 3)
            if (srv && M == 3) {
                auto syn_opt = NelderMead2D::minimize([&](Scalar a1, Scalar a2) {
                    Scalar q1 = syn_weights[0], q2 = syn_weights[1], q3 = syn_weights[2];
                    Scalar a3 = (1.0 - q1 * a1 - q2 * a2) / q3;
                    if (a1 < 0.01 || a2 < a1 || a3 < a2 || a3 > 50.0) return 1e20;
                    return -evaluate_log_l(s, omegas, weights, {a1, a2, a3}, syn_weights);
                }, syn_rates[0], syn_rates[1], 0.01, 5.0, 1e-3, 30);

                Scalar a1 = syn_opt.a;
                Scalar a2 = syn_opt.b;
                Scalar a3 = (1.0 - syn_weights[0] * a1 - syn_weights[1] * a2) / syn_weights[2];
                if (a1 >= 0.01 && a2 >= a1 && a3 >= a2) {
                    syn_rates[0] = a1;
                    syn_rates[1] = a2;
                    syn_rates[2] = a3;
                }
            }
        }

        // Branch length initialization and refinement
        std::vector<Scalar> bl(tree.num_nodes(), 0.0);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                bl[node.id] = s * branch_tau[node.id];
            }
        }

        if (refine_branches) {
            bl = refine_branch_lengths(bl, omegas, weights, syn_rates, syn_weights, 1);
            if (K >= 2) {
                weights = compute_expected_weights(1.0, omegas, weights, syn_rates, syn_weights, &bl);
            }
            if (srv && M >= 2) {
                syn_weights = compute_expected_syn_weights(1.0, omegas, weights, syn_rates, syn_weights, &bl);
            }
        }

        Scalar best_lnL = evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, &fit.site_log_likelihoods);
        fit.log_likelihood = best_lnL;
        fit.tree_scale = s;
        fit.branch_lengths = bl;
        fit.test_distribution.omegas = omegas;
        fit.test_distribution.weights = weights;
        if (srv) {
            fit.test_distribution.syn_rates = syn_rates;
            fit.test_distribution.syn_weights = syn_weights;
        }

        // Detect collapsed rate classes
        for (size_t i = 0; i < K; ++i) {
            if (weights[i] < 0.005) {
                fit.test_distribution.annotations.push_back("Rate class " + std::to_string(i) + " collapsed (weight < 0.5%)");
            }
        }
        for (size_t i = 0; i + 1 < K; ++i) {
            if (std::abs(omegas[i+1] - omegas[i]) < 0.015) {
                fit.test_distribution.annotations.push_back("Rate classes " + std::to_string(i) + " and " + std::to_string(i+1) + " collapsed into identical omega");
            }
        }

        size_t N_codons = aln.num_codons;
        size_t num_branches = tree.num_nodes() - 1;
        size_t num_mixture_params = (2 * K - 1) - (is_constrained ? 1 : 0);
        size_t num_syn_params = srv ? 2 * (M - 1) : 0;
        size_t K_params = 5 + (refine_branches ? num_branches : 1) + num_mixture_params + num_syn_params;
        fit.aicc = -2.0 * fit.log_likelihood + 2.0 * K_params * N_codons / (N_codons - K_params - 1.0);

        return fit;
    }

    // Step-up model selection across K in {1, ..., max_k} using AICc on the Null model
    std::pair<size_t, std::vector<BUSTEDFit>> select_optimal_k(
        size_t max_k = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3
    ) const {
        std::vector<BUSTEDFit> null_fits;
        null_fits.reserve(max_k);
        for (size_t k = 1; k <= max_k; ++k) {
            null_fits.push_back(fit_model(true, k, refine_branches, srv, num_syn_rates));
        }

        size_t optimal_k = 1;
        for (size_t k = 2; k <= max_k; ++k) {
            // Delta AICc threshold: improvement must be at least 2.0 to justify adding a rate class
            if (null_fits[k - 1].aicc < null_fits[optimal_k - 1].aicc - 2.0) {
                optimal_k = k;
            }
        }
        return {optimal_k, null_fits};
    }

    // Run complete BUSTED analysis with options
    BUSTEDResult run(BUSTEDSettings settings = {}) const {
        auto t0 = std::chrono::high_resolution_clock::now();
        BUSTEDResult res;

        if (settings.auto_select_k) {
            auto [opt_k, null_fits] = select_optimal_k(
                settings.max_k,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes
            );
            res.optimal_k = opt_k;
            res.k_null_fits = null_fits;
            res.constrained = null_fits[opt_k - 1];
            res.unconstrained = fit_model(
                false,
                opt_k,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes
            );
        } else {
            res.optimal_k = settings.num_rate_classes;
            res.unconstrained = fit_model(
                false,
                settings.num_rate_classes,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes
            );
            res.constrained = fit_model(
                true,
                settings.num_rate_classes,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes
            );
        }

        // Compute LRT and p-value
        if (res.unconstrained.log_likelihood < res.constrained.log_likelihood) {
            res.unconstrained.log_likelihood = res.constrained.log_likelihood;
            res.lrt = 0.0;
            res.p_value = 0.5;
        } else {
            res.lrt = 2.0 * (res.unconstrained.log_likelihood - res.constrained.log_likelihood);
            // Asymptotic null: 0.5 * chi^2_0 + 0.5 * chi^2_2
            res.p_value = 0.5 * std::exp(-res.lrt / 2.0);
        }

        // Compute per-site Evidence Ratios (Bayes Factors for positive selection)
        size_t num_sites = aln.num_codons;
        res.evidence_ratios.assign(num_sites, 1.0);
        Scalar p_pos = res.unconstrained.test_distribution.weights.back();
        Scalar prior_odds = (p_pos > 1e-6 && p_pos < 1.0 - 1e-6) ? (p_pos / (1.0 - p_pos)) : 1e-6;

        for (size_t s_idx = 0; s_idx < num_sites; ++s_idx) {
            size_t p_idx = aln.site_to_pattern[s_idx];
            Scalar null_ll = res.constrained.site_log_likelihoods[p_idx];
            Scalar alt_ll = res.unconstrained.site_log_likelihoods[p_idx];

            if (alt_ll > null_ll) {
                Scalar delta = alt_ll - null_ll;
                Scalar post_odds = std::exp(delta) * prior_odds;
                Scalar bf = post_odds / prior_odds;
                res.evidence_ratios[s_idx] = bf;
            } else {
                res.evidence_ratios[s_idx] = 1.0;
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    nlohmann::json to_json(const BUSTEDResult& res) const {
        nlohmann::json j;

        j["analysis"] = {
            {"authors", "Sergei L Kosakovsky Pond"},
            {"citation", "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71"},
            {"contact", "spond@temple.edu"},
            {"version", "4.7-next"},
            {"settings", {
                {"error-sink", 0},
                {"mss", "No"},
                {"multiple-hit", "None"},
                {"srv", (!res.unconstrained.test_distribution.syn_rates.empty() && res.unconstrained.test_distribution.syn_rates.size() > 1) ? "Yes" : "No"},
                {"optimal_k", res.optimal_k}
            }}
        };

        j["test results"] = {
            {"LRT", res.lrt},
            {"p-value", res.p_value},
            {"optimal_k", res.optimal_k}
        };

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

        fits["MG94xREV with separate rates for branch sets"] = {
            {"Log Likelihood", mg94_log_l},
            {"Rate Distributions", {
                {"non-synonymous/synonymous rate ratio for *test*", {{mg94_omega, 1.0}}}
            }}
        };

        auto format_distro = [](const BUSTEDRateDistribution& d) {
            nlohmann::json distro;
            for (size_t k = 0; k < d.omegas.size(); ++k) {
                nlohmann::json entry = {
                    {"omega", d.omegas[k]},
                    {"proportion", d.weights[k]}
                };
                if (!d.annotations.empty() && k < d.annotations.size()) {
                    entry["annotation"] = d.annotations[k];
                }
                distro[std::to_string(k)] = entry;
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

        nlohmann::json unc_rates = {{"Test", format_distro(res.unconstrained.test_distribution)}};
        if (!res.unconstrained.test_distribution.syn_rates.empty() && res.unconstrained.test_distribution.syn_rates.size() > 1) {
            unc_rates["Synonymous site-to-site rates"] = format_srv_distro(res.unconstrained.test_distribution);
        }
        fits["Unconstrained model"] = {
            {"AIC-c", res.unconstrained.aicc},
            {"Log Likelihood", res.unconstrained.log_likelihood},
            {"Rate Distributions", unc_rates}
        };

        nlohmann::json con_rates = {{"Test", format_distro(res.constrained.test_distribution)}};
        if (!res.constrained.test_distribution.syn_rates.empty() && res.constrained.test_distribution.syn_rates.size() > 1) {
            con_rates["Synonymous site-to-site rates"] = format_srv_distro(res.constrained.test_distribution);
        }
        fits["Constrained model"] = {
            {"AIC-c", res.constrained.aicc},
            {"Log Likelihood", res.constrained.log_likelihood},
            {"Rate Distributions", con_rates}
        };

        j["fits"] = fits;

        // Branch attributes
        nlohmann::json branch_attr;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                nlohmann::json b_entry;
                if (!res.unconstrained.branch_lengths.empty() && node.id >= 0 && static_cast<size_t>(node.id) < res.unconstrained.branch_lengths.size()) {
                    b_entry["unconstrained"] = res.unconstrained.branch_lengths[node.id];
                }
                if (!res.constrained.branch_lengths.empty() && node.id >= 0 && static_cast<size_t>(node.id) < res.constrained.branch_lengths.size()) {
                    b_entry["constrained"] = res.constrained.branch_lengths[node.id];
                }
                branch_attr[node.name] = b_entry;
            }
        }
        j["branch attributes"] = {{"0", branch_attr}};

        j["Evidence Ratios"] = {
            {"constrained", {res.evidence_ratios}}
        };

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

        j["runtime"] = res.runtime_seconds;
        return j;
    }
};

} // namespace hyphy::analyses
