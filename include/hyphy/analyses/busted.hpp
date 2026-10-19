#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/stats.hpp"
#include "hyphy/analyses/busted_types.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/squarem.hpp"
#include "nlohmann/json.hpp"
#include "hyphy/core/provenance.hpp"

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

class BUSTEDAnalyzer {
public:
    Tree tree;
    Alignment aln;
    std::shared_ptr<const GeneticCode> code;
    MG94Parameters base_params;

    // Relative branch lengths tau_b (normalized from MG94 baseline)
    std::vector<Scalar> branch_tau;
    mutable std::vector<bool> is_test_branch;
    // 0 = Test, 1 = Background/Comparison, 2 = Unclassified/Nuisance
    mutable std::vector<int> branch_partition;
    std::vector<size_t> leaf_to_taxon;
    mutable size_t num_test_branches = 0;
    mutable size_t num_background_branches = 0;
    mutable size_t num_unclassified_branches = 0;

    // Precomputed off-diagonal saturation replacement kernel P_sink(i, j) = pi_j / (1 - pi_i) for i != j, 0 for i == j
    Matrix P_sink;
    Vector inv_one_minus_pi;

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
        branch_partition.assign(num_nodes, 0);
        if (tree.root_id >= 0 && static_cast<size_t>(tree.root_id) < num_nodes) {
            is_test_branch[tree.root_id] = false;
        }
        num_test_branches = (num_nodes > 0) ? (num_nodes - 1) : 0; // Non-root branches
        num_background_branches = 0;
        num_unclassified_branches = 0;
        leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);
        build_sink_kernel();
        prepare_baseline_lengths();
    }

    void select_branch_partitions(
        const std::string& test_spec,
        const std::string& comp_spec = "",
        bool auto_detect_tags = false
    ) const {
        size_t num_nodes = tree.num_nodes();
        auto sel_test = tree.select_branches(test_spec, "", {}, auto_detect_tags);
        is_test_branch = sel_test.is_selected;
        num_test_branches = sel_test.num_selected;
        branch_partition.assign(num_nodes, 1);
        num_background_branches = 0;
        num_unclassified_branches = 0;

        if (!comp_spec.empty()) {
            auto sel_comp = tree.select_branches(comp_spec, "", {}, false);
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                int nid = node.id;
                if (is_test_branch[nid]) {
                    branch_partition[nid] = 0;
                } else if (sel_comp.is_selected[nid]) {
                    branch_partition[nid] = 1;
                    num_background_branches++;
                } else {
                    branch_partition[nid] = 2;
                    num_unclassified_branches++;
                }
            }
            if (num_background_branches == 0) {
                throw std::runtime_error("No non-test branches matched '--comparison " + comp_spec + "'.");
            }
        } else {
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                int nid = node.id;
                if (is_test_branch[nid]) {
                    branch_partition[nid] = 0;
                } else {
                    branch_partition[nid] = 1;
                    num_background_branches++;
                }
            }
        }
    }

    void select_test_branches(const std::string& spec) const {
        select_branch_partitions(spec, "", false);
    }

    void build_sink_kernel() {
        int S = code->num_sense_codons;
        P_sink = Matrix::Zero(S, S);
        inv_one_minus_pi = Vector::Zero(S);
        for (int i = 0; i < S; ++i) {
            Scalar denom = std::max(1e-12, 1.0 - aln.codon_frequencies_f3x4(i));
            inv_one_minus_pi(i) = 1.0 / denom;
            for (int j = 0; j < S; ++j) {
                if (i != j) {
                    P_sink(i, j) = aln.codon_frequencies_f3x4(j) * inv_one_minus_pi(i);
                }
            }
        }
    }

    // O(S) evaluation of V^T * P_sink * D using the rank-1 off-diagonal identity
    inline Scalar eval_sink_bilinear(const Vector& V, const Vector& D) const {
        Scalar pi_dot_D = aln.codon_frequencies_f3x4.dot(D);
        int S = static_cast<int>(V.size());
        Scalar acc = 0.0;
        for (int i = 0; i < S; ++i) {
            Scalar p_sink_D_i = (pi_dot_D - aln.codon_frequencies_f3x4(i) * D(i)) * inv_one_minus_pi(i);
            acc += V(i) * std::max(0.0, p_sink_D_i);
        }
        return acc;
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
        auto mg_res = mg_fitter.fit_full_model(1.0, base_p, nullptr, 35, 1e-4);
        base_p = mg_res.params;

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
        analyzer.mg94_omega = mg_res.params.beta;
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

    // Evaluates tree log-likelihood for arbitrary K rate categories, optional SRV (M classes), optional saturation error sink (eps), and custom branch lengths
    Scalar evaluate_log_l(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        std::vector<Scalar>* out_pattern_ll = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink = 0.0
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();
        size_t M = syn_rates.size();
        Scalar bio_scale = 1.0 - eps_sink;

        // Build K rate matrices
        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        bool has_bg = (num_test_branches + 1 < num_nodes);
        MG94Matrix mg_bg;
        if (has_bg) {
            MG94Parameters p_bg = base_params;
            p_bg.alpha = 1.0;
            p_bg.beta = mg94_omega;
            p_bg.delta = delta;
            p_bg.psi = psi;
            mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        // Precompute branch transition matrices P_mix[m][node_id]
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[node.id] : (s * branch_tau[node.id]);
                    if (is_test_branch[node.id]) {
                        P_mix[m][node.id] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink > 0.0) ? bio_scale : 1.0;
                        for (size_t k = 0; k < K; ++k) {
                            P_mix[m][node.id] += (use_bio * weights[k]) * mg[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink > 0.0) {
                            P_mix[m][node.id].noalias() += eps_sink * P_sink;
                        }
                    } else {
                        P_mix[m][node.id] = mg_bg.transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        if (out_pattern_ll) {
            out_pattern_ll->assign(num_patterns, 0.0);
        }

        Scalar total_log_l = 0.0;
        #pragma omp parallel
        {
            auto [node_L, node_log_scale] = LikelihoodEngine::make_pattern_scratch<Vector>(num_nodes, S);
            #pragma omp for reduction(+:total_log_l) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pll = -1e20;
                if (M == 1) {
                    LikelihoodEngine::eval_pattern_postorder(tree, pattern, leaf_to_taxon, P_mix[0], S, node_L, node_log_scale);
                    Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                    if (cond_L > 0.0) {
                        pll = std::log(cond_L) + node_log_scale[tree.root_id];
                    }
                } else {
                    std::vector<Scalar> terms(M, -std::numeric_limits<Scalar>::infinity());
                    for (size_t m = 0; m < M; ++m) {
                        LikelihoodEngine::eval_pattern_postorder(tree, pattern, leaf_to_taxon, P_mix[m], S, node_L, node_log_scale);
                        Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                        if (cond_L > 0.0) {
                            terms[m] = std::log(syn_weights[m]) + std::log(cond_L) + node_log_scale[tree.root_id];
                        }
                    }

                    Scalar lse = LikelihoodEngine::logsumexp(terms);
                    if (std::isfinite(lse)) pll = lse;
                }
                total_log_l += pattern.weight * pll;
                if (out_pattern_ll) {
                    (*out_pattern_ll)[p] = pll;
                }
            }
        }
        return total_log_l;
    }

    // Expectation Step: Computes new mixture weights (and optional saturation sink weight eps_sink) via inside-outside traversal
    std::vector<Scalar> compute_expected_weights(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink = 0.0,
        Scalar* out_new_eps_sink = nullptr,
        Scalar max_eps_sink = 0.01
    ) const {
        size_t K = omegas.size();
        if (K <= 1 && !out_new_eps_sink) return {1.0};

        size_t M = syn_rates.size();
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        Scalar bio_scale = (out_new_eps_sink != nullptr) ? (1.0 - eps_sink) : 1.0;

        std::vector<std::vector<std::vector<Matrix>>> P_comp(M, std::vector<std::vector<Matrix>>(K, std::vector<Matrix>(num_nodes)));
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
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
                    Scalar use_bio = (out_new_eps_sink != nullptr && is_test_branch[node.id]) ? bio_scale : 1.0;
                    for (size_t k = 0; k < K; ++k) {
                        P_mix[m][node.id] += (use_bio * weights[k]) * P_comp[m][k][node.id];
                    }
                    if (out_new_eps_sink != nullptr && is_test_branch[node.id]) {
                        P_mix[m][node.id].noalias() += eps_sink * P_sink;
                    }
                }
            }
        }

        Vector weight_accum = Vector::Zero(K);
        Scalar sink_accum = 0.0;
        Scalar total_pattern_weight = 0.0;

        #pragma omp parallel
        {
            Vector local_accum = Vector::Zero(K);
            Scalar local_sink_accum = 0.0;
            Scalar local_pat_wt = 0.0;
            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                std::vector<LikelihoodEngine::InsideOutsideResult> io_m(M);
                Scalar max_ll = -std::numeric_limits<Scalar>::infinity();
                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
                    io_m[m] = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S
                    );
                    if (io_m[m].likelihood > 0.0) {
                        Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                        if (term > max_ll) max_ll = term;
                    }
                }

                if (max_ll > -1e19) {
                    Scalar sum_exp = 0.0;
                    std::vector<Scalar> post_m(M, 0.0);
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                            Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood - max_ll);
                            post_m[m] = e;
                            sum_exp += e;
                        }
                    }
                    if (sum_exp > 0.0) {
                        Scalar inv_sum = 1.0 / sum_exp;
                        for (size_t m = 0; m < M; ++m) post_m[m] *= inv_sum;

                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id && is_test_branch[node.id]) {
                                int nid = node.id;
                                for (size_t k = 0; k < K; ++k) {
                                    Scalar ratio_k = 0.0;
                                    for (size_t m = 0; m < M; ++m) {
                                        if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                            ratio_k += post_m[m] * (io_m[m].V[nid].dot(P_comp[m][k][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                                        }
                                    }
                                    local_accum(k) += pattern.weight * (bio_scale * weights[k]) * ratio_k;
                                }
                                if (out_new_eps_sink != nullptr) {
                                    Scalar sink_ratio = 0.0;
                                    for (size_t m = 0; m < M; ++m) {
                                        if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                            sink_ratio += post_m[m] * (eval_sink_bilinear(io_m[m].V[nid], io_m[m].D[nid]) / io_m[m].likelihood);
                                        }
                                    }
                                    local_sink_accum += pattern.weight * eps_sink * sink_ratio;
                                }
                                local_pat_wt += pattern.weight;
                            }
                        }
                    }
                }
            }
            #pragma omp critical
            {
                weight_accum += local_accum;
                sink_accum += local_sink_accum;
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
            if (out_new_eps_sink != nullptr) {
                Scalar raw_eps = sink_accum / total_pattern_weight;
                *out_new_eps_sink = std::clamp(raw_eps, 1e-7, max_eps_sink);
            }
        } else {
            new_weights = weights;
            if (out_new_eps_sink != nullptr) {
                *out_new_eps_sink = eps_sink;
            }
        }
        return new_weights;
    }

    // Expectation Step: Computes new synonymous rate weights via site posteriors (with log-underflow scaling)
    std::vector<Scalar> compute_expected_syn_weights(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates,
        const std::vector<Scalar>& syn_weights,
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink = 0.0
    ) const {
        size_t M = syn_rates.size();
        if (M <= 1) return {1.0};

        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();
        Scalar bio_scale = 1.0 - eps_sink;

        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        bool has_bg = (num_test_branches + 1 < num_nodes);
        MG94Matrix mg_bg;
        if (has_bg) {
            MG94Parameters p_bg = base_params;
            p_bg.alpha = 1.0;
            p_bg.beta = mg94_omega;
            p_bg.delta = delta;
            p_bg.psi = psi;
            mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[node.id] : (s * branch_tau[node.id]);
                    if (is_test_branch[node.id]) {
                        P_mix[m][node.id] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink > 0.0) ? bio_scale : 1.0;
                        for (size_t k = 0; k < K; ++k) {
                            P_mix[m][node.id] += (use_bio * weights[k]) * mg[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink > 0.0) {
                            P_mix[m][node.id].noalias() += eps_sink * P_sink;
                        }
                    } else {
                        P_mix[m][node.id] = mg_bg.transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        Vector syn_accum = Vector::Zero(M);
        Scalar total_pattern_weight = 0.0;

        #pragma omp parallel
        {
            Vector local_accum = Vector::Zero(M);
            Scalar local_pat_wt = 0.0;
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));
            std::vector<double> node_log_scale(num_nodes, 0.0);

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar max_term = -std::numeric_limits<Scalar>::infinity();
                std::vector<Scalar> terms(M, -std::numeric_limits<Scalar>::infinity());
                std::vector<bool> valid(M, false);

                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
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
                            node_log_scale[node_id] = 0.0;
                        } else {
                            double child_scale_sum = 0.0;
                            tree.eval_internal_node(node_id, P_mix[m], node_L, node_log_scale, node_L[node_id], child_scale_sum);
                            node_log_scale[node_id] = child_scale_sum;
                        }
                    }
                    Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                    if (cond_L > 0.0) {
                        Scalar term = std::log(syn_weights[m]) + std::log(cond_L) + node_log_scale[tree.root_id];
                        terms[m] = term;
                        valid[m] = true;
                        if (term > max_term) max_term = term;
                    }
                }

                if (max_term > -1e19) {
                    Scalar sum_exp = 0.0;
                    std::vector<Scalar> rho(M, 0.0);
                    for (size_t m = 0; m < M; ++m) {
                        if (valid[m]) {
                            rho[m] = std::exp(terms[m] - max_term);
                            sum_exp += rho[m];
                        }
                    }
                    if (sum_exp > 0.0) {
                        for (size_t m = 0; m < M; ++m) {
                            local_accum(m) += pattern.weight * (rho[m] / sum_exp);
                        }
                        local_pat_wt += pattern.weight;
                    }
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
    // under the BUSTED mixture model using Inside-Outside adjoints (with optional SRV and error sink)
    std::pair<Scalar, std::vector<Scalar>> compute_busted_branch_gradients(
        const std::vector<Scalar>& branch_lengths,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink = 0.0
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();
        size_t M = syn_rates.size();
        Scalar bio_scale = 1.0 - eps_sink;

        // 1. Build K rate matrices
        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        bool has_bg = (num_test_branches + 1 < num_nodes);
        MG94Matrix mg_bg;
        if (has_bg) {
            MG94Parameters p_bg = base_params;
            p_bg.alpha = 1.0;
            p_bg.beta = mg94_omega;
            p_bg.delta = delta;
            p_bg.psi = psi;
            mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        // 2. Precompute P_mix[m] and dP_mix[m] for each synonymous rate category and branch
        // Note: d/dt_b (eps_sink * P_sink) == 0 because the saturation replacement kernel is branch-length invariant!
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<Matrix>> dP_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = branch_lengths[nid];
                    if (is_test_branch[nid]) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        dP_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink > 0.0) ? bio_scale : 1.0;
                        for (size_t k = 0; k < K; ++k) {
                            Matrix P_k = mg[k].transition_matrix(alpha_m * bl);
                            P_mix[m][nid] += (use_bio * weights[k]) * P_k;
                            dP_mix[m][nid] += (use_bio * weights[k]) * alpha_m * (mg[k].Q * P_k);
                        }
                        if (eps_sink > 0.0) {
                            P_mix[m][nid].noalias() += eps_sink * P_sink;
                        }
                    } else {
                        Matrix P_b = mg_bg.transition_matrix(alpha_m * bl);
                        P_mix[m][nid] = P_b;
                        dP_mix[m][nid] = alpha_m * (mg_bg.Q * P_b);
                    }
                }
            }
        }

        // 3. Parallel inside-outside traversal
        Scalar total_log_l = 0.0;
        std::vector<Scalar> grad_b(num_nodes, 0.0);

        #pragma omp parallel
        {
            std::vector<Scalar> local_grad(num_nodes, 0.0);
            Scalar local_ll = 0.0;

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                if (M == 1) {
                    auto io = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_mix[0], aln.codon_frequencies_f3x4, S
                    );
                    if (io.likelihood > 0.0) {
                        local_ll += pattern.weight * io.log_likelihood;
                        Scalar inv_L = pattern.weight / io.likelihood;
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int32_t vid = node.id;
                                Scalar dL_dt = io.V[vid].dot(dP_mix[0][vid] * io.D[vid]);
                                local_grad[vid] += inv_L * dL_dt;
                            }
                        }
                    }
                } else {
                    std::vector<LikelihoodEngine::InsideOutsideResult> io_m(M);
                    std::vector<Scalar> terms(M, -std::numeric_limits<Scalar>::infinity());
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] <= 1e-12) continue;
                        io_m[m] = LikelihoodEngine::compute_inside_outside(
                            tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S
                        );
                        if (io_m[m].likelihood > 0.0) {
                            terms[m] = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                        }
                    }

                    Scalar pat_ll = LikelihoodEngine::logsumexp(terms);
                    if (std::isfinite(pat_ll)) {
                        local_ll += pattern.weight * pat_ll;
                        std::vector<Scalar> post_w(M, 0.0);
                        for (size_t m = 0; m < M; ++m) {
                            if (std::isfinite(terms[m])) post_w[m] = std::exp(terms[m] - pat_ll);
                        }

                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int32_t vid = node.id;
                                Scalar d_lnL_dt = 0.0;
                                for (size_t m = 0; m < M; ++m) {
                                    if (post_w[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                        d_lnL_dt += post_w[m] * (io_m[m].V[vid].dot(dP_mix[m][vid] * io_m[m].D[vid]) / io_m[m].likelihood);
                                    }
                                }
                                local_grad[vid] += pattern.weight * d_lnL_dt;
                            }
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

    // High-performance joint branch length refinement via exact Inside-Outside adjoint gradients and L-BFGS in log-branch-length space
    std::vector<Scalar> refine_branch_lengths(
        const std::vector<Scalar>& init_bl,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        int max_iters = 25,
        Scalar eps_sink = 0.0
    ) const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Scalar> bl = init_bl;
        std::vector<Scalar> log_bl(num_nodes, 0.0);
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id)) {
                bl[i] = std::clamp(bl[i], 1e-8, 10.0);
                log_bl[i] = std::log(bl[i]);
            }
        }

        // L-BFGS state in log-branch-length space y_i = log(t_i), where d ln L / d y_i = t_i * (d ln L / d t_i)
        const size_t m_history = 8;
        std::vector<std::vector<Scalar>> s_hist;
        std::vector<std::vector<Scalar>> y_hist;
        std::vector<Scalar> rho_hist;

        auto eval_log_space = [&](const std::vector<Scalar>& y_vec) -> std::tuple<Scalar, std::vector<Scalar>, std::vector<Scalar>> {
            std::vector<Scalar> cur_bl(num_nodes, 0.0);
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    cur_bl[i] = std::exp(std::clamp(y_vec[i], -18.420680744, 2.302585093));
                }
            }
            auto [ll, grad_t] = compute_busted_branch_gradients(cur_bl, omegas, weights, syn_rates, syn_weights, delta, psi, eps_sink);
            std::vector<Scalar> grad_y(num_nodes, 0.0);
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    grad_y[i] = cur_bl[i] * grad_t[i];
                }
            }
            return {ll, grad_y, cur_bl};
        };

        auto [curr_ll, curr_grad, curr_bl] = eval_log_space(log_bl);
        bl = curr_bl;

        Scalar best_ll = curr_ll;
        std::vector<Scalar> best_bl = bl;

        for (int iter = 0; iter < max_iters; ++iter) {
            // Check convergence in log-space gradient
            Scalar max_g = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_g = std::max(max_g, std::abs(curr_grad[i]));
                }
            }
            if (max_g < 1e-4) break;

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
            Scalar gamma = 0.25;
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
                    gamma = std::clamp(s_dot_y / y_dot_y, 1e-3, 5.0);
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

            // Clamp per-branch log-step to [-1.2, +1.2] (at most ~3.3x multiplier per step) without throttling other branches!
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    d[i] = std::clamp(d[i], -1.2, 1.2);
                }
            }
            dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    dir_dot_grad += d[i] * curr_grad[i];
                }
            }

            Scalar step_size = 1.0;
            std::vector<Scalar> new_log_bl(num_nodes, 0.0);
            std::vector<Scalar> new_bl(num_nodes, 0.0);
            Scalar new_ll = curr_ll;
            std::vector<Scalar> new_grad;

            bool line_search_ok = false;
            for (int ls = 0; ls < 12; ++ls) {
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) {
                        new_log_bl[i] = std::clamp(log_bl[i] + step_size * d[i], -18.420680744, 2.302585093);
                    } else {
                        new_log_bl[i] = 0.0;
                    }
                }
                auto [test_ll, test_grad, test_bl] = eval_log_space(new_log_bl);
                if (test_ll > curr_ll + 1e-4 * step_size * dir_dot_grad || (test_ll > curr_ll + 1e-5 && ls >= 3)) {
                    new_ll = test_ll;
                    new_grad = test_grad;
                    new_bl = test_bl;
                    line_search_ok = true;
                    break;
                }
                step_size *= 0.5;
            }

            if (!line_search_ok) {
                break;
            }

            // Update L-BFGS history in log-branch-length space
            std::vector<Scalar> s_k(num_nodes);
            std::vector<Scalar> y_k(num_nodes);
            Scalar s_dot_y = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                s_k[i] = new_log_bl[i] - log_bl[i];
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

            log_bl = new_log_bl;
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

    // Fits BUSTED model for specified number of rate classes K in {1, 2, 3}, optional SRV, and optional saturation error sink
    BUSTEDFit fit_model(
        bool is_constrained,
        size_t K = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3,
        const std::string& multiple_hits = "None",
        ProgressBar* pbar = nullptr,
        const std::string& model_desc = "",
        bool use_error_sink = false,
        Scalar max_eps_sink = 0.01,
        const BUSTEDFit* warm_start_fit = nullptr
    ) const {
        BUSTEDFit fit;
        fit.num_rate_classes = K;
        Scalar s = 1.0;

        bool allow_double = (multiple_hits == "Double" || multiple_hits == "Double+Triple");
        bool allow_triple = (multiple_hits == "Double+Triple");
        Scalar delta = 0.0;
        Scalar psi = 0.0;
        Scalar eps_sink = use_error_sink ? std::min(0.002, max_eps_sink * 0.5) : 0.0;

        // Initialize SRV synonymous rate distribution
        size_t M = srv ? num_syn_rates : 1;
        std::vector<Scalar> syn_rates(M, 1.0);
        std::vector<Scalar> syn_weights(M, 1.0 / M);

        if (srv && M == 2) {
            syn_rates = {0.35, 1.65};
            syn_weights = {0.50, 0.50};
        } else if (srv && M == 3) {
            syn_rates = {0.30, 1.00, 3.10};
            syn_weights = {0.30, 0.60, 0.10};
        }
        if (warm_start_fit && srv && warm_start_fit->test_distribution.syn_rates.size() == M &&
            warm_start_fit->test_distribution.syn_weights.size() == M) {
            syn_rates = warm_start_fit->test_distribution.syn_rates;
            syn_weights = warm_start_fit->test_distribution.syn_weights;
        }

        const std::vector<Scalar>* ecm_bl_ptr = nullptr;
        if (warm_start_fit && !warm_start_fit->branch_lengths.empty()) {
            ecm_bl_ptr = &warm_start_fit->branch_lengths;
            s = 1.0;
        }

        // Initialize rates and weights
        std::vector<Scalar> omegas(K);
        std::vector<Scalar> weights(K);

        if (K == 1) {
            weights[0] = 1.0;
            omegas[0] = is_constrained ? std::min(0.5, mg94_omega) : mg94_omega;
        } else if (K == 2) {
            std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> candidates;
            if (is_constrained) {
                candidates = {
                    {{0.10, 1.0}, {0.75, 0.25}},
                    {{0.02, 1.0}, {0.85, 0.15}},
                    {{0.25, 1.0}, {0.60, 0.40}}
                };
            } else {
                candidates = {
                    {{0.10, 3.0}, {0.75, 0.25}},
                    {{0.02, 4.5}, {0.90, 0.10}},
                    {{0.05, 2.0}, {0.80, 0.20}}
                };
                if (warm_start_fit && warm_start_fit->test_distribution.omegas.size() == 2) {
                    Scalar w0 = warm_start_fit->test_distribution.omegas[0];
                    Scalar p0 = std::clamp(warm_start_fit->test_distribution.weights[0], 0.05, 0.95);
                    candidates.push_back({{w0, 2.5}, {p0, 1.0 - p0}});
                    candidates.push_back({{w0, 5.0}, {std::max(0.5, p0), 1.0 - std::max(0.5, p0)}});
                }
            }
            Scalar best_start_ll = -1e20;
            for (const auto& [cand_o, cand_w] : candidates) {
                Scalar ll = evaluate_log_l(s, cand_o, cand_w, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                if (ll > best_start_ll) {
                    best_start_ll = ll;
                    omegas = cand_o;
                    weights = cand_w;
                }
            }
        } else {
            // Default K = 3: evaluate grid of candidate distributions to seed ECM in best global basin
            if (is_constrained) {
                std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> candidates = {
                    {{0.05, 0.50, 1.0}, {0.45, 0.45, 0.10}},
                    {{0.01, 0.20, 1.0}, {0.30, 0.30, 0.40}},
                    {{0.10, 0.80, 1.0}, {0.20, 0.20, 0.60}},
                    {{0.001, 0.05, 1.0}, {0.15, 0.15, 0.70}}
                };
                Scalar best_start_ll = -1e20;
                for (const auto& [cand_o, cand_w] : candidates) {
                    Scalar ll = evaluate_log_l(s, cand_o, cand_w, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                    if (ll > best_start_ll) {
                        best_start_ll = ll;
                        omegas = cand_o;
                        weights = cand_w;
                    }
                }
            } else {
                std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> candidates = {
                    {{0.05, 0.50, 4.0}, {0.50, 0.40, 0.10}},
                    {{0.01, 0.50, 2.0}, {0.40, 0.10, 0.50}},
                    {{0.01, 1.00, 2.0}, {0.40, 0.05, 0.55}},
                    {{0.02, 0.20, 3.0}, {0.35, 0.35, 0.30}},
                    {{0.10, 0.90, 5.0}, {0.30, 0.40, 0.30}},
                    {{0.005, 0.10, 1.8}, {0.45, 0.15, 0.40}}
                };
                if (warm_start_fit && warm_start_fit->test_distribution.omegas.size() == 3) {
                    const auto& co = warm_start_fit->test_distribution.omegas;
                    const auto& cw = warm_start_fit->test_distribution.weights;
                    Scalar p2 = std::max(0.05, cw[2] * 0.5);
                    Scalar p0 = std::clamp(cw[0], 0.02, 0.90);
                    Scalar p1 = std::max(0.02, 1.0 - p0 - p2);
                    Scalar sum_p = p0 + p1 + p2;
                    candidates.push_back({{co[0], std::min(0.95, std::max(co[0] + 0.01, co[1])), 2.0}, {p0 / sum_p, p1 / sum_p, p2 / sum_p}});
                    candidates.push_back({{co[0], 1.0, 2.5}, {p0 / sum_p, p1 / sum_p, p2 / sum_p}});
                }
                Scalar best_start_ll = -1e20;
                for (const auto& [cand_o, cand_w] : candidates) {
                    Scalar ll = evaluate_log_l(s, cand_o, cand_w, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                    if (ll > best_start_ll) {
                        best_start_ll = ll;
                        omegas = cand_o;
                        weights = cand_w;
                    }
                }
            }
        }

        if (allow_double) {
            Scalar best_grid_ll = -1e20;
            std::vector<Scalar> d_cands = {0.0, 0.02, 0.08, 0.20};
            std::vector<Scalar> p_cands = allow_triple ? std::vector<Scalar>{0.0, 0.02, 0.05} : std::vector<Scalar>{0.0};
            for (Scalar d_c : d_cands) {
                for (Scalar p_c : p_cands) {
                    Scalar cand_ll = evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, d_c, p_c, eps_sink);
                    if (cand_ll > best_grid_ll) {
                        best_grid_ll = cand_ll;
                        delta = d_c;
                        psi = p_c;
                    }
                }
            }
        }

        // ECM Optimization Loop
        int max_ecm = (K == 1 && !srv && !allow_double && !use_error_sink) ? 1 : 6;
        for (int ecm_iter = 0; ecm_iter < max_ecm; ++ecm_iter) {
            Scalar* eps_ptr = use_error_sink ? &eps_sink : nullptr;

            // Block 1: SQUAREM acceleration for omega weights (and penalized EM update for eps_sink)
            if (K == 1 && use_error_sink) {
                weights = compute_expected_weights(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps_sink, eps_ptr, max_eps_sink);
            } else if (K == 2) {
                for (int sq = 0; sq < 2; ++sq) {
                    Scalar eps0 = eps_sink;
                    auto w0 = weights;
                    Scalar eps1 = eps0;
                    auto w1 = compute_expected_weights(s, omegas, w0, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps0, use_error_sink ? &eps1 : nullptr, max_eps_sink);
                    Scalar eps2 = eps1;
                    auto w2 = compute_expected_weights(s, omegas, w1, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps1, use_error_sink ? &eps2 : nullptr, max_eps_sink);
                    eps_sink = eps2;
                    Scalar r = w1[0] - w0[0];
                    Scalar v = (w2[0] - w1[0]) - r;
                    if (std::abs(v) > 1e-12) {
                        Scalar step = -std::abs(r) / std::abs(v);
                        Scalar p_acc = std::clamp(w0[0] - 2.0 * step * r + step * step * v, 1e-5, 1.0 - 1e-5);
                        std::vector<Scalar> w_acc = {p_acc, 1.0 - p_acc};
                        Scalar eps_stab = eps_sink;
                        auto w_stab = compute_expected_weights(s, omegas, w_acc, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps_sink, use_error_sink ? &eps_stab : nullptr, max_eps_sink);
                        if (evaluate_log_l(s, omegas, w_stab, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_stab) >=
                            evaluate_log_l(s, omegas, w2, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps2)) {
                            weights = w_stab;
                            eps_sink = eps_stab;
                        } else {
                            weights = w2;
                        }
                    } else {
                        weights = w2;
                    }
                }
            } else if (K == 3) {
                for (int sq = 0; sq < 2; ++sq) {
                    Scalar eps0 = eps_sink;
                    auto w0 = weights;
                    Scalar eps1 = eps0;
                    auto w1 = compute_expected_weights(s, omegas, w0, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps0, use_error_sink ? &eps1 : nullptr, max_eps_sink);
                    Scalar eps2 = eps1;
                    auto w2 = compute_expected_weights(s, omegas, w1, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps1, use_error_sink ? &eps2 : nullptr, max_eps_sink);
                    eps_sink = eps2;
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
                        Scalar eps_stab = eps_sink;
                        auto w_stab = compute_expected_weights(s, omegas, w_acc, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps_sink, use_error_sink ? &eps_stab : nullptr, max_eps_sink);
                        if (evaluate_log_l(s, omegas, w_stab, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_stab) >=
                            evaluate_log_l(s, omegas, w2, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps2)) {
                            weights = w_stab;
                            eps_sink = eps_stab;
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
                    auto q1 = compute_expected_syn_weights(s, omegas, weights, syn_rates, q0, ecm_bl_ptr, delta, psi, eps_sink);
                    auto q2 = compute_expected_syn_weights(s, omegas, weights, syn_rates, q1, ecm_bl_ptr, delta, psi, eps_sink);
                    if (M == 2) {
                        Scalar r = q1[0] - q0[0];
                        Scalar v = (q2[0] - q1[0]) - r;
                        if (std::abs(v) > 1e-12) {
                            Scalar step = -std::abs(r) / std::abs(v);
                            Scalar q_acc = std::clamp(q0[0] - 2.0 * step * r + step * step * v, 1e-5, 1.0 - 1e-5);
                            std::vector<Scalar> q_cand = {q_acc, 1.0 - q_acc};
                            auto q_stab = compute_expected_syn_weights(s, omegas, weights, syn_rates, q_cand, ecm_bl_ptr, delta, psi, eps_sink);
                            if (evaluate_log_l(s, omegas, weights, syn_rates, q_stab, ecm_bl_ptr, nullptr, delta, psi, eps_sink) >=
                                evaluate_log_l(s, omegas, weights, syn_rates, q2, ecm_bl_ptr, nullptr, delta, psi, eps_sink)) {
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
                            if (q1_acc + q2_acc > 1.0 - 1.e-5) {
                                Scalar sum_q = q1_acc + q2_acc;
                                q1_acc = (q1_acc / sum_q) * 0.99;
                                q2_acc = (q2_acc / sum_q) * 0.99;
                            }
                            Scalar q3_acc = 1.0 - q1_acc - q2_acc;
                            std::vector<Scalar> q_cand = {q1_acc, q2_acc, q3_acc};
                            auto q_stab = compute_expected_syn_weights(s, omegas, weights, syn_rates, q_cand, ecm_bl_ptr, delta, psi, eps_sink);
                            if (evaluate_log_l(s, omegas, weights, syn_rates, q_stab, ecm_bl_ptr, nullptr, delta, psi, eps_sink) >=
                                evaluate_log_l(s, omegas, weights, syn_rates, q2, ecm_bl_ptr, nullptr, delta, psi, eps_sink)) {
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
                        return -evaluate_log_l(s_c, {w1}, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                    }, omegas[0], s, 0.0, 1.0, 1e-4, 40);
                    omegas[0] = opt.a;
                    s = opt.b;
                } else if (K == 2) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_log_l(s_c, {w1, 1.0}, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                    }, omegas[0], s, 0.0, 1.0, 1e-4, 40);
                    omegas[0] = opt.a;
                    omegas[1] = 1.0;
                    s = opt.b;
                } else if (K == 3) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        if (x[1] > x[2]) return 1e20;
                        return -evaluate_log_l(x[0], {x[1], x[2], 1.0}, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                    }, {s, omegas[0], omegas[1]}, {0.01, 0.0, 0.0}, {50.0, 1.0, 1.0}, 1e-4, 50);
                    s = opt.x[0];
                    omegas[0] = opt.x[1];
                    omegas[1] = opt.x[2];
                    omegas[2] = 1.0;
                }
            } else {
                Scalar max_omega = use_error_sink ? 50.0 : 2000.0;
                if (K == 1) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_log_l(s_c, {w1}, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                    }, omegas[0], s, 0.0, max_omega, 1e-4, 40);
                    omegas[0] = opt.a;
                    s = opt.b;
                } else if (K == 2) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        return -evaluate_log_l(x[0], {x[1], x[2]}, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                    }, {s, omegas[0], omegas[1]}, {0.01, 0.0, 1.0}, {50.0, 1.0, max_omega}, 1e-4, 50);
                    s = opt.x[0];
                    omegas[0] = opt.x[1];
                    omegas[1] = opt.x[2];
                } else if (K == 3) {
                    auto opt = NelderMeadND<4>::minimize([&](const std::array<Scalar, 4>& x) {
                        if (x[1] > x[2] || x[2] > x[3]) return 1e20;
                        return -evaluate_log_l(x[0], {x[1], x[2], x[3]}, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                    }, {s, omegas[0], omegas[1], omegas[2]}, {0.01, 0.0, 0.0, 1.0}, {50.0, 1.0, max_omega, max_omega}, 1e-4, 80);
                    s = opt.x[0];
                    omegas[0] = opt.x[1];
                    omegas[1] = opt.x[2];
                    omegas[2] = opt.x[3];
                }
            }

            // Block 2b: Continuous search on synonymous rates (if SRV && M in {2, 3})
            if (srv && M == 2) {
                auto syn_opt = NelderMead2D::minimize([&](Scalar a1, Scalar /*unused*/) {
                    Scalar q1 = syn_weights[0], q2 = syn_weights[1];
                    Scalar a2 = (1.0 - q1 * a1) / q2;
                    if (a1 < 0.01 || a2 < a1 || a2 > 50.0) return 1e20;
                    return -evaluate_log_l(s, omegas, weights, {a1, a2}, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                }, syn_rates[0], 0.0, 0.01, 2.0, 1e-3, 25);
                Scalar a1 = syn_opt.a;
                Scalar a2 = (1.0 - syn_weights[0] * a1) / syn_weights[1];
                if (a1 >= 0.01 && a2 >= a1) {
                    syn_rates[0] = a1;
                    syn_rates[1] = a2;
                }
            } else if (srv && M == 3) {
                auto syn_opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                    Scalar a1 = x[0], a2 = x[1];
                    Scalar q1 = syn_weights[0], q2 = syn_weights[1], q3 = syn_weights[2];
                    Scalar a3 = (1.0 - q1 * a1 - q2 * a2) / q3;
                    if (a1 < 0.01 || a2 < a1 || a3 < a2 || a3 > 50.0) return 1e20;
                    return -evaluate_log_l(s, omegas, weights, {a1, a2, a3}, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                }, {syn_rates[0], syn_rates[1]}, {0.01, 0.01}, {5.0, 10.0}, 1e-4, 40);

                Scalar a1 = syn_opt.x[0];
                Scalar a2 = syn_opt.x[1];
                Scalar a3 = (1.0 - syn_weights[0] * a1 - syn_weights[1] * a2) / syn_weights[2];
                if (a1 >= 0.01 && a2 >= a1 && a3 >= a2) {
                    syn_rates[0] = a1;
                    syn_rates[1] = a2;
                    syn_rates[2] = a3;
                }
            }

            // Block 3: Multiple hits optimization
            if (allow_triple) {
                auto mh_opt = NelderMead2D::minimize([&](Scalar d, Scalar p) {
                    if (d < 0.0 || d > 10.0 || p < 0.0 || p > 10.0) return 1e20;
                    return -evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, d, p, eps_sink);
                }, delta, psi, 0.0, 10.0, 1e-4, 30);
                delta = std::max(0.0, mh_opt.a);
                psi = std::max(0.0, mh_opt.b);
            } else if (allow_double) {
                auto d_opt = NelderMead2D::minimize([&](Scalar d, Scalar /*unused*/) {
                    if (d < 0.0 || d > 10.0) return 1e20;
                    return -evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, d, 0.0, eps_sink);
                }, delta, 0.0, 0.0, 10.0, 1e-4, 25);
                delta = std::max(0.0, d_opt.a);
            }

            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << " (ECM " << (ecm_iter + 1) << "/" << max_ecm << "): lnL = "
                   << std::fixed << std::setprecision(2) << evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink);
                if (!is_constrained && !omegas.empty()) {
                    ss << ", ω" << K << " = " << std::setprecision(2) << omegas.back();
                }
                if (use_error_sink) {
                    ss << ", ε = " << std::setprecision(4) << eps_sink;
                }
                pbar->set_status(ss.str());
                pbar->tick();
            }
        }

        // Branch length initialization and refinement
        std::vector<Scalar> bl(tree.num_nodes(), 0.0);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                bl[node.id] = s * (ecm_bl_ptr ? (*ecm_bl_ptr)[node.id] : branch_tau[node.id]);
            }
        }

        if (refine_branches) {
            int num_cycles = (K == 1) ? 2 : 4;
            for (int cycle = 0; cycle < num_cycles; ++cycle) {
                bl = refine_branch_lengths(bl, omegas, weights, syn_rates, syn_weights, delta, psi, 25, eps_sink);

                if (K >= 2 || use_error_sink) {
                    for (int em = 0; em < 4; ++em) {
                        weights = compute_expected_weights(1.0, omegas, weights, syn_rates, syn_weights, &bl, delta, psi, eps_sink, use_error_sink ? &eps_sink : nullptr, max_eps_sink);
                        if (srv && M >= 2) {
                            syn_weights = compute_expected_syn_weights(1.0, omegas, weights, syn_rates, syn_weights, &bl, delta, psi, eps_sink);
                        }
                    }
                }

                // Polish omegas given refined branch lengths
                if (is_constrained) {
                    if (K == 2) {
                        auto [w1_opt, val] = Brent1D::minimize([&](Scalar w1) {
                            return -evaluate_log_l(1.0, {w1, 1.0}, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink);
                        }, 0.0, omegas[0], 1.0, 1e-4, 40);
                        omegas[0] = std::clamp(w1_opt, 0.0, 1.0);
                        omegas[1] = 1.0;
                    } else if (K == 3) {
                        auto opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                            if (x[0] > x[1]) return 1e20;
                            return -evaluate_log_l(1.0, {x[0], x[1], 1.0}, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink);
                        }, {omegas[0], omegas[1]}, {0.0, 0.0}, {1.0, 1.0}, 1e-4, 50);
                        omegas[0] = opt.x[0];
                        omegas[1] = opt.x[1];
                        omegas[2] = 1.0;
                    }
                } else {
                    Scalar max_omega = use_error_sink ? 50.0 : 2000.0;
                    if (K == 2) {
                        auto opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                            return -evaluate_log_l(1.0, {x[0], x[1]}, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink);
                        }, {omegas[0], omegas[1]}, {0.0, 1.0}, {1.0, max_omega}, 1e-4, 50);
                        omegas[0] = opt.x[0];
                        omegas[1] = opt.x[1];
                    } else if (K == 3) {
                        auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                            if (x[0] > x[1] || x[1] > x[2]) return 1e20;
                            return -evaluate_log_l(1.0, {x[0], x[1], x[2]}, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink);
                        }, {omegas[0], omegas[1], omegas[2]}, {0.0, 0.0, 1.0}, {1.0, max_omega, max_omega}, 1e-4, 80);
                        omegas[0] = opt.x[0];
                        omegas[1] = opt.x[1];
                        omegas[2] = opt.x[2];
                    }
                }

                // If SRV, polish syn_rates given bl
                if (srv && M == 3) {
                    auto syn_opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                        Scalar a1 = x[0], a2 = x[1];
                        Scalar q1 = syn_weights[0], q2 = syn_weights[1], q3 = syn_weights[2];
                        Scalar a3 = (1.0 - q1 * a1 - q2 * a2) / q3;
                        if (a1 < 0.01 || a2 < a1 || a3 < a2 || a3 > 50.0) return 1e20;
                        return -evaluate_log_l(1.0, omegas, weights, {a1, a2, a3}, syn_weights, &bl, nullptr, delta, psi, eps_sink);
                    }, {syn_rates[0], syn_rates[1]}, {0.01, 0.01}, {5.0, 10.0}, 1e-4, 40);

                    Scalar a1 = syn_opt.x[0];
                    Scalar a2 = syn_opt.x[1];
                    Scalar a3 = (1.0 - syn_weights[0] * a1 - syn_weights[1] * a2) / syn_weights[2];
                    if (a1 >= 0.01 && a2 >= a1 && a3 >= a2) {
                        syn_rates[0] = a1;
                        syn_rates[1] = a2;
                        syn_rates[2] = a3;
                    }
                }
            }

            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << ": Branch lengths refined";
                pbar->set_status(ss.str());
                pbar->tick();
            }

            if (allow_triple) {
                auto mh_polish = NelderMead2D::minimize([&](Scalar d, Scalar p) {
                    if (d < 0.0 || d > 10.0 || p < 0.0 || p > 10.0) return 1e20;
                    return -evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, nullptr, d, p, eps_sink);
                }, delta, psi, 0.0, 10.0, 1e-4, 25);
                delta = std::max(0.0, mh_polish.a);
                psi = std::max(0.0, mh_polish.b);
            } else if (allow_double) {
                auto d_polish = NelderMead2D::minimize([&](Scalar d, Scalar /*unused*/) {
                    if (d < 0.0 || d > 10.0) return 1e20;
                    return -evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, nullptr, d, 0.0, eps_sink);
                }, delta, 0.0, 0.0, 10.0, 1e-4, 25);
                delta = std::max(0.0, d_polish.a);
            }
            if (allow_double && pbar) {
                std::ostringstream ss;
                ss << model_desc << ": Multiple hits refined (δ=" << std::fixed << std::setprecision(3) << delta << ")";
                pbar->set_status(ss.str());
                pbar->tick();
            }
        }

        Scalar best_lnL = evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, &fit.site_log_likelihoods, delta, psi, eps_sink);
        fit.log_likelihood = best_lnL;
        fit.tree_scale = s;
        fit.branch_lengths = bl;
        fit.delta = delta;
        fit.psi = psi;
        fit.error_sink_weight = eps_sink;

        // Compute substitution rates and fractions
        Scalar total_rate = 0.0;
        Scalar double_rate = 0.0;
        Scalar triple_rate = 0.0;
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            MG94Matrix mat;
            mat.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            auto sr = mat.compute_substitution_rates(*code);
            total_rate += weights[k] * sr.total_rate;
            double_rate += weights[k] * sr.rate_double;
            triple_rate += weights[k] * sr.rate_triple;
        }
        fit.frac_delta = (total_rate > 1e-12) ? (double_rate / total_rate) : 0.0;
        fit.frac_psi = (total_rate > 1e-12) ? (triple_rate / total_rate) : 0.0;

        fit.test_distribution.omegas = omegas;
        fit.test_distribution.weights = weights;
        fit.test_distribution.error_sink_weight = eps_sink;
        if (srv) {
            fit.test_distribution.syn_rates = syn_rates;
            fit.test_distribution.syn_weights = syn_weights;
        }

        // Detect collapsed rate classes. annotations[k] annotates class k
        // directly (empty string = unannotated); a pairwise collapse note is
        // attached to the first class of the pair.
        fit.test_distribution.annotations.assign(K, "");
        for (size_t i = 0; i < K; ++i) {
            if (weights[i] < 0.005) {
                fit.test_distribution.annotations[i] = "Rate class " + std::to_string(i) + " collapsed (weight < 0.5%)";
            }
        }
        for (size_t i = 0; i + 1 < K; ++i) {
            if (std::abs(omegas[i+1] - omegas[i]) < 0.015 && fit.test_distribution.annotations[i].empty()) {
                fit.test_distribution.annotations[i] = "Rate classes " + std::to_string(i) + " and " + std::to_string(i+1) + " collapsed into identical omega";
            }
        }

        size_t N_codons = aln.num_codons;
        size_t num_branches = tree.num_nodes() - 1;
        size_t num_mixture_params = (2 * K - 1) - (is_constrained ? 1 : 0) + (use_error_sink ? 1 : 0);
        size_t num_syn_params = srv ? 2 * (M - 1) : 0;
        size_t num_mh_params = (allow_double ? 1 : 0) + (allow_triple ? 1 : 0);
        size_t K_params = 5 + (refine_branches ? num_branches : 1) + num_mixture_params + num_syn_params + num_mh_params;
        fit.aicc = -2.0 * fit.log_likelihood + 2.0 * K_params * N_codons / (N_codons - K_params - 1.0);

        return fit;
    }

    // Computes branch-site Error Sink Bayes Factors in a single O(patterns * branches) Inside-Outside pass
    void compute_error_sink_diagnostics(
        const BUSTEDFit& fit,
        Scalar bf_threshold,
        std::vector<Scalar>& out_site_max_bf,
        std::vector<BUSTEDErrorEvent>& out_flagged_events
    ) const {
        size_t num_sites = aln.num_codons;
        out_site_max_bf.assign(num_sites, 0.0);
        out_flagged_events.clear();

        Scalar eps = fit.error_sink_weight;
        if (eps <= 1e-9 || eps >= 1.0 - 1e-9) return;
        Scalar prior_odds = eps / (1.0 - eps);
        Scalar bio_scale = 1.0 - eps;

        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = fit.test_distribution.omegas.size();
        std::vector<Scalar> syn_rates = fit.test_distribution.syn_rates.empty() ? std::vector<Scalar>{1.0} : fit.test_distribution.syn_rates;
        std::vector<Scalar> syn_weights = fit.test_distribution.syn_weights.empty() ? std::vector<Scalar>{1.0} : fit.test_distribution.syn_weights;
        size_t M = syn_rates.size();

        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = fit.test_distribution.omegas[k];
            p.delta = fit.delta;
            p.psi = fit.psi;
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar bl = fit.branch_lengths[node.id];
                    P_mix[m][node.id] = Matrix::Zero(S, S);
                    Scalar use_bio = is_test_branch[node.id] ? bio_scale : 1.0;
                    for (size_t k = 0; k < K; ++k) {
                        P_mix[m][node.id] += (use_bio * fit.test_distribution.weights[k]) * mg[k].transition_matrix(alpha_m * bl);
                    }
                    if (is_test_branch[node.id]) {
                        P_mix[m][node.id].noalias() += eps * P_sink;
                    }
                }
            }
        }

        // Per-pattern, per-node posterior probability of the saturation sink (underflow-safe via log-likelihood)
        std::vector<std::vector<Scalar>> pat_node_post(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            const auto& pattern = aln.patterns[p];
            std::vector<LikelihoodEngine::InsideOutsideResult> io_m(M);
            Scalar max_ll = -std::numeric_limits<Scalar>::infinity();
            for (size_t m = 0; m < M; ++m) {
                if (syn_weights[m] <= 1e-12) continue;
                io_m[m] = LikelihoodEngine::compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S
                );
                if (io_m[m].likelihood > 0.0) {
                    Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                    if (term > max_ll) max_ll = term;
                }
            }
            if (max_ll > -1e19) {
                Scalar sum_exp = 0.0;
                std::vector<Scalar> post_m(M, 0.0);
                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                        Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood - max_ll);
                        post_m[m] = e;
                        sum_exp += e;
                    }
                }
                if (sum_exp > 0.0) {
                    for (size_t m = 0; m < M; ++m) post_m[m] /= sum_exp;
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id && is_test_branch[node.id]) {
                            int nid = node.id;
                            Scalar sink_ratio = 0.0;
                            for (size_t m = 0; m < M; ++m) {
                                if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                    sink_ratio += post_m[m] * (eval_sink_bilinear(io_m[m].V[nid], io_m[m].D[nid]) / io_m[m].likelihood);
                                }
                            }
                            Scalar post = std::clamp(eps * sink_ratio, 0.0, 1.0 - 1e-9);
                            pat_node_post[p][nid] = post;
                        }
                    }
                }
            }
        }

        for (size_t s_idx = 0; s_idx < num_sites; ++s_idx) {
            size_t p_idx = aln.site_to_pattern[s_idx];
            Scalar max_bf = 0.0;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && is_test_branch[node.id]) {
                    Scalar post = pat_node_post[p_idx][node.id];
                    Scalar post_odds = post / std::max(1e-12, 1.0 - post);
                    Scalar bf = post_odds / prior_odds;
                    if (bf > max_bf) max_bf = bf;
                    if (bf >= bf_threshold) {
                        out_flagged_events.push_back({s_idx + 1, node.name, post, bf});
                    }
                }
            }
            out_site_max_bf[s_idx] = max_bf;
        }
    }

    // Fast 4-Gate Surprisal & Score Automatic Model Selector (BUSTED-Auto)
    // Diagnoses Error Sink (O(1) saturation + spatial adjacency), Multi-Hits (O(t_b) quarantined tangent score),
    // Synonymous Rate Variation (Neyman-Scott overdispersion + SQUAREM probe), and Omega classes K in a single baseline pass.
    BUSTEDAutoDiagnostics diagnose_and_select_model(BUSTEDSettings& settings) const {
        BUSTEDAutoDiagnostics diag;
        diag.enabled = true;

        size_t num_sites = aln.num_codons;
        size_t num_patterns = aln.patterns.size();
        size_t num_nodes = tree.num_nodes();
        int S = code->num_sense_codons;

        // Step 0: Fast K=2 baseline warm-up (M0: M=1, delta=0, psi=0, eps=0)
        std::vector<Scalar> omegas0 = {std::min(0.05, 0.25 * mg94_omega), 1.0};
        std::vector<Scalar> weights0 = {0.80, 0.20};
        Scalar s0 = 1.0;
        for (int it = 0; it < 2; ++it) {
            weights0 = compute_expected_weights(s0, omegas0, weights0);
            auto opt0 = NelderMead2D::minimize([&](Scalar w1, Scalar sc) {
                return -evaluate_log_l(sc, {w1, 1.0}, weights0);
            }, omegas0[0], s0, 1e-4, 1.0, 1e-3, 20);
            omegas0[0] = opt0.a;
            s0 = opt0.b;
        }
        weights0 = compute_expected_weights(s0, omegas0, weights0);
        Scalar base_ll = evaluate_log_l(s0, omegas0, weights0);
        Scalar mean_omega0 = weights0[0] * omegas0[0] + weights0[1] * omegas0[1];

        // Build baseline branch transition matrices P0[node_id]
        std::vector<MG94Matrix> mg0(2);
        for (size_t k = 0; k < 2; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas0[k];
            mg0[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        std::vector<Matrix> P0(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar bl = s0 * branch_tau[node.id];
                P0[node.id] = weights0[0] * mg0[0].transition_matrix(bl) + weights0[1] * mg0[1].transition_matrix(bl);
            }
        }

        // Build pure 2-hit (Q2_off), 3-hit (Q3_off), and synonymous (Q_syn) tangent operators at mean_omega0
        Scalar r_nuc[4][4] = {
            {0.0, base_params.theta_AC, 1.0, base_params.theta_AT},
            {base_params.theta_AC, 0.0, base_params.theta_CG, base_params.theta_CT},
            {1.0, base_params.theta_CG, 0.0, base_params.theta_GT},
            {base_params.theta_AT, base_params.theta_CT, base_params.theta_GT, 0.0}
        };
        Scalar q_scale = std::max(1e-6, mg0[1].scale_factor);
        Matrix Q2_off = Matrix::Zero(S, S);
        Matrix Q3_off = Matrix::Zero(S, S);
        Matrix Q_syn = Matrix::Zero(S, S);
        for (int i = 0; i < S; ++i) {
            Scalar syn_row_sum = 0.0;
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                const auto& diff = code->diff_matrix[i][j];
                Scalar mod = diff.synonymous ? 1.0 : mean_omega0;
                if (diff.num_diffs == 1 && diff.synonymous) {
                    Scalar rate1 = r_nuc[diff.diff_from[0]][diff.diff_to[0]] * aln.pos_nuc_frequencies(diff.diff_pos[0], diff.diff_to[0]) / q_scale;
                    Q_syn(i, j) = rate1;
                    syn_row_sum += rate1;
                } else if (diff.num_diffs == 2) {
                    Scalar rate2 = mod * r_nuc[diff.diff_from[0]][diff.diff_to[0]] * r_nuc[diff.diff_from[1]][diff.diff_to[1]] *
                                   aln.pos_nuc_frequencies(diff.diff_pos[0], diff.diff_to[0]) * aln.pos_nuc_frequencies(diff.diff_pos[1], diff.diff_to[1]) / q_scale;
                    Q2_off(i, j) = rate2;
                } else if (diff.num_diffs == 3) {
                    Scalar rate3 = mod * r_nuc[diff.diff_from[0]][diff.diff_to[0]] * r_nuc[diff.diff_from[1]][diff.diff_to[1]] * r_nuc[diff.diff_from[2]][diff.diff_to[2]] *
                                   aln.pos_nuc_frequencies(diff.diff_pos[0], diff.diff_to[0]) * aln.pos_nuc_frequencies(diff.diff_pos[1], diff.diff_to[1]) * aln.pos_nuc_frequencies(diff.diff_pos[2], diff.diff_to[2]) / q_scale;
                    Q3_off(i, j) = rate3;
                }
            }
            Q_syn(i, i) = -syn_row_sum;
        }

        // Single parallel Inside-Outside pass on M0
        std::vector<std::vector<Scalar>> pat_R_sink(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        std::vector<std::vector<Scalar>> pat_T2(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        std::vector<std::vector<Scalar>> pat_T3(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        std::vector<Scalar> pat_syn_score(num_patterns, 0.0);
        std::vector<Scalar> pat_syn_info(num_patterns, 0.0);

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            const auto& pattern = aln.patterns[p];
            auto io = LikelihoodEngine::compute_inside_outside(
                tree, pattern, leaf_to_taxon, P0, aln.codon_frequencies_f3x4, S
            );
            Scalar L_p = io.likelihood;
            if (L_p > 0.0) {
                Scalar inv_L = 1.0 / L_p;
                Scalar u_syn = 0.0;
                Scalar v_syn = 0.0;
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) {
                        int nid = node.id;
                        Scalar bl = s0 * branch_tau[nid];
                        if (is_test_branch[nid]) {
                            Scalar sink_val = eval_sink_bilinear(io.V[nid], io.D[nid]);
                            pat_R_sink[p][nid] = sink_val * inv_L;
                            pat_T2[p][nid] = bl * io.V[nid].dot(Q2_off * io.D[nid]) * inv_L;
                            pat_T3[p][nid] = bl * io.V[nid].dot(Q3_off * io.D[nid]) * inv_L;
                        }
                        Scalar branch_u = bl * io.V[nid].dot(Q_syn * (P0[nid] * io.D[nid])) * inv_L;
                        u_syn += branch_u;
                        v_syn += branch_u * branch_u;
                    }
                }
                pat_syn_score[p] = u_syn;
                pat_syn_info[p] = std::max(1e-6, v_syn);
            }
        }

        // Gate 1: Evaluate Error Sink Surprisal & Spatial Adjacency
        Scalar max_R_sink = 0.0;
        int max_sink_nid = -1;
        size_t candidate_count = 0;
        size_t adjacent_clusters = 0;
        Scalar eps_init = 0.002;
        Scalar sink_resp_sum = 0.0;
        Scalar total_branch_sites = static_cast<Scalar>(num_sites * std::max<size_t>(1, num_test_branches));

        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && is_test_branch[node.id]) {
                int nid = node.id;
                int prev_high_site = -10;
                for (size_t s_idx = 0; s_idx < num_sites; ++s_idx) {
                    size_t p_idx = aln.site_to_pattern[s_idx];
                    Scalar r_val = pat_R_sink[p_idx][nid];
                    if (r_val > max_R_sink) {
                        max_R_sink = r_val;
                        max_sink_nid = nid;
                    }
                    Scalar bits = (r_val > 1.0) ? std::log2(r_val) : 0.0;
                    if (bits > 8.0) {
                        candidate_count++;
                    }
                    if (bits > 6.5) {
                        if (static_cast<int>(s_idx) - prev_high_site <= 2) {
                            adjacent_clusters++;
                        }
                        prev_high_site = static_cast<int>(s_idx);
                    }
                    Scalar resp = (eps_init * r_val) / ((1.0 - eps_init) + eps_init * r_val);
                    sink_resp_sum += resp;
                }
            }
        }

        diag.max_sink_bf = max_R_sink;
        diag.max_sink_surprisal_bits = (max_R_sink > 1.0) ? std::log2(max_R_sink) : 0.0;
        diag.sink_candidate_sites = candidate_count;
        diag.sink_adjacent_clusters = adjacent_clusters;

        Scalar eps_probe = std::clamp(sink_resp_sum / total_branch_sites, 1e-6, settings.error_sink_max_weight);
        Scalar sink_probe_ll = base_ll;
        if (diag.max_sink_surprisal_bits > 6.5) {
            sink_probe_ll = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, 0.0, 0.0, eps_probe);
            // If a single branch has extreme surprisal (which inflates its MG94 branch length), also check 1-step branch refinement under sink
            if (max_sink_nid >= 0 && diag.max_sink_surprisal_bits >= 10.0) {
                std::vector<Scalar> probe_bl(num_nodes, 0.0);
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) probe_bl[node.id] = s0 * branch_tau[node.id];
                }
                auto ref_bl = refine_branch_lengths(probe_bl, omegas0, weights0, {1.0}, {1.0}, 0.0, 0.0, 4, eps_probe);
                auto base_ref_bl = refine_branch_lengths(probe_bl, omegas0, weights0, {1.0}, {1.0}, 0.0, 0.0, 4, 0.0);
                Scalar ref_sink_ll = evaluate_log_l(1.0, omegas0, weights0, {1.0}, {1.0}, &ref_bl, nullptr, 0.0, 0.0, eps_probe);
                Scalar ref_base_ll = evaluate_log_l(1.0, omegas0, weights0, {1.0}, {1.0}, &base_ref_bl, nullptr, 0.0, 0.0, 0.0);
                if (ref_sink_ll - ref_base_ll > sink_probe_ll - base_ll) {
                    sink_probe_ll = base_ll + (ref_sink_ll - ref_base_ll);
                }
                // Recompute max_R_sink on the un-inflated branch length
                if (ref_bl[max_sink_nid] < probe_bl[max_sink_nid] * 0.5 && ref_bl[max_sink_nid] > 1e-8) {
                    Scalar ratio_gain = std::pow(probe_bl[max_sink_nid] / ref_bl[max_sink_nid], 2.0);
                    max_R_sink *= ratio_gain;
                    diag.max_sink_bf = max_R_sink;
                    diag.max_sink_surprisal_bits = (max_R_sink > 1.0) ? std::log2(max_R_sink) : 0.0;
                }
            }
        }
        diag.sink_delta_ll = std::max(0.0, sink_probe_ll - base_ll);

        // Activate Error Sink if Delta LL > 2.2 AND either extreme short-branch surprisal (>= 13.0 bits, BF >= 8192) or spatial adjacency cluster
        if (diag.sink_delta_ll > 2.2 && (diag.max_sink_surprisal_bits >= 13.0 || (diag.sink_adjacent_clusters >= 1 && diag.max_sink_surprisal_bits >= 8.0))) {
            diag.activated_error_sink = true;
        }

        // Gate 2: Evaluate Multi-Hit Tangent Score & Probe (Quarantining high-R_sink branch-sites if Error Sink is active)
        Scalar eps_active = diag.activated_error_sink ? eps_probe : 0.0;
        Scalar s_delta = 0.0;
        Scalar s_psi = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            Scalar w_p = aln.patterns[p].weight;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && is_test_branch[node.id]) {
                    int nid = node.id;
                    Scalar r_val = pat_R_sink[p][nid];
                    Scalar r_sink = (eps_active > 0.0) ? ((eps_active * r_val) / ((1.0 - eps_active) + eps_active * r_val)) : 0.0;
                    Scalar clean_wt = w_p * (1.0 - r_sink);
                    s_delta += clean_wt * pat_T2[p][nid];
                    s_psi += clean_wt * pat_T3[p][nid];
                }
            }
        }
        diag.double_hit_score = s_delta;
        diag.triple_hit_score = s_psi;

        Scalar mh_base_ll = diag.activated_error_sink ? sink_probe_ll : base_ll;
        if (s_delta > 0.15) {
            Scalar ll_d1 = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, 0.06, 0.0, eps_active);
            Scalar ll_d2 = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, 0.20, 0.0, eps_active);
            Scalar best_d = (ll_d2 > ll_d1) ? 0.20 : 0.06;
            Scalar best_d_ll = std::max(ll_d1, ll_d2);
            diag.double_hit_delta_ll = std::max(0.0, best_d_ll - mh_base_ll);

            // Check if diffuse 2-hit evolution (with eps=0) strictly outperforms an isolated non-clustered error sink
            if (diag.activated_error_sink && diag.sink_adjacent_clusters == 0 && diag.max_sink_surprisal_bits < 16.0) {
                Scalar pure_d_ll = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, best_d, 0.0, 0.0);
                if (pure_d_ll > sink_probe_ll + 0.5 && (pure_d_ll - base_ll) > 2.2) {
                    diag.activated_error_sink = false;
                    eps_active = 0.0;
                    mh_base_ll = base_ll;
                    best_d_ll = pure_d_ll;
                    diag.double_hit_delta_ll = pure_d_ll - base_ll;
                }
            }

            if (diag.double_hit_delta_ll > 2.2) {
                diag.selected_multiple_hits = "Double";
                if (s_psi > 0.05) {
                    Scalar ll_p = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, best_d, 0.06, eps_active);
                    diag.triple_hit_delta_ll = std::max(0.0, ll_p - best_d_ll);
                    if (diag.triple_hit_delta_ll > 2.2) {
                        diag.selected_multiple_hits = "Double+Triple";
                    }
                }
            }
        }

        // Gate 3: Evaluate Synonymous Rate Variation (SRV) Overdispersion & Probe
        Scalar num_disp = 0.0;
        Scalar den_disp = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            Scalar w_p = aln.patterns[p].weight;
            Scalar u2 = pat_syn_score[p] * pat_syn_score[p];
            Scalar v = pat_syn_info[p];
            num_disp += w_p * (u2 - v);
            den_disp += w_p * v * v;
        }
        diag.srv_overdispersion_z = (den_disp > 1e-12) ? (num_disp / std::sqrt(2.0 * den_disp)) : 0.0;

        Scalar active_d = (diag.selected_multiple_hits != "None") ? 0.10 : 0.0;
        Scalar active_p = (diag.selected_multiple_hits == "Double+Triple") ? 0.04 : 0.0;
        Scalar pre_srv_ll = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, active_d, active_p, eps_active);

        std::vector<Scalar> syn_r3 = {0.20, 0.75, 2.05};
        std::vector<Scalar> syn_w3 = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        for (int it = 0; it < 2; ++it) {
            syn_w3 = compute_expected_syn_weights(s0, omegas0, weights0, syn_r3, syn_w3, nullptr, active_d, active_p, eps_active);
        }
        Scalar srv3_ll = evaluate_log_l(s0, omegas0, weights0, syn_r3, syn_w3, nullptr, nullptr, active_d, active_p, eps_active);
        diag.srv_delta_ll = std::max(0.0, srv3_ll - pre_srv_ll);

        if (diag.srv_delta_ll > 4.5 && syn_w3[0] >= 0.02 && syn_w3[1] >= 0.02 && syn_w3[2] >= 0.02) {
            diag.activated_srv = true;
            diag.selected_syn_rates = 3;
        } else if (diag.srv_delta_ll > 2.2) {
            std::vector<Scalar> syn_r2 = {0.35, 1.65};
            std::vector<Scalar> syn_w2 = {0.50, 0.50};
            syn_w2 = compute_expected_syn_weights(s0, omegas0, weights0, syn_r2, syn_w2, nullptr, active_d, active_p, eps_active);
            Scalar srv2_ll = evaluate_log_l(s0, omegas0, weights0, syn_r2, syn_w2, nullptr, nullptr, active_d, active_p, eps_active);
            if (srv2_ll - pre_srv_ll > 2.2 && syn_w2[0] >= 0.02 && syn_w2[1] >= 0.02) {
                diag.activated_srv = true;
                diag.selected_syn_rates = 2;
                diag.srv_delta_ll = srv2_ll - pre_srv_ll;
            }
        }

        // Gate 4: Fast Unrefined Null Screen for Optimal Omega Classes K in {2, 3}
        BUSTEDFit null_k2 = fit_model(true, 2, false, diag.activated_srv, diag.selected_syn_rates, diag.selected_multiple_hits, nullptr, "", diag.activated_error_sink, settings.error_sink_max_weight);
        BUSTEDFit null_k3 = fit_model(true, 3, false, diag.activated_srv, diag.selected_syn_rates, diag.selected_multiple_hits, nullptr, "", diag.activated_error_sink, settings.error_sink_max_weight);

        bool k3_collapsed = (std::abs(null_k3.test_distribution.omegas[1] - null_k3.test_distribution.omegas[0]) < 0.025) ||
                            (null_k3.test_distribution.weights[0] < 0.015) ||
                            (null_k3.test_distribution.weights[1] < 0.015);
        if (!k3_collapsed && (null_k3.log_likelihood - null_k2.log_likelihood > 2.5)) {
            diag.selected_k = 3;
            diag.k_rationale = "K=3 retained (distinct two-rate purifying background, delta lnL = +" +
                               std::to_string(null_k3.log_likelihood - null_k2.log_likelihood) + ")";
        } else {
            diag.selected_k = 2;
            diag.k_rationale = k3_collapsed
                ? "K=2 selected (pruned collapsed/low-weight 3rd background rate class)"
                : "K=2 selected via AICc parsimony (K=3 gain < threshold)";
        }

        // Apply selected specification to settings
        settings.error_sink = diag.activated_error_sink;
        settings.multiple_hits = diag.selected_multiple_hits;
        settings.srv = diag.activated_srv;
        settings.num_syn_rate_classes = diag.selected_syn_rates;
        settings.num_rate_classes = diag.selected_k;
        settings.auto_select_k = false;

        return diag;
    }

    // Step-up model selection across K in {1, ..., max_k} using AICc on the Null model
    std::pair<size_t, std::vector<BUSTEDFit>> select_optimal_k(
        size_t max_k = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3,
        const std::string& multiple_hits = "None",
        ProgressBar* pbar = nullptr,
        bool use_error_sink = false,
        Scalar max_eps_sink = 0.01
    ) const {
        std::vector<BUSTEDFit> null_fits;
        null_fits.reserve(max_k);
        for (size_t k = 1; k <= max_k; ++k) {
            std::string desc = "Null K=" + std::to_string(k);
            null_fits.push_back(fit_model(true, k, refine_branches, srv, num_syn_rates, multiple_hits, pbar, desc, use_error_sink, max_eps_sink));
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

    // Two-set likelihood evaluation for BUSTED with separate Test (partition 0) and Background/Comparison (partition 1) K-class distributions
    Scalar evaluate_two_set_log_l(
        Scalar s,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        std::vector<Scalar>* out_pattern_ll = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K_t = omegas_test.size();
        size_t K_b = omegas_bg.size();
        size_t M = syn_rates.size();
        Scalar bio_scale_t = 1.0 - eps_sink_test;
        Scalar bio_scale_b = 1.0 - eps_sink_bg;

        std::vector<MG94Matrix> mg_t(K_t);
        for (size_t k = 0; k < K_t; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas_test[k];
            p.delta = delta;
            p.psi = psi;
            mg_t[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<MG94Matrix> mg_b(K_b);
        for (size_t k = 0; k < K_b; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas_bg[k];
            p.delta = delta;
            p.psi = psi;
            mg_b[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        MG94Matrix mg_unclass;
        if (num_unclassified_branches > 0) {
            MG94Parameters p_u = base_params;
            p_u.alpha = 1.0;
            p_u.beta = mg94_omega;
            p_u.delta = delta;
            p_u.psi = psi;
            mg_unclass.update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[nid] : (s * branch_tau[nid]);
                    int part = (static_cast<size_t>(nid) < branch_partition.size()) ? branch_partition[nid] : (is_test_branch[nid] ? 0 : 1);
                    if (part == 0) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_test > 0.0) ? bio_scale_t : 1.0;
                        for (size_t k = 0; k < K_t; ++k) {
                            P_mix[m][nid] += (use_bio * weights_test[k]) * mg_t[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_test > 0.0) {
                            P_mix[m][nid].noalias() += eps_sink_test * P_sink;
                        }
                    } else if (part == 1) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_bg > 0.0) ? bio_scale_b : 1.0;
                        for (size_t k = 0; k < K_b; ++k) {
                            P_mix[m][nid] += (use_bio * weights_bg[k]) * mg_b[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_bg > 0.0) {
                            P_mix[m][nid].noalias() += eps_sink_bg * P_sink;
                        }
                    } else {
                        P_mix[m][nid] = mg_unclass.transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        if (out_pattern_ll) {
            out_pattern_ll->assign(num_patterns, 0.0);
        }

        Scalar total_log_l = 0.0;
        #pragma omp parallel
        {
            auto [node_L, node_log_scale] = LikelihoodEngine::make_pattern_scratch<Vector>(num_nodes, S);
            #pragma omp for reduction(+:total_log_l) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pll = -1e20;
                if (M == 1) {
                    LikelihoodEngine::eval_pattern_postorder(tree, pattern, leaf_to_taxon, P_mix[0], S, node_L, node_log_scale);
                    Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                    if (cond_L > 0.0) {
                        pll = std::log(cond_L) + node_log_scale[tree.root_id];
                    }
                } else {
                    std::vector<Scalar> terms(M, -std::numeric_limits<Scalar>::infinity());
                    for (size_t m = 0; m < M; ++m) {
                        LikelihoodEngine::eval_pattern_postorder(tree, pattern, leaf_to_taxon, P_mix[m], S, node_L, node_log_scale);
                        Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                        if (cond_L > 0.0) {
                            terms[m] = std::log(syn_weights[m]) + std::log(cond_L) + node_log_scale[tree.root_id];
                        }
                    }

                    Scalar lse = LikelihoodEngine::logsumexp(terms);
                    if (std::isfinite(lse)) pll = lse;
                }
                total_log_l += pattern.weight * pll;
                if (out_pattern_ll) {
                    (*out_pattern_ll)[p] = pll;
                }
            }
        }
        return total_log_l;
    }

    // Two-set Expectation Step: simultaneously updates Test (partition 0) and Background (partition 1) mixture weights in one Inside-Outside pass
    std::pair<std::vector<Scalar>, std::vector<Scalar>> compute_expected_two_set_weights(
        Scalar s,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0,
        Scalar* out_new_eps_test = nullptr,
        Scalar* out_new_eps_bg = nullptr,
        Scalar max_eps_sink = 0.01
    ) const {
        size_t K_t = omegas_test.size();
        size_t K_b = omegas_bg.size();
        size_t M = syn_rates.size();
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        Scalar bio_scale_t = (out_new_eps_test != nullptr) ? (1.0 - eps_sink_test) : 1.0;
        Scalar bio_scale_b = (out_new_eps_bg != nullptr) ? (1.0 - eps_sink_bg) : 1.0;

        std::vector<std::vector<std::vector<Matrix>>> P_comp_t(M, std::vector<std::vector<Matrix>>(K_t, std::vector<Matrix>(num_nodes)));
        std::vector<std::vector<std::vector<Matrix>>> P_comp_b(M, std::vector<std::vector<Matrix>>(K_b, std::vector<Matrix>(num_nodes)));
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));

        for (size_t k = 0; k < K_t; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas_test[k];
            p.delta = delta;
            p.psi = psi;
            MG94Matrix mg;
            mg.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            for (size_t m = 0; m < M; ++m) {
                Scalar alpha_m = syn_rates[m];
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id && branch_partition[node.id] == 0) {
                        Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[node.id] : (s * branch_tau[node.id]);
                        P_comp_t[m][k][node.id] = mg.transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        for (size_t k = 0; k < K_b; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas_bg[k];
            p.delta = delta;
            p.psi = psi;
            MG94Matrix mg;
            mg.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            for (size_t m = 0; m < M; ++m) {
                Scalar alpha_m = syn_rates[m];
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id && branch_partition[node.id] == 1) {
                        Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[node.id] : (s * branch_tau[node.id]);
                        P_comp_b[m][k][node.id] = mg.transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        MG94Matrix mg_unclass;
        if (num_unclassified_branches > 0) {
            MG94Parameters p_u = base_params;
            p_u.alpha = 1.0;
            p_u.beta = mg94_omega;
            p_u.delta = delta;
            p_u.psi = psi;
            mg_unclass.update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    int part = branch_partition[nid];
                    if (part == 0) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_t; ++k) {
                            P_mix[m][nid] += (bio_scale_t * weights_test[k]) * P_comp_t[m][k][nid];
                        }
                        if (out_new_eps_test != nullptr) {
                            P_mix[m][nid].noalias() += eps_sink_test * P_sink;
                        }
                    } else if (part == 1) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_b; ++k) {
                            P_mix[m][nid] += (bio_scale_b * weights_bg[k]) * P_comp_b[m][k][nid];
                        }
                        if (out_new_eps_bg != nullptr) {
                            P_mix[m][nid].noalias() += eps_sink_bg * P_sink;
                        }
                    } else {
                        Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[nid] : (s * branch_tau[nid]);
                        P_mix[m][nid] = mg_unclass.transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        Vector accum_t = Vector::Zero(K_t);
        Vector accum_b = Vector::Zero(K_b);
        Scalar sink_accum_t = 0.0, sink_accum_b = 0.0;
        Scalar total_wt_t = 0.0, total_wt_b = 0.0;

        #pragma omp parallel
        {
            Vector local_t = Vector::Zero(K_t);
            Vector local_b = Vector::Zero(K_b);
            Scalar local_sink_t = 0.0, local_sink_b = 0.0;
            Scalar local_wt_t = 0.0, local_wt_b = 0.0;

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                std::vector<LikelihoodEngine::InsideOutsideResult> io_m(M);
                Scalar max_ll = -std::numeric_limits<Scalar>::infinity();

                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
                    io_m[m] = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S
                    );
                    if (io_m[m].likelihood > 0.0) {
                        Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                        if (term > max_ll) max_ll = term;
                    }
                }

                if (max_ll > -1e19) {
                    Scalar sum_exp = 0.0;
                    std::vector<Scalar> post_m(M, 0.0);
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                            Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood - max_ll);
                            post_m[m] = e;
                            sum_exp += e;
                        }
                    }
                    if (sum_exp > 0.0) {
                        for (size_t m = 0; m < M; ++m) post_m[m] /= sum_exp;
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int nid = node.id;
                                int part = branch_partition[nid];
                                if (part == 0) {
                                    for (size_t k = 0; k < K_t; ++k) {
                                        Scalar ratio_k = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                ratio_k += post_m[m] * (io_m[m].V[nid].dot(P_comp_t[m][k][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_t(k) += pattern.weight * (bio_scale_t * weights_test[k]) * ratio_k;
                                    }
                                    if (out_new_eps_test != nullptr) {
                                        Scalar sink_ratio = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                sink_ratio += post_m[m] * (eval_sink_bilinear(io_m[m].V[nid], io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_sink_t += pattern.weight * eps_sink_test * sink_ratio;
                                    }
                                    local_wt_t += pattern.weight;
                                } else if (part == 1) {
                                    for (size_t k = 0; k < K_b; ++k) {
                                        Scalar ratio_k = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                ratio_k += post_m[m] * (io_m[m].V[nid].dot(P_comp_b[m][k][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_b(k) += pattern.weight * (bio_scale_b * weights_bg[k]) * ratio_k;
                                    }
                                    if (out_new_eps_bg != nullptr) {
                                        Scalar sink_ratio = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                sink_ratio += post_m[m] * (eval_sink_bilinear(io_m[m].V[nid], io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_sink_b += pattern.weight * eps_sink_bg * sink_ratio;
                                    }
                                    local_wt_b += pattern.weight;
                                }
                            }
                        }
                    }
                }
            }
            #pragma omp critical
            {
                accum_t += local_t;
                accum_b += local_b;
                sink_accum_t += local_sink_t;
                sink_accum_b += local_sink_b;
                total_wt_t += local_wt_t;
                total_wt_b += local_wt_b;
            }
        }

        std::vector<Scalar> new_wt_t = weights_test;
        if (total_wt_t > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t k = 0; k < K_t; ++k) {
                new_wt_t[k] = std::clamp(accum_t(k) / total_wt_t, 1e-6, 1.0 - 1e-6);
                sum_w += new_wt_t[k];
            }
            for (size_t k = 0; k < K_t; ++k) new_wt_t[k] /= sum_w;
            if (out_new_eps_test != nullptr) {
                *out_new_eps_test = std::clamp(sink_accum_t / total_wt_t, 1e-7, max_eps_sink);
            }
        }

        std::vector<Scalar> new_wt_b = weights_bg;
        if (total_wt_b > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t k = 0; k < K_b; ++k) {
                new_wt_b[k] = std::clamp(accum_b(k) / total_wt_b, 1e-6, 1.0 - 1e-6);
                sum_w += new_wt_b[k];
            }
            for (size_t k = 0; k < K_b; ++k) new_wt_b[k] /= sum_w;
            if (out_new_eps_bg != nullptr) {
                *out_new_eps_bg = std::clamp(sink_accum_b / total_wt_b, 1e-7, max_eps_sink);
            }
        }

        return {new_wt_t, new_wt_b};
    }

    // Two-set synonymous rate E-step (with log-underflow scaling)
    std::vector<Scalar> compute_expected_two_set_syn_weights(
        Scalar s,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates,
        const std::vector<Scalar>& syn_weights,
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0
    ) const {
        size_t M = syn_rates.size();
        if (M <= 1) return {1.0};

        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K_t = omegas_test.size();
        size_t K_b = omegas_bg.size();
        Scalar bio_scale_t = 1.0 - eps_sink_test;
        Scalar bio_scale_b = 1.0 - eps_sink_bg;

        std::vector<MG94Matrix> mg_t(K_t), mg_b(K_b);
        for (size_t k = 0; k < K_t; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0; p.beta = omegas_test[k]; p.delta = delta; p.psi = psi;
            mg_t[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        for (size_t k = 0; k < K_b; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0; p.beta = omegas_bg[k]; p.delta = delta; p.psi = psi;
            mg_b[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        MG94Matrix mg_unclass;
        if (num_unclassified_branches > 0) {
            MG94Parameters p_u = base_params;
            p_u.alpha = 1.0; p_u.beta = mg94_omega; p_u.delta = delta; p_u.psi = psi;
            mg_unclass.update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[nid] : (s * branch_tau[nid]);
                    int part = branch_partition[nid];
                    if (part == 0) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_test > 0.0) ? bio_scale_t : 1.0;
                        for (size_t k = 0; k < K_t; ++k) {
                            P_mix[m][nid] += (use_bio * weights_test[k]) * mg_t[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_test > 0.0) P_mix[m][nid].noalias() += eps_sink_test * P_sink;
                    } else if (part == 1) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_bg > 0.0) ? bio_scale_b : 1.0;
                        for (size_t k = 0; k < K_b; ++k) {
                            P_mix[m][nid] += (use_bio * weights_bg[k]) * mg_b[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_bg > 0.0) P_mix[m][nid].noalias() += eps_sink_bg * P_sink;
                    } else {
                        P_mix[m][nid] = mg_unclass.transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        Vector syn_accum = Vector::Zero(M);
        Scalar total_pattern_weight = 0.0;

        #pragma omp parallel
        {
            Vector local_accum = Vector::Zero(M);
            Scalar local_pat_wt = 0.0;
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));
            std::vector<double> node_log_scale(num_nodes, 0.0);

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar max_term = -std::numeric_limits<Scalar>::infinity();
                std::vector<Scalar> terms(M, -std::numeric_limits<Scalar>::infinity());
                std::vector<bool> valid(M, false);

                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
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
                            node_log_scale[node_id] = 0.0;
                        } else {
                            double child_scale_sum = 0.0;
                            tree.eval_internal_node(node_id, P_mix[m], node_L, node_log_scale, node_L[node_id], child_scale_sum);
                            node_log_scale[node_id] = child_scale_sum;
                        }
                    }
                    Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                    if (cond_L > 0.0) {
                        Scalar term = std::log(syn_weights[m]) + std::log(cond_L) + node_log_scale[tree.root_id];
                        terms[m] = term;
                        valid[m] = true;
                        if (term > max_term) max_term = term;
                    }
                }

                if (max_term > -1e19) {
                    Scalar sum_exp = 0.0;
                    std::vector<Scalar> rho(M, 0.0);
                    for (size_t m = 0; m < M; ++m) {
                        if (valid[m]) {
                            rho[m] = std::exp(terms[m] - max_term);
                            sum_exp += rho[m];
                        }
                    }
                    if (sum_exp > 0.0) {
                        for (size_t m = 0; m < M; ++m) {
                            local_accum(m) += pattern.weight * (rho[m] / sum_exp);
                        }
                        local_pat_wt += pattern.weight;
                    }
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
            for (size_t m = 0; m < M; ++m) new_syn_weights[m] /= sum_w;
        } else {
            new_syn_weights = syn_weights;
        }
        return new_syn_weights;
    }

    // Two-set analytical branch gradients
    std::pair<Scalar, std::vector<Scalar>> compute_two_set_branch_gradients(
        const std::vector<Scalar>& branch_lengths,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K_t = omegas_test.size();
        size_t K_b = omegas_bg.size();
        size_t M = syn_rates.size();
        Scalar bio_scale_t = 1.0 - eps_sink_test;
        Scalar bio_scale_b = 1.0 - eps_sink_bg;

        std::vector<MG94Matrix> mg_t(K_t), mg_b(K_b);
        for (size_t k = 0; k < K_t; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0; p.beta = omegas_test[k]; p.delta = delta; p.psi = psi;
            mg_t[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        for (size_t k = 0; k < K_b; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0; p.beta = omegas_bg[k]; p.delta = delta; p.psi = psi;
            mg_b[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        MG94Matrix mg_unclass;
        if (num_unclassified_branches > 0) {
            MG94Parameters p_u = base_params;
            p_u.alpha = 1.0; p_u.beta = mg94_omega; p_u.delta = delta; p_u.psi = psi;
            mg_unclass.update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<Matrix>> dP_mix(M, std::vector<Matrix>(num_nodes));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = branch_lengths[nid];
                    int part = branch_partition[nid];
                    if (part == 0) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        dP_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_test > 0.0) ? bio_scale_t : 1.0;
                        for (size_t k = 0; k < K_t; ++k) {
                            Matrix P_k = mg_t[k].transition_matrix(alpha_m * bl);
                            P_mix[m][nid] += (use_bio * weights_test[k]) * P_k;
                            dP_mix[m][nid] += (use_bio * weights_test[k]) * alpha_m * (mg_t[k].Q * P_k);
                        }
                        if (eps_sink_test > 0.0) P_mix[m][nid].noalias() += eps_sink_test * P_sink;
                    } else if (part == 1) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        dP_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_bg > 0.0) ? bio_scale_b : 1.0;
                        for (size_t k = 0; k < K_b; ++k) {
                            Matrix P_k = mg_b[k].transition_matrix(alpha_m * bl);
                            P_mix[m][nid] += (use_bio * weights_bg[k]) * P_k;
                            dP_mix[m][nid] += (use_bio * weights_bg[k]) * alpha_m * (mg_b[k].Q * P_k);
                        }
                        if (eps_sink_bg > 0.0) P_mix[m][nid].noalias() += eps_sink_bg * P_sink;
                    } else {
                        Matrix P_u = mg_unclass.transition_matrix(alpha_m * bl);
                        P_mix[m][nid] = P_u;
                        dP_mix[m][nid] = alpha_m * (mg_unclass.Q * P_u);
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
                if (M == 1) {
                    auto io = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_mix[0], aln.codon_frequencies_f3x4, S
                    );
                    if (io.likelihood > 0.0) {
                        local_ll += pattern.weight * io.log_likelihood;
                        Scalar inv_L = pattern.weight / io.likelihood;
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int32_t vid = node.id;
                                Scalar dL_dt = io.V[vid].dot(dP_mix[0][vid] * io.D[vid]);
                                local_grad[vid] += inv_L * dL_dt;
                            }
                        }
                    }
                } else {
                    std::vector<LikelihoodEngine::InsideOutsideResult> io_m(M);
                    std::vector<Scalar> terms(M, -std::numeric_limits<Scalar>::infinity());
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] <= 1e-12) continue;
                        io_m[m] = LikelihoodEngine::compute_inside_outside(
                            tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S
                        );
                        if (io_m[m].likelihood > 0.0) {
                            terms[m] = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                        }
                    }

                    Scalar pat_ll = LikelihoodEngine::logsumexp(terms);
                    if (std::isfinite(pat_ll)) {
                        local_ll += pattern.weight * pat_ll;
                        std::vector<Scalar> post_w(M, 0.0);
                        for (size_t m = 0; m < M; ++m) {
                            if (std::isfinite(terms[m])) post_w[m] = std::exp(terms[m] - pat_ll);
                        }

                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int32_t vid = node.id;
                                Scalar d_lnL_dt = 0.0;
                                for (size_t m = 0; m < M; ++m) {
                                    if (post_w[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                        d_lnL_dt += post_w[m] * (io_m[m].V[vid].dot(dP_mix[m][vid] * io_m[m].D[vid]) / io_m[m].likelihood);
                                    }
                                }
                                local_grad[vid] += pattern.weight * d_lnL_dt;
                            }
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

    std::vector<Scalar> refine_two_set_branch_lengths(
        const std::vector<Scalar>& init_bl,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        int max_iters = 10,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0
    ) const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Scalar> bl = init_bl;
        std::vector<Scalar> log_bl(num_nodes, 0.0);
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id)) {
                bl[i] = std::clamp(bl[i], 1e-8, 10.0);
                log_bl[i] = std::log(bl[i]);
            }
        }

        const size_t m_history = 8;
        std::vector<std::vector<Scalar>> s_hist, y_hist;
        std::vector<Scalar> rho_hist;

        auto eval_log_space = [&](const std::vector<Scalar>& y_vec) -> std::tuple<Scalar, std::vector<Scalar>, std::vector<Scalar>> {
            std::vector<Scalar> cur_bl(num_nodes, 0.0);
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    cur_bl[i] = std::exp(std::clamp(y_vec[i], -18.420680744, 2.302585093));
                }
            }
            auto [ll, grad_t] = compute_two_set_branch_gradients(
                cur_bl, omegas_test, weights_test, omegas_bg, weights_bg, syn_rates, syn_weights, delta, psi, eps_sink_test, eps_sink_bg
            );
            std::vector<Scalar> grad_y(num_nodes, 0.0);
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    grad_y[i] = cur_bl[i] * grad_t[i];
                }
            }
            return {ll, grad_y, cur_bl};
        };

        auto [curr_ll, curr_grad, curr_bl] = eval_log_space(log_bl);
        bl = curr_bl;
        Scalar best_ll = curr_ll;
        std::vector<Scalar> best_bl = bl;

        for (int iter = 0; iter < max_iters; ++iter) {
            Scalar max_g = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_g = std::max(max_g, std::abs(curr_grad[i]));
                }
            }
            if (max_g < 1e-4) break;

            std::vector<Scalar> q = curr_grad;
            std::vector<Scalar> alpha(s_hist.size());
            for (int j = static_cast<int>(s_hist.size()) - 1; j >= 0; --j) {
                Scalar s_dot_q = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) s_dot_q += s_hist[j][i] * q[i];
                alpha[j] = rho_hist[j] * s_dot_q;
                for (size_t i = 0; i < num_nodes; ++i) q[i] -= alpha[j] * y_hist[j][i];
            }

            Scalar gamma = 0.25;
            if (!s_hist.empty()) {
                Scalar s_dot_y = 0.0, y_dot_y = 0.0;
                const auto& s_last = s_hist.back();
                const auto& y_last = y_hist.back();
                for (size_t i = 0; i < num_nodes; ++i) {
                    s_dot_y += s_last[i] * y_last[i];
                    y_dot_y += y_last[i] * y_last[i];
                }
                if (y_dot_y > 1e-12) gamma = std::clamp(s_dot_y / y_dot_y, 1e-3, 5.0);
            }

            std::vector<Scalar> r(num_nodes);
            for (size_t i = 0; i < num_nodes; ++i) r[i] = gamma * q[i];
            for (size_t j = 0; j < s_hist.size(); ++j) {
                Scalar y_dot_r = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) y_dot_r += y_hist[j][i] * r[i];
                Scalar beta = rho_hist[j] * y_dot_r;
                for (size_t i = 0; i < num_nodes; ++i) r[i] += s_hist[j][i] * (alpha[j] - beta);
            }

            std::vector<Scalar> d = r;
            d[tree.root_id] = 0.0;
            Scalar dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) dir_dot_grad += d[i] * curr_grad[i];
            }
            if (dir_dot_grad <= 0.0) {
                d = curr_grad;
                d[tree.root_id] = 0.0;
                dir_dot_grad = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) dir_dot_grad += d[i] * curr_grad[i];
                }
                s_hist.clear(); y_hist.clear(); rho_hist.clear();
            }

            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    d[i] = std::clamp(d[i], -1.2, 1.2);
                }
            }
            dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) dir_dot_grad += d[i] * curr_grad[i];
            }

            Scalar step_size = 1.0;
            std::vector<Scalar> new_log_bl(num_nodes, 0.0);
            std::vector<Scalar> new_bl(num_nodes, 0.0);
            Scalar new_ll = curr_ll;
            std::vector<Scalar> new_grad;
            bool line_search_ok = false;
            for (int ls = 0; ls < 10; ++ls) {
                for (size_t i = 0; i < num_nodes; ++i) {
                    new_log_bl[i] = (i != static_cast<size_t>(tree.root_id))
                        ? std::clamp(log_bl[i] + step_size * d[i], -18.420680744, 2.302585093)
                        : 0.0;
                }
                auto [test_ll, test_grad, test_bl] = eval_log_space(new_log_bl);
                if (test_ll > curr_ll + 1e-4 * step_size * dir_dot_grad || (test_ll > curr_ll + 1e-5 && ls >= 3)) {
                    new_ll = test_ll;
                    new_grad = test_grad;
                    new_bl = test_bl;
                    line_search_ok = true;
                    break;
                }
                step_size *= 0.5;
            }
            if (!line_search_ok) break;

            std::vector<Scalar> s_k(num_nodes), y_k(num_nodes);
            Scalar s_dot_y = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                s_k[i] = new_log_bl[i] - log_bl[i];
                y_k[i] = curr_grad[i] - new_grad[i];
                s_dot_y += s_k[i] * y_k[i];
            }
            if (s_dot_y > 1e-8) {
                if (s_hist.size() >= m_history) {
                    s_hist.erase(s_hist.begin()); y_hist.erase(y_hist.begin()); rho_hist.erase(rho_hist.begin());
                }
                s_hist.push_back(s_k); y_hist.push_back(y_k); rho_hist.push_back(1.0 / s_dot_y);
            }

            log_bl = new_log_bl;
            bl = new_bl; curr_ll = new_ll; curr_grad = new_grad;
            if (curr_ll > best_ll) {
                best_ll = curr_ll; best_bl = bl;
            }
        }
        return best_bl;
    }

    // Fits two-set BUSTED / BUSTED-PH model with independent K-class omega distributions on Test and Background branches
    BUSTEDFit fit_two_set_model(
        bool constrain_test,
        bool constrain_bg,
        size_t K = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3,
        const std::string& multiple_hits = "None",
        ProgressBar* pbar = nullptr,
        const std::string& model_desc = "",
        bool use_error_sink = false,
        Scalar max_eps_sink = 0.01,
        const BUSTEDRateDistribution* init_test_dist = nullptr,
        const BUSTEDRateDistribution* init_bg_dist = nullptr
    ) const {
        BUSTEDFit fit;
        fit.num_rate_classes = K;
        fit.has_background_distribution = true;
        Scalar s = 1.0;

        bool allow_double = (multiple_hits == "Double" || multiple_hits == "Double+Triple");
        bool allow_triple = (multiple_hits == "Double+Triple");
        Scalar delta = 0.0;
        Scalar psi = 0.0;
        Scalar eps_sink_t = use_error_sink ? std::min(0.002, max_eps_sink * 0.5) : 0.0;
        Scalar eps_sink_b = use_error_sink ? std::min(0.002, max_eps_sink * 0.5) : 0.0;

        auto init_dist = [&](bool is_con, const BUSTEDRateDistribution* warm, std::vector<Scalar>& omegas, std::vector<Scalar>& weights) {
            if (warm && warm->omegas.size() == K && warm->weights.size() == K) {
                omegas = warm->omegas;
                weights = warm->weights;
                if (is_con) {
                    omegas.back() = (K == 1) ? std::min(1.0, omegas.back()) : 1.0;
                    for (size_t k = 0; k + 1 < K; ++k) omegas[k] = std::min(1.0, omegas[k]);
                } else if (omegas.back() <= 1.001) {
                    omegas.back() = 2.5;
                }
                return;
            }
            omegas.resize(K);
            weights.resize(K);
            if (K == 1) {
                weights[0] = 1.0;
                omegas[0] = is_con ? std::min(0.5, mg94_omega) : mg94_omega;
            } else if (K == 2) {
                weights = {0.75, 0.25};
                omegas = is_con ? std::vector<Scalar>{0.10, 1.0} : std::vector<Scalar>{0.10, 3.0};
            } else {
                if (is_con) {
                    omegas = {0.05, 0.50, 1.0};
                    weights = {0.45, 0.45, 0.10};
                } else {
                    omegas = {0.05, 0.50, 4.0};
                    weights = {0.50, 0.40, 0.10};
                }
            }
        };

        std::vector<Scalar> omegas_t, weights_t, omegas_b, weights_b;
        init_dist(constrain_test, init_test_dist, omegas_t, weights_t);
        init_dist(constrain_bg, init_bg_dist, omegas_b, weights_b);

        size_t M = srv ? num_syn_rates : 1;
        std::vector<Scalar> syn_rates(M, 1.0);
        std::vector<Scalar> syn_weights(M, 1.0 / M);
        if (srv && M == 2) {
            syn_rates = {0.35, 1.65};
            syn_weights = {0.50, 0.50};
        } else if (srv && M == 3) {
            syn_rates = {0.20, 0.60, 2.20};
            syn_weights = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        }

        if (allow_double) {
            Scalar best_grid_ll = -1e20;
            std::vector<Scalar> d_cands = {0.0, 0.02, 0.08, 0.20};
            std::vector<Scalar> p_cands = allow_triple ? std::vector<Scalar>{0.0, 0.02, 0.05} : std::vector<Scalar>{0.0};
            for (Scalar d_c : d_cands) {
                for (Scalar p_c : p_cands) {
                    Scalar cand_ll = evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, d_c, p_c, eps_sink_t, eps_sink_b);
                    if (cand_ll > best_grid_ll) {
                        best_grid_ll = cand_ll;
                        delta = d_c;
                        psi = p_c;
                    }
                }
            }
        }

        int max_ecm = (K == 1 && !srv && !allow_double && !use_error_sink) ? 1 : 6;
        for (int ecm_iter = 0; ecm_iter < max_ecm; ++ecm_iter) {
            for (int em = 0; em < 2; ++em) {
                auto [wt_t, wt_b] = compute_expected_two_set_weights(
                    s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, delta, psi,
                    eps_sink_t, eps_sink_b,
                    use_error_sink ? &eps_sink_t : nullptr,
                    use_error_sink ? &eps_sink_b : nullptr,
                    max_eps_sink
                );
                weights_t = wt_t;
                weights_b = wt_b;
            }

            if (srv && M >= 2) {
                for (int em = 0; em < 2; ++em) {
                    syn_weights = compute_expected_two_set_syn_weights(
                        s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, delta, psi, eps_sink_t, eps_sink_b
                    );
                }
            }

            Scalar max_omega = use_error_sink ? 50.0 : 2000.0;

            // Optimize Test omegas + scale s
            if (constrain_test) {
                if (K == 1) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_two_set_log_l(s_c, {w1}, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, omegas_t[0], s, 1e-4, 1.0, 1e-4, 35);
                    omegas_t[0] = opt.a; s = opt.b;
                } else if (K == 2) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_two_set_log_l(s_c, {w1, 1.0}, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, omegas_t[0], s, 1e-4, 1.0, 1e-4, 35);
                    omegas_t[0] = opt.a; omegas_t[1] = 1.0; s = opt.b;
                } else if (K == 3) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        if (x[1] > x[2]) return 1e20;
                        return -evaluate_two_set_log_l(x[0], {x[1], x[2], 1.0}, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, {s, omegas_t[0], omegas_t[1]}, {0.01, 1e-4, 1e-4}, {50.0, 1.0, 1.0}, 1e-4, 35);
                    s = opt.x[0]; omegas_t[0] = opt.x[1]; omegas_t[1] = opt.x[2]; omegas_t[2] = 1.0;
                }
            } else {
                if (K == 1) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_two_set_log_l(s_c, {w1}, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, omegas_t[0], s, 1e-4, max_omega, 1e-4, 35);
                    omegas_t[0] = opt.a; s = opt.b;
                } else if (K == 2) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        return -evaluate_two_set_log_l(x[0], {x[1], x[2]}, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, {s, omegas_t[0], omegas_t[1]}, {0.01, 1e-4, 1.0}, {50.0, 1.0, max_omega}, 1e-4, 45);
                    s = opt.x[0]; omegas_t[0] = opt.x[1]; omegas_t[1] = opt.x[2];
                } else if (K == 3) {
                    auto opt = NelderMeadND<4>::minimize([&](const std::array<Scalar, 4>& x) {
                        if (x[1] > x[2]) return 1e20;
                        return -evaluate_two_set_log_l(x[0], {x[1], x[2], x[3]}, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, {s, omegas_t[0], omegas_t[1], omegas_t[2]}, {0.01, 1e-4, 1e-4, 1.0}, {50.0, 1.0, 1.0, max_omega}, 1e-4, 45);
                    s = opt.x[0]; omegas_t[0] = opt.x[1]; omegas_t[1] = opt.x[2]; omegas_t[2] = opt.x[3];
                }
            }

            // Optimize Background omegas
            if (constrain_bg) {
                if (K == 1) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_two_set_log_l(s_c, omegas_t, weights_t, {w1}, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, omegas_b[0], s, 1e-4, 1.0, 1e-4, 35);
                    omegas_b[0] = opt.a; s = opt.b;
                } else if (K == 2) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_two_set_log_l(s_c, omegas_t, weights_t, {w1, 1.0}, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, omegas_b[0], s, 1e-4, 1.0, 1e-4, 35);
                    omegas_b[0] = opt.a; omegas_b[1] = 1.0; s = opt.b;
                } else if (K == 3) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        if (x[1] > x[2]) return 1e20;
                        return -evaluate_two_set_log_l(x[0], omegas_t, weights_t, {x[1], x[2], 1.0}, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, {s, omegas_b[0], omegas_b[1]}, {0.01, 1e-4, 1e-4}, {50.0, 1.0, 1.0}, 1e-4, 35);
                    s = opt.x[0]; omegas_b[0] = opt.x[1]; omegas_b[1] = opt.x[2]; omegas_b[2] = 1.0;
                }
            } else {
                if (K == 1) {
                    auto opt = NelderMead2D::minimize([&](Scalar w1, Scalar s_c) {
                        return -evaluate_two_set_log_l(s_c, omegas_t, weights_t, {w1}, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, omegas_b[0], s, 1e-4, max_omega, 1e-4, 35);
                    omegas_b[0] = opt.a; s = opt.b;
                } else if (K == 2) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        return -evaluate_two_set_log_l(x[0], omegas_t, weights_t, {x[1], x[2]}, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, {s, omegas_b[0], omegas_b[1]}, {0.01, 1e-4, 1.0}, {50.0, 1.0, max_omega}, 1e-4, 45);
                    s = opt.x[0]; omegas_b[0] = opt.x[1]; omegas_b[1] = opt.x[2];
                } else if (K == 3) {
                    auto opt = NelderMeadND<4>::minimize([&](const std::array<Scalar, 4>& x) {
                        if (x[1] > x[2]) return 1e20;
                        return -evaluate_two_set_log_l(x[0], omegas_t, weights_t, {x[1], x[2], x[3]}, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                    }, {s, omegas_b[0], omegas_b[1], omegas_b[2]}, {0.01, 1e-4, 1e-4, 1.0}, {50.0, 1.0, 1.0, max_omega}, 1e-4, 45);
                    s = opt.x[0]; omegas_b[0] = opt.x[1]; omegas_b[1] = opt.x[2]; omegas_b[2] = opt.x[3];
                }
            }

            if (srv && M == 2) {
                auto syn_opt = NelderMead2D::minimize([&](Scalar a1, Scalar /*unused*/) {
                    Scalar q1 = syn_weights[0], q2 = syn_weights[1];
                    Scalar a2 = (1.0 - q1 * a1) / q2;
                    if (a1 < 0.01 || a2 < a1 || a2 > 50.0) return 1e20;
                    return -evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, {a1, a2}, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                }, syn_rates[0], 0.0, 0.01, 2.0, 1e-3, 25);
                Scalar a1 = syn_opt.a;
                Scalar a2 = (1.0 - syn_weights[0] * a1) / syn_weights[1];
                if (a1 >= 0.01 && a2 >= a1) {
                    syn_rates[0] = a1; syn_rates[1] = a2;
                }
            } else if (srv && M == 3) {
                auto syn_opt = NelderMead2D::minimize([&](Scalar a1, Scalar a2) {
                    Scalar q1 = syn_weights[0], q2 = syn_weights[1], q3 = syn_weights[2];
                    Scalar a3 = (1.0 - q1 * a1 - q2 * a2) / q3;
                    if (a1 < 0.01 || a2 < a1 || a3 < a2 || a3 > 50.0) return 1e20;
                    return -evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, {a1, a2, a3}, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                }, syn_rates[0], syn_rates[1], 0.01, 5.0, 1e-3, 30);
                Scalar a1 = syn_opt.a, a2 = syn_opt.b;
                Scalar a3 = (1.0 - syn_weights[0] * a1 - syn_weights[1] * a2) / syn_weights[2];
                if (a1 >= 0.01 && a2 >= a1 && a3 >= a2) {
                    syn_rates[0] = a1; syn_rates[1] = a2; syn_rates[2] = a3;
                }
            }

            if (allow_triple) {
                auto mh_opt = NelderMead2D::minimize([&](Scalar d, Scalar p) {
                    if (d < 0.0 || d > 10.0 || p < 0.0 || p > 10.0) return 1e20;
                    return -evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, d, p, eps_sink_t, eps_sink_b);
                }, delta, psi, 0.0, 10.0, 1e-4, 25);
                delta = std::max(0.0, mh_opt.a); psi = std::max(0.0, mh_opt.b);
            } else if (allow_double) {
                auto d_opt = NelderMead2D::minimize([&](Scalar d, Scalar /*unused*/) {
                    if (d < 0.0 || d > 10.0) return 1e20;
                    return -evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, d, 0.0, eps_sink_t, eps_sink_b);
                }, delta, 0.0, 0.0, 10.0, 1e-4, 25);
                delta = std::max(0.0, d_opt.a);
            }

            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << " (ECM " << (ecm_iter + 1) << "/" << max_ecm << "): lnL = "
                   << std::fixed << std::setprecision(2)
                   << evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b);
                pbar->set_status(ss.str());
                pbar->tick();
            }
        }

        std::vector<Scalar> bl(tree.num_nodes(), 0.0);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                bl[node.id] = s * branch_tau[node.id];
            }
        }

        if (refine_branches) {
            bl = refine_two_set_branch_lengths(bl, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, delta, psi, 1, eps_sink_t, eps_sink_b);
            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << ": Branch lengths refined";
                pbar->set_status(ss.str());
                pbar->tick();
            }
            if (K >= 2 || use_error_sink) {
                auto [wt_t, wt_b] = compute_expected_two_set_weights(
                    1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, delta, psi,
                    eps_sink_t, eps_sink_b,
                    use_error_sink ? &eps_sink_t : nullptr,
                    use_error_sink ? &eps_sink_b : nullptr,
                    max_eps_sink
                );
                weights_t = wt_t;
                weights_b = wt_b;
            }
            if (srv && M >= 2) {
                syn_weights = compute_expected_two_set_syn_weights(1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, delta, psi, eps_sink_t, eps_sink_b);
            }
            if (allow_triple) {
                auto mh_polish = NelderMead2D::minimize([&](Scalar d, Scalar p) {
                    if (d < 0.0 || d > 10.0 || p < 0.0 || p > 10.0) return 1e20;
                    return -evaluate_two_set_log_l(1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, nullptr, d, p, eps_sink_t, eps_sink_b);
                }, delta, psi, 0.0, 10.0, 1e-4, 25);
                delta = std::max(0.0, mh_polish.a); psi = std::max(0.0, mh_polish.b);
            } else if (allow_double) {
                auto d_polish = NelderMead2D::minimize([&](Scalar d, Scalar /*unused*/) {
                    if (d < 0.0 || d > 10.0) return 1e20;
                    return -evaluate_two_set_log_l(1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, nullptr, d, 0.0, eps_sink_t, eps_sink_b);
                }, delta, 0.0, 0.0, 10.0, 1e-4, 25);
                delta = std::max(0.0, d_polish.a);
            }
            if (allow_double && pbar) {
                pbar->tick();
            }
        }

        Scalar best_lnL = evaluate_two_set_log_l(
            1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, &fit.site_log_likelihoods, delta, psi, eps_sink_t, eps_sink_b
        );
        fit.log_likelihood = best_lnL;
        fit.tree_scale = s;
        fit.branch_lengths = bl;
        fit.delta = delta;
        fit.psi = psi;
        fit.error_sink_weight = eps_sink_t;

        fit.test_distribution.omegas = omegas_t;
        fit.test_distribution.weights = weights_t;
        fit.test_distribution.error_sink_weight = eps_sink_t;
        fit.background_distribution.omegas = omegas_b;
        fit.background_distribution.weights = weights_b;
        fit.background_distribution.error_sink_weight = eps_sink_b;
        if (srv) {
            fit.test_distribution.syn_rates = syn_rates;
            fit.test_distribution.syn_weights = syn_weights;
            fit.background_distribution.syn_rates = syn_rates;
            fit.background_distribution.syn_weights = syn_weights;
        }

        size_t N_codons = aln.num_codons;
        size_t num_branches = tree.num_nodes() - 1;
        size_t num_mixture_params = (2 * K - 1 - (constrain_test ? 1 : 0)) + (2 * K - 1 - (constrain_bg ? 1 : 0)) + (use_error_sink ? 2 : 0);
        size_t num_syn_params = srv ? 2 * (M - 1) : 0;
        size_t num_mh_params = (allow_double ? 1 : 0) + (allow_triple ? 1 : 0);
        size_t K_params = 5 + (refine_branches ? num_branches : 1) + num_mixture_params + num_syn_params + num_mh_params;
        fit.aicc = -2.0 * fit.log_likelihood + 2.0 * K_params * N_codons / std::max(1.0, static_cast<double>(N_codons) - static_cast<double>(K_params) - 1.0);

        return fit;
    }

    // Fits coupled-k RELAX two-set model:
    //   Reference (Background, partition 1): K-class distribution (omegas_b[c], weights[c])
    //   Test (Foreground, partition 0):      K-class distribution (omegas_t[c] = (omegas_b[c])^k, weights[c])
    // Warm-started from shared_null_fit (the exact k = 1 RELAX Null model), guaranteeing lnL_RELAX >= lnL_shared.
    BUSTEDFit fit_relax_coupled_model(
        const BUSTEDFit& shared_null_fit,
        size_t K = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3,
        const std::string& multiple_hits = "None",
        ProgressBar* pbar = nullptr,
        const std::string& model_desc = "",
        bool use_error_sink = false,
        Scalar max_eps_sink = 0.01,
        Scalar* out_relax_k = nullptr
    ) const {
        BUSTEDFit fit = shared_null_fit;
        fit.num_rate_classes = K;
        fit.has_background_distribution = true;

        bool allow_double = (multiple_hits == "Double" || multiple_hits == "Double+Triple");
        bool allow_triple = (multiple_hits == "Double+Triple");
        Scalar max_omega = use_error_sink ? 50.0 : 2000.0;

        std::vector<Scalar> omegas_b = shared_null_fit.test_distribution.omegas;
        std::vector<Scalar> weights = shared_null_fit.test_distribution.weights;
        if (omegas_b.size() != K || weights.size() != K) {
            omegas_b.resize(K, 1.0);
            weights.assign(K, 1.0 / static_cast<Scalar>(K));
        }

        size_t M = srv ? num_syn_rates : 1;
        std::vector<Scalar> syn_rates = (srv && shared_null_fit.test_distribution.syn_rates.size() == M)
            ? shared_null_fit.test_distribution.syn_rates : std::vector<Scalar>(M, 1.0);
        std::vector<Scalar> syn_weights = (srv && shared_null_fit.test_distribution.syn_weights.size() == M)
            ? shared_null_fit.test_distribution.syn_weights : std::vector<Scalar>(M, 1.0 / static_cast<Scalar>(M));

        Scalar delta = shared_null_fit.delta;
        Scalar psi = shared_null_fit.psi;
        Scalar eps_sink_t = shared_null_fit.error_sink_weight;
        Scalar eps_sink_b = shared_null_fit.error_sink_weight;

        std::vector<Scalar> base_bl = shared_null_fit.branch_lengths;
        if (base_bl.size() != tree.num_nodes()) {
            base_bl.assign(tree.num_nodes(), 0.0);
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    base_bl[node.id] = shared_null_fit.tree_scale * branch_tau[node.id];
                }
            }
        }
        std::vector<Scalar> cur_bl = base_bl;
        auto get_scaled_bl = [&](Scalar scale_s) -> const std::vector<Scalar>* {
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    cur_bl[node.id] = scale_s * base_bl[node.id];
                }
            }
            return &cur_bl;
        };

        auto couple_omegas = [&](const std::vector<Scalar>& wb, Scalar kv) {
            std::vector<Scalar> wt(wb.size());
            for (size_t c = 0; c < wb.size(); ++c) {
                if (std::abs(wb[c] - 1.0) < 1e-12) {
                    wt[c] = 1.0;
                } else if (wb[c] <= 1e-12) {
                    wt[c] = (kv > 0.0) ? 1e-8 : 1.0;
                } else {
                    wt[c] = std::clamp(std::pow(wb[c], kv), 1e-8, max_omega);
                }
            }
            return wt;
        };

        Scalar s = 1.0;
        Scalar k_relax = 1.0;
        std::vector<Scalar> omegas_t = couple_omegas(omegas_b, k_relax);

        Scalar best_tracked_ll = evaluate_two_set_log_l(
            1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b
        );

        // Initial 1D search on k_relax holding shared distribution fixed
        {
            std::vector<Scalar> k_grid = {0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0};
            Scalar best_k = k_relax;
            for (Scalar kc : k_grid) {
                auto wt_c = couple_omegas(omegas_b, kc);
                Scalar ll_c = evaluate_two_set_log_l(
                    1.0, wt_c, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b
                );
                if (ll_c > best_tracked_ll) {
                    best_tracked_ll = ll_c;
                    best_k = kc;
                }
            }
            auto [brent_k, brent_nll] = Brent1D::minimize([&](Scalar kc) {
                Scalar kc_clamped = std::clamp(kc, 0.0, 50.0);
                auto wt_c = couple_omegas(omegas_b, kc_clamped);
                return -evaluate_two_set_log_l(
                    1.0, wt_c, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b
                );
            }, std::max(0.0, best_k * 0.25), best_k, std::min(50.0, std::max(2.0, best_k * 2.5)), 1e-4, 35);
            if (-brent_nll >= best_tracked_ll) {
                best_tracked_ll = -brent_nll;
                k_relax = std::clamp(brent_k, 0.0, 50.0);
            } else {
                k_relax = best_k;
            }
            omegas_t = couple_omegas(omegas_b, k_relax);
        }

        size_t B_t = std::max<size_t>(1, num_test_branches);
        size_t B_b = std::max<size_t>(1, num_background_branches);
        Scalar frac_t = static_cast<Scalar>(B_t) / static_cast<Scalar>(B_t + B_b);
        Scalar frac_b = static_cast<Scalar>(B_b) / static_cast<Scalar>(B_t + B_b);

        int max_ecm = (K == 1 && !srv && !allow_double && !use_error_sink) ? 1 : 3;
        for (int ecm_iter = 0; ecm_iter < max_ecm; ++ecm_iter) {
            // 1. Shared weights E-step
            if (K >= 2 || use_error_sink) {
                Scalar cand_eps_t = eps_sink_t, cand_eps_b = eps_sink_b;
                auto [wt_t, wt_b] = compute_expected_two_set_weights(
                    1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), delta, psi,
                    eps_sink_t, eps_sink_b,
                    use_error_sink ? &cand_eps_t : nullptr,
                    use_error_sink ? &cand_eps_b : nullptr,
                    max_eps_sink
                );
                std::vector<Scalar> w_cand(K);
                Scalar sum_w = 0.0;
                for (size_t c = 0; c < K; ++c) {
                    w_cand[c] = std::clamp(frac_t * wt_t[c] + frac_b * wt_b[c], 1e-6, 1.0 - 1e-6);
                    sum_w += w_cand[c];
                }
                for (size_t c = 0; c < K; ++c) w_cand[c] /= sum_w;
                Scalar cand_eps = frac_t * cand_eps_t + frac_b * cand_eps_b;
                Scalar cand_ll = evaluate_two_set_log_l(
                    1.0, omegas_t, w_cand, omegas_b, w_cand, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, cand_eps, cand_eps
                );
                if (cand_ll >= best_tracked_ll - 1e-9) {
                    weights = w_cand;
                    eps_sink_t = cand_eps;
                    eps_sink_b = cand_eps;
                    best_tracked_ll = std::max(best_tracked_ll, cand_ll);
                }
            }

            // 2. SRV weights E-step
            if (srv && M >= 2) {
                auto q_cand = compute_expected_two_set_syn_weights(
                    1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), delta, psi, eps_sink_t, eps_sink_b
                );
                Scalar cand_ll = evaluate_two_set_log_l(
                    1.0, omegas_t, weights, omegas_b, weights, syn_rates, q_cand, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b
                );
                if (cand_ll >= best_tracked_ll - 1e-9) {
                    syn_weights = q_cand;
                    best_tracked_ll = std::max(best_tracked_ll, cand_ll);
                }
            }

            // 3. Joint M-step on (s, omegas_b, k_relax)
            if (K == 1) {
                auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                    std::vector<Scalar> wb = {x[1]};
                    auto wt = couple_omegas(wb, x[2]);
                    return -evaluate_two_set_log_l(1.0, wt, weights, wb, weights, syn_rates, syn_weights, get_scaled_bl(x[0]), nullptr, delta, psi, eps_sink_t, eps_sink_b);
                }, {s, omegas_b[0], k_relax}, {0.05, 1e-4, 0.0}, {20.0, max_omega, 50.0}, 1e-4, 40);
                if (-opt.f >= best_tracked_ll) {
                    s = opt.x[0]; omegas_b[0] = opt.x[1]; k_relax = opt.x[2];
                    best_tracked_ll = -opt.f;
                }
            } else if (K == 2) {
                std::vector<std::array<Scalar, 4>> extra_grid;
                if (ecm_iter == 0) {
                    extra_grid = {
                        {s, omegas_b[0], std::max(2.0, omegas_b[1]), 0.2},
                        {s, omegas_b[0], std::max(2.0, omegas_b[1]), 1.0},
                        {s, omegas_b[0], std::max(2.0, omegas_b[1]), 2.0}
                    };
                }
                auto opt = NelderMeadND<4>::minimize([&](const std::array<Scalar, 4>& x) {
                    std::vector<Scalar> wb = {x[1], x[2]};
                    auto wt = couple_omegas(wb, x[3]);
                    return -evaluate_two_set_log_l(1.0, wt, weights, wb, weights, syn_rates, syn_weights, get_scaled_bl(x[0]), nullptr, delta, psi, eps_sink_t, eps_sink_b);
                }, {s, omegas_b[0], omegas_b[1], k_relax}, {0.05, 1e-4, 1.0, 0.0}, {20.0, 1.0, max_omega, 50.0}, 1e-4, 45, extra_grid);
                if (-opt.f >= best_tracked_ll) {
                    s = opt.x[0]; omegas_b[0] = opt.x[1]; omegas_b[1] = opt.x[2]; k_relax = opt.x[3];
                    best_tracked_ll = -opt.f;
                }
            } else if (K == 3) {
                std::vector<std::array<Scalar, 5>> extra_grid;
                if (ecm_iter == 0) {
                    extra_grid = {
                        {s, omegas_b[0], omegas_b[1], std::max(2.0, omegas_b[2]), 0.2},
                        {s, omegas_b[0], omegas_b[1], std::max(2.0, omegas_b[2]), 1.0},
                        {s, omegas_b[0], omegas_b[1], std::max(2.0, omegas_b[2]), 2.0}
                    };
                }
                auto opt = NelderMeadND<5>::minimize([&](const std::array<Scalar, 5>& x) {
                    if (x[1] > x[2]) return 1e20;
                    std::vector<Scalar> wb = {x[1], x[2], x[3]};
                    auto wt = couple_omegas(wb, x[4]);
                    return -evaluate_two_set_log_l(1.0, wt, weights, wb, weights, syn_rates, syn_weights, get_scaled_bl(x[0]), nullptr, delta, psi, eps_sink_t, eps_sink_b);
                }, {s, omegas_b[0], omegas_b[1], omegas_b[2], k_relax}, {0.05, 1e-4, 1e-4, 1.0, 0.0}, {20.0, 1.0, 1.0, max_omega, 50.0}, 1e-4, 50, extra_grid);
                if (-opt.f >= best_tracked_ll) {
                    s = opt.x[0]; omegas_b[0] = opt.x[1]; omegas_b[1] = opt.x[2]; omegas_b[2] = opt.x[3]; k_relax = opt.x[4];
                    best_tracked_ll = -opt.f;
                }
            }
            omegas_t = couple_omegas(omegas_b, k_relax);

            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << " (ECM " << (ecm_iter + 1) << "/" << max_ecm << ", k="
                   << std::fixed << std::setprecision(2) << k_relax << "): lnL = " << best_tracked_ll;
                pbar->set_status(ss.str());
                pbar->tick();
            }
        }

        std::vector<Scalar> bl = *get_scaled_bl(s);

        if (refine_branches) {
            auto cand_bl = refine_two_set_branch_lengths(
                bl, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, delta, psi, 1, eps_sink_t, eps_sink_b
            );
            Scalar cand_bl_ll = evaluate_two_set_log_l(
                1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, &cand_bl, nullptr, delta, psi, eps_sink_t, eps_sink_b
            );
            if (cand_bl_ll >= best_tracked_ll) {
                bl = cand_bl;
                best_tracked_ll = cand_bl_ll;
            }
            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << ": Branch lengths refined";
                pbar->set_status(ss.str());
                pbar->tick();
            }

            // Final 1D polish on k_relax with refined branch lengths
            auto [brent_k, brent_nll] = Brent1D::minimize([&](Scalar kc) {
                Scalar kc_clamped = std::clamp(kc, 0.0, 50.0);
                auto wt_c = couple_omegas(omegas_b, kc_clamped);
                return -evaluate_two_set_log_l(
                    1.0, wt_c, weights, omegas_b, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink_t, eps_sink_b
                );
            }, std::max(0.0, k_relax * 0.5), k_relax, std::min(50.0, std::max(1.5, k_relax * 1.5)), 1e-4, 30);
            if (-brent_nll >= best_tracked_ll) {
                best_tracked_ll = -brent_nll;
                k_relax = std::clamp(brent_k, 0.0, 50.0);
                omegas_t = couple_omegas(omegas_b, k_relax);
            }
            if (allow_double && pbar) {
                pbar->tick();
            }
        }

        Scalar best_lnL = evaluate_two_set_log_l(
            1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, &bl, &fit.site_log_likelihoods, delta, psi, eps_sink_t, eps_sink_b
        );
        if (best_lnL < shared_null_fit.log_likelihood) {
            fit = shared_null_fit;
            k_relax = 1.0;
        } else {
            fit.log_likelihood = best_lnL;
            fit.tree_scale = s * shared_null_fit.tree_scale;
            fit.branch_lengths = bl;
            fit.delta = delta;
            fit.psi = psi;
            fit.error_sink_weight = eps_sink_t;

            fit.test_distribution.omegas = omegas_t;
            fit.test_distribution.weights = weights;
            fit.test_distribution.error_sink_weight = eps_sink_t;
            fit.background_distribution.omegas = omegas_b;
            fit.background_distribution.weights = weights;
            fit.background_distribution.error_sink_weight = eps_sink_b;
            if (srv) {
                fit.test_distribution.syn_rates = syn_rates;
                fit.test_distribution.syn_weights = syn_weights;
                fit.background_distribution.syn_rates = syn_rates;
                fit.background_distribution.syn_weights = syn_weights;
            }
        }

        if (out_relax_k) *out_relax_k = k_relax;

        size_t N_codons = aln.num_codons;
        size_t num_branches = tree.num_nodes() - 1;
        size_t num_mixture_params = (2 * K - 1) + 1 + (use_error_sink ? 1 : 0); // shared weights + BG omegas + k_relax
        size_t num_syn_params = srv ? 2 * (M - 1) : 0;
        size_t num_mh_params = (allow_double ? 1 : 0) + (allow_triple ? 1 : 0);
        size_t K_params = 5 + (refine_branches ? num_branches : 1) + num_mixture_params + num_syn_params + num_mh_params;
        fit.aicc = -2.0 * fit.log_likelihood + 2.0 * K_params * N_codons / std::max(1.0, static_cast<double>(N_codons) - static_cast<double>(K_params) - 1.0);

        return fit;
    }

    // Run complete BUSTED or BUSTED-PH analysis with options
    BUSTEDResult run(
        BUSTEDSettings settings = {},
        bool show_progress = false,
        bool force_progress = false
    ) const {
        auto t0 = std::chrono::high_resolution_clock::now();
        select_branch_partitions(settings.test_branches, settings.comparison_branches, settings.phenotype_mode || settings.relax_mode);
        if ((settings.phenotype_mode || settings.relax_mode) && num_background_branches == 0) {
            throw std::runtime_error(
                "BUSTED-PH / RELAX mode requires a proper subset of foreground/test branches and a non-empty background/reference set. "
                "Specify foreground branches via --branches <spec> or label branches with {TEST}/{Foreground} in the tree."
            );
        }

        BUSTEDResult res;
        if (settings.full_auto) {
            res.auto_diagnostics = diagnose_and_select_model(settings);
        }
        res.settings = settings;

        bool allow_double = (settings.multiple_hits == "Double" || settings.multiple_hits == "Double+Triple");
        auto count_steps = [&](size_t K) -> size_t {
            size_t max_ecm = (K == 1 && !settings.srv && !allow_double && !settings.error_sink) ? 1 : 6;
            return max_ecm + (settings.refine_branch_lengths ? 1 : 0) + ((settings.refine_branch_lengths && allow_double) ? 1 : 0);
        };
        auto count_relax_steps = [&](size_t K) -> size_t {
            size_t max_ecm = (K == 1 && !settings.srv && !allow_double && !settings.error_sink) ? 1 : 3;
            return max_ecm + (settings.refine_branch_lengths ? 1 : 0) + ((settings.refine_branch_lengths && allow_double) ? 1 : 0);
        };

        bool run_relax_in_ph = settings.phenotype_mode && (settings.include_relax_in_ph || settings.relax_mode);
        size_t total_steps = 0;
        if (settings.phenotype_mode) {
            total_steps = count_steps(settings.num_rate_classes) * 4
                        + (run_relax_in_ph ? count_relax_steps(settings.num_rate_classes) : 0);
        } else if (settings.auto_select_k) {
            for (size_t k = 1; k <= settings.max_k; ++k) {
                total_steps += count_steps(k);
            }
            total_steps += count_steps(settings.max_k);
            if (settings.relax_mode) {
                total_steps += count_steps(settings.max_k) + count_relax_steps(settings.max_k);
            }
        } else {
            total_steps = count_steps(settings.num_rate_classes) * 2
                        + (settings.relax_mode ? (count_steps(settings.num_rate_classes) + count_relax_steps(settings.num_rate_classes)) : 0);
        }

        std::unique_ptr<ProgressBar> pbar;
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar = std::make_unique<ProgressBar>(
                total_steps,
                settings.phenotype_mode ? "[BUSTED-PH] Mixture Models" : "[BUSTED] Mixture Models",
                "steps",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        if (settings.phenotype_mode) {
            size_t K = settings.num_rate_classes;
            res.optimal_k = K;
            res.has_background_test = true;
            res.shared_df = static_cast<int>(2 * K - 1);

            // 1. Fit Shared Distributions Model (M_shared: D_T == D_B on Test U Background; also exact RELAX Null k = 1)
            {
                auto saved_is_test = is_test_branch;
                auto saved_num_test = num_test_branches;
                num_test_branches = 0;
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) {
                        bool in_shared = (branch_partition[node.id] <= 1);
                        is_test_branch[node.id] = in_shared;
                        if (in_shared) num_test_branches++;
                    }
                }
                res.shared_distributions = fit_model(
                    false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, pbar.get(), "Shared Dist (RELAX Null k=1) K=" + std::to_string(K),
                    settings.error_sink, settings.error_sink_max_weight
                );
                res.shared_distributions.has_background_distribution = true;
                res.shared_distributions.background_distribution = res.shared_distributions.test_distribution;
                is_test_branch = saved_is_test;
                num_test_branches = saved_num_test;
            }

            // 1b. Fit Coupled-k RELAX Alternative Model (M_RELAX: omega_FG = omega_BG^k, shared weights)
            if (run_relax_in_ph) {
                res.has_relax_test = true;
                res.relax_alternative = fit_relax_coupled_model(
                    res.shared_distributions, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, pbar.get(), "RELAX Coupled-k K=" + std::to_string(K),
                    settings.error_sink, settings.error_sink_max_weight, &res.relax_k
                );
                res.relax_lrt = std::max(0.0, 2.0 * (res.relax_alternative.log_likelihood - res.shared_distributions.log_likelihood));
                res.relax_p_value = busted_chi2_p_value(res.relax_lrt, 1.0);
                if (res.relax_p_value <= settings.p_value_threshold) {
                    res.relax_direction = (res.relax_k > 1.0) ? "Intensification" : "Relaxation";
                } else {
                    res.relax_direction = "None";
                }
            }

            // 2. Fit Constrained Test Model (M_null_T: omega_K^(T) = 1, omega_K^(B) >= 1)
            res.constrained = fit_two_set_model(
                true, false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar.get(), "Constrained Test K=" + std::to_string(K),
                settings.error_sink, settings.error_sink_max_weight,
                &res.shared_distributions.test_distribution, &res.shared_distributions.test_distribution
            );

            // 3. Fit Constrained Background Model (M_null_B: omega_K^(T) >= 1, omega_K^(B) = 1)
            res.constrained_background = fit_two_set_model(
                false, true, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar.get(), "Constrained Background K=" + std::to_string(K),
                settings.error_sink, settings.error_sink_max_weight,
                &res.shared_distributions.test_distribution, &res.shared_distributions.test_distribution
            );

            // 4. Fit Unconstrained Two-Set Model (M_unconstrained: omega_K^(T) >= 1, omega_K^(B) >= 1)
            // Warm-start from the best nested model's Test/Background distributions
            const BUSTEDRateDistribution* warm_t = &res.constrained_background.test_distribution;
            const BUSTEDRateDistribution* warm_b = &res.constrained.background_distribution;
            if (res.has_relax_test && res.relax_alternative.log_likelihood > res.constrained_background.log_likelihood &&
                res.relax_alternative.log_likelihood > res.constrained.log_likelihood) {
                warm_t = &res.relax_alternative.test_distribution;
                warm_b = &res.relax_alternative.background_distribution;
            }
            res.unconstrained = fit_two_set_model(
                false, false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar.get(), "Unconstrained Two-Set K=" + std::to_string(K),
                settings.error_sink, settings.error_sink_max_weight,
                warm_t, warm_b
            );

            // Ensure proper nesting of log-likelihoods (M_shared c M_RELAX c M_unconstrained)
            Scalar max_nested_ll = std::max({
                res.constrained.log_likelihood,
                res.constrained_background.log_likelihood,
                res.shared_distributions.log_likelihood,
                res.has_relax_test ? res.relax_alternative.log_likelihood : res.shared_distributions.log_likelihood
            });
            if (res.unconstrained.log_likelihood < max_nested_ll) {
                res.unconstrained.log_likelihood = max_nested_ll;
            }

            bool unc_failed = (res.unconstrained.log_likelihood <= -1e19 || !std::isfinite(res.unconstrained.log_likelihood));
            bool con_failed = (res.constrained.log_likelihood <= -1e19 || !std::isfinite(res.constrained.log_likelihood));
            bool bg_failed = (res.constrained_background.log_likelihood <= -1e19 || !std::isfinite(res.constrained_background.log_likelihood));
            bool shared_failed = (res.shared_distributions.log_likelihood <= -1e19 || !std::isfinite(res.shared_distributions.log_likelihood));

            if (unc_failed || con_failed || bg_failed || shared_failed) {
                res.is_failed = true;
                res.error_message = "Optimization failed: model log-likelihood collapsed to penalty threshold (-1e20). The fit did not converge.";
                res.lrt = 0.0;
                res.p_value = 1.0;
                res.background_lrt = 0.0;
                res.background_p_value = 1.0;
                res.shared_lrt = 0.0;
                res.shared_p_value = 1.0;
                res.union_intersection_p_value = 1.0;
                res.phenotype_associated = false;
                res.phenotype_verdict = "Analysis failed: likelihood collapsed or did not converge.";
            } else {
                // Test 1: Foreground / Test selection (LRT_T)
                res.lrt = std::max(0.0, 2.0 * (res.unconstrained.log_likelihood - res.constrained.log_likelihood));
                res.p_value = (res.lrt > 0.0) ? (0.5 * std::exp(-res.lrt / 2.0)) : 0.5;

                // Test 2: Background / Comparison selection (LRT_B)
                res.background_lrt = std::max(0.0, 2.0 * (res.unconstrained.log_likelihood - res.constrained_background.log_likelihood));
                res.background_p_value = (res.background_lrt > 0.0) ? (0.5 * std::exp(-res.background_lrt / 2.0)) : 0.5;

                // Test 3: Shared vs Separate Distributions (LRT_shared ~ chi^2_{2K-1})
                res.shared_lrt = std::max(0.0, 2.0 * (res.unconstrained.log_likelihood - res.shared_distributions.log_likelihood));
                res.shared_p_value = busted_chi2_p_value(res.shared_lrt, static_cast<Scalar>(res.shared_df));
            }

            // Canonical BUSTED-PH Testing Procedure:
            // Step 1: Union-Intersection Test (UIT) on (p_FG, p_diff) -> p_UIT = max(p_FG, p_diff) <= alpha
            // Step 2: Lack of evidence for background selection via AIC-c:
            //         Delta AICc = AICc(Constrained BG) - AICc(Unconstrained) <= 2.0
            //         (asymptotically equivalent to LRT_B <= 4.0 <=> p_BG > 0.5 * exp(-2) ≈ 0.068)
            Scalar alpha = settings.p_value_threshold;
            res.union_intersection_p_value = std::max(res.p_value, res.shared_p_value);
            res.background_delta_aicc = res.constrained_background.aicc - res.unconstrained.aicc;
            res.background_neutral_supported = (res.background_delta_aicc <= 2.0) || (res.background_p_value > 0.068);
            res.ranking_score = -std::log10(std::max(1e-8, res.p_value))
                                - std::log10(std::max(1e-8, res.shared_p_value))
                                + std::log10(std::max(1e-8, res.background_p_value));

            if (res.union_intersection_p_value <= alpha) {
                if (res.background_neutral_supported) {
                    res.phenotype_associated = true;
                    res.phenotype_verdict =
                        "The composite null hypothesis (no FG selection or no FG/BG difference) is rejected (p_UIT <= alpha), "
                        "and the neutral model for background branches is sufficiently supported by AIC-c: "
                        "statistical evidence that selection is associated with the trait.";
                } else {
                    res.phenotype_associated = false;
                    res.phenotype_verdict =
                        "The composite null hypothesis (no FG selection or no FG/BG difference) is rejected (p_UIT <= alpha), "
                        "but the neutral model for background branches is not supported by AIC-c: "
                        "selection is acting broadly on the tree, not just on branches with the trait.";
                }
            } else {
                res.phenotype_associated = false;
                res.phenotype_verdict =
                    "The composite null hypothesis for no selection on foreground or no difference between foreground and background "
                    "could not be rejected (p_UIT = max(p_FG, p_diff) > alpha): no statistical evidence that selection is associated with the trait.";
            }
        } else if (settings.auto_select_k) {
            auto [opt_k, null_fits] = select_optimal_k(
                settings.max_k,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes,
                settings.multiple_hits,
                pbar.get(),
                settings.error_sink,
                settings.error_sink_max_weight
            );
            res.optimal_k = opt_k;
            res.k_null_fits = null_fits;
            res.constrained = null_fits[opt_k - 1];
            res.unconstrained = fit_model(
                false,
                opt_k,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes,
                settings.multiple_hits,
                pbar.get(),
                "Unconstrained K=" + std::to_string(opt_k),
                settings.error_sink,
                settings.error_sink_max_weight,
                &res.constrained
            );
        } else {
            res.optimal_k = settings.num_rate_classes;
            res.constrained = fit_model(
                true,
                settings.num_rate_classes,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes,
                settings.multiple_hits,
                pbar.get(),
                "Constrained K=" + std::to_string(settings.num_rate_classes),
                settings.error_sink,
                settings.error_sink_max_weight
            );
            res.unconstrained = fit_model(
                false,
                settings.num_rate_classes,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes,
                settings.multiple_hits,
                pbar.get(),
                "Unconstrained K=" + std::to_string(settings.num_rate_classes),
                settings.error_sink,
                settings.error_sink_max_weight,
                &res.constrained
            );
        }

        if (!settings.phenotype_mode) {
            bool unc_failed = (res.unconstrained.log_likelihood <= -1e19 || !std::isfinite(res.unconstrained.log_likelihood));
            bool con_failed = (res.constrained.log_likelihood <= -1e19 || !std::isfinite(res.constrained.log_likelihood));
            if (unc_failed || con_failed) {
                res.is_failed = true;
                res.error_message = "Optimization failed: model log-likelihood collapsed to penalty threshold (-1e20). The fit did not converge.";
                res.lrt = 0.0;
                res.p_value = 1.0;
            } else if (res.unconstrained.log_likelihood < res.constrained.log_likelihood) {
                res.unconstrained.log_likelihood = res.constrained.log_likelihood;
                res.lrt = 0.0;
                res.p_value = 0.5;
            } else {
                res.lrt = 2.0 * (res.unconstrained.log_likelihood - res.constrained.log_likelihood);
                // Asymptotic null: 0.5 * chi^2_0 + 0.5 * chi^2_2
                res.p_value = 0.5 * std::exp(-res.lrt / 2.0);
            }

            // If --relax was requested on standard BUSTED, fit Shared (k=1) and RELAX (k free)
            if (settings.relax_mode && !res.is_failed) {
                size_t K = res.optimal_k;
                res.has_relax_test = true;
                {
                    auto saved_is_test = is_test_branch;
                    auto saved_num_test = num_test_branches;
                    num_test_branches = 0;
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            bool in_shared = (branch_partition[node.id] <= 1);
                            is_test_branch[node.id] = in_shared;
                            if (in_shared) num_test_branches++;
                        }
                    }
                    res.shared_distributions = fit_model(
                        false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                        settings.multiple_hits, pbar.get(), "RELAX Null (k=1) K=" + std::to_string(K),
                        settings.error_sink, settings.error_sink_max_weight
                    );
                    res.shared_distributions.has_background_distribution = true;
                    res.shared_distributions.background_distribution = res.shared_distributions.test_distribution;
                    is_test_branch = saved_is_test;
                    num_test_branches = saved_num_test;
                }
                res.relax_alternative = fit_relax_coupled_model(
                    res.shared_distributions, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, pbar.get(), "RELAX Coupled-k K=" + std::to_string(K),
                    settings.error_sink, settings.error_sink_max_weight, &res.relax_k
                );
                res.relax_lrt = std::max(0.0, 2.0 * (res.relax_alternative.log_likelihood - res.shared_distributions.log_likelihood));
                res.relax_p_value = busted_chi2_p_value(res.relax_lrt, 1.0);
                if (res.relax_p_value <= settings.p_value_threshold) {
                    res.relax_direction = (res.relax_k > 1.0) ? "Intensification" : "Relaxation";
                } else {
                    res.relax_direction = "None";
                }
            }
        }

        if (pbar) {
            std::ostringstream summary;
            if (settings.phenotype_mode) {
                summary << (res.phenotype_associated ? "\033[1;32m" : "")
                        << (res.phenotype_associated ? "Phenotype-associated selection detected" : "No phenotype-specific selection")
                        << (res.phenotype_associated ? "\033[0m" : "")
                        << " (p_UIT=" << std::fixed << std::setprecision(4) << res.union_intersection_p_value
                        << ", ΔAICc_BG=" << std::setprecision(2) << res.background_delta_aicc
                        << ", p_BG=" << std::setprecision(4) << res.background_p_value;
                if (res.has_relax_test) {
                    summary << ", RELAX k=" << std::setprecision(2) << res.relax_k
                            << " [p=" << std::setprecision(4) << res.relax_p_value << "]";
                }
                summary << ")";
            } else if (res.p_value <= settings.p_value_threshold) {
                summary << "\033[1;32mEvidence of episodic diversifying selection\033[0m: LRT = "
                        << std::fixed << std::setprecision(2) << res.lrt
                        << ", p = " << std::setprecision(4) << res.p_value;
            } else {
                summary << "No evidence of diversifying selection: LRT = "
                        << std::fixed << std::setprecision(2) << res.lrt
                        << ", p = " << std::setprecision(4) << res.p_value;
            }
            pbar->finish(summary.str());
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

        if (settings.error_sink) {
            compute_error_sink_diagnostics(
                res.unconstrained,
                settings.error_sink_bf_threshold,
                res.site_error_bfs,
                res.flagged_error_events
            );
            if (!settings.output_filtered_alignment.empty()) {
                write_filtered_alignment(settings.output_filtered_alignment, res);
            }
        } else {
            res.site_error_bfs.assign(num_sites, 0.0);
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    // Masks branch-site error events flagged by BUSTED-E (BF >= error_sink_bf_threshold)
    // as "---" in all descendant leaf sequences of the flagged branch at that codon site,
    // matching HyPhy 2.5's error-filter.bf workflow.
    Alignment filter_alignment_by_error_sink(const BUSTEDResult& res) const {
        Alignment filtered;
        filtered.code = aln.code;
        filtered.embedded_tree_newick = aln.embedded_tree_newick;
        for (const auto& seq : aln.sequences) {
            filtered.add_sequence(seq.name, seq.raw_seq);
        }

        if (res.flagged_error_events.empty()) {
            filtered.finalize();
            return filtered;
        }

        // Precompute descendant leaf taxon indices for each tree node
        size_t num_nodes = tree.num_nodes();
        std::vector<std::vector<size_t>> node_descendant_taxa(num_nodes);
        for (int32_t nid : tree.post_order) {
            const auto& node = tree.nodes[nid];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                if (t_idx != static_cast<size_t>(-1) && t_idx < filtered.sequences.size()) {
                    node_descendant_taxa[nid].push_back(t_idx);
                }
            } else {
                for (int32_t cid : node.children) {
                    const auto& child_taxa = node_descendant_taxa[cid];
                    node_descendant_taxa[nid].insert(node_descendant_taxa[nid].end(), child_taxa.begin(), child_taxa.end());
                }
            }
        }

        std::unordered_map<std::string, int32_t> name_to_nid;
        for (const auto& node : tree.nodes) {
            if (!node.name.empty()) {
                name_to_nid[node.name] = node.id;
            }
        }

        for (const auto& ev : res.flagged_error_events) {
            if (ev.site == 0 || ev.site > aln.num_codons) continue;
            size_t codon_idx = ev.site - 1;
            size_t nt_offset = codon_idx * 3;
            auto it = name_to_nid.find(ev.branch_name);
            if (it != name_to_nid.end()) {
                int32_t nid = it->second;
                for (size_t t_idx : node_descendant_taxa[nid]) {
                    if (nt_offset + 3 <= filtered.sequences[t_idx].raw_seq.size()) {
                        filtered.sequences[t_idx].raw_seq.replace(nt_offset, 3, "---");
                    }
                }
            }
        }

        filtered.finalize();
        return filtered;
    }

    // Writes the BUSTED-E error-filtered codon alignment to a FASTA file
    void write_filtered_alignment(const std::string& filepath, const BUSTEDResult& res) const {
        Alignment filtered = filter_alignment_by_error_sink(res);
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open filtered alignment file for writing: " + filepath);
        }
        for (const auto& seq : filtered.sequences) {
            out << ">" << seq.name << "\n" << seq.raw_seq << "\n";
        }
    }

nlohmann::json to_legacy_json(const BUSTEDResult& res) const;

nlohmann::json to_modern_json(const BUSTEDResult& res, const Provenance& prov = {}) const;
    nlohmann::json to_json(
        const BUSTEDResult& res,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        if (format == JSONFormat::Legacy) {
            return to_legacy_json(res);
        }
        return to_modern_json(res, prov);
    }

    void save_json(
        const std::string& filepath,
        const BUSTEDResult& res,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open file for writing: " + filepath);
        }
        nlohmann::json j = to_json(res, format, prov);
        out << j.dump(2) << "\n";
    }

};

} // namespace hyphy::analyses
