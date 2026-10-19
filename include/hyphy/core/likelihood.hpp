#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <cctype>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::core {

class LikelihoodEngine {
public:
    static size_t find_taxon_index(const Alignment& aln, const std::string& name) {
        std::string clean_name = name;
        if (clean_name.size() >= 2 && ((clean_name.front() == '\'' && clean_name.back() == '\'') || (clean_name.front() == '"' && clean_name.back() == '"'))) {
            clean_name = clean_name.substr(1, clean_name.size() - 2);
        }
        auto it = aln.taxon_to_index.find(clean_name);
        if (it != aln.taxon_to_index.end()) return it->second;

        // Fast O(1) case-insensitive lookup
        std::string lower = clean_name;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
        auto it_lower = aln.lower_taxon_to_index.find(lower);
        if (it_lower != aln.lower_taxon_to_index.end()) return it_lower->second;

        return static_cast<size_t>(-1);
    }

    static std::vector<size_t> build_leaf_to_taxon_map(const Tree& tree, const Alignment& aln) {
        size_t num_nodes = tree.num_nodes();
        std::vector<size_t> leaf_to_taxon(num_nodes, static_cast<size_t>(-1));
        size_t leaves = 0, matched = 0;
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaves++;
                size_t idx = find_taxon_index(aln, node.name);
                leaf_to_taxon[node.id] = idx;
                if (idx != static_cast<size_t>(-1)) matched++;
            }
        }
        if (leaves > 0 && matched == 0) {
            throw std::runtime_error("Phylogenetic tree leaves do not match any sequences in alignment! (0/" + std::to_string(leaves) + " matched)");
        }
        return leaf_to_taxon;
    }

    // Nucleotide GTR Tree Likelihood
    static Scalar compute_gtr_log_likelihood(
        const Tree& tree,
        const Alignment& aln,
        const GTRMatrix& gtr_model,
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr
    ) {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();

        // Precompute transition matrices for all branches
        std::vector<Matrix4> P_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P_branches[node.id] = gtr_model.transition_matrix(node.branch_length);
            }
        }

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        Scalar total_log_l = 0.0;

        #pragma omp parallel
        {
            std::vector<Vector4> node_L(num_nodes);
            std::vector<double> node_log_scale(num_nodes, 0.0);
            #pragma omp for reduction(+:total_log_l) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];

                for (int32_t node_id : tree.post_order) {
                    const auto& node = tree.nodes[node_id];

                    if (node.is_leaf) {
                        size_t t_idx = leaf_to_taxon[node_id];
                        if (t_idx != static_cast<size_t>(-1)) {
                            int8_t state = pattern.states[t_idx];
                            if (state >= 0 && state < 4) {
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
                        tree.eval_internal_node(node_id, P_branches, node_L, node_log_scale, node_L[node_id], child_scale_sum);
                        node_log_scale[node_id] = child_scale_sum;
                    }
                }

                Scalar pattern_likelihood = gtr_model.pi.dot(node_L[tree.root_id]);
                if (pattern_likelihood > 0.0) {
                    total_log_l += pattern.weight * (std::log(pattern_likelihood) + node_log_scale[tree.root_id]);
                } else {
                    total_log_l += pattern.weight * (-1e20);
                }
            }
        }

        return total_log_l;
    }

    // Codon MG94 Tree Likelihood
    static Scalar compute_mg94_log_likelihood(
        const Tree& tree,
        const Alignment& aln,
        const MG94Matrix& mg94_model,
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr,
        Scalar branch_scale = 1.0
    ) {
        size_t num_nodes = tree.num_nodes();
        int S = mg94_model.pi.size();

        // Precompute transition matrices for all branches
        std::vector<Matrix> P_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P_branches[node.id] = mg94_model.transition_matrix(node.branch_length * branch_scale);
            }
        }

        return compute_tree_log_likelihood_from_P(tree, aln, P_branches, mg94_model.pi, S, custom_leaf_to_taxon);
    }

    // Fast single-pass codon tree likelihood from precomputed transition matrices P_branches
    // Eliminates pre-order traversal and intermediate heap allocations
    static Scalar compute_tree_log_likelihood_from_P(
        const Tree& tree,
        const Alignment& aln,
        const std::vector<Matrix>& P_branches,
        const Vector& pi,
        int S,
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr
    ) {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double LOG_SCALER = std::log(SCALER_UP);

        Scalar total_log_l = 0.0;

        #pragma omp parallel
        {
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));
            std::vector<double> node_log_scale(num_nodes, 0.0);

            #pragma omp for reduction(+:total_log_l) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];

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
                        tree.eval_internal_node(node_id, P_branches, node_L, node_log_scale, node_L[node_id], child_scale_sum);

                        double max_val = node_L[node_id].maxCoeff();
                        int boost_count = 0;
                        while (max_val < SCALER_THRESH && max_val > 0.0) {
                            node_L[node_id] *= SCALER_UP;
                            max_val *= SCALER_UP;
                            boost_count++;
                        }
                        node_log_scale[node_id] = child_scale_sum - boost_count * LOG_SCALER;
                    }
                }

                Scalar pattern_likelihood = pi.dot(node_L[tree.root_id]);
                if (pattern_likelihood > 0.0) {
                    total_log_l += pattern.weight * (std::log(pattern_likelihood) + node_log_scale[tree.root_id]);
                } else {
                    total_log_l += pattern.weight * (-1e20);
                }
            }
        }

        return total_log_l;
    }

    // Inside-Outside two-pass traversal for a single site pattern
    // Returns downward vectors D, parent-sibling context vectors V, and total likelihood L
    struct InsideOutsideResult {
        std::vector<Vector> D; // Downward subtree likelihoods D[node]
        std::vector<Vector> V; // Upward parent-sibling messages V[node]
        Scalar likelihood = 0.0;
        Scalar log_likelihood = 0.0;
    };

    static InsideOutsideResult compute_inside_outside(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<Matrix>& P_branches,
        const Vector& pi,
        int S
    ) {
        size_t num_nodes = tree.num_nodes();
        InsideOutsideResult res;
        res.D.assign(num_nodes, Vector::Zero(S));
        res.V.assign(num_nodes, Vector::Zero(S));
        std::vector<Vector> P_D(num_nodes, Vector::Zero(S));
        std::vector<double> B(num_nodes, 1.0);
        std::vector<double> node_log_scale(num_nodes, 0.0);

        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double LOG_SCALER = std::log(SCALER_UP);

        // 1. Post-order (bottom-up) pass to compute D and cache P_D
        for (int32_t nid : tree.post_order) {
            const auto& node = tree.nodes[nid];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                if (t_idx != static_cast<size_t>(-1)) {
                    int8_t state = pattern.states[t_idx];
                    if (state >= 0 && state < S) {
                        res.D[nid].setZero();
                        res.D[nid](state) = 1.0;
                    } else {
                        res.D[nid].setOnes();
                    }
                } else {
                    res.D[nid].setOnes();
                }
                node_log_scale[nid] = 0.0;
                B[nid] = 1.0;
            } else {
                double child_scale_sum = 0.0;
                for (int32_t cid : node.children) {
                    P_D[cid] = P_branches[cid] * res.D[cid];
                }
                tree.eval_internal_node(nid, P_branches, res.D, node_log_scale, res.D[nid], child_scale_sum);

                double max_val = res.D[nid].maxCoeff();
                int boost_count = 0;
                double boost_mult = 1.0;
                while (max_val < SCALER_THRESH && max_val > 0.0) {
                    res.D[nid] *= SCALER_UP;
                    max_val *= SCALER_UP;
                    boost_mult *= SCALER_UP;
                    boost_count++;
                }
                B[nid] = boost_mult;
                node_log_scale[nid] = child_scale_sum - boost_count * LOG_SCALER;
            }
        }

        res.likelihood = res.D[tree.root_id].dot(pi);
        if (res.likelihood > 0.0) {
            res.log_likelihood = std::log(res.likelihood) + node_log_scale[tree.root_id];
        } else {
            res.log_likelihood = -1e20;
        }

        // 2. Pre-order (top-down) pass to compute U and V reusing cached P_D
        std::vector<Vector> U(num_nodes, Vector::Zero(S));
        U[tree.root_id] = pi;

        for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
            int32_t uid = *it;
            const auto& node = tree.nodes[uid];
            if (node.is_leaf) continue;

            size_t num_children = node.children.size();
            for (size_t i = 0; i < num_children; ++i) {
                int32_t vid = node.children[i];
                Vector v_msg = U[uid];
                for (size_t j = 0; j < num_children; ++j) {
                    if (i != j) {
                        v_msg.array() *= P_D[node.children[j]].array();
                    }
                }
                v_msg *= B[uid];
                res.V[vid] = v_msg;
                U[vid] = P_branches[vid].transpose() * v_msg;
            }
        }

        return res;
    }

    // Computes posterior branch responsibilities gamma_{b, k} for K mixture components
    // for a single pattern in a single inside-outside pass.
    // Returns: responsibilities (num_nodes x K) and pattern likelihood L.
    static std::pair<Matrix, Scalar> compute_pattern_mixture_responsibilities(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<std::vector<Matrix>>& P_components, // [K][num_nodes]
        const std::vector<Scalar>& weights,                    // [K]
        const Vector& pi,
        int S
    ) {
        size_t num_nodes = tree.num_nodes();
        size_t K = weights.size();

        // Build mixture transition matrices P_mix[b] = sum_k w_k P_k[b]
        std::vector<Matrix> P_mix(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P_mix[node.id] = Matrix::Zero(S, S);
                for (size_t k = 0; k < K; ++k) {
                    P_mix[node.id] += weights[k] * P_components[k][node.id];
                }
            }
        }

        // Run Inside-Outside pass
        auto io_res = compute_inside_outside(tree, pattern, leaf_to_taxon, P_mix, pi, S);

        Matrix gamma = Matrix::Zero(num_nodes, K);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                int32_t nid = node.id;
                Vector comp_L(K);
                Scalar denom = 0.0;
                for (size_t k = 0; k < K; ++k) {
                    Scalar L_k = io_res.V[nid].dot(P_components[k][nid] * io_res.D[nid]);
                    comp_L(k) = weights[k] * L_k;
                    denom += comp_L(k);
                }
                if (denom > 1e-300) {
                    for (size_t k = 0; k < K; ++k) {
                        gamma(nid, k) = comp_L(k) / denom;
                    }
                } else {
                    for (size_t k = 0; k < K; ++k) {
                        gamma(nid, k) = weights[k];
                    }
                }
            }
        }

        return {gamma, io_res.likelihood};
    }

    // Evaluates log-likelihood and exact analytical branch gradients d ln L / d t_b
    // using the 2-pass Inside-Outside adjoint state formulation
    static std::pair<Scalar, std::vector<Scalar>> compute_branch_length_gradients(
        const Tree& tree,
        const Alignment& aln,
        const MG94Parameters& params,
        const std::vector<Scalar>& branch_lengths = {},
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr
    ) {
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();

        MG94Matrix mg;
        mg.update(params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

        bool use_custom_bl = !branch_lengths.empty();
        std::vector<Matrix> P_branches(num_nodes);
        std::vector<Matrix> dP_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar bl = (use_custom_bl && static_cast<size_t>(node.id) < branch_lengths.size())
                    ? branch_lengths[node.id] : node.branch_length;
                P_branches[node.id] = mg.transition_matrix(bl);
                dP_branches[node.id] = mg.Q * P_branches[node.id];
            }
        }

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        Scalar total_log_l = 0.0;
        std::vector<Scalar> grad_b(num_nodes, 0.0);

        #pragma omp parallel
        {
            std::vector<Scalar> local_grad(num_nodes, 0.0);
            Scalar local_ll = 0.0;

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                auto io = compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_branches, aln.codon_frequencies_f3x4, S
                );
                Scalar L = io.likelihood;
                if (L > 0.0) {
                    local_ll += pattern.weight * io.log_likelihood;
                    Scalar inv_L = pattern.weight / L;
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t vid = node.id;
                            Scalar dL_dt = io.V[vid].dot(dP_branches[vid] * io.D[vid]);
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

    struct JointMG94Gradients {
        Scalar log_likelihood = 0.0;
        std::vector<Scalar> grad_branches; // d ln L / d t_b (in substitutions/site)
        Scalar grad_omega = 0.0;           // d ln L / d omega
        Scalar grad_theta_AC = 0.0;        // d ln L / d theta_AC
        Scalar grad_theta_AT = 0.0;        // d ln L / d theta_AT
        Scalar grad_theta_CG = 0.0;        // d ln L / d theta_CG
        Scalar grad_theta_CT = 0.0;        // d ln L / d theta_CT
        Scalar grad_theta_GT = 0.0;        // d ln L / d theta_GT
    };

    // Evaluates log-likelihood and exact joint analytical gradients for all branch lengths
    // and all 6 global parameters (omega, theta_AC, theta_AT, theta_CG, theta_CT, theta_GT)
    // using the 2-pass Inside-Outside adjoint state and spectral Fréchet derivatives of the matrix exponential.
    static JointMG94Gradients compute_joint_frechet_gradients(
        const Tree& tree,
        const Alignment& aln,
        const MG94Parameters& params,
        const std::vector<Scalar>& branch_lengths = {},
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr
    ) {
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();

        MG94Matrix mg;
        mg.update(params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

        Scalar mu = mg.scale_factor;
        Scalar conv = (mu > 1e-12) ? (3.0 / mu) : 1.0;

        std::vector<Scalar> tau(num_nodes, 0.0);
        std::vector<Matrix> P_branches(num_nodes);
        std::vector<Matrix> dP_branches(num_nodes);

        bool use_custom_bl = !branch_lengths.empty();
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                int32_t bid = node.id;
                Scalar bl = (use_custom_bl && static_cast<size_t>(bid) < branch_lengths.size())
                    ? branch_lengths[bid] : node.branch_length;
                tau[bid] = bl * conv;
                P_branches[bid] = mg.transition_matrix(tau[bid]);
                dP_branches[bid] = mg.Q * P_branches[bid];
            }
        }

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        Scalar total_log_l = 0.0;
        std::vector<Scalar> dL_dtau(num_nodes, 0.0);
        std::vector<Matrix> A_branches(num_nodes, Matrix::Zero(S, S));

        #pragma omp parallel
        {
            Scalar local_ll = 0.0;
            std::vector<Scalar> local_dL(num_nodes, 0.0);
            std::vector<Matrix> local_A(num_nodes, Matrix::Zero(S, S));

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                auto io = compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_branches, aln.codon_frequencies_f3x4, S
                );
                Scalar L = io.likelihood;
                if (L > 0.0) {
                    local_ll += pattern.weight * io.log_likelihood;
                    Scalar inv_L = pattern.weight / L;

                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t bid = node.id;
                            Scalar dL_dt = io.V[bid].dot(dP_branches[bid] * io.D[bid]);
                            local_dL[bid] += inv_L * dL_dt;
                            local_A[bid].noalias() += inv_L * (io.V[bid] * io.D[bid].transpose());
                        }
                    }
                }
            }

            #pragma omp critical
            {
                total_log_l += local_ll;
                for (size_t i = 0; i < num_nodes; ++i) {
                    dL_dtau[i] += local_dL[i];
                    A_branches[i] += local_A[i];
                }
            }
        }

        Matrix H = Matrix::Zero(S, S);
        Scalar sum_tau_dL_dtau = 0.0;

        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                int32_t bid = node.id;
                Scalar tb = tau[bid];
                sum_tau_dL_dtau += tb * dL_dtau[bid];

                Matrix A_tilde = mg.U_pi.transpose() * A_branches[bid] * mg.V_pi.transpose();

                Matrix E(S, S);
                for (int i = 0; i < S; ++i) {
                    for (int j = 0; j < S; ++j) {
                        Scalar diff = mg.eigenvalues(i) - mg.eigenvalues(j);
                        if (std::abs(diff) < 1e-9) {
                            E(i, j) = tb * std::exp(mg.eigenvalues(i) * tb);
                        } else {
                            E(i, j) = (std::exp(mg.eigenvalues(i) * tb) - std::exp(mg.eigenvalues(j) * tb)) / diff;
                        }
                    }
                }

                H.noalias() += A_tilde.cwiseProduct(E);
            }
        }

        auto compute_global_grad = [&](const Matrix& dQ_k) -> Scalar {
            Matrix M_k = mg.V_pi * dQ_k * mg.U_pi;
            Scalar part1 = (M_k.cwiseProduct(H)).sum();
            Scalar dmu_k = 0.0;
            for (int i = 0; i < S; ++i) {
                dmu_k -= mg.pi(i) * dQ_k(i, i);
            }
            Scalar part2 = -(dmu_k / mu) * sum_tau_dL_dtau;
            return part1 + part2;
        };

        Scalar r_nuc[4][4] = {
            {0.0, params.theta_AC, 1.0, params.theta_AT},
            {params.theta_AC, 0.0, params.theta_CG, params.theta_CT},
            {1.0, params.theta_CG, 0.0, params.theta_GT},
            {params.theta_AT, params.theta_CT, params.theta_GT, 0.0}
        };

        Matrix dQ_omega = Matrix::Zero(S, S);
        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                const auto& diff = gcode.diff_matrix[i][j];
                if (diff.num_diffs == 1 && !diff.synonymous) {
                    Scalar nuc_rate = r_nuc[diff.diff_from[0]][diff.diff_to[0]];
                    Scalar target_nuc_freq = aln.pos_nuc_frequencies(diff.diff_pos[0], diff.diff_to[0]);
                    dQ_omega(i, j) = 1.0 * nuc_rate * target_nuc_freq;
                }
            }
            dQ_omega(i, i) = -dQ_omega.row(i).sum();
        }

        auto build_dQ_nuc = [&](int u, int v) -> Matrix {
            Matrix dQ = Matrix::Zero(S, S);
            for (int i = 0; i < S; ++i) {
                for (int j = 0; j < S; ++j) {
                    if (i == j) continue;
                    const auto& diff = gcode.diff_matrix[i][j];
                    if (diff.num_diffs == 1) {
                        int from = diff.diff_from[0];
                        int to = diff.diff_to[0];
                        if ((from == u && to == v) || (from == v && to == u)) {
                            Scalar rate_mod = diff.synonymous ? params.alpha : params.beta;
                            Scalar target_nuc_freq = aln.pos_nuc_frequencies(diff.diff_pos[0], to);
                            dQ(i, j) = rate_mod * 1.0 * target_nuc_freq;
                        }
                    }
                }
                dQ(i, i) = -dQ.row(i).sum();
            }
            return dQ;
        };

        JointMG94Gradients res;
        res.log_likelihood = total_log_l;
        res.grad_branches.resize(num_nodes, 0.0);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                res.grad_branches[node.id] = conv * dL_dtau[node.id];
            }
        }

        res.grad_omega = compute_global_grad(dQ_omega);
        res.grad_theta_AC = compute_global_grad(build_dQ_nuc(0, 1));
        res.grad_theta_AT = compute_global_grad(build_dQ_nuc(0, 3));
        res.grad_theta_CG = compute_global_grad(build_dQ_nuc(1, 2));
        res.grad_theta_CT = compute_global_grad(build_dQ_nuc(1, 3));
        res.grad_theta_GT = compute_global_grad(build_dQ_nuc(2, 3));

        return res;
    }

    /**
     * Compute stationary distribution for ExpCM at a single site.
     *
     * p_{r,x} = (π_{r,A(x)})^β · φ_{x1} φ_{x2} φ_{x3} / Z_r
     *
     * @param gcode        Genetic code
     * @param prefs        Site preferences (20 AAs)
     * @param site         Site index (0-based)
     * @param beta         Stringency parameter
     * @param phi          Nucleotide frequencies (4-element vector)
     * @return             Stationary distribution (61 codons)
     */
    static Eigen::Matrix<Scalar, 61, 1> compute_expcm_stationary(
        const GeneticCode& gcode,
        const Preferences& prefs,
        size_t site,
        Scalar beta,
        const Vector4& phi
    ) {
        Eigen::Matrix<Scalar, 61, 1> p = Eigen::Matrix<Scalar, 61, 1>::Zero();
        Scalar Z = 0.0;

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            char aa = gcode.sense_codon_aa[x];
            Scalar pref_beta = std::pow(prefs.at(site, aa), beta);

            // Nucleotide product: φ_{x1} · φ_{x2} · φ_{x3}
            const std::string& codon = gcode.sense_codons[x];
            int8_t n1 = char_to_nuc(codon[0]);
            int8_t n2 = char_to_nuc(codon[1]);
            int8_t n3 = char_to_nuc(codon[2]);
            Scalar nuc_freq = phi(n1) * phi(n2) * phi(n3);

            p(x) = pref_beta * nuc_freq;
            Z += p(x);
        }

        if (Z > 0.0) {
            p /= Z;  // Normalize
        }

        return p;
    }

    /**
     * Compute ExpCM log-likelihood at fixed parameters.
     *
     * Iterates over sites, building per-site P matrices and computing likelihoods
     * via Felsenstein pruning. Uses scaled matrix exponential for transition matrices.
     *
     * @param tree         Phylogenetic tree (with branch lengths)
     * @param aln          Codon alignment
     * @param gcode        Genetic code
     * @param prefs        Site-specific preferences
     * @param params       ExpCM parameters (beta, omega, kappa, phi)
     * @param branch_scale Factor to scale all branch lengths (default 1.0)
     * @return             Log-likelihood
     */
    static Scalar compute_expcm_log_likelihood(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0
    ) {
        params.validate();

        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        size_t num_sites = prefs.num_sites();

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = build_leaf_to_taxon_map(tree, aln);

        // Build Q matrix (shared across all sites)
        ExpCMQMatrix Q_matrix;
        Q_matrix.build(gcode, params);
        const Matrix61Codon& Q = Q_matrix.Q;

        Scalar total_log_l = 0.0;
        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double LOG_SCALER = std::log(SCALER_UP);

        // Process each alignment site (codon position)
        // Each site has a specific pattern and requires its own P matrix
        #pragma omp parallel for reduction(+:total_log_l) schedule(dynamic)
        for (size_t site = 0; site < num_sites; ++site) {
            // 1. Get the pattern for this site
            size_t pattern_idx = aln.site_to_pattern[site];
            const auto& pattern = aln.patterns[pattern_idx];

            // 2. Build F matrix for this site
            ExpCMFMatrix F_matrix;
            F_matrix.build(gcode, prefs, site, params);

            // 3. Compute P = F ⊙ Q (element-wise product)
            Matrix61Codon P = F_matrix.F.cwiseProduct(Q);

            // 4. Set diagonals so rows sum to 0
            for (int x = 0; x < 61; ++x) {
                P(x, x) = -P.row(x).sum();
            }

            // 5. Compute stationary distribution for this site
            Eigen::Matrix<Scalar, 61, 1> pi_site = compute_expcm_stationary(
                gcode, prefs, site, params.beta, params.phi
            );

            // 6. Compute transition matrices for each branch using scaled matrix exponential
            std::vector<Matrix61Codon> P_branches(num_nodes);
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar scaled_bl = node.branch_length * branch_scale;
                    // Compute exp(P · scaled_bl) via scaling and squaring
                    P_branches[node.id] = matrix_exp_scaling_squaring(P * scaled_bl);
                }
            }

            // 7. Compute site likelihood via Felsenstein pruning
            std::vector<Eigen::Matrix<Scalar, 61, 1>> node_L(num_nodes);
            std::vector<double> node_log_scale(num_nodes, 0.0);

            for (int32_t node_id : tree.post_order) {
                const auto& node = tree.nodes[node_id];

                if (node.is_leaf) {
                    size_t t_idx = leaf_to_taxon[node_id];
                    if (t_idx != static_cast<size_t>(-1)) {
                        int8_t state = pattern.states[t_idx];
                        if (state >= 0 && state < 61) {
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
                    // Internal node: use tree's eval_internal_node function
                    double child_scale_sum = 0.0;
                    tree.eval_internal_node(node_id, P_branches, node_L, node_log_scale, node_L[node_id], child_scale_sum);

                    // Rescale if underflow
                    double max_val = node_L[node_id].maxCoeff();
                    int boost_count = 0;
                    while (max_val < SCALER_THRESH && max_val > 0.0) {
                        node_L[node_id] *= SCALER_UP;
                        max_val *= SCALER_UP;
                        boost_count++;
                    }
                    node_log_scale[node_id] = child_scale_sum - boost_count * LOG_SCALER;
                }
            }

            // Compute site likelihood: pi_site · node_L[root]
            Scalar site_likelihood = pi_site.dot(node_L[tree.root_id]);
            if (site_likelihood > 0.0) {
                total_log_l += std::log(site_likelihood) + node_log_scale[tree.root_id];
            } else {
                total_log_l += (-1e20);
            }
        }

        return total_log_l;
    }

    /**
     * Compute matrix exponential using scaling and squaring method.
     * Robust for general matrices, including those with complex eigenvalues.
     */
    static Matrix61Codon matrix_exp_scaling_squaring(const Matrix61Codon& M) {
        // Scaling and squaring method from EXPOKIT
        // For simplicity, use repeated squaring with scaling

        const int max_scale = 30;
        Scalar norm_inf = M.rowwise().lpNorm<1>().maxCoeff();

        if (norm_inf < 1e-10) {
            // For very small matrices, use I + M
            return Matrix61Codon::Identity() + M;
        }

        // Find s such that ||M/2^s|| < 1
        int s = 0;
        Scalar one_half_power = 1.0;
        while (norm_inf > 1.0 && s < max_scale) {
            norm_inf *= 0.5;
            one_half_power *= 0.5;
            s++;
        }

        // Compute exp(M/2^s) using Padé approximation
        Matrix61Codon M_scaled = M * one_half_power;
        Matrix61Codon result = pade_exp(M_scaled);

        // Square result s times to get exp(M)
        for (int i = 0; i < s; ++i) {
            result = result * result;
        }

        return result;
    }

    /**
     * Padé (5,5) approximation for matrix exponential.
     * Higher-order approximation for better accuracy with scaled matrices.
     */
    static Matrix61Codon pade_exp(const Matrix61Codon& X) {
        Matrix61Codon I = Matrix61Codon::Identity();
        Matrix61Codon X2 = X * X;
        Matrix61Codon X3 = X2 * X;
        Matrix61Codon X4 = X3 * X;
        Matrix61Codon X5 = X4 * X;

        // (5,5) Padé approximation coefficients
        // exp(X) ≈ N(X) / D(X)
        // N(X) = I + c1*X + c2*X² + c3*X³ + c4*X⁴ + c5*X⁵
        // D(X) = I - c1*X + c2*X² - c3*X³ + c4*X⁴ - c5*X⁵
        // where: c1=1/2, c2=1/9, c3=1/72, c4=1/1008, c5=1/30240

        // Using integer arithmetic to avoid floating point errors:
        // Multiply everything by 30240 (LCM of denominators)
        // N = 30240*I + 15120*X + 3360*X² + 420*X³ + 30*X⁴ + X⁵
        // D = 30240*I - 15120*X + 3360*X² - 420*X³ + 30*X⁴ - X⁵

        Matrix61Codon num = 30240.0 * I + 15120.0 * X + 3360.0 * X2 + 420.0 * X3 + 30.0 * X4 + X5;
        Matrix61Codon den = 30240.0 * I - 15120.0 * X + 3360.0 * X2 - 420.0 * X3 + 30.0 * X4 - X5;

        return den.inverse() * num;
    }
};

} // namespace hyphy::core
