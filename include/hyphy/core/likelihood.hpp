#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <utility>
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

    // Per-thread traversal scratch: node likelihood vectors + per-node log
    // scales, destructured at each call site to keep the usual local names:
    //   auto [node_L, node_log_scale] = make_pattern_scratch<Vector>(n, S);
    template <typename VecType>
    static std::pair<std::vector<VecType>, std::vector<double>>
    make_pattern_scratch(size_t num_nodes, int S) {
        return {std::vector<VecType>(num_nodes, VecType::Zero(S)),
                std::vector<double>(num_nodes, 0.0)};
    }

    // Max-shifted log-sum-exp over a container of log-domain terms (each term
    // already folds in log class weight and node log-scale). Entries equal to
    // -infinity mark absent or degenerate classes; -infinity is returned when
    // no term is finite so callers keep their own sentinel.
    template <typename Terms>
    static Scalar logsumexp(const Terms& terms) {
        Scalar max_t = -std::numeric_limits<Scalar>::infinity();
        for (const Scalar t : terms) {
            if (t > max_t) max_t = t;
        }
        if (!std::isfinite(max_t)) return max_t;
        Scalar sum_exp = 0.0;
        for (const Scalar t : terms) {
            if (std::isfinite(t)) sum_exp += std::exp(t - max_t);
        }
        return (sum_exp > 0.0)
            ? max_t + std::log(sum_exp)
            : -std::numeric_limits<Scalar>::infinity();
    }

    // Standard per-branch transition set: P[node.id] =
    // model.transition_matrix(node.branch_length * branch_scale) for every
    // non-root node. The root entry stays default-constructed; traversals
    // never read it.
    template <typename MatType, typename Model>
    static std::vector<MatType> build_P_branches(
        const Tree& tree, const Model& model, Scalar branch_scale = 1.0
    ) {
        std::vector<MatType> P(tree.nodes.size());
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P[node.id] = model.transition_matrix(node.branch_length * branch_scale);
            }
        }
        return P;
    }

    // One post-order pass over the tree for a single site pattern: leaf
    // likelihoods are initialized from pattern.states and internal-node
    // partials are computed via eval_internal_node (which rescales underflowing
    // partials internally and folds the result into node_log_scale).
    // Pass init_leaves=false on repeat traversals of the same pattern when leaf
    // likelihoods are already populated — e.g. rate-class/mixture loops that
    // re-evaluate internals under a different P_branches set.
    template <typename VecType, typename MatType>
    static void eval_pattern_postorder(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<MatType>& P_branches,
        int S,
        std::vector<VecType>& node_L,
        std::vector<double>& node_log_scale,
        bool init_leaves = true
    ) {
        for (int32_t node_id : tree.post_order) {
            const auto& node = tree.nodes[node_id];
            if (node.is_leaf) {
                if (!init_leaves) continue;
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
                node_log_scale[node_id] = child_scale_sum;
            }
        }
    }

    // Per-pattern log-likelihood under an MG94 model: postorder pass, root
    // mixture over equilibrium frequencies, -1e20 sentinel when the root
    // conditional likelihood is not positive.
    static Scalar compute_pattern_log_likelihood(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const MG94Matrix& mg94_model,
        const std::vector<Matrix>& P_branches,
        std::vector<Vector>& node_L,
        std::vector<double>& node_log_scale
    ) {
        eval_pattern_postorder(tree, pattern, leaf_to_taxon, P_branches,
                               static_cast<int>(mg94_model.pi.size()), node_L, node_log_scale);
        Scalar pat_l = mg94_model.pi.dot(node_L[tree.root_id]);
        return (pat_l > 0.0) ? std::log(pat_l) + node_log_scale[tree.root_id] : -1e20;
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
        auto P_branches = build_P_branches<Matrix4>(tree, gtr_model);

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        Scalar total_log_l = 0.0;

        #pragma omp parallel
        {
            auto [node_L, node_log_scale] = make_pattern_scratch<Vector4>(num_nodes, 4);
            #pragma omp for reduction(+:total_log_l) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];

                eval_pattern_postorder(tree, pattern, leaf_to_taxon, P_branches, 4, node_L, node_log_scale);

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
        int S = mg94_model.pi.size();

        // Precompute transition matrices for all branches
        auto P_branches = build_P_branches<Matrix>(tree, mg94_model, branch_scale);

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

        Scalar total_log_l = 0.0;

        #pragma omp parallel
        {
            auto [node_L, node_log_scale] = make_pattern_scratch<Vector>(num_nodes, S);

            #pragma omp for reduction(+:total_log_l) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];

                eval_pattern_postorder(tree, pattern, leaf_to_taxon, P_branches, S, node_L, node_log_scale);

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
    template <typename VecType>
    struct InsideOutsideResultT {
        std::vector<VecType> D; // Downward subtree likelihoods D[node]
        std::vector<VecType> V; // Upward parent-sibling messages V[node]
        Scalar likelihood = 0.0;
        Scalar log_likelihood = 0.0;
    };
    using InsideOutsideResult = InsideOutsideResultT<Vector>;

    // Inside-outside (sum-product) pass for one pattern. When out_P_D is
    // non-null it is filled with the per-child messages P_b[c] * D[c] computed
    // during the inside pass — callers needing pairwise child products (e.g.
    // polytomy resolution) reuse them instead of recomputing.
    template <typename VecType, typename MatType>
    static InsideOutsideResultT<VecType> compute_inside_outside(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<MatType>& P_branches,
        const VecType& pi,
        int S,
        std::vector<VecType>* out_P_D = nullptr
    ) {
        size_t num_nodes = tree.num_nodes();
        InsideOutsideResultT<VecType> res;
        res.D.assign(num_nodes, VecType::Zero(S));
        res.V.assign(num_nodes, VecType::Zero(S));
        std::vector<VecType> local_P_D;
        std::vector<VecType>& P_D = out_P_D ? *out_P_D : local_P_D;
        P_D.assign(num_nodes, VecType::Zero(S));
        std::vector<double> B(num_nodes, 1.0);
        std::vector<double> node_log_scale(num_nodes, 0.0);

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
                double boost_mult = 1.0;
                for (int32_t cid : node.children) {
                    P_D[cid] = P_branches[cid] * res.D[cid];
                }
                tree.eval_internal_node(nid, P_branches, res.D, node_log_scale, res.D[nid], child_scale_sum, &boost_mult);
                // B[nid] must record the node's own boost so the pre-order pass
                // re-applies the same scaling to outside messages (v_msg *= B[uid]).
                B[nid] = boost_mult;
                node_log_scale[nid] = child_scale_sum;
            }
        }

        res.likelihood = res.D[tree.root_id].dot(pi);
        if (res.likelihood > 0.0) {
            res.log_likelihood = std::log(res.likelihood) + node_log_scale[tree.root_id];
        } else {
            res.log_likelihood = -1e20;
        }

        // 2. Pre-order (top-down) pass to compute U and V reusing cached P_D
        std::vector<VecType> U(num_nodes, VecType::Zero(S));
        U[tree.root_id] = pi;

        for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
            int32_t uid = *it;
            const auto& node = tree.nodes[uid];
            if (node.is_leaf) continue;

            size_t num_children = node.children.size();
            for (size_t i = 0; i < num_children; ++i) {
                int32_t vid = node.children[i];
                VecType v_msg = U[uid];
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
};

} // namespace hyphy::core
