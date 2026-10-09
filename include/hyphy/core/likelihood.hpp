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
        auto it = aln.taxon_to_index.find(name);
        if (it != aln.taxon_to_index.end()) return it->second;
        // Fallback: case-insensitive match
        for (const auto& [tname, tidx] : aln.taxon_to_index) {
            if (tname.size() == name.size() &&
                std::equal(tname.begin(), tname.end(), name.begin(), [](char a, char b) {
                    return std::toupper(static_cast<unsigned char>(a)) == std::toupper(static_cast<unsigned char>(b));
                })) {
                return tidx;
            }
        }
        return static_cast<size_t>(-1);
    }

    // Nucleotide GTR Tree Likelihood
    static Scalar compute_gtr_log_likelihood(
        const Tree& tree,
        const Alignment& aln,
        const GTRMatrix& gtr_model
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

        // Precompute leaf-to-taxon map for fast direct indexing
        std::vector<size_t> leaf_to_taxon(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = find_taxon_index(aln, node.name);
            }
        }

        Scalar total_log_l = 0.0;

        #pragma omp parallel
        {
            std::vector<Vector4> node_L(num_nodes);
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
                    } else {
                        node_L[node_id].setOnes();
                        for (int32_t child_id : node.children) {
                            Vector4 child_msg = P_branches[child_id] * node_L[child_id];
                            node_L[node_id] = node_L[node_id].cwiseProduct(child_msg);
                        }
                    }
                }

                Scalar pattern_likelihood = gtr_model.pi.dot(node_L[tree.root_id]);
                if (pattern_likelihood > 0.0) {
                    total_log_l += pattern.weight * std::log(pattern_likelihood);
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
        const MG94Matrix& mg94_model
    ) {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = mg94_model.pi.size();

        // Precompute transition matrices for all branches
        std::vector<Matrix> P_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P_branches[node.id] = mg94_model.transition_matrix(node.branch_length);
            }
        }

        // Precompute leaf-to-taxon map
        std::vector<size_t> leaf_to_taxon(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = find_taxon_index(aln, node.name);
            }
        }

        Scalar total_log_l = 0.0;

        #pragma omp parallel
        {
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));
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
                    } else {
                        node_L[node_id].setOnes();
                        for (int32_t child_id : node.children) {
                            Vector child_msg = P_branches[child_id] * node_L[child_id];
                            node_L[node_id] = node_L[node_id].cwiseProduct(child_msg);
                        }
                    }
                }

                Scalar pattern_likelihood = mg94_model.pi.dot(node_L[tree.root_id]);
                if (pattern_likelihood > 0.0) {
                    total_log_l += pattern.weight * std::log(pattern_likelihood);
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

        // 1. Post-order (bottom-up) pass to compute D
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
            } else {
                res.D[nid].setOnes();
                for (int32_t cid : node.children) {
                    Vector child_msg = P_branches[cid] * res.D[cid];
                    res.D[nid] = res.D[nid].cwiseProduct(child_msg);
                }
            }
        }

        res.likelihood = res.D[tree.root_id].dot(pi);

        // 2. Pre-order (top-down) pass to compute U and V
        std::vector<Vector> U(num_nodes, Vector::Zero(S));
        U[tree.root_id] = pi;

        for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
            int32_t uid = *it;
            const auto& node = tree.nodes[uid];
            if (node.is_leaf) continue;

            size_t num_children = node.children.size();
            std::vector<Vector> child_branch_L(num_children);
            for (size_t i = 0; i < num_children; ++i) {
                int32_t cid = node.children[i];
                child_branch_L[i] = P_branches[cid] * res.D[cid];
            }

            for (size_t i = 0; i < num_children; ++i) {
                int32_t vid = node.children[i];
                Vector v_msg = U[uid];
                for (size_t j = 0; j < num_children; ++j) {
                    if (i != j) {
                        v_msg = v_msg.cwiseProduct(child_branch_L[j]);
                    }
                }
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
        const std::vector<Scalar>& branch_lengths = {}
    ) {
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();

        MG94Matrix mg;
        mg.update(params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

        bool use_custom_bl = !branch_lengths.empty();
        std::vector<Matrix> P_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar bl = (use_custom_bl && static_cast<size_t>(node.id) < branch_lengths.size())
                    ? branch_lengths[node.id] : node.branch_length;
                P_branches[node.id] = mg.transition_matrix(bl);
            }
        }

        std::vector<size_t> leaf_to_taxon(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = find_taxon_index(aln, node.name);
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
                auto io = compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_branches, aln.codon_frequencies_f3x4, S
                );
                Scalar L = io.likelihood;
                if (L > 0.0) {
                    local_ll += pattern.weight * std::log(L);
                    Scalar inv_L = pattern.weight / L;
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t vid = node.id;
                            Matrix dP = mg.Q * P_branches[vid];
                            Scalar dL_dt = io.V[vid].dot(dP * io.D[vid]);
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
};

} // namespace hyphy::core
