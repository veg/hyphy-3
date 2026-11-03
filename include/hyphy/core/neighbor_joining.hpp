#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace hyphy::core {

class NeighborJoining {
public:
    // Computes pairwise Tamura-Nei 93 (TN93) distances between sequences in an alignment
    // Falls back to Jukes-Cantor if logarithms are undefined.
    static Matrix compute_tn93_distances(const Alignment& aln) {
        size_t N = aln.num_taxa;
        Matrix D = Matrix::Zero(N, N);
        if (N <= 1) return D;

        // Base frequencies
        Vector4 freqs = aln.nuc_frequencies;
        Scalar f_sum = freqs.sum();
        if (f_sum <= 0.0) {
            freqs.setConstant(0.25);
        } else {
            freqs /= f_sum;
        }

        Scalar fA = freqs(0);
        Scalar fC = freqs(1);
        Scalar fG = freqs(2);
        Scalar fT = freqs(3);

        Scalar fY = fC + fT;
        Scalar fR = fA + fG;

        bool use_k2p = (freqs.minCoeff() <= 1e-6 || fY <= 1e-6 || fR <= 1e-6);
        Scalar K1 = 0.0, K2 = 0.0, K3 = 0.0;
        if (!use_k2p) {
            K1 = 2.0 * fA * fG / fR;
            K2 = 2.0 * fC * fT / fY;
            K3 = 2.0 * (fR * fY - (fA * fG * fY / fR) - (fC * fT * fR / fY));
        }

        size_t num_sites = aln.num_nucleotides;

        for (size_t s1 = 0; s1 < N; ++s1) {
            const std::string& seq1 = aln.sequences[s1].raw_seq;
            for (size_t s2 = s1 + 1; s2 < N; ++s2) {
                const std::string& seq2 = aln.sequences[s2].raw_seq;

                size_t sites_compared = 0;
                size_t n_AG = 0;
                size_t n_CT = 0;
                size_t n_tv = 0;

                for (size_t k = 0; k < num_sites; ++k) {
                    int8_t n1 = char_to_nuc(seq1[k]);
                    int8_t n2 = char_to_nuc(seq2[k]);
                    if (n1 >= 0 && n1 < 4 && n2 >= 0 && n2 < 4) {
                        sites_compared++;
                        if (n1 != n2) {
                            if ((n1 == 0 && n2 == 2) || (n1 == 2 && n2 == 0)) {
                                n_AG++;
                            } else if ((n1 == 1 && n2 == 3) || (n1 == 3 && n2 == 1)) {
                                n_CT++;
                            } else {
                                n_tv++;
                            }
                        }
                    }
                }

                Scalar d = 1.0;
                if (sites_compared > 0) {
                    Scalar AG = static_cast<Scalar>(n_AG) / sites_compared;
                    Scalar CT = static_cast<Scalar>(n_CT) / sites_compared;
                    Scalar TV = static_cast<Scalar>(n_tv) / sites_compared;

                    bool valid_tn93 = false;
                    if (!use_k2p) {
                        Scalar d1 = std::max<Scalar>(0.02, 1.0 - (K1 > 1e-12 ? (AG / K1) : 0.0) - (0.5 * TV / fR));
                        Scalar d2 = std::max<Scalar>(0.02, 1.0 - (K2 > 1e-12 ? (CT / K2) : 0.0) - (0.5 * TV / fY));
                        Scalar d3 = std::max<Scalar>(0.02, 1.0 - (0.5 * TV / (fR * fY)));

                        d = -K1 * std::log(d1) - K2 * std::log(d2) - K3 * std::log(d3);
                        if (d >= 0.0 && !std::isnan(d) && !std::isinf(d)) {
                            valid_tn93 = true;
                        }
                    }

                    if (!valid_tn93) {
                        // Fall back to continuous clamped Jukes-Cantor
                        Scalar p = AG + CT + TV;
                        Scalar jc_arg = std::max<Scalar>(0.02, 1.0 - 4.0 / 3.0 * p);
                        d = -0.75 * std::log(jc_arg);
                    }
                }

                d = std::max(d, static_cast<Scalar>(1e-6));
                D(s1, s2) = d;
                D(s2, s1) = d;
            }
        }

        return D;
    }

    // Infers a Neighbor-Joining tree from an alignment using Saitou-Nei algorithm
    static Tree infer(const Alignment& aln) {
        Matrix D = compute_tn93_distances(aln);
        return infer_from_distances(D, aln.taxon_names);
    }

    // Infers a Neighbor-Joining tree from a precomputed distance matrix and taxon names
    static Tree infer_from_distances(const Matrix& D, const std::vector<std::string>& taxon_names) {
        size_t N = taxon_names.size();
        if (N == 0) {
            throw std::runtime_error("Cannot build tree with 0 taxa");
        }

        Tree tree;

        if (N == 1) {
            TreeNode leaf;
            leaf.id = 0;
            leaf.name = taxon_names[0];
            leaf.is_leaf = true;
            leaf.branch_length = 0.0;
            tree.nodes.push_back(std::move(leaf));
            tree.root_id = 0;
            tree.leaf_name_to_id[taxon_names[0]] = 0;
            tree.compute_post_order();
            return tree;
        }

        if (N == 2) {
            TreeNode l0;
            l0.id = 0;
            l0.name = taxon_names[0];
            l0.is_leaf = true;
            l0.parent_id = 2;
            l0.branch_length = std::max(D(0, 1) * 0.5, static_cast<Scalar>(1e-6));

            TreeNode l1;
            l1.id = 1;
            l1.name = taxon_names[1];
            l1.is_leaf = true;
            l1.parent_id = 2;
            l1.branch_length = std::max(D(0, 1) * 0.5, static_cast<Scalar>(1e-6));

            TreeNode root;
            root.id = 2;
            root.name = "Node1";
            root.is_leaf = false;
            root.children = {0, 1};

            tree.nodes.push_back(std::move(l0));
            tree.nodes.push_back(std::move(l1));
            tree.nodes.push_back(std::move(root));
            tree.root_id = 2;
            tree.leaf_name_to_id[taxon_names[0]] = 0;
            tree.leaf_name_to_id[taxon_names[1]] = 1;
            tree.compute_post_order();
            return tree;
        }

        // Initialize leaf nodes
        tree.nodes.reserve(2 * N);
        for (size_t i = 0; i < N; ++i) {
            TreeNode leaf;
            leaf.id = static_cast<int32_t>(i);
            leaf.name = taxon_names[i];
            leaf.is_leaf = true;
            tree.nodes.push_back(std::move(leaf));
            tree.leaf_name_to_id[taxon_names[i]] = static_cast<int32_t>(i);
        }

        std::vector<int32_t> active(N);
        for (size_t i = 0; i < N; ++i) active[i] = static_cast<int32_t>(i);
        Matrix D_curr = D;
        Matrix V_curr = D;

        int node_counter = 1;

        // Reduce active clusters down to 3 using BioNJ variance-weighted updates (Gascuel 1997)
        while (active.size() > 3) {
            size_t m = active.size();
            Vector R = Vector::Zero(m);
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    if (i != j) R(i) += D_curr(i, j);
                }
            }

            Scalar min_q = std::numeric_limits<Scalar>::infinity();
            size_t best_i = 0;
            size_t best_j = 1;

            Scalar scale = static_cast<Scalar>(m - 2);
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = i + 1; j < m; ++j) {
                    Scalar q = scale * D_curr(i, j) - R(i) - R(j);
                    if (q < min_q) {
                        min_q = q;
                        best_i = i;
                        best_j = j;
                    }
                }
            }

            int32_t id_i = active[best_i];
            int32_t id_j = active[best_j];
            int32_t new_id = static_cast<int32_t>(tree.nodes.size());

            TreeNode internal_node;
            internal_node.id = new_id;
            internal_node.name = "Node" + std::to_string(node_counter++);
            internal_node.is_leaf = false;
            internal_node.children = {id_i, id_j};

            Scalar dist_ij = D_curr(best_i, best_j);
            Scalar v_i = 0.5 * dist_ij + (R(best_i) - R(best_j)) / (2.0 * scale);
            Scalar v_j = dist_ij - v_i;

            tree.nodes[id_i].parent_id = new_id;
            tree.nodes[id_i].branch_length = std::max(v_i, static_cast<Scalar>(1e-6));

            tree.nodes[id_j].parent_id = new_id;
            tree.nodes[id_j].branch_length = std::max(v_j, static_cast<Scalar>(1e-6));

            tree.nodes.push_back(std::move(internal_node));

            // Construct new distance and variance matrices of size m - 1
            std::vector<int32_t> next_active;
            next_active.reserve(m - 1);
            std::vector<size_t> old_indices;
            old_indices.reserve(m - 2);

            for (size_t k = 0; k < m; ++k) {
                if (k != best_i && k != best_j) {
                    next_active.push_back(active[k]);
                    old_indices.push_back(k);
                }
            }
            next_active.push_back(new_id);

            Scalar var_ij = V_curr(best_i, best_j);
            Scalar lam = 0.5;
            if (var_ij > 1e-12 && scale > 0.0) {
                Scalar sum_v_diff = 0.0;
                for (size_t k : old_indices) {
                    sum_v_diff += (V_curr(best_j, k) - V_curr(best_i, k));
                }
                lam = std::clamp<Scalar>(0.5 + sum_v_diff / (2.0 * scale * var_ij), 0.0, 1.0);
            }

            size_t new_m = m - 1;
            Matrix D_next = Matrix::Zero(new_m, new_m);
            Matrix V_next = Matrix::Zero(new_m, new_m);
            for (size_t a = 0; a < old_indices.size(); ++a) {
                size_t ka = old_indices[a];
                for (size_t b = a + 1; b < old_indices.size(); ++b) {
                    size_t kb = old_indices[b];
                    D_next(a, b) = D_curr(ka, kb);
                    D_next(b, a) = D_curr(ka, kb);
                    V_next(a, b) = V_curr(ka, kb);
                    V_next(b, a) = V_curr(ka, kb);
                }
                Scalar d_uk = lam * (D_curr(best_i, ka) - v_i) + (1.0 - lam) * (D_curr(best_j, ka) - v_j);
                d_uk = std::max(d_uk, static_cast<Scalar>(1e-6));
                Scalar v_uk = std::max<Scalar>(
                    1e-12,
                    lam * V_curr(best_i, ka) + (1.0 - lam) * V_curr(best_j, ka) - lam * (1.0 - lam) * var_ij
                );
                D_next(a, new_m - 1) = d_uk;
                D_next(new_m - 1, a) = d_uk;
                V_next(a, new_m - 1) = v_uk;
                V_next(new_m - 1, a) = v_uk;
            }

            active = std::move(next_active);
            D_curr = std::move(D_next);
            V_curr = std::move(V_next);
        }

        // Final 3 nodes connect to root (trifurcating unrooted tree)
        int32_t root_id = static_cast<int32_t>(tree.nodes.size());
        TreeNode root_node;
        root_node.id = root_id;
        root_node.name = "Node" + std::to_string(node_counter++);
        root_node.is_leaf = false;
        root_node.children = {active[0], active[1], active[2]};

        Scalar d01 = D_curr(0, 1);
        Scalar d02 = D_curr(0, 2);
        Scalar d12 = D_curr(1, 2);

        Scalar v0 = 0.5 * (d01 + d02 - d12);
        Scalar v1 = 0.5 * (d01 + d12 - d02);
        Scalar v2 = 0.5 * (d02 + d12 - d01);

        tree.nodes[active[0]].parent_id = root_id;
        tree.nodes[active[0]].branch_length = std::max(v0, static_cast<Scalar>(1e-6));

        tree.nodes[active[1]].parent_id = root_id;
        tree.nodes[active[1]].branch_length = std::max(v1, static_cast<Scalar>(1e-6));

        tree.nodes[active[2]].parent_id = root_id;
        tree.nodes[active[2]].branch_length = std::max(v2, static_cast<Scalar>(1e-6));

        tree.nodes.push_back(std::move(root_node));
        tree.root_id = root_id;

        tree.compute_post_order();
        return tree;
    }
};

} // namespace hyphy::core
