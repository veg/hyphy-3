#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/neighbor_joining.hpp"
#include "hyphy/core/ohmic_network.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "LBFGSB.h"
#include <nlohmann/json.hpp>
#include <Eigen/Eigenvalues>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <iomanip>
#include <chrono>

namespace hyphy::analyses {

using namespace hyphy::core;

enum class TreeSearchModel {
    GTR,   // Fast nucleotide GTR
    MG94   // Codon MG94 3x4
};

struct TreeEngineSettings {
    TreeSearchModel model = TreeSearchModel::GTR;
    bool allow_polytomies = true;
    Scalar polytomy_collapse_threshold = 1e-4; // Branch lengths <= threshold collapsed into polytomies
    Scalar spectral_star_threshold = 0.05;      // If lambda_2 <= threshold, consider cluster star-like
    Scalar spectral_gap_ratio = 0.15;           // If lambda_2 / lambda_3 <= ratio, don't force binary split
    Scalar topological_grad_threshold = 2.0;    // Minimum positive gradient to resolve a polytomy
    bool perform_topological_gradients = true;  // Test resolving polytomies via directional derivatives
    bool perform_spr = true;                   // Perform Ohmic resistance-bounded SPR
    size_t spr_max_rounds = 3;                 // Maximum SPR improvement rounds
    size_t spr_max_candidates = 8;             // Nearest electrical neighbors to evaluate per pruning
    bool perform_model_escalation = false;      // Escalate GTR -> MG94 on codon alignments (default: false, Option C)
    size_t lbfgs_max_iterations = 40;           // Max L-BFGS iterations per branch length sweep
    Scalar lbfgs_epsilon = 1e-4;                // L-BFGS gradient convergence tolerance

    // Ohmic Divide-and-Conquer settings (enabled by default)
    bool enable_dnc = true;                    // Enable Ohmic Divide-and-Conquer by default
    size_t dnc_min_taxa = 16;                  // Minimum taxa to trigger D&C partitioning (default: 16)
    size_t dnc_target_island_size = 14;        // Target maximum taxa per Ohmic island (default: 14)
    size_t dnc_min_island_size = 4;            // Minimum taxa per Ohmic island (default: 4)

    bool quiet = false;
};

struct TreeEngineResult {
    Tree initial_spectral_tree;
    Tree optimized_tree;
    Scalar initial_gtr_log_l = 0.0;
    Scalar optimized_gtr_log_l = 0.0;
    Scalar final_log_l = 0.0;
    std::string active_model = "GTR";
    size_t num_taxa = 0;
    size_t num_patterns = 0;
    size_t num_initial_nodes = 0;
    size_t num_final_nodes = 0;
    size_t num_polytomies_retained = 0;
    size_t num_branches_collapsed = 0;
    size_t num_polytomies_resolved = 0;
    size_t num_spr_moves_accepted = 0;
    double runtime_seconds = 0.0;

    // Ohmic Divide-and-Conquer statistics
    bool dnc_used = false;
    size_t num_islands = 0;
    std::vector<size_t> island_sizes;
    size_t dnc_spr_moves = 0;

    // GTR parameter estimates
    Vector4 gtr_frequencies = Vector4::Constant(0.25);
    std::array<Scalar, 5> gtr_rates = {1.0, 1.0, 1.0, 1.0, 1.0}; // AC, AT, CG, CT, GT (AG=1)

    // MG94 parameter estimates (if escalated)
    Scalar mg94_omega = 1.0;
};

class TreeEngine {
public:
    TreeEngineSettings settings;

    explicit TreeEngine(TreeEngineSettings s = TreeEngineSettings{}) : settings(s) {}

    // =========================================================================
    // 1. Spectral Tree Builder: Native Polytomous Scaffold from Ohmic Conductance
    // =========================================================================
    static Tree build_spectral_scaffold(
        const Alignment& aln,
        bool allow_polytomies = true,
        Scalar star_threshold = 0.05,
        Scalar gap_ratio = 0.15
    ) {
        size_t N = aln.num_taxa;
        if (N < 2) {
            throw std::runtime_error("Spectral tree scaffolding requires at least 2 taxa.");
        }
        Matrix D = NeighborJoining::compute_tn93_distances(aln);
        return build_spectral_scaffold_from_distances(D, aln.taxon_names, allow_polytomies, star_threshold, gap_ratio);
    }

    static Tree build_spectral_scaffold_from_distances(
        const Matrix& D,
        const std::vector<std::string>& taxon_names,
        bool allow_polytomies = true,
        Scalar star_threshold = 0.05,
        Scalar gap_ratio = 0.15
    ) {
        (void)gap_ratio;
        size_t N = taxon_names.size();
        if (N < 2) {
            throw std::runtime_error("Spectral tree scaffolding requires at least 2 taxa.");
        }

        Tree tree;
        int32_t node_counter = 0;

        auto create_leaf = [&](size_t taxon_idx) -> int32_t {
            int32_t id = static_cast<int32_t>(tree.nodes.size());
            TreeNode node;
            node.id = id;
            node.name = taxon_names[taxon_idx];
            node.is_leaf = true;
            node.branch_length = 0.01;
            tree.nodes.push_back(std::move(node));
            tree.leaf_name_to_id[taxon_names[taxon_idx]] = id;
            return id;
        };

        auto create_internal = [&]() -> int32_t {
            int32_t id = static_cast<int32_t>(tree.nodes.size());
            TreeNode node;
            node.id = id;
            node.name = "Node" + std::to_string(++node_counter);
            node.is_leaf = false;
            node.branch_length = 0.01;
            tree.nodes.push_back(std::move(node));
            return id;
        };

        // Recursive Grounded Dirichlet Green's Function (Transfer Impedance) Spectral Partitioning:
        // Grounds each sub-circuit U at its external boundary port A (or internal bulk reference potential r_i)
        // to form the transfer-impedance matrix G_{ij}^{(A)} = 0.5 * (r_i + r_j - D_{ij}), which cancels
        // terminal branch-length rate heterogeneity identically and isolates shared ancestral stem paths.
        std::function<int32_t(const std::vector<size_t>&, const std::vector<size_t>&)> recurse_partition;
        recurse_partition = [&](const std::vector<size_t>& indices, const std::vector<size_t>& anchor_indices) -> int32_t {
            size_t m = indices.size();
            if (m == 1) {
                return create_leaf(indices[0]);
            }
            if (m == 2) {
                int32_t u = create_internal();
                int32_t c1 = create_leaf(indices[0]);
                int32_t c2 = create_leaf(indices[1]);
                tree.nodes[u].children = {c1, c2};
                tree.nodes[c1].parent_id = u;
                tree.nodes[c2].parent_id = u;
                Scalar d = D(static_cast<Eigen::Index>(indices[0]), static_cast<Eigen::Index>(indices[1]));
                tree.nodes[c1].branch_length = std::max<Scalar>(1e-4, d * 0.5);
                tree.nodes[c2].branch_length = std::max<Scalar>(1e-4, d * 0.5);
                return u;
            }

            Matrix D_loc = Matrix::Zero(m, m);
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = i + 1; j < m; ++j) {
                    Scalar d = D(static_cast<Eigen::Index>(indices[i]), static_cast<Eigen::Index>(indices[j]));
                    D_loc(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = d;
                    D_loc(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = d;
                }
            }

            // Internal bulk reference potential r_int(i) and proximal boundary-port potential r_anc(i)
            Vector r_int = Vector::Zero(m);
            for (size_t i = 0; i < m; ++i) {
                r_int(static_cast<Eigen::Index>(i)) =
                    D_loc.row(static_cast<Eigen::Index>(i)).sum() / static_cast<Scalar>(std::max<size_t>(1, m - 2));
            }
            Vector r_anc = r_int;
            if (!anchor_indices.empty()) {
                std::vector<std::pair<Scalar, size_t>> anc_dists;
                anc_dists.reserve(anchor_indices.size());
                for (size_t a : anchor_indices) {
                    Scalar min_d = 1e18;
                    for (size_t i = 0; i < m; ++i) {
                        min_d = std::min(min_d, D(static_cast<Eigen::Index>(indices[i]), static_cast<Eigen::Index>(a)));
                    }
                    anc_dists.push_back({min_d, a});
                }
                std::sort(anc_dists.begin(), anc_dists.end());
                size_t k_port = std::min<size_t>(4, anc_dists.size());
                for (size_t i = 0; i < m; ++i) {
                    Scalar s = 0.0;
                    for (size_t k = 0; k < k_port; ++k) {
                        s += D(static_cast<Eigen::Index>(indices[i]), static_cast<Eigen::Index>(anc_dists[k].second));
                    }
                    r_anc(static_cast<Eigen::Index>(i)) = s / static_cast<Scalar>(k_port);
                }
            }

            // Grounded Dirichlet Green's Matrices (Transfer Impedance):
            // G_anc(i, j) = 0.5 * (r_anc(i) + r_anc(j) - D_loc(i, j))
            // G_int(i, j) = 0.5 * (r_int(i) + r_int(j) - D_loc(i, j))
            Matrix G_anc = Matrix::Zero(m, m);
            Matrix G_int = Matrix::Zero(m, m);
            Scalar min_ga = 1e18, min_gi = 1e18;
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = i + 1; j < m; ++j) {
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    Scalar ga = 0.5 * (r_anc(ii) + r_anc(jj) - D_loc(ii, jj));
                    Scalar gi = 0.5 * (r_int(ii) + r_int(jj) - D_loc(ii, jj));
                    G_anc(ii, jj) = G_anc(jj, ii) = ga;
                    G_int(ii, jj) = G_int(jj, ii) = gi;
                    min_ga = std::min(min_ga, ga);
                    min_gi = std::min(min_gi, gi);
                }
            }
            Matrix G = Matrix::Zero(m, m);
            std::vector<Scalar> pos_g;
            pos_g.reserve(m * (m - 1) / 2);
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = i + 1; j < m; ++j) {
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    G_anc(ii, jj) = G_anc(jj, ii) = std::max<Scalar>(0.0, G_anc(ii, jj) - min_ga);
                    G_int(ii, jj) = G_int(jj, ii) = std::max<Scalar>(0.0, G_int(ii, jj) - min_gi);
                    G(ii, jj) = G(jj, ii) = anchor_indices.empty()
                        ? G_int(ii, jj)
                        : (0.60 * G_anc(ii, jj) + 0.40 * G_int(ii, jj));
                    if (G(ii, jj) > 1e-6) pos_g.push_back(G(ii, jj));
                }
            }

            if (m == 3) {
                // For m == 3, the pair (i, j) with the largest grounded transfer impedance shares the internal stem
                Scalar best_g = -1e18;
                size_t best_solo = 0;
                for (size_t solo = 0; solo < 3; ++solo) {
                    Eigen::Index ii = static_cast<Eigen::Index>((solo + 1) % 3);
                    Eigen::Index jj = static_cast<Eigen::Index>((solo + 2) % 3);
                    Scalar g_val = anchor_indices.empty() ? G(ii, jj) : G_anc(ii, jj);
                    if (g_val > best_g) {
                        best_g = g_val;
                        best_solo = solo;
                    }
                }
                if (allow_polytomies && best_g < star_threshold * 0.05) {
                    int32_t poly_node = create_internal();
                    for (size_t idx : indices) {
                        int32_t c = create_leaf(idx);
                        tree.nodes[poly_node].children.push_back(c);
                        tree.nodes[c].parent_id = poly_node;
                        Scalar avg_d = 0.0;
                        for (size_t o : indices) {
                            if (o != idx) avg_d += D(static_cast<Eigen::Index>(idx), static_cast<Eigen::Index>(o));
                        }
                        avg_d /= static_cast<Scalar>(m - 1);
                        tree.nodes[c].branch_length = std::max<Scalar>(1e-4, avg_d * 0.5);
                    }
                    return poly_node;
                }
                std::vector<size_t> grp1 = {indices[best_solo]};
                std::vector<size_t> grp2 = {indices[(best_solo + 1) % 3], indices[(best_solo + 2) % 3]};
                std::vector<size_t> anc1 = grp2;
                anc1.insert(anc1.end(), anchor_indices.begin(), anchor_indices.end());
                std::vector<size_t> anc2 = grp1;
                anc2.insert(anc2.end(), anchor_indices.begin(), anchor_indices.end());
                int32_t parent = create_internal();
                int32_t child1 = recurse_partition(grp1, anc1);
                int32_t child2 = recurse_partition(grp2, anc2);
                tree.nodes[parent].children = {child1, child2};
                tree.nodes[child1].parent_id = parent;
                tree.nodes[child2].parent_id = parent;
                Scalar d_solo = 0.5 * (D_loc(best_solo, (best_solo + 1) % 3) + D_loc(best_solo, (best_solo + 2) % 3));
                tree.nodes[child1].branch_length = std::max<Scalar>(1e-4, d_solo * 0.5);
                tree.nodes[child2].branch_length = std::max<Scalar>(1e-4, best_g);
                return parent;
            }

            Scalar tau_g = 0.05;
            Scalar max_g = 0.0;
            if (!pos_g.empty()) {
                std::sort(pos_g.begin(), pos_g.end());
                tau_g = std::max<Scalar>(1e-4, pos_g[pos_g.size() / 2]);
                max_g = pos_g.back();
            }

            // Build Grounded Ohmic Coupling Matrix W from shared ancestral stem impedance G(i, j)
            Matrix W = Matrix::Zero(m, m);
            Vector deg = Vector::Zero(m);
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = i + 1; j < m; ++j) {
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    Scalar w = G(ii, jj) * std::exp(std::min<Scalar>(5.0, G(ii, jj) / tau_g)) + 1e-6;
                    W(ii, jj) = W(jj, ii) = w;
                    deg(ii) += w;
                    deg(jj) += w;
                }
            }

            Vector inv_sqrt_deg(m);
            for (size_t i = 0; i < m; ++i) {
                inv_sqrt_deg(static_cast<Eigen::Index>(i)) =
                    (deg(static_cast<Eigen::Index>(i)) > 1e-12)
                        ? (1.0 / std::sqrt(deg(static_cast<Eigen::Index>(i))))
                        : 1.0;
            }
            Matrix L_norm = Matrix::Identity(m, m) - inv_sqrt_deg.asDiagonal() * W * inv_sqrt_deg.asDiagonal();
            L_norm = 0.5 * (L_norm + L_norm.transpose());
            Eigen::SelfAdjointEigenSolver<Matrix> es_nm(L_norm);
            if (es_nm.info() != Eigen::Success) {
                int32_t u = create_internal();
                for (size_t idx : indices) {
                    int32_t c = create_leaf(idx);
                    tree.nodes[u].children.push_back(c);
                    tree.nodes[c].parent_id = u;
                }
                return u;
            }

            // In a Grounded Green's transfer impedance graph G(i, j), a star radiation has near-zero
            // shared internal stem impedance across all pairs (max_g < 0.05 * star_threshold)
            bool is_star = allow_polytomies && (m <= 8) && (max_g < star_threshold * 0.05);
            if (is_star) {
                int32_t poly_node = create_internal();
                for (size_t idx : indices) {
                    int32_t c = create_leaf(idx);
                    tree.nodes[poly_node].children.push_back(c);
                    tree.nodes[c].parent_id = poly_node;
                    Scalar avg_d = 0.0;
                    for (size_t o : indices) {
                        if (o != idx) avg_d += D(static_cast<Eigen::Index>(idx), static_cast<Eigen::Index>(o));
                    }
                    avg_d /= static_cast<Scalar>(m - 1);
                    tree.nodes[c].branch_length = std::max<Scalar>(1e-4, avg_d * 0.5);
                }
                return poly_node;
            }

            Matrix L_unnorm = Matrix(deg.asDiagonal()) - W;
            L_unnorm = 0.5 * (L_unnorm + L_unnorm.transpose());
            Eigen::SelfAdjointEigenSolver<Matrix> es_un(L_unnorm);
            Eigen::SelfAdjointEigenSolver<Matrix> es_G(G);

            std::vector<std::vector<size_t>> orders;
            auto add_vec_order = [&](const Vector& v) {
                std::vector<size_t> ord(m);
                std::iota(ord.begin(), ord.end(), 0);
                std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
                    return v(static_cast<Eigen::Index>(a)) < v(static_cast<Eigen::Index>(b));
                });
                orders.push_back(std::move(ord));
            };
            add_vec_order(deg);
            add_vec_order(G_anc.rowwise().sum());
            if (es_G.info() == Eigen::Success) {
                add_vec_order(es_G.eigenvectors().col(static_cast<Eigen::Index>(m - 1)));
                add_vec_order(es_G.eigenvectors().col(static_cast<Eigen::Index>(m - 2)));
            }
            if (es_un.info() == Eigen::Success) {
                add_vec_order(es_un.eigenvectors().col(1));
                if (m >= 4) add_vec_order(es_un.eigenvectors().col(2));
            }
            add_vec_order(inv_sqrt_deg.asDiagonal() * es_nm.eigenvectors().col(1));
            if (m >= 4) {
                add_vec_order(inv_sqrt_deg.asDiagonal() * es_nm.eigenvectors().col(2));
            }

            // Evaluate candidate bipartitions (L_idx, R_idx) along the spectral modes using
            // Grounded Transfer-Impedance Monophyly Margin & Cross-Cut Stem Drop:
            Scalar best_score = -1e18;
            std::vector<size_t> grp1, grp2;

            for (const auto& ord : orders) {
                for (size_t k = 1; k < m; ++k) {
                    if (anchor_indices.empty() && (k == 1 || k == m - 1) && m >= 4) {
                        continue;
                    }
                    std::vector<size_t> L_idx(ord.begin(), ord.begin() + k);
                    std::vector<size_t> R_idx(ord.begin() + k, ord.end());

                    Scalar sum_cross = 0.0, max_cross = -1e18;
                    for (size_t l : L_idx) {
                        for (size_t r : R_idx) {
                            Scalar g = G_anc(static_cast<Eigen::Index>(l), static_cast<Eigen::Index>(r));
                            sum_cross += g;
                            if (g > max_cross) max_cross = g;
                        }
                    }
                    Scalar mean_cross = sum_cross / static_cast<Scalar>(L_idx.size() * R_idx.size());

                    Scalar min_node_margin = 1e18;
                    Scalar sum_node_margin = 0.0;
                    if (L_idx.size() >= 2) {
                        for (size_t l : L_idx) {
                            Scalar max_in = -1e18, max_out = -1e18;
                            for (size_t l2 : L_idx) {
                                if (l2 != l) max_in = std::max(max_in, G(static_cast<Eigen::Index>(l), static_cast<Eigen::Index>(l2)));
                            }
                            for (size_t r : R_idx) {
                                max_out = std::max(max_out, G(static_cast<Eigen::Index>(l), static_cast<Eigen::Index>(r)));
                            }
                            Scalar marg = max_in - max_out;
                            min_node_margin = std::min(min_node_margin, marg);
                            sum_node_margin += marg;
                        }
                    }
                    if (R_idx.size() >= 2) {
                        for (size_t r : R_idx) {
                            Scalar max_in = -1e18, max_out = -1e18;
                            for (size_t r2 : R_idx) {
                                if (r2 != r) max_in = std::max(max_in, G(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(r2)));
                            }
                            for (size_t l : L_idx) {
                                max_out = std::max(max_out, G(static_cast<Eigen::Index>(l), static_cast<Eigen::Index>(r)));
                            }
                            Scalar marg = max_in - max_out;
                            min_node_margin = std::min(min_node_margin, marg);
                            sum_node_margin += marg;
                        }
                    }
                    Scalar mean_node_margin = sum_node_margin / static_cast<Scalar>(m);

                    Scalar sum_in = 0.0, cnt_in = 0.0;
                    for (size_t i = 0; i < L_idx.size(); ++i) {
                        for (size_t j = i + 1; j < L_idx.size(); ++j) {
                            sum_in += G_anc(static_cast<Eigen::Index>(L_idx[i]), static_cast<Eigen::Index>(L_idx[j]));
                            cnt_in += 1.0;
                        }
                    }
                    for (size_t i = 0; i < R_idx.size(); ++i) {
                        for (size_t j = i + 1; j < R_idx.size(); ++j) {
                            sum_in += G_anc(static_cast<Eigen::Index>(R_idx[i]), static_cast<Eigen::Index>(R_idx[j]));
                            cnt_in += 1.0;
                        }
                    }
                    Scalar mean_in = (cnt_in > 0.0) ? (sum_in / cnt_in) : 0.0;
                    Scalar stem_gap = mean_in - mean_cross;

                    Scalar score = 0.50 * min_node_margin + 0.30 * mean_node_margin + 0.20 * (stem_gap - 0.5 * max_cross);
                    if (score > best_score) {
                        best_score = score;
                        grp1.clear();
                        grp2.clear();
                        for (size_t idx : L_idx) grp1.push_back(indices[idx]);
                        for (size_t idx : R_idx) grp2.push_back(indices[idx]);
                    }
                }
            }

            if (grp1.empty() || grp2.empty()) {
                grp1.clear();
                grp2.clear();
                for (size_t i = 0; i < m; ++i) {
                    if (i < m / 2) grp1.push_back(indices[i]);
                    else grp2.push_back(indices[i]);
                }
            }

            std::vector<size_t> anc1 = grp2;
            anc1.insert(anc1.end(), anchor_indices.begin(), anchor_indices.end());
            std::vector<size_t> anc2 = grp1;
            anc2.insert(anc2.end(), anchor_indices.begin(), anchor_indices.end());

            int32_t parent = create_internal();
            int32_t child1 = recurse_partition(grp1, anc1);
            int32_t child2 = recurse_partition(grp2, anc2);

            tree.nodes[parent].children = {child1, child2};
            tree.nodes[child1].parent_id = parent;
            tree.nodes[child2].parent_id = parent;

            Scalar cross_dist = 0.0;
            for (size_t i : grp1) {
                for (size_t j : grp2) {
                    cross_dist += D(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                }
            }
            cross_dist /= static_cast<Scalar>(grp1.size() * grp2.size());
            tree.nodes[child1].branch_length = std::max<Scalar>(1e-4, cross_dist * 0.25);
            tree.nodes[child2].branch_length = std::max<Scalar>(1e-4, cross_dist * 0.25);

            return parent;
        };

        std::vector<size_t> all_taxa(N);
        std::iota(all_taxa.begin(), all_taxa.end(), 0);
        std::vector<size_t> empty_anchor;

        int32_t root = recurse_partition(all_taxa, empty_anchor);
        tree.root_id = root;
        tree.compute_post_order();
        return tree;
    }

    // =========================================================================
    // 2. Simultaneous All-Branch Quasi-Newton (L-BFGS) Optimization under GTR
    // =========================================================================
    static Scalar optimize_branches_gtr(
        Tree& tree,
        const Alignment& aln,
        const GTRMatrix& gtr,
        size_t max_iterations = 40,
        Scalar epsilon = 1e-4,
        const std::vector<int32_t>* target_branches = nullptr
    ) {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

        // Identify branch indices to optimize
        std::vector<int32_t> branches;
        if (target_branches != nullptr) {
            branches = *target_branches;
        } else {
            branches.reserve(num_nodes);
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                    branches.push_back(node.id);
                }
            }
        }
        size_t n_b = branches.size();
        if (n_b == 0) return 0.0;

        // Bounded L-BFGS Objective over log-branch lengths y_i = ln(t_i)
        struct GTRBranchObjective {
            const Tree& tree_ref;
            const Alignment& aln_ref;
            const GTRMatrix& gtr_ref;
            const std::vector<int32_t>& brs;
            const std::vector<size_t>& l2t;
            size_t n_nodes;
            size_t n_pat;
            std::vector<Matrix4> P_b_base;
            std::vector<Scalar> bl_base;

            GTRBranchObjective(
                const Tree& t, const Alignment& a, const GTRMatrix& g,
                const std::vector<int32_t>& b, const std::vector<size_t>& l,
                size_t nn, size_t np
            ) : tree_ref(t), aln_ref(a), gtr_ref(g), brs(b), l2t(l), n_nodes(nn), n_pat(np),
                P_b_base(nn), bl_base(nn, 0.0)
            {
                for (const auto& node : tree_ref.nodes) {
                    if (node.id != tree_ref.root_id && node.parent_id != INVALID_INDEX) {
                        int32_t bid = node.id;
                        bl_base[bid] = std::max<Scalar>(1e-6, node.branch_length);
                        P_b_base[bid] = gtr_ref.transition_matrix(bl_base[bid]);
                    }
                }
            }

            Scalar operator()(const Vector& y, Vector& grad) {
                std::vector<Matrix4> P_b = P_b_base;
                std::vector<Matrix4> dP_b(n_nodes);
                std::vector<Scalar> bl = bl_base;

                for (size_t i = 0; i < brs.size(); ++i) {
                    int32_t bid = brs[i];
                    Scalar t = std::exp(y(static_cast<Eigen::Index>(i)));
                    bl[bid] = t;
                    P_b[bid] = gtr_ref.transition_matrix(t);
                    dP_b[bid] = gtr_ref.Q * P_b[bid];
                }

                Scalar total_ll = 0.0;
                std::vector<Scalar> dL_dt(n_nodes, 0.0);
                std::vector<Scalar> pat_ll(n_pat, 0.0);
                std::vector<std::vector<Scalar>> pat_dL(n_nodes, std::vector<Scalar>(n_pat, 0.0));

                #ifdef _OPENMP
                #pragma omp parallel if(!omp_in_parallel())
                #endif
                {
                    std::vector<Vector4> node_L(n_nodes);
                    std::vector<double> node_log_scale(n_nodes, 0.0);
                    std::vector<Vector4> P_D(n_nodes);
                    std::vector<int> boost_counts(n_nodes, 0);
                    std::vector<Vector4> U(n_nodes, Vector4::Zero());
                    std::vector<Vector4> V(n_nodes, Vector4::Zero());

                    #pragma omp for schedule(static)
                    for (size_t p = 0; p < n_pat; ++p) {
                        const auto& pattern = aln_ref.patterns[p];

                        // 1. Post-order traversal (Inside pass)
                        for (int32_t nid : tree_ref.post_order) {
                            const auto& node = tree_ref.nodes[nid];
                            if (node.is_leaf) {
                                size_t t_idx = l2t[nid];
                                if (t_idx != static_cast<size_t>(-1)) {
                                    int8_t s = pattern.states[t_idx];
                                    if (s >= 0 && s < 4) {
                                        node_L[nid].setZero();
                                        node_L[nid](s) = 1.0;
                                    } else {
                                        node_L[nid].setOnes();
                                    }
                                } else {
                                    node_L[nid].setOnes();
                                }
                                node_log_scale[nid] = 0.0;
                                boost_counts[nid] = 0;
                            } else {
                                double child_scale_sum = 0.0;
                                for (int32_t cid : node.children) {
                                    P_D[cid] = P_b[cid] * node_L[cid];
                                }
                                int boost_cnt = 0;
                                const std::vector<std::array<Matrix4, 3>>* no_bridges = nullptr;
                                tree_ref.eval_internal_node_with_model(nid, P_b, gtr_ref, node_L, node_log_scale, node_L[nid], child_scale_sum, 1.0, nullptr, no_bridges, &boost_cnt);
                                boost_counts[nid] = boost_cnt;
                                node_log_scale[nid] = child_scale_sum;
                            }
                        }

                        Scalar pat_l = gtr_ref.pi.dot(node_L[tree_ref.root_id]);
                        if (pat_l > 0.0) {
                            pat_ll[p] = pattern.weight * (std::log(pat_l) + node_log_scale[tree_ref.root_id]);
                            Scalar inv_L = pattern.weight / pat_l;

                            // 2. Pre-order traversal (Outside pass)
                            U[tree_ref.root_id] = gtr_ref.pi;
                            for (auto it = tree_ref.post_order.rbegin(); it != tree_ref.post_order.rend(); ++it) {
                                int32_t uid = *it;
                                const auto& node = tree_ref.nodes[uid];
                                if (node.is_leaf) continue;

                                const double SCALER_THRESH = 1e-150;
                                const double SCALER_UP = 1e150;
                                const double SCALER_DOWN = 1e-150;
                                size_t num_children = node.children.size();
                                for (size_t i = 0; i < num_children; ++i) {
                                    int32_t vid = node.children[i];
                                    Vector4 v_msg = U[uid];
                                    int sib_boost = 0;
                                    for (size_t j = 0; j < num_children; ++j) {
                                        if (i != j) {
                                            v_msg.array() *= P_D[node.children[j]].array();
                                            double cur_max = v_msg.maxCoeff();
                                            while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                                                v_msg *= SCALER_UP;
                                                cur_max *= SCALER_UP;
                                                sib_boost++;
                                            }
                                        }
                                    }
                                    int net_boost = boost_counts[uid] - sib_boost;
                                    while (net_boost > 0) {
                                        v_msg *= SCALER_UP;
                                        net_boost--;
                                    }
                                    while (net_boost < 0) {
                                        v_msg *= SCALER_DOWN;
                                        net_boost++;
                                    }
                                    V[vid] = v_msg;
                                    U[vid] = P_b[vid].transpose() * v_msg;

                                    Scalar dL = v_msg.dot(dP_b[vid] * node_L[vid]);
                                    pat_dL[vid][p] = inv_L * dL;
                                }
                            }
                        }
                    }
                }

                total_ll = LikelihoodEngine::neumaier_sum(pat_ll);
                for (int32_t bid : brs) {
                    dL_dt[bid] = LikelihoodEngine::neumaier_sum(pat_dL[bid]);
                }

                grad.resize(static_cast<Eigen::Index>(brs.size()));
                for (size_t i = 0; i < brs.size(); ++i) {
                    int32_t bid = brs[i];
                    Scalar t = bl[bid];
                    // Chain rule: d(-ln L) / dy = -t * (d ln L / dt)
                    grad(static_cast<Eigen::Index>(i)) = -t * dL_dt[bid];
                }

                return -total_ll;
            }
        };

        GTRBranchObjective obj{tree, aln, gtr, branches, leaf_to_taxon, num_nodes, num_patterns};

        Vector y(n_b), lb(n_b), ub(n_b);
        for (size_t i = 0; i < n_b; ++i) {
            Scalar t = std::clamp(tree.nodes[branches[i]].branch_length, 1e-6, 50.0);
            y(static_cast<Eigen::Index>(i)) = std::log(t);
            lb(static_cast<Eigen::Index>(i)) = -14.0; // t >= ~8e-7
            ub(static_cast<Eigen::Index>(i)) =   4.5; // t <= ~90
        }

        LBFGSpp::LBFGSBParam<Scalar> opt_param;
        opt_param.m = std::max<int>(3, std::min<int>(6, static_cast<int>(n_b)));
        opt_param.epsilon = epsilon;
        opt_param.max_iterations = static_cast<int>(max_iterations);

        LBFGSpp::LBFGSBSolver<Scalar> solver(opt_param);
        Scalar fx = 0.0;
        try {
            solver.minimize(obj, y, fx, lb, ub);
        } catch (...) {}

        // Store back optimized branch lengths
        for (size_t i = 0; i < n_b; ++i) {
            tree.nodes[branches[i]].branch_length = std::exp(y(static_cast<Eigen::Index>(i)));
        }

        return -fx;
    }

    // =========================================================================
    // 2b. Optimize 5 GTR Exchangeability Rates (AC, AT, CG, CT, GT)
    // =========================================================================
    static Scalar optimize_gtr_rates(
        const Tree& tree,
        const Alignment& aln,
        GTRMatrix& gtr_model,
        GTRParameters& gtr_params
    ) {
        auto opt_rate = [&](Scalar& r_ref) {
            auto rate_obj = [&](Scalar val) -> Scalar {
                Scalar old_val = r_ref;
                r_ref = val;
                GTRMatrix temp_mat;
                temp_mat.update(gtr_params, aln.nuc_frequencies);
                Scalar lnl = LikelihoodEngine::compute_gtr_log_likelihood(tree, aln, temp_mat);
                r_ref = old_val;
                return -lnl;
            };
            auto [best_r, _] = opt::Brent1D::minimize(rate_obj, 1e-4, r_ref, 30.0, 1e-3, 20);
            r_ref = best_r;
            gtr_model.update(gtr_params, aln.nuc_frequencies);
        };

        for (int pass = 0; pass < 2; ++pass) {
            opt_rate(gtr_params.theta_AC);
            opt_rate(gtr_params.theta_AT);
            opt_rate(gtr_params.theta_CG);
            opt_rate(gtr_params.theta_CT);
            opt_rate(gtr_params.theta_GT);
        }

        return LikelihoodEngine::compute_gtr_log_likelihood(tree, aln, gtr_model);
    }

    // =========================================================================
    // 3. Dynamic Polytomy Collapse: Collapses Zero-Signal Internal Branches
    // =========================================================================
    static size_t collapse_short_branches(Tree& tree, Scalar threshold = 1e-4) {
        size_t collapsed_count = 0;
        bool changed = true;

        while (changed) {
            changed = false;
            for (auto& node : tree.nodes) {
                // Only collapse internal non-root nodes whose branch length falls below threshold
                if (!node.is_leaf && node.id != tree.root_id && node.branch_length <= threshold) {
                    int32_t u = node.id;
                    int32_t p = node.parent_id;
                    if (p == INVALID_INDEX) continue;

                    // Remove u from parent's children list
                    auto& p_children = tree.nodes[p].children;
                    p_children.erase(std::remove(p_children.begin(), p_children.end(), u), p_children.end());

                    // Adopt all of u's children under p
                    for (int32_t cid : node.children) {
                        tree.nodes[cid].parent_id = p;
                        p_children.push_back(cid);
                    }

                    // Invalidate u
                    node.children.clear();
                    node.parent_id = INVALID_INDEX;
                    node.branch_length = 0.0;

                    collapsed_count++;
                    changed = true;
                    break; // Restart loop after structural mutation
                }
            }
        }

        if (collapsed_count > 0) {
            tree.compute_post_order();
        }
        return collapsed_count;
    }

    // =========================================================================
    // 4. Directional Topological Gradient at Polytomies: "Pay Your Rent" Testing
    // =========================================================================
    static size_t resolve_polytomies_by_gradient(
        Tree& tree,
        const Alignment& aln,
        const GTRMatrix& gtr,
        Scalar gradient_threshold = 2.0
    ) {
        size_t resolved_count = 0;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

        // Precompute P and dP for current branch lengths
        std::vector<Matrix4> P_b(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                P_b[node.id] = gtr.transition_matrix(node.branch_length);
            }
        }

        // Find internal nodes with >= 3 children (polytomies)
        std::vector<int32_t> polytomies;
        for (const auto& node : tree.nodes) {
            if (!node.is_leaf && node.children.size() >= 3) {
                polytomies.push_back(node.id);
            }
        }

        if (polytomies.empty()) return 0;

        // Inside-Outside pass to get partials at all nodes
        std::vector<Vector4> D_all(num_nodes);
        std::vector<Vector4> V_all(num_nodes);
        std::vector<double> node_log_scale(num_nodes, 0.0);
        std::vector<Vector4> P_D(num_nodes);
        std::vector<int> boost_counts(num_nodes, 0);

        for (int32_t poly_id : polytomies) {
            const std::vector<int32_t> p_children = tree.nodes[poly_id].children;
            size_t k = p_children.size();
            if (k < 3) continue;

            // Evaluate directional derivatives for all pairs of children (ci, cj)
            Scalar best_grad = -1e9;
            size_t best_i = 0, best_j = 1;

            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];

                // Inside pass for pattern
                for (int32_t nid : tree.post_order) {
                    const auto& node = tree.nodes[nid];
                    if (node.is_leaf) {
                        size_t t_idx = leaf_to_taxon[nid];
                        if (t_idx != static_cast<size_t>(-1)) {
                            int8_t s = pattern.states[t_idx];
                            if (s >= 0 && s < 4) {
                                D_all[nid].setZero();
                                D_all[nid](s) = 1.0;
                            } else {
                                D_all[nid].setOnes();
                            }
                        } else {
                            D_all[nid].setOnes();
                        }
                        node_log_scale[nid] = 0.0;
                        boost_counts[nid] = 0;
                    } else {
                        double child_scale_sum = 0.0;
                        for (int32_t cid : node.children) {
                            P_D[cid] = P_b[cid] * D_all[cid];
                        }
                        int boost_cnt = 0;
                        const std::vector<std::array<Matrix4, 3>>* no_bridges = nullptr;
                        tree.eval_internal_node_with_model(nid, P_b, gtr, D_all, node_log_scale, D_all[nid], child_scale_sum, 1.0, nullptr, no_bridges, &boost_cnt);
                        boost_counts[nid] = boost_cnt;
                        node_log_scale[nid] = child_scale_sum;
                    }
                }

                Scalar pat_l = gtr.pi.dot(D_all[tree.root_id]);
                if (pat_l <= 0.0) continue;
                Scalar inv_L = pattern.weight / pat_l;

                // Outside pass
                const double SCALER_THRESH = 1e-150;
                const double SCALER_UP = 1e150;
                const double SCALER_DOWN = 1e-150;
                std::vector<Vector4> U(num_nodes, Vector4::Zero());
                U[tree.root_id] = gtr.pi;
                for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
                    int32_t uid = *it;
                    const auto& node = tree.nodes[uid];
                    if (node.is_leaf) continue;

                    size_t num_children = node.children.size();
                    for (size_t i = 0; i < num_children; ++i) {
                        int32_t vid = node.children[i];
                        Vector4 v_msg = U[uid];
                        int sib_boost = 0;
                        for (size_t j = 0; j < num_children; ++j) {
                            if (i != j) {
                                v_msg.array() *= P_D[node.children[j]].array();
                                double cur_max = v_msg.maxCoeff();
                                while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                                    v_msg *= SCALER_UP;
                                    cur_max *= SCALER_UP;
                                    sib_boost++;
                                }
                            }
                        }
                        int net_boost = boost_counts[uid] - sib_boost;
                        while (net_boost > 0) {
                            v_msg *= SCALER_UP;
                            net_boost--;
                        }
                        while (net_boost < 0) {
                            v_msg *= SCALER_DOWN;
                            net_boost++;
                        }
                        V_all[vid] = v_msg;
                        U[vid] = P_b[vid].transpose() * v_msg;
                    }
                }

                Vector4 V_poly = (poly_id == tree.root_id) ? gtr.pi : V_all[poly_id];

                // Screen all pairs of children
                for (size_t i = 0; i < k; ++i) {
                    for (size_t j = i + 1; j < k; ++j) {
                        int32_t ci = p_children[i];
                        int32_t cj = p_children[j];

                        Vector4 D_pair = P_D[ci].cwiseProduct(P_D[cj]);
                        Vector4 V_pair = V_poly;
                        for (size_t m = 0; m < k; ++m) {
                            if (m != i && m != j) {
                                V_pair.array() *= P_D[p_children[m]].array();
                            }
                        }
                        Scalar dL = V_pair.dot(gtr.Q * D_pair);
                        Scalar grad = inv_L * dL;

                        if (grad > best_grad) {
                            best_grad = grad;
                            best_i = i;
                            best_j = j;
                        }
                    }
                }
            }

            // If the best directional gradient exceeds threshold, resolve the pair!
            if (best_grad > gradient_threshold) {
                int32_t ci = p_children[best_i];
                int32_t cj = p_children[best_j];

                // Create new internal node grouping (ci, cj)
                int32_t new_id = static_cast<int32_t>(tree.nodes.size());
                TreeNode new_node;
                new_node.id = new_id;
                new_node.name = "Resolved_" + std::to_string(new_id);
                new_node.parent_id = poly_id;
                new_node.children = {ci, cj};
                new_node.branch_length = 0.01;
                tree.nodes.push_back(std::move(new_node));

                tree.nodes[ci].parent_id = new_id;
                tree.nodes[cj].parent_id = new_id;

                // Replace ci and cj in poly_id's children with new_id
                auto& p_ch = tree.nodes[poly_id].children;
                p_ch.erase(std::remove(p_ch.begin(), p_ch.end(), ci), p_ch.end());
                p_ch.erase(std::remove(p_ch.begin(), p_ch.end(), cj), p_ch.end());
                p_ch.push_back(new_id);

                // Keep working buffers and post-order synchronized for subsequent polytomies in this pass
                num_nodes = tree.num_nodes();
                P_b.push_back(gtr.transition_matrix(0.01));
                D_all.resize(num_nodes, Vector4::Zero());
                V_all.resize(num_nodes, Vector4::Zero());
                node_log_scale.resize(num_nodes, 0.0);
                P_D.resize(num_nodes, Vector4::Zero());
                boost_counts.resize(num_nodes, 0);
                leaf_to_taxon.resize(num_nodes, static_cast<size_t>(-1));
                tree.compute_post_order();

                resolved_count++;
            }
        }

        return resolved_count;
    }

    // =========================================================================
    // 4b. Topological Rearrangement via Ohmic Resistance-Bounded SPR & Port Grounding
    // =========================================================================
    static bool is_in_subtree(const Tree& tree, int32_t query, int32_t root) {
        if (query == root) return true;
        int32_t curr = query;
        size_t depth = 0;
        size_t max_depth = tree.nodes.size();
        while (curr != INVALID_INDEX && depth < max_depth) {
            if (curr == root) return true;
            curr = tree.nodes[curr].parent_id;
            depth++;
        }
        return false;
    }

    struct PortRegroundResult {
        Tree tree;
        std::vector<int32_t> affected_branches;
    };

    // Native Ohmic Operator: Re-grounds a subtree C_u (currently attached to p_u = parent(u))
    // by sliding its external attachment port from u to an internal node x strictly inside C_u,
    // reversing the parent-child orientation along the internal path u <-> ... <-> x.
    static std::optional<PortRegroundResult> reground_subtree_at_internal_port(
        const Tree& orig,
        int32_t u,
        int32_t x
    ) {
        if (u == orig.root_id || u == INVALID_INDEX || x == INVALID_INDEX || u == x) return std::nullopt;
        if (orig.nodes[x].is_leaf) return std::nullopt;
        if (!is_in_subtree(orig, x, u)) return std::nullopt;

        int32_t p_u = orig.nodes[u].parent_id;
        if (p_u == INVALID_INDEX) return std::nullopt;

        // Collect path from u down to x: path[0] = u, path[1] = ..., path[k] = x
        std::vector<int32_t> path;
        int32_t curr = x;
        while (curr != INVALID_INDEX) {
            path.push_back(curr);
            if (curr == u) break;
            curr = orig.nodes[curr].parent_id;
        }
        if (path.empty() || path.back() != u || path.size() < 2) return std::nullopt;
        std::reverse(path.begin(), path.end()); // Now path[0] == u, path.back() == x

        Tree tree = orig;
        Scalar old_tu = std::max<Scalar>(1e-4, orig.nodes[u].branch_length);

        // Save original edge lengths along the path: edge_len[i] is branch length of path[i] (edge path[i-1] -> path[i])
        std::vector<Scalar> edge_len(path.size(), 0.0);
        for (size_t i = 1; i < path.size(); ++i) {
            edge_len[i] = std::max<Scalar>(1e-4, orig.nodes[path[i]].branch_length);
        }

        // Reverse edges along the path: for i = 0 .. k-1, path[i] becomes a child of path[i+1]
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            int32_t parent_old = path[i];
            int32_t child_old = path[i + 1];

            // Remove child_old from parent_old.children
            auto& ch_list = tree.nodes[parent_old].children;
            ch_list.erase(std::remove(ch_list.begin(), ch_list.end(), child_old), ch_list.end());

            // Add parent_old to child_old.children
            tree.nodes[child_old].children.push_back(parent_old);
            tree.nodes[parent_old].parent_id = child_old;
            tree.nodes[parent_old].branch_length = edge_len[i + 1];
        }

        // Attach x directly to p_u in place of u
        auto& pu_ch = tree.nodes[p_u].children;
        for (size_t i = 0; i < pu_ch.size(); ++i) {
            if (pu_ch[i] == u) {
                pu_ch[i] = x;
                break;
            }
        }
        tree.nodes[x].parent_id = p_u;
        tree.nodes[x].branch_length = std::max<Scalar>(1e-4, 0.5 * old_tu);

        std::vector<int32_t> raw_affected;
        raw_affected.push_back(x);
        for (int32_t c : tree.nodes[x].children) raw_affected.push_back(c);

        if (tree.nodes[u].children.size() == 1) {
            int32_t s = tree.nodes[u].children[0];
            int32_t p_new = tree.nodes[u].parent_id; // which is path[1]
            auto& pnew_ch = tree.nodes[p_new].children;
            for (size_t i = 0; i < pnew_ch.size(); ++i) {
                if (pnew_ch[i] == u) {
                    pnew_ch[i] = s;
                    break;
                }
            }
            tree.nodes[s].parent_id = p_new;
            tree.nodes[s].branch_length += tree.nodes[u].branch_length;
            tree.nodes[u].children.clear();
            tree.nodes[u].parent_id = INVALID_INDEX;
            tree.nodes[u].branch_length = 0.0;
            raw_affected.push_back(s);
        }

        for (size_t i = 0; i < path.size(); ++i) {
            raw_affected.push_back(path[i]);
        }

        std::vector<int32_t> affected;
        for (int32_t nid : raw_affected) {
            if (nid >= 0 && static_cast<size_t>(nid) < tree.nodes.size() &&
                nid != tree.root_id && tree.nodes[nid].parent_id != INVALID_INDEX) {
                affected.push_back(nid);
            }
        }
        std::sort(affected.begin(), affected.end());
        affected.erase(std::unique(affected.begin(), affected.end()), affected.end());

        tree.compute_post_order();
        if (tree.post_order.empty()) return std::nullopt;
        for (int32_t nid : tree.post_order) {
            if (nid < 0 || static_cast<size_t>(nid) >= tree.nodes.size()) return std::nullopt;
        }
        return PortRegroundResult{std::move(tree), std::move(affected)};
    }

    // Native Ohmic Operator: Sweeps internal clades C_u and slides their entry port to the
    // internal node x in C_u that minimizes effective resistance / maximizes GTR log-likelihood.
    static size_t run_ohmic_port_regrounding(
        Tree& tree,
        const Alignment& aln,
        const GTRMatrix& gtr,
        const Matrix& Omega,
        size_t max_ports_per_clade = 6
    ) {
        size_t accepted_regroundings = 0;
        size_t N = aln.num_taxa;
        if (N < 5) return 0;

        Scalar current_best_ll = LikelihoodEngine::compute_gtr_log_likelihood(tree, aln, gtr);

        auto compute_node_leaves = [&]() {
            std::vector<std::vector<size_t>> leaves(tree.nodes.size());
            for (int32_t nid : tree.post_order) {
                if (nid == INVALID_INDEX || static_cast<size_t>(nid) >= tree.nodes.size()) continue;
                if (tree.nodes[nid].is_leaf) {
                    size_t tid = LikelihoodEngine::find_taxon_index(aln, tree.nodes[nid].name);
                    if (tid != static_cast<size_t>(-1) && tid < N) {
                        leaves[nid].push_back(tid);
                    }
                } else {
                    for (int32_t c : tree.nodes[nid].children) {
                        if (c != INVALID_INDEX && static_cast<size_t>(c) < leaves.size()) {
                            leaves[nid].insert(leaves[nid].end(), leaves[c].begin(), leaves[c].end());
                        }
                    }
                }
            }
            return leaves;
        };

        std::vector<std::vector<size_t>> node_leaves = compute_node_leaves();

        std::vector<int32_t> internal_clades;
        for (const auto& node : tree.nodes) {
            if (!node.is_leaf && node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                if (static_cast<size_t>(node.id) < node_leaves.size() && node_leaves[node.id].size() >= 4) {
                    internal_clades.push_back(node.id);
                }
            }
        }

        // Prioritize clades with higher Ohmic tension (longer subtending branch or larger clade)
        std::sort(internal_clades.begin(), internal_clades.end(), [&](int32_t a, int32_t b) {
            return tree.nodes[a].branch_length > tree.nodes[b].branch_length;
        });

        for (int32_t u : internal_clades) {
            if (u < 0 || static_cast<size_t>(u) >= tree.nodes.size() || static_cast<size_t>(u) >= node_leaves.size()) continue;
            if (tree.nodes[u].is_leaf || tree.nodes[u].parent_id == INVALID_INDEX || u == tree.root_id) continue;
            // Verify u is still reachable from root_id
            if (!is_in_subtree(tree, u, tree.root_id)) continue;
            if (node_leaves[u].size() < 4) continue;

            const std::vector<size_t> leaves_u = node_leaves[u];

            // Identify external reference leaves (prefer sibling leaves of u, else all outside leaves)
            std::vector<uint8_t> in_u(N, 0);
            for (size_t lu : leaves_u) {
                if (lu < N) in_u[lu] = 1;
            }

            std::vector<size_t> ext_leaves;
            int32_t p_u = tree.nodes[u].parent_id;
            if (p_u != INVALID_INDEX && p_u >= 0 && static_cast<size_t>(p_u) < node_leaves.size()) {
                for (size_t lp : node_leaves[p_u]) {
                    if (lp < N && !in_u[lp]) ext_leaves.push_back(lp);
                }
            }
            if (ext_leaves.empty()) {
                for (size_t i = 0; i < N; ++i) {
                    if (!in_u[i]) ext_leaves.push_back(i);
                }
            }
            if (ext_leaves.empty()) continue;

            // Collect all internal nodes strictly inside C_u and rank by minimum resistance to ext_leaves
            struct PortCandidate {
                int32_t x;
                Scalar score;
            };
            std::vector<PortCandidate> ports;
            auto collect_internal = [&](int32_t nid, auto& self) -> void {
                for (int32_t c : tree.nodes[nid].children) {
                    if (c == INVALID_INDEX || static_cast<size_t>(c) >= tree.nodes.size()) continue;
                    if (!tree.nodes[c].is_leaf && static_cast<size_t>(c) < node_leaves.size() && !node_leaves[c].empty()) {
                        Scalar min_r = 1e9;
                        Scalar sum_r = 0.0;
                        size_t cnt = 0;
                        size_t step_c = std::max<size_t>(1, node_leaves[c].size() / 4);
                        size_t step_e = std::max<size_t>(1, ext_leaves.size() / 6);
                        for (size_t ic = 0; ic < node_leaves[c].size(); ic += step_c) {
                            size_t lc = node_leaves[c][ic];
                            for (size_t ie = 0; ie < ext_leaves.size(); ie += step_e) {
                                size_t le = ext_leaves[ie];
                                Scalar r = Omega(static_cast<Eigen::Index>(lc), static_cast<Eigen::Index>(le));
                                if (r < min_r) min_r = r;
                                sum_r += r;
                                cnt++;
                            }
                        }
                        Scalar mean_r = (cnt > 0) ? (sum_r / static_cast<Scalar>(cnt)) : min_r;
                        ports.push_back({c, 0.8 * min_r + 0.2 * mean_r});
                        self(c, self);
                    }
                }
            };
            collect_internal(u, collect_internal);
            if (ports.empty()) continue;

            std::sort(ports.begin(), ports.end(), [](const auto& a, const auto& b) {
                return a.score < b.score;
            });

            size_t limit = std::min(ports.size(), max_ports_per_clade);
            for (size_t i = 0; i < limit; ++i) {
                int32_t x = ports[i].x;
                auto maybe_reg = reground_subtree_at_internal_port(tree, u, x);
                if (!maybe_reg.has_value()) continue;

                Tree cand_tree = std::move(maybe_reg->tree);
                auto& aff = maybe_reg->affected_branches;
                Scalar cand_ll = optimize_branches_gtr(cand_tree, aln, gtr, 4, 1e-2, &aff);

                if (cand_ll > current_best_ll + 0.05) {
                    current_best_ll = cand_ll;
                    tree = std::move(cand_tree);
                    node_leaves = compute_node_leaves();
                    accepted_regroundings++;
                    break;
                }
            }
        }

        return accepted_regroundings;
    }

    struct SPRMoveResult {
        Tree tree;
        int32_t u;
        int32_t v;
        int32_t w;
        int32_t s;
    };

    static std::optional<SPRMoveResult> apply_spr_move(
        const Tree& orig,
        int32_t u,
        int32_t v
    ) {
        if (u == orig.root_id || v == orig.root_id) return std::nullopt;
        if (is_in_subtree(orig, v, u)) return std::nullopt;

        int32_t p = orig.nodes[u].parent_id;
        if (p == INVALID_INDEX) return std::nullopt;
        if (v == p) return std::nullopt;

        for (int32_t sib : orig.nodes[p].children) {
            if (sib == v) return std::nullopt; // Sibling regraft is trivial / NNI
        }

        Tree tree = orig;
        int32_t s = INVALID_INDEX;

        // Step 1: Prune u from p
        auto& p_ch = tree.nodes[p].children;
        p_ch.erase(std::remove(p_ch.begin(), p_ch.end(), u), p_ch.end());

        // Handle parent collapse if p has only 1 child left
        if (p_ch.size() == 1) {
            s = p_ch[0];
            int32_t gp = tree.nodes[p].parent_id;
            if (p == tree.root_id) {
                tree.root_id = s;
                tree.nodes[s].parent_id = INVALID_INDEX;
            } else {
                auto& gp_ch = tree.nodes[gp].children;
                for (size_t i = 0; i < gp_ch.size(); ++i) {
                    if (gp_ch[i] == p) {
                        gp_ch[i] = s;
                        break;
                    }
                }
                tree.nodes[s].parent_id = gp;
                tree.nodes[s].branch_length += tree.nodes[p].branch_length;
            }
            tree.nodes[p].children.clear();
            tree.nodes[p].parent_id = INVALID_INDEX;
        }

        // Step 2: Regraft u onto branch leading to v (bisecting edge (pv, v))
        int32_t pv = tree.nodes[v].parent_id;
        if (pv == INVALID_INDEX) return std::nullopt;

        int32_t w = -1;
        if (tree.nodes[p].parent_id == INVALID_INDEX && tree.nodes[p].children.empty() && p != tree.root_id) {
            w = p; // Reuse detached node
        } else {
            w = static_cast<int32_t>(tree.nodes.size());
            TreeNode new_w;
            new_w.id = w;
            tree.nodes.push_back(std::move(new_w));
        }

        std::string cand_name = "SPR_Node_" + std::to_string(w);
        int suffix = 1;
        auto name_exists = [&](const std::string& nm) {
            for (const auto& nd : tree.nodes) {
                if (nd.id != w && nd.name == nm) return true;
            }
            return false;
        };
        while (name_exists(cand_name)) {
            cand_name = "SPR_Node_" + std::to_string(w) + "_" + std::to_string(++suffix);
        }
        tree.nodes[w].name = cand_name;
        tree.nodes[w].is_leaf = false;
        tree.nodes[w].parent_id = pv;
        tree.nodes[w].children = {v, u};

        auto& pv_ch = tree.nodes[pv].children;
        for (size_t i = 0; i < pv_ch.size(); ++i) {
            if (pv_ch[i] == v) {
                pv_ch[i] = w;
                break;
            }
        }

        Scalar old_tv = std::max<Scalar>(1e-4, tree.nodes[v].branch_length);
        tree.nodes[w].branch_length = 0.5 * old_tv;
        tree.nodes[v].branch_length = 0.5 * old_tv;
        tree.nodes[v].parent_id = w;
        tree.nodes[u].parent_id = w;

        tree.compute_post_order();
        return SPRMoveResult{std::move(tree), u, v, w, s};
    }

    static size_t run_resistance_bounded_spr(
        Tree& tree,
        const Alignment& aln,
        GTRMatrix& gtr,
        size_t max_rounds = 4,
        size_t max_candidates = 16
    ) {
        size_t accepted_moves = 0;
        size_t N = aln.num_taxa;
        if (N < 4) return 0;

        // Compute N x N Foster effective resistance matrix Omega
        Matrix D = NeighborJoining::compute_tn93_distances(aln);
        Scalar sigma = 0.05;
        std::vector<Scalar> p_d;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                if (D(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) > 1e-5) {
                    p_d.push_back(D(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
                }
            }
        }
        if (!p_d.empty()) {
            std::sort(p_d.begin(), p_d.end());
            sigma = std::clamp(p_d[p_d.size() / 2], 1e-4, 5.0);
        }

        Matrix W = Matrix::Zero(N, N);
        Vector deg = Vector::Zero(N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                Scalar w = std::exp(-D(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) / sigma);
                W(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = w;
                W(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = w;
                deg(static_cast<Eigen::Index>(i)) += w;
                deg(static_cast<Eigen::Index>(j)) += w;
            }
        }
        Vector inv_sqrt_deg(N);
        for (size_t i = 0; i < N; ++i) {
            inv_sqrt_deg(static_cast<Eigen::Index>(i)) = (deg(static_cast<Eigen::Index>(i)) > 1e-12)
                ? (1.0 / std::sqrt(deg(static_cast<Eigen::Index>(i)))) : 1.0;
        }

        Matrix L_norm = Matrix::Identity(N, N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                if (i != j) {
                    L_norm(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
                        -W(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) *
                         inv_sqrt_deg(static_cast<Eigen::Index>(i)) * inv_sqrt_deg(static_cast<Eigen::Index>(j));
                }
            }
        }

        Eigen::SelfAdjointEigenSolver<Matrix> es(L_norm);
        Matrix L_pinv = Matrix::Zero(N, N);
        if (es.info() == Eigen::Success) {
            for (int k = 1; k < static_cast<int>(N); ++k) {
                Scalar lambda = es.eigenvalues()(k);
                if (lambda > 1e-8) {
                    Vector vk = es.eigenvectors().col(k);
                    L_pinv.noalias() += (1.0 / lambda) * (vk * vk.transpose());
                }
            }
        }

        Matrix Omega = Matrix::Zero(N, N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                Omega(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
                    (i == j) ? 0.0 : std::max<Scalar>(0.0,
                        L_pinv(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) +
                        L_pinv(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(j)) -
                        2.0 * L_pinv(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
            }
        }

        Scalar current_best_ll = LikelihoodEngine::compute_gtr_log_likelihood(tree, aln, gtr);

        for (size_t round = 0; round < max_rounds; ++round) {
            bool improved_in_round = false;

            // Phase A: Native Ohmic Port-Regrounding (relaxes inverted clades by sliding entry ports)
            size_t regrounded = run_ohmic_port_regrounding(tree, aln, gtr, Omega, 6);
            if (regrounded > 0) {
                accepted_moves += regrounded;
                improved_in_round = true;
                current_best_ll = optimize_branches_gtr(tree, aln, gtr, 10, 1e-3);
            }

            auto compute_node_leaves = [&]() {
                std::vector<std::vector<size_t>> leaves(tree.nodes.size());
                for (int32_t nid : tree.post_order) {
                    if (nid == INVALID_INDEX || static_cast<size_t>(nid) >= tree.nodes.size()) continue;
                    if (tree.nodes[nid].is_leaf) {
                        size_t tid = LikelihoodEngine::find_taxon_index(aln, tree.nodes[nid].name);
                        if (tid != static_cast<size_t>(-1) && tid < N) {
                            leaves[nid].push_back(tid);
                        }
                    } else {
                        for (int32_t c : tree.nodes[nid].children) {
                            if (c != INVALID_INDEX && static_cast<size_t>(c) < leaves.size()) {
                                leaves[nid].insert(leaves[nid].end(),
                                                   leaves[c].begin(),
                                                   leaves[c].end());
                            }
                        }
                    }
                }
                return leaves;
            };

            // Single linear pass over post_order to compute descendant leaves for every node
            std::vector<std::vector<size_t>> node_leaves = compute_node_leaves();

            // Snapshot candidate nodes to avoid iterator invalidation during structural SPR updates
            std::vector<int32_t> candidate_nodes;
            candidate_nodes.reserve(tree.nodes.size());
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && node.id != INVALID_INDEX) {
                    candidate_nodes.push_back(node.id);
                }
            }

            for (int32_t u : candidate_nodes) {
                if (u < 0 || static_cast<size_t>(u) >= tree.nodes.size() || static_cast<size_t>(u) >= node_leaves.size()) continue;
                if (tree.nodes[u].parent_id == INVALID_INDEX || node_leaves[u].empty()) continue;

                const std::vector<size_t> leaves_u = node_leaves[u];

                struct TargetCandidate {
                    int32_t v;
                    Scalar score;
                    int32_t parent_id;
                };
                std::vector<TargetCandidate> all_candidates;

                for (const auto& v_node : tree.nodes) {
                    int32_t v = v_node.id;
                    if (v == tree.root_id || v == INVALID_INDEX || v == u) continue;
                    if (v < 0 || static_cast<size_t>(v) >= node_leaves.size()) continue;
                    if (is_in_subtree(tree, v, u)) continue;
                    if (node_leaves[v].empty()) continue;

                    Scalar min_r = 1e9;
                    Scalar sum_r = 0.0;
                    size_t count_r = 0;
                    size_t step_u = std::max<size_t>(1, leaves_u.size() / 3);
                    size_t step_v = std::max<size_t>(1, node_leaves[v].size() / 3);
                    for (size_t iu = 0; iu < leaves_u.size(); iu += step_u) {
                        size_t lu = leaves_u[iu];
                        for (size_t iv = 0; iv < node_leaves[v].size(); iv += step_v) {
                            size_t lv = node_leaves[v][iv];
                            if (lu < N && lv < N) {
                                Scalar r = Omega(static_cast<Eigen::Index>(lu), static_cast<Eigen::Index>(lv));
                                if (r < min_r) min_r = r;
                                sum_r += r;
                                count_r++;
                            }
                        }
                    }
                    Scalar mean_r = (count_r > 0) ? (sum_r / static_cast<Scalar>(count_r)) : min_r;
                    Scalar composite_r = 0.7 * min_r + 0.3 * mean_r;
                    all_candidates.push_back({v, composite_r, v_node.parent_id});
                }

                std::sort(all_candidates.begin(), all_candidates.end(), [](const auto& a, const auto& b) {
                    return a.score < b.score;
                });

                // Multi-level Ohmic Radius Selection:
                // Tier 1: immediate low-resistance neighbors
                // Tier 2: diverse clade/parent bridges
                std::vector<int32_t> selected_targets;
                std::unordered_set<int32_t> seen_targets;
                std::unordered_map<int32_t, size_t> parent_counts;

                for (const auto& cand : all_candidates) {
                    if (selected_targets.size() >= max_candidates) break;
                    if (parent_counts[cand.parent_id] < 2 || selected_targets.size() < 6) {
                        selected_targets.push_back(cand.v);
                        seen_targets.insert(cand.v);
                        parent_counts[cand.parent_id]++;
                    }
                }

                for (const auto& cand : all_candidates) {
                    if (selected_targets.size() >= max_candidates) break;
                    if (seen_targets.find(cand.v) == seen_targets.end()) {
                        selected_targets.push_back(cand.v);
                        seen_targets.insert(cand.v);
                    }
                }

                for (int32_t v : selected_targets) {
                    auto maybe_res = apply_spr_move(tree, u, v);
                    if (!maybe_res.has_value()) continue;

                    Tree cand_tree = std::move(maybe_res->tree);
                    int32_t w = maybe_res->w;
                    int32_t s = maybe_res->s;

                    // Incident branches to optimize: {u, v, w} + optionally {s}
                    std::vector<int32_t> local_brs = {u, v, w};
                    if (s != INVALID_INDEX && s != cand_tree.root_id && cand_tree.nodes[s].parent_id != INVALID_INDEX) {
                        local_brs.push_back(s);
                    }

                    // Ultra-fast 3/4-branch polish (2 iterations of L-BFGS with exact analytic gradients)
                    Scalar cand_ll = optimize_branches_gtr(cand_tree, aln, gtr, 2, 1e-2, &local_brs);

                    if (cand_ll > current_best_ll + 0.05) {
                        current_best_ll = cand_ll;
                        tree = std::move(cand_tree);
                        node_leaves = compute_node_leaves();
                        accepted_moves++;
                        improved_in_round = true;
                        break; // Move accepted for node u; advance to next u
                    }
                }
            }

            if (improved_in_round) {
                // Synchronize all tree branches after an active SPR round
                current_best_ll = optimize_branches_gtr(tree, aln, gtr, 15, 1e-3);
            } else {
                break; // Converged
            }
        }

        return accepted_moves;
    }

    // =========================================================================
    // 4c. Ohmic Divide-and-Conquer: Subtree Extraction, Parallel Search & Stitch
    // =========================================================================
    static std::vector<std::string> get_subtree_leaf_names(const Tree& tree, int32_t nid) {
        std::vector<std::string> names;
        auto dfs = [&](int32_t u, auto& self) -> void {
            if (tree.nodes[u].is_leaf) {
                names.push_back(tree.nodes[u].name);
            } else {
                for (int32_t c : tree.nodes[u].children) {
                    self(c, self);
                }
            }
        };
        dfs(nid, dfs);
        return names;
    }

    static std::string subtree_to_newick(
        const Tree& tree,
        int32_t nid,
        bool include_branch_lengths = true
    ) {
        auto write_node = [&](auto& self, int32_t u) -> std::string {
            const auto& node = tree.nodes[u];
            std::string res;
            if (!node.children.empty()) {
                res += "(";
                for (size_t i = 0; i < node.children.size(); ++i) {
                    if (i > 0) res += ",";
                    res += self(self, node.children[i]);
                }
                res += ")";
            }
            res += node.name;
            if (include_branch_lengths && u != nid) {
                std::ostringstream ss;
                ss << ":" << std::setprecision(8) << node.branch_length;
                res += ss.str();
            }
            return res;
        };
        return write_node(write_node, nid) + ";";
    }

    static std::vector<int32_t> identify_ohmic_islands(
        const Tree& tree,
        size_t min_island_size = 4,
        size_t max_island_size = 14
    ) {
        std::vector<int32_t> island_roots;

        auto traverse = [&](int32_t u, auto& self) -> void {
            if (tree.nodes[u].is_leaf) return;

            auto leaves = get_subtree_leaf_names(tree, u);
            size_t m = leaves.size();

            // If this subtree is within the target island size, and not the entire root tree
            if (m <= max_island_size && u != tree.root_id) {
                if (m >= min_island_size) {
                    island_roots.push_back(u);
                }
                return; // Do not recurse further into this island
            }

            // Otherwise, recurse on children
            for (int32_t c : tree.nodes[u].children) {
                self(c, self);
            }
        };

        traverse(tree.root_id, traverse);
        return island_roots;
    }

    static Alignment create_taxon_sub_alignment(
        const Alignment& full_aln,
        const std::vector<std::string>& taxon_names
    ) {
        Alignment sub;
        sub.code = full_aln.code;
        sub.num_codons = 0;
        sub.num_nucleotides = full_aln.num_nucleotides;
        sub.sequences.reserve(taxon_names.size());
        sub.taxon_names.reserve(taxon_names.size());

        for (const auto& name : taxon_names) {
            auto it = full_aln.taxon_to_index.find(name);
            if (it != full_aln.taxon_to_index.end()) {
                sub.sequences.push_back(full_aln.sequences[it->second]);
                sub.taxon_names.push_back(name);
                sub.taxon_to_index[name] = sub.sequences.size() - 1;
                std::string lower = name;
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                sub.lower_taxon_to_index[lower] = sub.sequences.size() - 1;
            }
        }
        sub.num_taxa = sub.sequences.size();
        sub.finalize_nucleotide();
        return sub;
    }

    static std::string stitch_tree_newick(
        const Tree& orig_tree,
        const std::unordered_map<int32_t, std::string>& island_newicks
    ) {
        auto write_node = [&](auto& self, int32_t u) -> std::string {
            auto it = island_newicks.find(u);
            if (it != island_newicks.end()) {
                std::string s = it->second;
                while (!s.empty() && (s.back() == ';' || s.back() == ' ' || s.back() == '\n' || s.back() == '\r')) {
                    s.pop_back();
                }
                if (u != orig_tree.root_id) {
                    std::ostringstream ss;
                    ss << ":" << std::setprecision(8) << orig_tree.nodes[u].branch_length;
                    s += ss.str();
                }
                return s;
            }
            const auto& node = orig_tree.nodes[u];
            std::string res;
            if (!node.children.empty()) {
                res += "(";
                for (size_t i = 0; i < node.children.size(); ++i) {
                    if (i > 0) res += ",";
                    res += self(self, node.children[i]);
                }
                res += ")";
            }
            res += node.name;
            if (!node.model_tag.empty()) {
                res += "{" + node.model_tag + "}";
            }
            if (u != orig_tree.root_id) {
                std::ostringstream ss;
                ss << ":" << std::setprecision(8) << node.branch_length;
                res += ss.str();
            }
            return res;
        };
        return write_node(write_node, orig_tree.root_id) + ";";
    }

    static size_t run_ohmic_divide_and_conquer(
        Tree& tree,
        const Alignment& aln,
        const GTRMatrix& gtr,
        const TreeEngineSettings& s,
        std::vector<size_t>& out_island_sizes
    ) {
        auto island_roots = identify_ohmic_islands(tree, s.dnc_min_island_size, s.dnc_target_island_size);
        if (island_roots.size() < 2) {
            return 0; // Not enough islands to justify D&C
        }

        size_t num_islands = island_roots.size();
        out_island_sizes.resize(num_islands);
        std::vector<std::string> island_newicks_res(num_islands);
        std::vector<size_t> island_spr_res(num_islands, 0);

        Matrix D = NeighborJoining::compute_tn93_distances(aln);

        // Prepare Kron-Reduced / Boundary-Anchored island sub-alignments and initial trees
        // Each island k includes 1 external boundary anchor taxon a* (closest non-island leaf)
        // pinned at the root of the island tree so local SPR & port-regrounding preserve root polarity!
        std::vector<Alignment> sub_alns(num_islands);
        std::vector<Tree> island_trees(num_islands);
        std::vector<std::string> anchor_names(num_islands);

        for (size_t k = 0; k < num_islands; ++k) {
            int32_t r = island_roots[k];
            auto taxa = get_subtree_leaf_names(tree, r);
            out_island_sizes[k] = taxa.size();

            std::unordered_set<std::string> in_island(taxa.begin(), taxa.end());
            std::string best_anchor;
            Scalar best_d = 1e9;

            for (size_t i = 0; i < aln.num_taxa; ++i) {
                const std::string& cand_name = aln.sequences[i].name;
                if (in_island.count(cand_name)) continue;
                for (const auto& tname : taxa) {
                    auto it = aln.taxon_to_index.find(tname);
                    if (it != aln.taxon_to_index.end()) {
                        Scalar d = D(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(it->second));
                        if (d < best_d) {
                            best_d = d;
                            best_anchor = cand_name;
                        }
                    }
                }
            }

            anchor_names[k] = best_anchor;
            std::vector<std::string> taxa_with_anchor = taxa;
            std::string sub_nwk = subtree_to_newick(tree, r, true);
            if (!best_anchor.empty()) {
                taxa_with_anchor.push_back(best_anchor);
                while (!sub_nwk.empty() && sub_nwk.back() == ';') sub_nwk.pop_back();
                Scalar anchor_bl = std::max<Scalar>(1e-3, best_d * 0.5);
                std::ostringstream wrapped;
                wrapped << "(" << sub_nwk << ":" << std::setprecision(8)
                        << std::max<Scalar>(1e-4, tree.nodes[r].branch_length)
                        << "," << best_anchor << ":" << anchor_bl << ");";
                sub_nwk = wrapped.str();
            }

            sub_alns[k] = create_taxon_sub_alignment(aln, taxa_with_anchor);
            island_trees[k] = Tree::from_newick(sub_nwk);
        }

        // Parallel optimization of disjoint Boundary-Anchored Ohmic islands
        #pragma omp parallel for schedule(dynamic)
        for (size_t k = 0; k < num_islands; ++k) {
            auto& i_tree = island_trees[k];
            const auto& i_aln = sub_alns[k];
            GTRMatrix local_gtr = gtr;

            // 1. Local branch length optimization under GTR
            optimize_branches_gtr(i_tree, i_aln, local_gtr, 20, 1e-3);

            // 2. Local Ohmic resistance-bounded SPR + Port Regrounding
            size_t moves = run_resistance_bounded_spr(
                i_tree, i_aln, local_gtr, 2, std::min<size_t>(6, i_aln.num_taxa)
            );
            island_spr_res[k] = moves;

            // 3. Local polytomy handling
            if (s.allow_polytomies) {
                collapse_short_branches(i_tree, s.polytomy_collapse_threshold);
                if (s.perform_topological_gradients) {
                    resolve_polytomies_by_gradient(i_tree, i_aln, local_gtr, s.topological_grad_threshold);
                }
            }

            // 4. Final local polish
            optimize_branches_gtr(i_tree, i_aln, local_gtr, 15, 1e-3);

            // 5. Extract rooted ingroup Newick by stripping the boundary anchor taxon
            const std::string& anc = anchor_names[k];
            if (!anc.empty()) {
                int32_t anc_id = INVALID_INDEX;
                for (const auto& nd : i_tree.nodes) {
                    if (nd.is_leaf && nd.name == anc) {
                        anc_id = nd.id;
                        break;
                    }
                }
                // Only extract if anchor remained a direct child of root and root has 1 ingroup child
                if (anc_id != INVALID_INDEX &&
                    i_tree.nodes[anc_id].parent_id == i_tree.root_id &&
                    i_tree.nodes[i_tree.root_id].children.size() == 2) {
                    int32_t ingroup_root = (i_tree.nodes[i_tree.root_id].children[0] == anc_id)
                        ? i_tree.nodes[i_tree.root_id].children[1]
                        : i_tree.nodes[i_tree.root_id].children[0];
                    island_newicks_res[k] = subtree_to_newick(i_tree, ingroup_root, true);
                } else {
                    // Fallback: re-extract original ingroup if anchor migrated inside
                    island_newicks_res[k] = subtree_to_newick(tree, island_roots[k], true);
                }
            } else {
                island_newicks_res[k] = i_tree.to_newick(true);
            }
        }

        size_t total_spr_moves = 0;
        std::unordered_map<int32_t, std::string> map_newicks;
        for (size_t k = 0; k < num_islands; ++k) {
            map_newicks[island_roots[k]] = island_newicks_res[k];
            total_spr_moves += island_spr_res[k];
        }

        // Stitch optimized subtrees back into global tree
        std::string stitched_nwk = stitch_tree_newick(tree, map_newicks);
        Tree stitched_tree = Tree::from_newick(stitched_nwk);

        // Global harmonization: all-branch L-BFGS on stitched tree
        Scalar pre_ll = LikelihoodEngine::compute_gtr_log_likelihood(tree, aln, gtr);
        Scalar post_ll = optimize_branches_gtr(stitched_tree, aln, gtr, s.lbfgs_max_iterations, s.lbfgs_epsilon);

        if (post_ll >= pre_ll - 1.0) { // Accept if better or within numerical noise
            tree = std::move(stitched_tree);
        }

        return total_spr_moves;
    }

    // =========================================================================
    // 5. Unified ML Tree Search Engine Execution Pipeline
    // =========================================================================
    TreeEngineResult run(const Alignment& aln) {
        auto t_start = std::chrono::high_resolution_clock::now();
        TreeEngineResult res;
        res.num_taxa = aln.num_taxa;
        res.num_patterns = aln.patterns.size();

        if (!settings.quiet) {
            Panel::print_banner(
                "HYPHY-3 ML/OHMIC TREE ENGINE",
                "Fast Spectral Scaffolding • Native Polytomy Detection • All-Branch Quasi-Newton",
                "Next-generation topology search and branch optimization",
                "Authors: Sergei L. Kosakovsky Pond & HyPhy Team • v3.0.0"
            );
            std::cout << Console::bold("Input Data:") << " " << res.num_taxa << " taxa, "
                      << res.num_patterns << " unique site patterns ("
                      << aln.num_nucleotides << " nucleotides)\n";
        }

        // Extract nucleotide alignment for GTR phases
        Alignment nuc_aln = aln.to_nucleotide_alignment();
        res.num_patterns = nuc_aln.patterns.size();

        // Initialize Empirical GTR Model
        GTRMatrix gtr_model;
        GTRParameters gtr_params;
        gtr_model.update(gtr_params, nuc_aln.nuc_frequencies);
        res.gtr_frequencies = nuc_aln.nuc_frequencies;

        // Stage 1: Spectral Scaffolding
        if (!settings.quiet) {
            std::cout << Console::brand("[1/4] Constructing Spectral Tree Scaffold from Ohmic Conductance...\n");
        }
        Tree tree = build_spectral_scaffold(
            nuc_aln,
            settings.allow_polytomies,
            settings.spectral_star_threshold,
            settings.spectral_gap_ratio
        );
        res.initial_spectral_tree = tree;
        res.num_initial_nodes = tree.num_nodes();

        res.initial_gtr_log_l = LikelihoodEngine::compute_gtr_log_likelihood(tree, nuc_aln, gtr_model);

        if (!settings.quiet) {
            std::cout << "      Scaffold constructed (" << res.num_initial_nodes << " nodes). Initial GTR ln L: "
                      << std::fixed << std::setprecision(2) << res.initial_gtr_log_l << "\n";
        }

        // Stage 2: Simultaneous All-Branch L-BFGS & GTR Rate Optimization
        if (!settings.quiet) {
            std::cout << Console::brand("[2/4] Simultaneous All-Branch Quasi-Newton (L-BFGS) & GTR Rate Optimization...\n");
        }
        Scalar opt_ll = optimize_branches_gtr(
            tree, nuc_aln, gtr_model, settings.lbfgs_max_iterations, settings.lbfgs_epsilon
        );
        opt_ll = optimize_gtr_rates(tree, nuc_aln, gtr_model, gtr_params);
        opt_ll = optimize_branches_gtr(
            tree, nuc_aln, gtr_model, settings.lbfgs_max_iterations / 2, settings.lbfgs_epsilon
        );
        res.optimized_gtr_log_l = opt_ll;
        res.gtr_rates = {gtr_params.theta_AC, gtr_params.theta_AT, gtr_params.theta_CG, gtr_params.theta_CT, gtr_params.theta_GT};

        if (!settings.quiet) {
            std::cout << "      Optimized GTR ln L: " << std::fixed << std::setprecision(2) << opt_ll
                      << " (Δln L = +" << (opt_ll - res.initial_gtr_log_l) << ")\n"
                      << "      GTR Rates (AG=1): AC=" << std::setprecision(3) << gtr_params.theta_AC
                      << ", AT=" << gtr_params.theta_AT
                      << ", CG=" << gtr_params.theta_CG
                      << ", CT=" << gtr_params.theta_CT
                      << ", GT=" << gtr_params.theta_GT << "\n";
        }

        // Stage 2.5: Topological Rearrangement via Ohmic Divide-and-Conquer + Resistance-Bounded SPR
        if (settings.perform_spr) {
            size_t dnc_moves = 0;
            if (settings.enable_dnc && nuc_aln.num_taxa >= settings.dnc_min_taxa) {
                if (!settings.quiet) {
                    std::cout << Console::brand("[2.2/4] Ohmic Divide-and-Conquer Subtree Optimization (Parallel Islands)...\n");
                }
                dnc_moves = run_ohmic_divide_and_conquer(
                    tree, nuc_aln, gtr_model, settings, res.island_sizes
                );
                if (!res.island_sizes.empty()) {
                    res.dnc_used = true;
                    res.num_islands = res.island_sizes.size();
                    res.dnc_spr_moves = dnc_moves;
                    opt_ll = LikelihoodEngine::compute_gtr_log_likelihood(tree, nuc_aln, gtr_model);
                    res.optimized_gtr_log_l = opt_ll;
                    if (!settings.quiet) {
                        std::cout << "      Partitioned into " << res.num_islands << " Ohmic islands (taxa:";
                        for (size_t sz : res.island_sizes) std::cout << " " << sz;
                        std::cout << "). Accepted " << dnc_moves << " island SPR moves (GTR ln L = "
                                  << std::fixed << std::setprecision(2) << opt_ll << ")\n";
                    }
                }
            }

            if (!settings.quiet) {
                std::cout << Console::brand("[2.5/4] Global Backbone Polish via Ohmic Resistance-Bounded SPR...\n");
            }
            size_t spr_rounds = settings.spr_max_rounds;
            size_t spr_moves = run_resistance_bounded_spr(
                tree, nuc_aln, gtr_model, spr_rounds, settings.spr_max_candidates
            );
            res.num_spr_moves_accepted = dnc_moves + spr_moves;
            if (spr_moves > 0) {
                opt_ll = optimize_branches_gtr(
                    tree, nuc_aln, gtr_model, settings.lbfgs_max_iterations, settings.lbfgs_epsilon
                );
                opt_ll = optimize_gtr_rates(tree, nuc_aln, gtr_model, gtr_params);
                opt_ll = optimize_branches_gtr(
                    tree, nuc_aln, gtr_model, settings.lbfgs_max_iterations / 2, settings.lbfgs_epsilon
                );
                res.optimized_gtr_log_l = opt_ll;
                res.gtr_rates = {gtr_params.theta_AC, gtr_params.theta_AT, gtr_params.theta_CG, gtr_params.theta_CT, gtr_params.theta_GT};
                if (!settings.quiet) {
                    std::cout << "      Accepted " << spr_moves << " global SPR moves (New GTR ln L = "
                              << std::fixed << std::setprecision(2) << opt_ll << ")\n";
                }
            } else {
                if (!settings.quiet) {
                    std::cout << "      Backbone topology is optimal; 0 global SPR moves needed.\n";
                }
            }
        }

        // Stage 3: Polytomy Collapse (Zero-Signal Branches)
        if (settings.allow_polytomies) {
            if (!settings.quiet) {
                std::cout << Console::brand("[3/4] Screening Zero-Signal Branches for Polytomy Collapse...\n");
            }
            size_t collapsed = collapse_short_branches(tree, settings.polytomy_collapse_threshold);
            res.num_branches_collapsed = collapsed;

            if (collapsed > 0) {
                // Re-optimize remaining branches after collapse
                Scalar post_col_ll = optimize_branches_gtr(
                    tree, nuc_aln, gtr_model, settings.lbfgs_max_iterations / 2, settings.lbfgs_epsilon
                );
                res.optimized_gtr_log_l = post_col_ll;
                if (!settings.quiet) {
                    std::cout << "      Collapsed " << collapsed << " near-zero branches into native polytomies.\n";
                }
            } else {
                if (!settings.quiet) {
                    std::cout << "      No branches collapsed; all branches exceed threshold t > "
                              << settings.polytomy_collapse_threshold << ".\n";
                }
            }

            // Directional Topological Gradient testing
            if (settings.perform_topological_gradients) {
                size_t resolved = resolve_polytomies_by_gradient(
                    tree, nuc_aln, gtr_model, settings.topological_grad_threshold
                );
                res.num_polytomies_resolved = resolved;
                if (resolved > 0) {
                    Scalar post_res_ll = optimize_branches_gtr(
                        tree, nuc_aln, gtr_model, settings.lbfgs_max_iterations / 2, settings.lbfgs_epsilon
                    );
                    res.optimized_gtr_log_l = post_res_ll;
                    if (!settings.quiet) {
                        std::cout << "      Resolved " << resolved
                                  << " polytomies with positive topological gradients.\n";
                    }
                }
            }
        }

        // Count final retained polytomies
        size_t final_polys = 0;
        for (const auto& node : tree.nodes) {
            if (!node.is_leaf && node.children.size() >= 3) {
                final_polys++;
            }
        }
        res.num_polytomies_retained = final_polys;

        // Stage 3.5: Final Global GTR Polish (All Branches & Rates)
        if (!settings.quiet) {
            std::cout << Console::brand("[3.5/4] Final Global GTR Polish (All Branches & Rates)...\n");
        }
        opt_ll = optimize_branches_gtr(
            tree, nuc_aln, gtr_model, settings.lbfgs_max_iterations, settings.lbfgs_epsilon
        );
        opt_ll = optimize_gtr_rates(tree, nuc_aln, gtr_model, gtr_params);
        opt_ll = optimize_branches_gtr(
            tree, nuc_aln, gtr_model, settings.lbfgs_max_iterations / 2, settings.lbfgs_epsilon
        );
        res.optimized_gtr_log_l = opt_ll;
        res.gtr_rates = {gtr_params.theta_AC, gtr_params.theta_AT, gtr_params.theta_CG, gtr_params.theta_CT, gtr_params.theta_GT};

        if (!settings.quiet) {
            std::cout << "      Polished GTR ln L: " << std::fixed << std::setprecision(2) << opt_ll
                      << " (Rates: AC=" << std::setprecision(3) << gtr_params.theta_AC
                      << ", AT=" << gtr_params.theta_AT
                      << ", CG=" << gtr_params.theta_CG
                      << ", CT=" << gtr_params.theta_CT
                      << ", GT=" << gtr_params.theta_GT << ")\n";
        }

        // Stage 4: Model Escalation (GTR -> MG94 Codon Polish if explicitly requested)
        res.final_log_l = res.optimized_gtr_log_l;
        res.active_model = "GTR";

        if (settings.model == TreeSearchModel::MG94 || settings.perform_model_escalation) {
            if (aln.num_codons == 0) {
                if (!settings.quiet) {
                    std::cout << Console::warning("Warning: Codon MG94 model requested, but alignment is not an in-frame coding sequence. Retaining GTR model.\n");
                }
            } else {
                if (!settings.quiet) {
                    std::cout << Console::brand("[4/4] Escalating to Codon MG94 3x4 Full Model Optimization...\n");
                }
                try {
                    MG94Parameters base_params;
                    base_params.alpha = 1.0;
                    base_params.beta = 0.5;
                    base_params.theta_AC = gtr_params.theta_AC;
                    base_params.theta_AT = gtr_params.theta_AT;
                    base_params.theta_CG = gtr_params.theta_CG;
                    base_params.theta_CT = gtr_params.theta_CT;
                    base_params.theta_GT = gtr_params.theta_GT;

                    opt::MG94Fitter fitter(tree, aln);
                    auto fit_res = fitter.fit_full_model(
                        0.5,
                        base_params,
                        nullptr,
                        settings.lbfgs_max_iterations,
                        settings.lbfgs_epsilon
                    );

                    tree = fit_res.tree;
                    res.final_log_l = fit_res.log_likelihood;
                    res.active_model = "MG94";
                    res.mg94_omega = fit_res.x_opt(0);

                    if (!settings.quiet) {
                        std::cout << "      Codon MG94 ln L: " << std::fixed << std::setprecision(2) << res.final_log_l
                                  << ", Global omega: " << std::setprecision(4) << res.mg94_omega
                                  << " (iterations: " << fit_res.iterations << ")\n";
                    }
                } catch (const std::exception& e) {
                    if (!settings.quiet) {
                        std::cout << "      Codon escalation failed: " << e.what() << "\n";
                    }
                }
            }
        } else {
            if (!settings.quiet) {
                std::cout << Console::brand("[4/4] Finalizing ML Topology...\n");
            }
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t_end - t_start).count();
        res.optimized_tree = tree;
        res.num_final_nodes = tree.num_nodes();

        if (!settings.quiet) {
            std::cout << "\n" << Console::bold("=== ML Tree Search Summary ===") << "\n"
                      << "  Final Log-Likelihood   : " << std::fixed << std::setprecision(2) << res.final_log_l << "\n"
                      << "  Active Model           : " << res.active_model << "\n"
                      << "  Total Nodes            : " << res.num_final_nodes << " (Leaves: " << res.num_taxa << ")\n"
                      << "  Polytomies Retained    : " << res.num_polytomies_retained << "\n"
                      << "  Collapsed Zero-Branches: " << res.num_branches_collapsed << "\n"
                      << "  Runtime                : " << std::setprecision(2) << res.runtime_seconds << " seconds\n"
                      << "  Newick Topology        : " << tree.to_newick(true) << "\n\n";
        }

        return res;
    }

    nlohmann::json to_json(const TreeEngineResult& res, const Provenance& prov = {}) const;

    // -------------------------------------------------------------
    // Unified CLI Tree Resolution Helper
    // -------------------------------------------------------------
    static inline bool is_auto_tree_spec(const std::string& path) {
        if (path.empty()) return false;
        std::string lower = path;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        return lower == "auto";
    }

    static inline bool is_ohmic_tree_spec(const std::string& path) {
        if (path.empty()) return false;
        std::string lower = path;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        return lower == "ohmic";
    }

    static Tree resolve_tree(
        const Alignment& aln,
        const std::string& tree_file,
        bool force_ohmic = false,
        bool use_hyphaeon = false,
        const std::string& embeddings_file = "",
        hyphy::core::OhmicNetwork* out_net = nullptr,
        std::string* out_source_desc = nullptr,
        bool quiet = false
    ) {
        bool is_auto = is_auto_tree_spec(tree_file);
        bool is_ohmic = is_ohmic_tree_spec(tree_file);

        // 1. Explicit Ohmic network requested (--ohmic, --tree ohmic, or embeddings provided)
        if (force_ohmic || is_ohmic || !embeddings_file.empty()) {
            return hyphy::core::OhmicNetwork::resolve_tree(
                aln, "", true, use_hyphaeon, embeddings_file, out_net, out_source_desc
            );
        }

        // 2. Explicit "auto" requested (--tree auto)
        if (is_auto) {
            if (!quiet) {
                std::cout << Console::brand("[Phylogeny] Inferring ML tree via TreeEngine (Spectral + GTR + Foster SPR)...\n");
            }
            TreeEngineSettings settings;
            settings.quiet = true;
            TreeEngine engine(settings);
            auto res = engine.run(aln);
            if (out_source_desc) {
                *out_source_desc = "TreeEngine (auto, lnL=" +
                                   std::to_string(static_cast<int>(res.final_log_l)) +
                                   ", nodes=" + std::to_string(res.num_final_nodes) +
                                   ", poly=" + std::to_string(res.num_polytomies_retained) + ")";
            }
            if (!quiet) {
                std::cout << Console::brand("[Phylogeny] Tree inferred: ")
                          << res.num_taxa << " taxa, lnL = "
                          << std::fixed << std::setprecision(2) << res.final_log_l
                          << ", " << res.num_final_nodes << " nodes ("
                          << res.num_polytomies_retained << " polytomies retained)\n";
            }
            return res.optimized_tree;
        }

        // 3. Explicit tree file specified (and not "auto" / "ohmic")
        if (!tree_file.empty()) {
            if (out_source_desc) *out_source_desc = tree_file;
            return Tree::from_newick_file(tree_file);
        }

        // 4. Embedded Newick tree in alignment
        if (!aln.embedded_tree_newick.empty()) {
            if (out_source_desc) *out_source_desc = "Embedded in alignment";
            return Tree::from_newick(aln.embedded_tree_newick);
        }

        // 5. Default fallback when no tree, no embedded tree, and no --ohmic:
        // Automatically infer with TreeEngine!
        if (!quiet) {
            std::cout << Console::brand("[Phylogeny] No tree provided; auto-inferring ML tree via TreeEngine...\n");
        }
        TreeEngineSettings settings;
        settings.quiet = true;
        TreeEngine engine(settings);
        auto res = engine.run(aln);
        if (out_source_desc) {
            *out_source_desc = "TreeEngine (auto, lnL=" +
                               std::to_string(static_cast<int>(res.final_log_l)) +
                               ", nodes=" + std::to_string(res.num_final_nodes) +
                               ", poly=" + std::to_string(res.num_polytomies_retained) + ")";
        }
        if (!quiet) {
            std::cout << Console::brand("[Phylogeny] Tree inferred: ")
                      << res.num_taxa << " taxa, lnL = "
                      << std::fixed << std::setprecision(2) << res.final_log_l
                      << ", " << res.num_final_nodes << " nodes ("
                      << res.num_polytomies_retained << " polytomies retained)\n";
        }
        return res.optimized_tree;
    }

    // Overload without out_net pointer
    static Tree resolve_tree(
        const Alignment& aln,
        const std::string& tree_file,
        bool force_ohmic,
        bool use_hyphaeon,
        const std::string& embeddings_file,
        std::string* out_source_desc,
        bool quiet = false
    ) {
        return resolve_tree(aln, tree_file, force_ohmic, use_hyphaeon, embeddings_file, nullptr, out_source_desc, quiet);
    }
};

} // namespace hyphy::analyses
