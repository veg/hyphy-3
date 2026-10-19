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

        // Compute median positive distance for RBF kernel bandwidth sigma
        std::vector<Scalar> pos_dists;
        pos_dists.reserve(N * (N - 1) / 2);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                Scalar dij = D(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                if (dij > 1e-6 && std::isfinite(dij)) {
                    pos_dists.push_back(dij);
                }
            }
        }
        Scalar sigma = 0.05;
        if (!pos_dists.empty()) {
            std::sort(pos_dists.begin(), pos_dists.end());
            sigma = pos_dists[pos_dists.size() / 2];
            sigma = std::clamp(sigma, 1e-4, 5.0);
        }

        Tree tree;
        int32_t node_counter = 0;

        auto create_leaf = [&](size_t taxon_idx) -> int32_t {
            int32_t id = static_cast<int32_t>(tree.nodes.size());
            TreeNode node;
            node.id = id;
            node.name = aln.sequences[taxon_idx].name;
            node.is_leaf = true;
            node.branch_length = 0.01;
            tree.nodes.push_back(std::move(node));
            tree.leaf_name_to_id[aln.sequences[taxon_idx].name] = id;
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

        // Recursive Spectral Cheeger Partition
        std::function<int32_t(const std::vector<size_t>&)> recurse_partition;
        recurse_partition = [&](const std::vector<size_t>& indices) -> int32_t {
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

            // Build local conductance matrix W and normalized Laplacian L_norm
            Matrix W = Matrix::Zero(m, m);
            Vector deg = Vector::Zero(m);
            for (size_t i = 0; i < m; ++i) {
                size_t ti = indices[i];
                for (size_t j = i + 1; j < m; ++j) {
                    size_t tj = indices[j];
                    Scalar dij = D(static_cast<Eigen::Index>(ti), static_cast<Eigen::Index>(tj));
                    Scalar w = std::exp(-dij / sigma);
                    W(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = w;
                    W(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = w;
                    deg(static_cast<Eigen::Index>(i)) += w;
                    deg(static_cast<Eigen::Index>(j)) += w;
                }
            }

            Vector inv_sqrt_deg(m);
            for (size_t i = 0; i < m; ++i) {
                inv_sqrt_deg(static_cast<Eigen::Index>(i)) = (deg(static_cast<Eigen::Index>(i)) > 1e-12)
                    ? (1.0 / std::sqrt(deg(static_cast<Eigen::Index>(i)))) : 1.0;
            }

            Matrix L_norm = Matrix::Identity(m, m);
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    if (i != j) {
                        L_norm(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
                            -W(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) *
                             inv_sqrt_deg(static_cast<Eigen::Index>(i)) * inv_sqrt_deg(static_cast<Eigen::Index>(j));
                    }
                }
            }

            // Compute Eigendecomposition
            Eigen::SelfAdjointEigenSolver<Matrix> es(L_norm);
            if (es.info() != Eigen::Success) {
                // Fallback: direct star polytomy
                int32_t u = create_internal();
                for (size_t idx : indices) {
                    int32_t c = create_leaf(idx);
                    tree.nodes[u].children.push_back(c);
                    tree.nodes[c].parent_id = u;
                }
                return u;
            }

            Vector evals = es.eigenvalues();
            Scalar lambda2 = evals(1);
            Scalar lambda3 = (m > 2) ? evals(2) : 1.0;

            // Check if this cluster is a star-like radiation (polytomy)
            bool is_star = allow_polytomies && (m <= 8) &&
                           (lambda2 < star_threshold || (lambda3 > 1e-4 && (lambda2 / lambda3) < gap_ratio));

            if (is_star) {
                int32_t poly_node = create_internal();
                for (size_t idx : indices) {
                    int32_t c = create_leaf(idx);
                    tree.nodes[poly_node].children.push_back(c);
                    tree.nodes[c].parent_id = poly_node;
                    // Assign branch length from average distance to other members
                    Scalar avg_d = 0.0;
                    for (size_t o : indices) {
                        if (o != idx) avg_d += D(static_cast<Eigen::Index>(idx), static_cast<Eigen::Index>(o));
                    }
                    avg_d /= static_cast<Scalar>(m - 1);
                    tree.nodes[c].branch_length = std::max<Scalar>(1e-4, avg_d * 0.5);
                }
                return poly_node;
            }

            // Bipartition along the Fiedler vector v2 using maximum eigengap (Cheeger sweep cut)
            Vector v2 = es.eigenvectors().col(1);
            std::vector<size_t> grp1, grp2;

            std::vector<size_t> sort_order(m);
            std::iota(sort_order.begin(), sort_order.end(), 0);
            std::sort(sort_order.begin(), sort_order.end(), [&](size_t a, size_t b) {
                return v2(static_cast<Eigen::Index>(a)) < v2(static_cast<Eigen::Index>(b));
            });

            // Find the maximum eigengap in sorted v2
            Scalar max_gap = -1.0;
            size_t best_cut = m / 2; // fallback

            for (size_t k = 0; k < m - 1; ++k) {
                Scalar gap = v2(static_cast<Eigen::Index>(sort_order[k + 1])) - 
                             v2(static_cast<Eigen::Index>(sort_order[k]));
                if (gap > max_gap) {
                    max_gap = gap;
                    best_cut = k + 1;
                }
            }

            for (size_t k = 0; k < m; ++k) {
                if (k < best_cut) {
                    grp1.push_back(indices[sort_order[k]]);
                } else {
                    grp2.push_back(indices[sort_order[k]]);
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

            int32_t parent = create_internal();
            int32_t child1 = recurse_partition(grp1);
            int32_t child2 = recurse_partition(grp2);

            tree.nodes[parent].children = {child1, child2};
            tree.nodes[child1].parent_id = parent;
            tree.nodes[child2].parent_id = parent;

            // Set branch lengths between clusters
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

        int32_t root = recurse_partition(all_taxa);
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

                #ifdef _OPENMP
                #pragma omp parallel if(!omp_in_parallel())
                #endif
                {
                    Scalar local_ll = 0.0;
                    std::vector<Scalar> local_dL(n_nodes, 0.0);

                    #pragma omp for schedule(dynamic)
                    for (size_t p = 0; p < n_pat; ++p) {
                        const auto& pattern = aln_ref.patterns[p];

                        auto io = LikelihoodEngine::compute_inside_outside(
                            tree_ref, pattern, l2t, P_b, gtr_ref.pi, 4);
                        if (io.likelihood > 0.0) {
                            local_ll += pattern.weight * io.log_likelihood;
                            Scalar inv_L = pattern.weight / io.likelihood;
                            // Only optimized branches consume gradients; dP_b
                            // is uninitialized on the others.
                            for (int32_t vid : brs) {
                                local_dL[vid] += inv_L * io.V[vid].dot(dP_b[vid] * io.D[vid]);
                            }
                        }
                    }

                    #pragma omp critical
                    {
                        total_ll += local_ll;
                        for (size_t i = 0; i < n_nodes; ++i) {
                            dL_dt[i] += local_dL[i];
                        }
                    }
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
        opt_param.m = std::min<int>(6, std::max<int>(2, static_cast<int>(n_b)));
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
        auto P_b = LikelihoodEngine::build_P_branches<Matrix4>(tree, gtr);

        // Find internal nodes with >= 3 children (polytomies)
        std::vector<int32_t> polytomies;
        for (const auto& node : tree.nodes) {
            if (!node.is_leaf && node.children.size() >= 3) {
                polytomies.push_back(node.id);
            }
        }

        if (polytomies.empty()) return 0;

        // Inside-outside partials at all nodes; P_D keeps the per-child
        // messages P_b[c] * D[c] for the pairwise child screening below.
        std::vector<Vector4> P_D(num_nodes);

        for (int32_t poly_id : polytomies) {
            const auto& pnode = tree.nodes[poly_id];
            size_t k = pnode.children.size();
            if (k < 3) continue;

            // Evaluate directional derivatives for all pairs of children (ci, cj)
            Scalar best_grad = -1e9;
            size_t best_i = 0, best_j = 1;

            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];

                auto io = LikelihoodEngine::compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_b, gtr.pi, 4, &P_D);
                if (io.likelihood <= 0.0) continue;
                Scalar inv_L = pattern.weight / io.likelihood;

                Vector4 V_poly = (poly_id == tree.root_id) ? gtr.pi : io.V[poly_id];

                // Screen all pairs of children
                for (size_t i = 0; i < k; ++i) {
                    for (size_t j = i + 1; j < k; ++j) {
                        int32_t ci = pnode.children[i];
                        int32_t cj = pnode.children[j];

                        Vector4 D_pair = P_D[ci].cwiseProduct(P_D[cj]);
                        Vector4 V_pair = V_poly;
                        for (size_t m = 0; m < k; ++m) {
                            if (m != i && m != j) {
                                V_pair.array() *= P_D[pnode.children[m]].array();
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
                int32_t ci = pnode.children[best_i];
                int32_t cj = pnode.children[best_j];

                // Create new internal node grouping (ci, cj)
                int32_t new_id = static_cast<int32_t>(tree.nodes.size());
                TreeNode new_node;
                new_node.id = new_id;
                new_node.name = "Resolved_" + std::to_string(resolved_count + 1);
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

                resolved_count++;
            }
        }

        if (resolved_count > 0) {
            tree.compute_post_order();
        }
        return resolved_count;
    }

    // =========================================================================
    // 4b. Topological Rearrangement via Ohmic Resistance-Bounded SPR
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

        tree.nodes[w].name = "SPR_Node_" + std::to_string(w);
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

                const auto& leaves_u = node_leaves[u];

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

        // Prepare island sub-alignments and initial trees
        std::vector<Alignment> sub_alns(num_islands);
        std::vector<Tree> island_trees(num_islands);

        for (size_t k = 0; k < num_islands; ++k) {
            int32_t r = island_roots[k];
            auto taxa = get_subtree_leaf_names(tree, r);
            out_island_sizes[k] = taxa.size();
            sub_alns[k] = create_taxon_sub_alignment(aln, taxa);
            std::string sub_nwk = subtree_to_newick(tree, r, true);
            island_trees[k] = Tree::from_newick(sub_nwk);
        }

        // Parallel optimization of disjoint Ohmic islands
        #pragma omp parallel for schedule(dynamic)
        for (size_t k = 0; k < num_islands; ++k) {
            auto& i_tree = island_trees[k];
            const auto& i_aln = sub_alns[k];
            GTRMatrix local_gtr = gtr;

            // 1. Local branch length optimization under GTR
            optimize_branches_gtr(i_tree, i_aln, local_gtr, 20, 1e-3);

            // 2. Local Ohmic resistance-bounded SPR
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

            // Export local Newick representation
            island_newicks_res[k] = i_tree.to_newick(true);
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
