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
    static inline int& user_requested_threads_storage() {
        static int val = 0;
        return val;
    }

    static void set_user_requested_threads(int n) {
        user_requested_threads_storage() = n;
    }

    static int get_user_requested_threads() {
        return user_requested_threads_storage();
    }

    // Neumaier compensated summation for deterministic, high-precision reduction across site patterns
    static inline void neumaier_acc(Scalar& sum, Scalar& c, Scalar x) {
        Scalar t = sum + x;
        if (std::abs(sum) >= std::abs(x)) {
            c += (sum - t) + x;
        } else {
            c += (x - t) + sum;
        }
        sum = t;
    }

    static Scalar neumaier_sum(const std::vector<Scalar>& terms) {
        Scalar sum = 0.0;
        Scalar c = 0.0;
        for (Scalar x : terms) {
            neumaier_acc(sum, c, x);
        }
        return sum + c;
    }

    template <typename VecType>
    static inline void init_leaf_likelihood(
        const SitePattern& pattern,
        size_t t_idx,
        int S,
        VecType& out_L
    ) {
        if (t_idx == static_cast<size_t>(-1)) {
            out_L.setOnes();
            return;
        }
        int8_t state = pattern.states[t_idx];
        if (state >= 0 && state < S) {
            out_L.setZero();
            out_L(state) = 1.0;
        } else if (t_idx < pattern.ambiguous_masks.size() && pattern.ambiguous_masks[t_idx] != 0) {
            uint64_t mask = pattern.ambiguous_masks[t_idx];
            uint64_t all_sense = (S >= 64) ? ~uint64_t(0) : ((uint64_t(1) << S) - 1);
            if ((mask & all_sense) == all_sense) {
                out_L.setOnes();
            } else {
                out_L.setZero();
                for (int s = 0; s < S; ++s) {
                    if (mask & (uint64_t(1) << s)) {
                        out_L(s) = 1.0;
                    }
                }
            }
        } else {
            out_L.setOnes();
        }
    }

    // Workload-aware OpenMP thread cap to prevent thread oversubscription and allocator/barrier contention
    // on small alignments (e.g. 135 patterns on 48- or 64-core servers).
    static int optimal_omp_threads(size_t num_items, size_t min_items_per_thread = 16, int default_cap = 8) {
#ifdef _OPENMP
        int max_t = omp_get_max_threads();
        int user_t = get_user_requested_threads();
        int cap = (user_t > 0) ? std::min(max_t, user_t) : std::min(max_t, std::max(default_cap, 1));
        int work_cap = static_cast<int>(std::max<size_t>(1, (num_items + min_items_per_thread - 1) / min_items_per_thread));
        return std::max(1, std::min({max_t, cap, work_cap}));
#else
        (void)num_items;
        (void)min_items_per_thread;
        (void)default_cap;
        return 1;
#endif
    }

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
        std::vector<Matrix4> P_branches(num_nodes, Matrix4::Identity());
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P_branches[node.id] = gtr_model.transition_matrix(node.branch_length);
            }
        }

        std::vector<std::array<Matrix4, 3>> P_bridges;
        const std::vector<std::array<Matrix4, 3>>* p_bridges_ptr = nullptr;
        if (!tree.bridges.empty()) {
            P_bridges.resize(tree.bridges.size());
            for (size_t k = 0; k < tree.bridges.size(); ++k) {
                for (int r = 0; r < 3; ++r) {
                    if (tree.bridges[k].weights[r] > 0.0 || r == 0) {
                        P_bridges[k][r] = gtr_model.transition_matrix(tree.bridges[k].branch_lengths[r]);
                    } else {
                        P_bridges[k][r] = Matrix4::Identity();
                    }
                }
            }
            p_bridges_ptr = &P_bridges;
        }

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        std::vector<Scalar> pattern_ll(num_patterns, 0.0);
        int n_threads = optimal_omp_threads(num_patterns, 32);

        #pragma omp parallel num_threads(n_threads)
        {
            std::vector<Vector4> node_L(num_nodes);
            std::vector<double> node_log_scale(num_nodes, 0.0);
            #pragma omp for schedule(static)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];

                for (int32_t node_id : tree.post_order) {
                    const auto& node = tree.nodes[node_id];

                    if (node.is_leaf) {
                        size_t t_idx = leaf_to_taxon[node_id];
                        init_leaf_likelihood(pattern, t_idx, 4, node_L[node_id]);
                        node_log_scale[node_id] = 0.0;
                    } else {
                        double child_scale_sum = 0.0;
                        tree.eval_internal_node(node_id, P_branches, node_L, node_log_scale, node_L[node_id], child_scale_sum, p_bridges_ptr);
                        node_log_scale[node_id] = child_scale_sum;
                    }
                }

                Scalar pattern_likelihood = gtr_model.pi.dot(node_L[tree.root_id]);
                if (pattern_likelihood > 0.0) {
                    pattern_ll[p] = pattern.weight * (std::log(pattern_likelihood) + node_log_scale[tree.root_id]);
                } else {
                    pattern_ll[p] = pattern.weight * (-1e20);
                }
            }
        }

        return neumaier_sum(pattern_ll);
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
        std::vector<Matrix> P_branches(num_nodes, Matrix::Identity(S, S));
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P_branches[node.id] = mg94_model.transition_matrix(node.branch_length * branch_scale);
            }
        }

        std::vector<std::array<Matrix, 3>> P_bridges;
        const std::vector<std::array<Matrix, 3>>* p_bridges_ptr = nullptr;
        if (!tree.bridges.empty()) {
            P_bridges.resize(tree.bridges.size());
            for (size_t k = 0; k < tree.bridges.size(); ++k) {
                for (int r = 0; r < 3; ++r) {
                    if (tree.bridges[k].weights[r] > 0.0 || r == 0) {
                        P_bridges[k][r] = mg94_model.transition_matrix(tree.bridges[k].branch_lengths[r] * branch_scale);
                    } else {
                        P_bridges[k][r] = Matrix::Identity(S, S);
                    }
                }
            }
            p_bridges_ptr = &P_bridges;
        }

        return compute_tree_log_likelihood_from_P(tree, aln, P_branches, mg94_model.pi, S, custom_leaf_to_taxon, p_bridges_ptr);
    }

    // Reusable thread-local workspace for batched 61 x B_tile Inside-only (Felsenstein pruning) GEMM traversals
    struct BatchedInsideWorkspace {
        std::vector<Matrix> D;              // [num_nodes] each S x B_cap
        Matrix M_tmp;                       // S x B_cap scratch for P_branches[cid] * D[cid]
        Eigen::MatrixXd node_log_scale;     // num_nodes x B_cap
        Vector root_L;                      // B_cap
        Vector root_log_L;                  // B_cap

        void init(size_t num_nodes, int S, size_t B_cap) {
            if (D.size() != num_nodes || (num_nodes > 0 && (D[0].rows() != S || static_cast<size_t>(D[0].cols()) < B_cap))) {
                D.assign(num_nodes, Matrix::Zero(S, B_cap));
                M_tmp.resize(S, B_cap);
                node_log_scale.resize(num_nodes, B_cap);
                root_L.resize(B_cap);
                root_log_L.resize(B_cap);
            }
        }
    };

    // Reusable thread-local workspace for batched 61 x B_tile Inside-Outside GEMM traversals
    struct BatchedInsideOutsideWorkspace {
        std::vector<Matrix> D;              // [num_nodes] each S x B_cap
        std::vector<Matrix> V;              // [num_nodes] each S x B_cap
        std::vector<Matrix> P_D;            // [num_nodes] each S x B_cap
        std::vector<Matrix> U;              // [num_nodes] each S x B_cap
        Eigen::MatrixXi boost_counts;       // num_nodes x B_cap
        Eigen::MatrixXd node_log_scale;     // num_nodes x B_cap
        std::vector<Matrix> pref;           // [max_children] each S x B_cap
        std::vector<Matrix> suff;           // [max_children] each S x B_cap
        Eigen::MatrixXi pref_boost;         // max_children x B_cap
        Eigen::MatrixXi suff_boost;         // max_children x B_cap
        Vector likelihood;                  // B_cap
        Vector log_likelihood;              // B_cap

        void init(size_t num_nodes, int S, size_t B_cap) {
            if (D.size() != num_nodes || (num_nodes > 0 && (D[0].rows() != S || static_cast<size_t>(D[0].cols()) < B_cap))) {
                D.assign(num_nodes, Matrix::Zero(S, B_cap));
                V.assign(num_nodes, Matrix::Zero(S, B_cap));
                P_D.assign(num_nodes, Matrix::Zero(S, B_cap));
                U.assign(num_nodes, Matrix::Zero(S, B_cap));
                boost_counts.resize(num_nodes, B_cap);
                node_log_scale.resize(num_nodes, B_cap);
                likelihood.resize(B_cap);
                log_likelihood.resize(B_cap);
            }
        }
    };

    // Batched Inside-only (Felsenstein pruning) pass over pattern range [p_start, p_end) using S x B GEMM
    static void compute_inside_batch_in_place(
        const Tree& tree,
        const Alignment& aln,
        size_t p_start,
        size_t p_end,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<Matrix>& P_branches,
        const Vector& pi,
        int S,
        BatchedInsideWorkspace& ws,
        const std::vector<std::array<Matrix, 3>>* P_bridges = nullptr
    ) {
        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double LOG_SCALER = 345.38776394910684;

        const Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);
        if (B <= 0) return;
        size_t num_nodes = tree.num_nodes();
        ws.init(num_nodes, S, static_cast<size_t>(B));

        for (int32_t nid : tree.post_order) {
            const auto& node = tree.nodes[nid];
            auto D_nid = ws.D[nid].leftCols(B);

            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                for (Eigen::Index b = 0; b < B; ++b) {
                    auto col = D_nid.col(b);
                    init_leaf_likelihood(aln.patterns[p_start + b], t_idx, S, col);
                    ws.node_log_scale(nid, b) = 0.0;
                }
            } else {
                int32_t b_idx = tree.bridge_index(nid);
                if (b_idx >= 0) {
                    const auto& wb = tree.bridges[static_cast<size_t>(b_idx)];
                    int32_t id_A = wb.children[0];
                    int32_t id_B = wb.children[1];
                    int32_t id_C = wb.children[2];

                    Matrix M_A(S, B), M_B(S, B), M_C(S, B);
                    M_A.noalias() = P_branches[id_A] * ws.D[id_A].leftCols(B);
                    M_B.noalias() = P_branches[id_B] * ws.D[id_B].leftCols(B);
                    M_C.noalias() = P_branches[id_C] * ws.D[id_C].leftCols(B);
                    ws.node_log_scale.row(nid).head(B) =
                        ws.node_log_scale.row(id_A).head(B) +
                        ws.node_log_scale.row(id_B).head(B) +
                        ws.node_log_scale.row(id_C).head(B);

                    auto get_P_br = [&](int r) -> const Matrix& {
                        if (P_bridges == nullptr || static_cast<size_t>(b_idx) >= P_bridges->size()) {
                            throw std::logic_error("Wheatstone bridge transition matrices P_bridges required when tree has bridges");
                        }
                        return (*P_bridges)[static_cast<size_t>(b_idx)][r];
                    };

                    D_nid.setZero();
                    Matrix Y_tmp(S, B), PB_Y(S, B);
                    if (wb.weights[0] > 0.0) {
                        Y_tmp.array() = M_A.array() * M_B.array();
                        PB_Y.noalias() = get_P_br(0) * Y_tmp;
                        D_nid.array() += wb.weights[0] * (PB_Y.array() * M_C.array());
                    }
                    if (wb.weights[1] > 0.0) {
                        Y_tmp.array() = M_A.array() * M_C.array();
                        PB_Y.noalias() = get_P_br(1) * Y_tmp;
                        D_nid.array() += wb.weights[1] * (PB_Y.array() * M_B.array());
                    }
                    if (wb.weights[2] > 0.0) {
                        Y_tmp.array() = M_B.array() * M_C.array();
                        PB_Y.noalias() = get_P_br(2) * Y_tmp;
                        D_nid.array() += wb.weights[2] * (PB_Y.array() * M_A.array());
                    }
                    for (Eigen::Index b = 0; b < B; ++b) {
                        double max_val = D_nid.col(b).maxCoeff();
                        int boost_count = 0;
                        while (max_val < SCALER_THRESH && max_val > 0.0) {
                            D_nid.col(b) *= SCALER_UP;
                            max_val *= SCALER_UP;
                            boost_count++;
                        }
                        ws.node_log_scale(nid, b) -= boost_count * LOG_SCALER;
                    }
                } else {
                    D_nid.setOnes();
                    ws.node_log_scale.row(nid).head(B).setZero();
                    auto M_sub = ws.M_tmp.leftCols(B);

                    for (int32_t cid : node.children) {
                        // Single BLAS GEMM: (S x S) * (S x B) -> (S x B)
                        M_sub.noalias() = P_branches[cid] * ws.D[cid].leftCols(B);
                        D_nid.array() *= M_sub.array();
                        ws.node_log_scale.row(nid).head(B) += ws.node_log_scale.row(cid).head(B);

                        for (Eigen::Index b = 0; b < B; ++b) {
                            double cur_max = D_nid.col(b).maxCoeff();
                            int boost_count = 0;
                            while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                                D_nid.col(b) *= SCALER_UP;
                                cur_max *= SCALER_UP;
                                boost_count++;
                            }
                            if (boost_count > 0) {
                                ws.node_log_scale(nid, b) -= boost_count * LOG_SCALER;
                            }
                        }
                    }

                    for (Eigen::Index b = 0; b < B; ++b) {
                        double max_val = D_nid.col(b).maxCoeff();
                        int boost_count = 0;
                        while (max_val < SCALER_THRESH && max_val > 0.0) {
                            D_nid.col(b) *= SCALER_UP;
                            max_val *= SCALER_UP;
                            boost_count++;
                        }
                        if (boost_count > 0) {
                            ws.node_log_scale(nid, b) -= boost_count * LOG_SCALER;
                        }
                    }
                }
            }
        }

        ws.root_L.head(B).noalias() = ws.D[tree.root_id].leftCols(B).transpose() * pi;
        for (Eigen::Index b = 0; b < B; ++b) {
            Scalar L_b = ws.root_L(b);
            if (L_b > 0.0) {
                ws.root_log_L(b) = std::log(L_b) + ws.node_log_scale(tree.root_id, b);
            } else {
                ws.root_log_L(b) = -1e20;
            }
        }
    }

    // Batched Inside-Outside two-pass traversal over pattern range [p_start, p_end) using S x B GEMM
    static void compute_inside_outside_batch_in_place(
        const Tree& tree,
        const Alignment& aln,
        size_t p_start,
        size_t p_end,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<Matrix>& P_branches,
        const Vector& pi,
        int S,
        BatchedInsideOutsideWorkspace& ws,
        const std::vector<std::array<Matrix, 3>>* P_bridges = nullptr
    ) {
        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double SCALER_DOWN = 1e-150;
        const double LOG_SCALER = 345.38776394910684;

        const Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);
        if (B <= 0) return;
        size_t num_nodes = tree.num_nodes();
        ws.init(num_nodes, S, static_cast<size_t>(B));

        // 1. Post-order (bottom-up) pass: compute D[nid] (S x B) and cache P_D[cid] = P_branches[cid] * D[cid] via GEMM
        for (int32_t nid : tree.post_order) {
            const auto& node = tree.nodes[nid];
            auto D_nid = ws.D[nid].leftCols(B);

            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                for (Eigen::Index b = 0; b < B; ++b) {
                    auto col = D_nid.col(b);
                    init_leaf_likelihood(aln.patterns[p_start + b], t_idx, S, col);
                    ws.node_log_scale(nid, b) = 0.0;
                    ws.boost_counts(nid, b) = 0;
                }
            } else {
                for (int32_t cid : node.children) {
                    // Single BLAS GEMM: (S x S) * (S x B) -> (S x B)
                    ws.P_D[cid].leftCols(B).noalias() = P_branches[cid] * ws.D[cid].leftCols(B);
                }

                int32_t b_idx = tree.bridge_index(nid);
                if (b_idx >= 0) {
                    const auto& wb = tree.bridges[static_cast<size_t>(b_idx)];
                    int32_t id_A = wb.children[0];
                    int32_t id_B = wb.children[1];
                    int32_t id_C = wb.children[2];

                    auto M_A = ws.P_D[id_A].leftCols(B);
                    auto M_B = ws.P_D[id_B].leftCols(B);
                    auto M_C = ws.P_D[id_C].leftCols(B);
                    ws.node_log_scale.row(nid).head(B) =
                        ws.node_log_scale.row(id_A).head(B) +
                        ws.node_log_scale.row(id_B).head(B) +
                        ws.node_log_scale.row(id_C).head(B);

                    auto get_P_br = [&](int r) -> const Matrix& {
                        if (P_bridges == nullptr || static_cast<size_t>(b_idx) >= P_bridges->size()) {
                            throw std::logic_error("Wheatstone bridge transition matrices P_bridges required when tree has bridges");
                        }
                        return (*P_bridges)[static_cast<size_t>(b_idx)][r];
                    };

                    D_nid.setZero();
                    Matrix Y_tmp(S, B), PB_Y(S, B);
                    if (wb.weights[0] > 0.0) {
                        Y_tmp.array() = M_A.array() * M_B.array();
                        PB_Y.noalias() = get_P_br(0) * Y_tmp;
                        D_nid.array() += wb.weights[0] * (PB_Y.array() * M_C.array());
                    }
                    if (wb.weights[1] > 0.0) {
                        Y_tmp.array() = M_A.array() * M_C.array();
                        PB_Y.noalias() = get_P_br(1) * Y_tmp;
                        D_nid.array() += wb.weights[1] * (PB_Y.array() * M_B.array());
                    }
                    if (wb.weights[2] > 0.0) {
                        Y_tmp.array() = M_B.array() * M_C.array();
                        PB_Y.noalias() = get_P_br(2) * Y_tmp;
                        D_nid.array() += wb.weights[2] * (PB_Y.array() * M_A.array());
                    }
                    for (Eigen::Index b = 0; b < B; ++b) {
                        double max_val = D_nid.col(b).maxCoeff();
                        int boost_count = 0;
                        while (max_val < SCALER_THRESH && max_val > 0.0) {
                            D_nid.col(b) *= SCALER_UP;
                            max_val *= SCALER_UP;
                            boost_count++;
                        }
                        ws.boost_counts(nid, b) = boost_count;
                        ws.node_log_scale(nid, b) -= boost_count * LOG_SCALER;
                    }
                } else {
                    D_nid.setOnes();
                    ws.node_log_scale.row(nid).head(B).setZero();
                    ws.boost_counts.row(nid).head(B).setZero();

                    for (int32_t cid : node.children) {
                        D_nid.array() *= ws.P_D[cid].leftCols(B).array();
                        ws.node_log_scale.row(nid).head(B) += ws.node_log_scale.row(cid).head(B);

                        for (Eigen::Index b = 0; b < B; ++b) {
                            double cur_max = D_nid.col(b).maxCoeff();
                            while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                                D_nid.col(b) *= SCALER_UP;
                                cur_max *= SCALER_UP;
                                ws.boost_counts(nid, b) += 1;
                            }
                        }
                    }

                    for (Eigen::Index b = 0; b < B; ++b) {
                        double max_val = D_nid.col(b).maxCoeff();
                        while (max_val < SCALER_THRESH && max_val > 0.0) {
                            D_nid.col(b) *= SCALER_UP;
                            max_val *= SCALER_UP;
                            ws.boost_counts(nid, b) += 1;
                        }
                        ws.node_log_scale(nid, b) -= ws.boost_counts(nid, b) * LOG_SCALER;
                    }
                }
            }
        }

        ws.likelihood.head(B).noalias() = ws.D[tree.root_id].leftCols(B).transpose() * pi;
        for (Eigen::Index b = 0; b < B; ++b) {
            Scalar L_b = ws.likelihood(b);
            if (L_b > 0.0) {
                ws.log_likelihood(b) = std::log(L_b) + ws.node_log_scale(tree.root_id, b);
            } else {
                ws.log_likelihood(b) = -1e20;
            }
        }

        // 2. Pre-order (top-down) pass: compute V[vid] and U[vid] = P_branches[vid]^T * V[vid] via GEMM
        ws.U[tree.root_id].leftCols(B).colwise() = pi;

        for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
            int32_t uid = *it;
            const auto& node = tree.nodes[uid];
            if (node.is_leaf) continue;

            size_t num_children = node.children.size();
            if (num_children == 0) continue;

            int32_t b_idx = tree.bridge_index(uid);
            if (b_idx >= 0) {
                const auto& wb = tree.bridges[static_cast<size_t>(b_idx)];
                int32_t id_A = wb.children[0];
                int32_t id_B = wb.children[1];
                int32_t id_C = wb.children[2];

                auto M_A = ws.P_D[id_A].leftCols(B);
                auto M_B = ws.P_D[id_B].leftCols(B);
                auto M_C = ws.P_D[id_C].leftCols(B);
                auto U_uid = ws.U[uid].leftCols(B);

                auto get_P_br = [&](int r) -> const Matrix& {
                    if (P_bridges == nullptr) {
                        throw std::logic_error(
                            "LikelihoodEngine::compute_inside_outside_batch_in_place requires non-null P_bridges when tree.bridges is non-empty");
                    }
                    return (*P_bridges)[static_cast<size_t>(b_idx)][r];
                };

                auto V_A = ws.V[id_A].leftCols(B);
                auto V_B = ws.V[id_B].leftCols(B);
                auto V_C = ws.V[id_C].leftCols(B);
                V_A.setZero();
                V_B.setZero();
                V_C.setZero();

                Matrix Y_tmp(S, B), PB_Y(S, B), QT_Y(S, B);
                if (wb.weights[0] > 0.0) {
                    const Matrix& P0 = get_P_br(0);
                    Y_tmp.array() = M_A.array() * M_B.array();
                    PB_Y.noalias() = P0 * Y_tmp;
                    Y_tmp.array() = U_uid.array() * M_C.array();
                    QT_Y.noalias() = P0.transpose() * Y_tmp;
                    V_A.array() += wb.weights[0] * (QT_Y.array() * M_B.array());
                    V_B.array() += wb.weights[0] * (QT_Y.array() * M_A.array());
                    V_C.array() += wb.weights[0] * (U_uid.array() * PB_Y.array());
                }
                if (wb.weights[1] > 0.0) {
                    const Matrix& P1 = get_P_br(1);
                    Y_tmp.array() = M_A.array() * M_C.array();
                    PB_Y.noalias() = P1 * Y_tmp;
                    Y_tmp.array() = U_uid.array() * M_B.array();
                    QT_Y.noalias() = P1.transpose() * Y_tmp;
                    V_A.array() += wb.weights[1] * (QT_Y.array() * M_C.array());
                    V_B.array() += wb.weights[1] * (U_uid.array() * PB_Y.array());
                    V_C.array() += wb.weights[1] * (QT_Y.array() * M_A.array());
                }
                if (wb.weights[2] > 0.0) {
                    const Matrix& P2 = get_P_br(2);
                    Y_tmp.array() = M_B.array() * M_C.array();
                    PB_Y.noalias() = P2 * Y_tmp;
                    Y_tmp.array() = U_uid.array() * M_A.array();
                    QT_Y.noalias() = P2.transpose() * Y_tmp;
                    V_A.array() += wb.weights[2] * (U_uid.array() * PB_Y.array());
                    V_B.array() += wb.weights[2] * (QT_Y.array() * M_C.array());
                    V_C.array() += wb.weights[2] * (QT_Y.array() * M_B.array());
                }

                for (int32_t vid : wb.children) {
                    auto V_vid = ws.V[vid].leftCols(B);
                    for (Eigen::Index b = 0; b < B; ++b) {
                        int net_boost = ws.boost_counts(uid, b);
                        while (net_boost > 0) {
                            V_vid.col(b) *= SCALER_UP;
                            net_boost--;
                        }
                        while (net_boost < 0) {
                            V_vid.col(b) *= SCALER_DOWN;
                            net_boost++;
                        }
                    }
                    ws.U[vid].leftCols(B).noalias() = P_branches[vid].transpose() * V_vid;
                }
                continue;
            }

            if (ws.pref.size() < num_children || (num_children > 0 && static_cast<size_t>(ws.pref[0].cols()) < static_cast<size_t>(B))) {
                ws.pref.assign(num_children, Matrix::Ones(S, ws.D[0].cols()));
                ws.suff.assign(num_children, Matrix::Ones(S, ws.D[0].cols()));
                ws.pref_boost.resize(num_children, ws.D[0].cols());
                ws.suff_boost.resize(num_children, ws.D[0].cols());
            }

            // Prefix products across children (with per-column dynamic rescaling)
            ws.pref[0].leftCols(B).setOnes();
            ws.pref_boost.row(0).head(B).setZero();
            for (size_t i = 1; i < num_children; ++i) {
                int32_t prev_cid = node.children[i - 1];
                ws.pref[i].leftCols(B).array() = ws.pref[i - 1].leftCols(B).array() * ws.P_D[prev_cid].leftCols(B).array();
                for (Eigen::Index b = 0; b < B; ++b) {
                    int b_cnt = ws.pref_boost(i - 1, b);
                    double cur_max = ws.pref[i].col(b).maxCoeff();
                    while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                        ws.pref[i].col(b) *= SCALER_UP;
                        cur_max *= SCALER_UP;
                        b_cnt++;
                    }
                    ws.pref_boost(i, b) = b_cnt;
                }
            }

            // Suffix products across children (with per-column dynamic rescaling)
            ws.suff[num_children - 1].leftCols(B).setOnes();
            ws.suff_boost.row(num_children - 1).head(B).setZero();
            for (size_t i = num_children - 1; i > 0; --i) {
                int32_t next_cid = node.children[i];
                ws.suff[i - 1].leftCols(B).array() = ws.suff[i].leftCols(B).array() * ws.P_D[next_cid].leftCols(B).array();
                for (Eigen::Index b = 0; b < B; ++b) {
                    int b_cnt = ws.suff_boost(i, b);
                    double cur_max = ws.suff[i - 1].col(b).maxCoeff();
                    while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                        ws.suff[i - 1].col(b) *= SCALER_UP;
                        cur_max *= SCALER_UP;
                        b_cnt++;
                    }
                    ws.suff_boost(i - 1, b) = b_cnt;
                }
            }

            for (size_t i = 0; i < num_children; ++i) {
                int32_t vid = node.children[i];
                auto V_vid = ws.V[vid].leftCols(B);
                V_vid.array() = ws.U[uid].leftCols(B).array() * ws.pref[i].leftCols(B).array() * ws.suff[i].leftCols(B).array();
                for (Eigen::Index b = 0; b < B; ++b) {
                    int net_boost = ws.boost_counts(uid, b) - ws.pref_boost(i, b) - ws.suff_boost(i, b);
                    while (net_boost > 0) {
                        V_vid.col(b) *= SCALER_UP;
                        net_boost--;
                    }
                    while (net_boost < 0) {
                        V_vid.col(b) *= SCALER_DOWN;
                        net_boost++;
                    }
                }
                // Single BLAS GEMM: (S x S)^T * (S x B) -> (S x B)
                ws.U[vid].leftCols(B).noalias() = P_branches[vid].transpose() * V_vid;
            }
        }
    }

    // Fast single-pass codon tree likelihood from precomputed transition matrices P_branches
    // Uses batched 61 x B_tile GEMM across site patterns while preserving exact deterministic Neumaier summation
    static Scalar compute_tree_log_likelihood_from_P(
        const Tree& tree,
        const Alignment& aln,
        const std::vector<Matrix>& P_branches,
        const Vector& pi,
        int S,
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr,
        const std::vector<std::array<Matrix, 3>>* P_bridges = nullptr
    ) {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        if (num_patterns == 0) return 0.0;

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        std::vector<Scalar> pattern_ll(num_patterns, 0.0);
        int n_threads = optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            BatchedInsideWorkspace ws;
            ws.init(num_nodes, S, tile_size);

            #pragma omp for schedule(static)
            for (size_t t = 0; t < num_tiles; ++t) {
                size_t p_start = t * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                compute_inside_batch_in_place(tree, aln, p_start, p_end, leaf_to_taxon, P_branches, pi, S, ws, P_bridges);
                for (size_t p = p_start; p < p_end; ++p) {
                    pattern_ll[p] = aln.patterns[p].weight * ws.root_log_L(static_cast<Eigen::Index>(p - p_start));
                }
            }
        }

        return neumaier_sum(pattern_ll);
    }

    // Inside-Outside two-pass traversal for a single site pattern
    // Returns downward vectors D, parent-sibling context vectors V, and total likelihood L
    struct InsideOutsideResult {
        std::vector<Vector> D; // Downward subtree likelihoods D[node]
        std::vector<Vector> V; // Upward parent-sibling messages V[node]
        Scalar likelihood = 0.0;
        Scalar log_likelihood = 0.0;
    };

    // Reusable thread-local workspace for zero-allocation Inside-Outside traversals
    struct InsideOutsideWorkspace {
        std::vector<Vector> D;
        std::vector<Vector> V;
        std::vector<Vector> P_D;
        std::vector<Vector> U;
        std::vector<int> boost_counts;
        std::vector<double> node_log_scale;
        std::vector<Vector> pref;
        std::vector<Vector> suff;
        std::vector<int> pref_boost;
        std::vector<int> suff_boost;
        Vector v_msg;
        Scalar likelihood = 0.0;
        Scalar log_likelihood = 0.0;

        void init(size_t num_nodes, int S) {
            if (D.size() != num_nodes || (num_nodes > 0 && D[0].size() != S)) {
                D.assign(num_nodes, Vector::Zero(S));
                V.assign(num_nodes, Vector::Zero(S));
                P_D.assign(num_nodes, Vector::Zero(S));
                U.assign(num_nodes, Vector::Zero(S));
                boost_counts.assign(num_nodes, 0);
                node_log_scale.assign(num_nodes, 0.0);
                v_msg.resize(S);
            }
        }
    };

    static void compute_inside_outside_in_place(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<Matrix>& P_branches,
        const Vector& pi,
        int S,
        InsideOutsideWorkspace& ws,
        const std::vector<std::array<Matrix, 3>>* P_bridges = nullptr
    ) {
        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double SCALER_DOWN = 1e-150;

        size_t num_nodes = tree.num_nodes();
        ws.init(num_nodes, S);

        // 1. Post-order (bottom-up) pass to compute D and cache P_D
        for (int32_t nid : tree.post_order) {
            const auto& node = tree.nodes[nid];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                init_leaf_likelihood(pattern, t_idx, S, ws.D[nid]);
                ws.node_log_scale[nid] = 0.0;
                ws.boost_counts[nid] = 0;
            } else {
                double child_scale_sum = 0.0;
                for (int32_t cid : node.children) {
                    ws.P_D[cid].noalias() = P_branches[cid] * ws.D[cid];
                }
                int boost_cnt = 0;
                tree.eval_internal_node(nid, P_branches, ws.D, ws.node_log_scale, ws.D[nid], child_scale_sum, P_bridges, nullptr, &boost_cnt);
                ws.boost_counts[nid] = boost_cnt;
                ws.node_log_scale[nid] = child_scale_sum;
            }
        }

        ws.likelihood = ws.D[tree.root_id].dot(pi);
        if (ws.likelihood > 0.0) {
            ws.log_likelihood = std::log(ws.likelihood) + ws.node_log_scale[tree.root_id];
        } else {
            ws.log_likelihood = -1e20;
        }

        // 2. Pre-order (top-down) pass to compute U and V reusing cached P_D
        // Uses scaled prefix/suffix sibling products + integer boost counts to prevent
        // both sibling underflow and B[uid] = (1e150)^k overflow on high-degree polytomies.
        ws.U[tree.root_id] = pi;

        for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
            int32_t uid = *it;
            const auto& node = tree.nodes[uid];
            if (node.is_leaf) continue;

            size_t num_children = node.children.size();
            if (num_children == 0) continue;

            int32_t b_idx = tree.bridge_index(uid);
            if (b_idx >= 0) {
                const auto& wb = tree.bridges[static_cast<size_t>(b_idx)];
                int32_t id_A = wb.children[0];
                int32_t id_B = wb.children[1];
                int32_t id_C = wb.children[2];

                const Vector& M_A = ws.P_D[id_A];
                const Vector& M_B = ws.P_D[id_B];
                const Vector& M_C = ws.P_D[id_C];
                const Vector& U_uid = ws.U[uid];

                auto get_P_br = [&](int r) -> const Matrix& {
                    if (P_bridges == nullptr) {
                        throw std::logic_error(
                            "LikelihoodEngine::compute_inside_outside_in_place requires non-null P_bridges when tree.bridges is non-empty");
                    }
                    return (*P_bridges)[static_cast<size_t>(b_idx)][r];
                };

                ws.V[id_A].setZero();
                ws.V[id_B].setZero();
                ws.V[id_C].setZero();

                if (wb.weights[0] > 0.0) {
                    const Matrix& P0 = get_P_br(0);
                    Vector PB0 = P0 * M_A.cwiseProduct(M_B);
                    Vector QT0 = P0.transpose() * U_uid.cwiseProduct(M_C);
                    ws.V[id_A].noalias() += wb.weights[0] * QT0.cwiseProduct(M_B);
                    ws.V[id_B].noalias() += wb.weights[0] * QT0.cwiseProduct(M_A);
                    ws.V[id_C].noalias() += wb.weights[0] * U_uid.cwiseProduct(PB0);
                }
                if (wb.weights[1] > 0.0) {
                    const Matrix& P1 = get_P_br(1);
                    Vector PB1 = P1 * M_A.cwiseProduct(M_C);
                    Vector QT1 = P1.transpose() * U_uid.cwiseProduct(M_B);
                    ws.V[id_A].noalias() += wb.weights[1] * QT1.cwiseProduct(M_C);
                    ws.V[id_B].noalias() += wb.weights[1] * U_uid.cwiseProduct(PB1);
                    ws.V[id_C].noalias() += wb.weights[1] * QT1.cwiseProduct(M_A);
                }
                if (wb.weights[2] > 0.0) {
                    const Matrix& P2 = get_P_br(2);
                    Vector PB2 = P2 * M_B.cwiseProduct(M_C);
                    Vector QT2 = P2.transpose() * U_uid.cwiseProduct(M_A);
                    ws.V[id_A].noalias() += wb.weights[2] * U_uid.cwiseProduct(PB2);
                    ws.V[id_B].noalias() += wb.weights[2] * QT2.cwiseProduct(M_C);
                    ws.V[id_C].noalias() += wb.weights[2] * QT2.cwiseProduct(M_B);
                }

                int total_boost = ws.boost_counts[uid];
                for (int32_t vid : wb.children) {
                    int net_boost = total_boost;
                    while (net_boost > 0) {
                        ws.V[vid] *= SCALER_UP;
                        net_boost--;
                    }
                    while (net_boost < 0) {
                        ws.V[vid] *= SCALER_DOWN;
                        net_boost++;
                    }
                    ws.U[vid].noalias() = P_branches[vid].transpose() * ws.V[vid];
                }
                continue;
            }

            if (ws.pref.size() < num_children) {
                ws.pref.resize(num_children, Vector::Ones(S));
                ws.suff.resize(num_children, Vector::Ones(S));
                ws.pref_boost.resize(num_children, 0);
                ws.suff_boost.resize(num_children, 0);
            }

            // Prefix products: pref[i] = prod_{j < i} P_D[children[j]] (with dynamic rescaling)
            ws.pref[0].setOnes(S);
            ws.pref_boost[0] = 0;
            for (size_t i = 1; i < num_children; ++i) {
                int32_t prev_cid = node.children[i - 1];
                ws.pref[i] = ws.pref[i - 1].cwiseProduct(ws.P_D[prev_cid]);
                int b_cnt = ws.pref_boost[i - 1];
                double cur_max = ws.pref[i].maxCoeff();
                while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                    ws.pref[i] *= SCALER_UP;
                    cur_max *= SCALER_UP;
                    b_cnt++;
                }
                ws.pref_boost[i] = b_cnt;
            }

            // Suffix products: suff[i] = prod_{j > i} P_D[children[j]] (with dynamic rescaling)
            ws.suff[num_children - 1].setOnes(S);
            ws.suff_boost[num_children - 1] = 0;
            for (size_t i = num_children - 1; i > 0; --i) {
                int32_t next_cid = node.children[i];
                ws.suff[i - 1] = ws.suff[i].cwiseProduct(ws.P_D[next_cid]);
                int b_cnt = ws.suff_boost[i];
                double cur_max = ws.suff[i - 1].maxCoeff();
                while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                    ws.suff[i - 1] *= SCALER_UP;
                    cur_max *= SCALER_UP;
                    b_cnt++;
                }
                ws.suff_boost[i - 1] = b_cnt;
            }

            int total_boost = ws.boost_counts[uid];
            for (size_t i = 0; i < num_children; ++i) {
                int32_t vid = node.children[i];
                ws.v_msg = ws.U[uid].cwiseProduct(ws.pref[i]).cwiseProduct(ws.suff[i]);
                int net_boost = total_boost - ws.pref_boost[i] - ws.suff_boost[i];
                while (net_boost > 0) {
                    ws.v_msg *= SCALER_UP;
                    net_boost--;
                }
                while (net_boost < 0) {
                    ws.v_msg *= SCALER_DOWN;
                    net_boost++;
                }
                ws.V[vid] = ws.v_msg;
                ws.U[vid].noalias() = P_branches[vid].transpose() * ws.v_msg;
            }
        }
    }

    static InsideOutsideResult compute_inside_outside(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<Matrix>& P_branches,
        const Vector& pi,
        int S,
        const std::vector<std::array<Matrix, 3>>* P_bridges = nullptr
    ) {
        InsideOutsideWorkspace ws;
        compute_inside_outside_in_place(tree, pattern, leaf_to_taxon, P_branches, pi, S, ws, P_bridges);
        InsideOutsideResult res;
        res.D = std::move(ws.D);
        res.V = std::move(ws.V);
        res.likelihood = ws.likelihood;
        res.log_likelihood = ws.log_likelihood;
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
    // using the batched 2-pass Inside-Outside GEMM formulation
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
        std::vector<Matrix> P_branches(num_nodes, Matrix::Identity(S, S));
        std::vector<Matrix> dP_branches(num_nodes, Matrix::Zero(S, S));
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar bl = (use_custom_bl && static_cast<size_t>(node.id) < branch_lengths.size())
                    ? branch_lengths[node.id] : node.branch_length;
                P_branches[node.id] = mg.transition_matrix(bl);
                dP_branches[node.id] = mg.Q * P_branches[node.id];
            }
        }

        std::vector<std::array<Matrix, 3>> P_bridges;
        const std::vector<std::array<Matrix, 3>>* p_bridges_ptr = nullptr;
        if (!tree.bridges.empty()) {
            P_bridges.resize(tree.bridges.size());
            for (size_t k = 0; k < tree.bridges.size(); ++k) {
                for (int r = 0; r < 3; ++r) {
                    if (tree.bridges[k].weights[r] > 0.0 || r == 0) {
                        P_bridges[k][r] = mg.transition_matrix(tree.bridges[k].branch_lengths[r]);
                    } else {
                        P_bridges[k][r] = Matrix::Identity(S, S);
                    }
                }
            }
            p_bridges_ptr = &P_bridges;
        }

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        std::vector<Scalar> pattern_ll(num_patterns, 0.0);
        Matrix pattern_grad = Matrix::Zero(num_nodes, num_patterns);
        int n_threads = optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            BatchedInsideOutsideWorkspace ws;
            ws.init(num_nodes, S, tile_size);
            Matrix dP_D(S, tile_size);
            Vector inv_L_vec(tile_size);

            #pragma omp for schedule(static)
            for (size_t t = 0; t < num_tiles; ++t) {
                size_t p_start = t * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);

                compute_inside_outside_batch_in_place(
                    tree, aln, p_start, p_end, leaf_to_taxon, P_branches, aln.codon_frequencies_f3x4, S, ws, p_bridges_ptr
                );

                for (Eigen::Index b = 0; b < B; ++b) {
                    size_t p = p_start + b;
                    Scalar L = ws.likelihood(b);
                    if (L > 0.0) {
                        pattern_ll[p] = aln.patterns[p].weight * ws.log_likelihood(b);
                        inv_L_vec(b) = aln.patterns[p].weight / L;
                    } else {
                        inv_L_vec(b) = 0.0;
                    }
                }

                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) {
                        int32_t vid = node.id;
                        auto dP_D_sub = dP_D.leftCols(B);
                        // Batched GEMM: (S x S) * (S x B) -> (S x B)
                        dP_D_sub.noalias() = dP_branches[vid] * ws.D[vid].leftCols(B);
                        for (Eigen::Index b = 0; b < B; ++b) {
                            if (inv_L_vec(b) > 0.0) {
                                pattern_grad(vid, p_start + b) = inv_L_vec(b) * ws.V[vid].col(b).dot(dP_D_sub.col(b));
                            }
                        }
                    }
                }
            }
        }

        Scalar total_log_l = neumaier_sum(pattern_ll);
        std::vector<Scalar> grad_b(num_nodes, 0.0);
        std::vector<Scalar> row_terms(num_patterns, 0.0);
        for (size_t i = 0; i < num_nodes; ++i) {
            for (size_t p = 0; p < num_patterns; ++p) {
                row_terms[p] = pattern_grad(i, p);
            }
            grad_b[i] = neumaier_sum(row_terms);
        }

        return {total_log_l, grad_b};
    }

    struct JointMG94Gradients {
        Scalar log_likelihood = 0.0;
        std::vector<Scalar> grad_branches; // d ln L / d t_b (in substitutions/site)
        Scalar grad_omega = 0.0;           // d ln L / d omega (or omega_test when branch partition is active)
        Scalar grad_omega_bg = 0.0;        // d ln L / d omega_bg (when branch partition is active)
        Scalar grad_theta_AC = 0.0;        // d ln L / d theta_AC
        Scalar grad_theta_AT = 0.0;        // d ln L / d theta_AT
        Scalar grad_theta_CG = 0.0;        // d ln L / d theta_CG
        Scalar grad_theta_CT = 0.0;        // d ln L / d theta_CT
        Scalar grad_theta_GT = 0.0;        // d ln L / d theta_GT
    };

    // Evaluates log-likelihood and exact joint analytical gradients for all branch lengths
    // and global parameters (omega [or omega_test + omega_bg], theta_AC, theta_AT, theta_CG, theta_CT, theta_GT)
    // using the batched 2-pass Inside-Outside GEMM and spectral Fréchet derivatives of the matrix exponential.
    static JointMG94Gradients compute_joint_frechet_gradients(
        const Tree& tree,
        const Alignment& aln,
        const MG94Parameters& params,
        const std::vector<Scalar>& branch_lengths = {},
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr,
        const std::vector<bool>* is_test_branch = nullptr,
        Scalar omega_bg = -1.0
    ) {
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();

        bool has_partition = (is_test_branch != nullptr && is_test_branch->size() == num_nodes && omega_bg >= 0.0);

        MG94Matrix mg_test;
        mg_test.update(params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
        Scalar mu_test = std::max(mg_test.scale_factor, Scalar(1e-12));
        Scalar conv_test = 3.0 / mu_test;

        MG94Parameters params_bg = params;
        params_bg.beta = has_partition ? omega_bg : params.beta;
        MG94Matrix mg_bg;
        if (has_partition) {
            mg_bg.update(params_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
        }
        Scalar mu_bg = has_partition ? std::max(mg_bg.scale_factor, Scalar(1e-12)) : mu_test;
        Scalar conv_bg = 3.0 / mu_bg;

        std::vector<Scalar> tau(num_nodes, 0.0);
        std::vector<Matrix> P_branches(num_nodes, Matrix::Identity(S, S));
        std::vector<Matrix> dP_branches(num_nodes, Matrix::Zero(S, S));

        bool use_custom_bl = !branch_lengths.empty();
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                int32_t bid = node.id;
                bool is_test = (!has_partition || (*is_test_branch)[bid]);
                const MG94Matrix& mg_b = is_test ? mg_test : mg_bg;
                Scalar conv_b = is_test ? conv_test : conv_bg;
                Scalar bl = (use_custom_bl && static_cast<size_t>(bid) < branch_lengths.size())
                    ? branch_lengths[bid] : node.branch_length;
                tau[bid] = bl * conv_b;
                P_branches[bid] = mg_b.transition_matrix(tau[bid]);
                dP_branches[bid] = mg_b.Q * P_branches[bid];
            }
        }

        std::vector<std::array<Matrix, 3>> P_bridges;
        std::vector<std::array<Matrix, 3>> dP_bridges;
        const std::vector<std::array<Matrix, 3>>* p_bridges_ptr = nullptr;
        if (!tree.bridges.empty()) {
            P_bridges.resize(tree.bridges.size());
            dP_bridges.resize(tree.bridges.size());
            for (size_t k = 0; k < tree.bridges.size(); ++k) {
                const auto& wb = tree.bridges[k];
                int32_t bid = (wb.bridge_node_id >= 0 && wb.bridge_node_id != wb.hub_node_id) ? wb.bridge_node_id : wb.hub_node_id;
                bool is_test = (!has_partition || (bid >= 0 && static_cast<size_t>(bid) < is_test_branch->size() && (*is_test_branch)[bid]));
                const MG94Matrix& mg_b = is_test ? mg_test : mg_bg;
                Scalar conv_b = is_test ? conv_test : conv_bg;
                for (int r = 0; r < 3; ++r) {
                    if (wb.weights[r] > 0.0 || r == 0) {
                        P_bridges[k][r] = mg_b.transition_matrix(wb.branch_lengths[r] * conv_b);
                        dP_bridges[k][r] = mg_b.Q * P_bridges[k][r];
                    } else {
                        P_bridges[k][r] = Matrix::Identity(S, S);
                        dP_bridges[k][r] = Matrix::Zero(S, S);
                    }
                }
            }
            p_bridges_ptr = &P_bridges;
        }

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        std::vector<Scalar> pattern_ll(num_patterns, 0.0);
        Matrix pattern_dL = Matrix::Zero(num_nodes, num_patterns);
        // Pack V and D across patterns for each branch so A_branches[bid] = V_mat[bid] * D_mat[bid]^T is a single deterministic BLAS GEMM!
        std::vector<Matrix> V_mat(num_nodes, Matrix::Zero(S, num_patterns));
        std::vector<Matrix> D_mat(num_nodes, Matrix::Zero(S, num_patterns));
        size_t num_bridge_res = 3 * tree.bridges.size();
        Matrix pattern_dL_bridge = Matrix::Zero(num_bridge_res, num_patterns);
        std::vector<Matrix> V_bridge_mat(num_bridge_res, Matrix::Zero(S, num_patterns));
        std::vector<Matrix> D_bridge_mat(num_bridge_res, Matrix::Zero(S, num_patterns));
        int n_threads = optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            const double SCALER_UP = 1e150;
            const double SCALER_DOWN = 1e-150;
            BatchedInsideOutsideWorkspace ws;
            ws.init(num_nodes, S, tile_size);
            Matrix dP_D(S, tile_size);
            Matrix D_br(S, tile_size);
            Matrix V_br(S, tile_size);
            Vector inv_L_vec(tile_size);

            #pragma omp for schedule(static)
            for (size_t t = 0; t < num_tiles; ++t) {
                size_t p_start = t * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);

                compute_inside_outside_batch_in_place(
                    tree, aln, p_start, p_end, leaf_to_taxon, P_branches, aln.codon_frequencies_f3x4, S, ws, p_bridges_ptr
                );

                for (Eigen::Index b = 0; b < B; ++b) {
                    size_t p = p_start + b;
                    Scalar L = ws.likelihood(b);
                    if (L > 0.0) {
                        pattern_ll[p] = aln.patterns[p].weight * ws.log_likelihood(b);
                        inv_L_vec(b) = aln.patterns[p].weight / L;
                    } else {
                        inv_L_vec(b) = 0.0;
                    }
                }

                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                        int32_t bid = node.id;
                        auto D_sub = ws.D[bid].leftCols(B);
                        auto V_sub = ws.V[bid].leftCols(B);
                        auto dP_D_sub = dP_D.leftCols(B);
                        dP_D_sub.noalias() = dP_branches[bid] * D_sub;
                        for (Eigen::Index b = 0; b < B; ++b) {
                            Scalar w_inv_L = inv_L_vec(b);
                            if (w_inv_L > 0.0) {
                                size_t p = p_start + b;
                                pattern_dL(bid, p) = w_inv_L * V_sub.col(b).dot(dP_D_sub.col(b));
                                V_mat[bid].col(p) = w_inv_L * V_sub.col(b);
                                D_mat[bid].col(p) = D_sub.col(b);
                            }
                        }
                    }
                }

                for (size_t k = 0; k < tree.bridges.size(); ++k) {
                    const auto& wb = tree.bridges[k];
                    int32_t uid = wb.hub_node_id;
                    int32_t id_A = wb.children[0];
                    int32_t id_B = wb.children[1];
                    int32_t id_C = wb.children[2];
                    auto M_A = ws.P_D[id_A].leftCols(B);
                    auto M_B = ws.P_D[id_B].leftCols(B);
                    auto M_C = ws.P_D[id_C].leftCols(B);
                    auto U_uid = ws.U[uid].leftCols(B);
                    auto D_br_sub = D_br.leftCols(B);
                    auto V_br_sub = V_br.leftCols(B);
                    auto dP_D_sub = dP_D.leftCols(B);

                    for (int r = 0; r < 3; ++r) {
                        if (wb.weights[r] <= 0.0) continue;
                        size_t kr = 3 * k + r;
                        if (r == 0) {
                            D_br_sub.array() = M_A.array() * M_B.array();
                            V_br_sub.array() = wb.weights[0] * (U_uid.array() * M_C.array());
                        } else if (r == 1) {
                            D_br_sub.array() = M_A.array() * M_C.array();
                            V_br_sub.array() = wb.weights[1] * (U_uid.array() * M_B.array());
                        } else {
                            D_br_sub.array() = M_B.array() * M_C.array();
                            V_br_sub.array() = wb.weights[2] * (U_uid.array() * M_A.array());
                        }
                        for (Eigen::Index b = 0; b < B; ++b) {
                            int net_boost = ws.boost_counts(uid, b);
                            while (net_boost > 0) {
                                V_br_sub.col(b) *= SCALER_UP;
                                net_boost--;
                            }
                            while (net_boost < 0) {
                                V_br_sub.col(b) *= SCALER_DOWN;
                                net_boost++;
                            }
                        }
                        dP_D_sub.noalias() = dP_bridges[k][r] * D_br_sub;
                        for (Eigen::Index b = 0; b < B; ++b) {
                            Scalar w_inv_L = inv_L_vec(b);
                            if (w_inv_L > 0.0) {
                                size_t p = p_start + b;
                                pattern_dL_bridge(kr, p) = w_inv_L * V_br_sub.col(b).dot(dP_D_sub.col(b));
                                V_bridge_mat[kr].col(p) = w_inv_L * V_br_sub.col(b);
                                D_bridge_mat[kr].col(p) = D_br_sub.col(b);
                            }
                        }
                    }
                }
            }
        }

        Scalar total_log_l = neumaier_sum(pattern_ll);
        std::vector<Scalar> dL_dtau(num_nodes, 0.0);
        std::vector<Matrix> A_branches(num_nodes);
        std::vector<Scalar> row_terms(num_patterns, 0.0);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                int32_t bid = node.id;
                for (size_t p = 0; p < num_patterns; ++p) {
                    row_terms[p] = pattern_dL(bid, p);
                }
                dL_dtau[bid] = neumaier_sum(row_terms);
                A_branches[bid].noalias() = V_mat[bid] * D_mat[bid].transpose();
            }
        }

        Matrix H_test = Matrix::Zero(S, S);
        Scalar sum_tau_dL_dtau_test = 0.0;
        Matrix H_bg = Matrix::Zero(S, S);
        Scalar sum_tau_dL_dtau_bg = 0.0;

        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                int32_t bid = node.id;
                bool is_test = (!has_partition || (*is_test_branch)[bid]);
                const MG94Matrix& mg_b = is_test ? mg_test : mg_bg;
                Matrix& H_target = is_test ? H_test : H_bg;
                Scalar& sum_target = is_test ? sum_tau_dL_dtau_test : sum_tau_dL_dtau_bg;

                Scalar tb = tau[bid];
                sum_target += tb * dL_dtau[bid];

                Matrix A_tilde = mg_b.U_pi.transpose() * A_branches[bid] * mg_b.V_pi.transpose();

                Matrix E(S, S);
                for (int i = 0; i < S; ++i) {
                    for (int j = 0; j < S; ++j) {
                        Scalar diff = mg_b.eigenvalues(i) - mg_b.eigenvalues(j);
                        if (std::abs(diff) < 1e-9) {
                            E(i, j) = tb * std::exp(mg_b.eigenvalues(i) * tb);
                        } else {
                            E(i, j) = (std::exp(mg_b.eigenvalues(i) * tb) - std::exp(mg_b.eigenvalues(j) * tb)) / diff;
                        }
                    }
                }

                H_target.noalias() += A_tilde.cwiseProduct(E);
            }
        }

        for (size_t k = 0; k < tree.bridges.size(); ++k) {
            const auto& wb = tree.bridges[k];
            int32_t bid = (wb.bridge_node_id >= 0 && wb.bridge_node_id != wb.hub_node_id) ? wb.bridge_node_id : wb.hub_node_id;
            bool is_test = (!has_partition || (bid >= 0 && static_cast<size_t>(bid) < is_test_branch->size() && (*is_test_branch)[bid]));
            const MG94Matrix& mg_b = is_test ? mg_test : mg_bg;
            Scalar conv_b = is_test ? conv_test : conv_bg;
            Matrix& H_target = is_test ? H_test : H_bg;
            Scalar& sum_target = is_test ? sum_tau_dL_dtau_test : sum_tau_dL_dtau_bg;

            for (int r = 0; r < 3; ++r) {
                if (wb.weights[r] <= 0.0) continue;
                size_t kr = 3 * k + r;
                for (size_t p = 0; p < num_patterns; ++p) {
                    row_terms[p] = pattern_dL_bridge(kr, p);
                }
                Scalar dL_dtau_kr = neumaier_sum(row_terms);
                Scalar tb = wb.branch_lengths[r] * conv_b;
                sum_target += tb * dL_dtau_kr;

                Matrix A_kr = V_bridge_mat[kr] * D_bridge_mat[kr].transpose();
                Matrix A_tilde = mg_b.U_pi.transpose() * A_kr * mg_b.V_pi.transpose();

                Matrix E(S, S);
                for (int i = 0; i < S; ++i) {
                    for (int j = 0; j < S; ++j) {
                        Scalar diff = mg_b.eigenvalues(i) - mg_b.eigenvalues(j);
                        if (std::abs(diff) < 1e-9) {
                            E(i, j) = tb * std::exp(mg_b.eigenvalues(i) * tb);
                        } else {
                            E(i, j) = (std::exp(mg_b.eigenvalues(i) * tb) - std::exp(mg_b.eigenvalues(j) * tb)) / diff;
                        }
                    }
                }

                H_target.noalias() += A_tilde.cwiseProduct(E);
            }
        }

        auto compute_part_grad = [&](const MG94Matrix& mg_part, Scalar mu_part, const Matrix& H_part, Scalar sum_part, const Matrix& dQ_k) -> Scalar {
            Matrix M_k = mg_part.V_pi * dQ_k * mg_part.U_pi;
            Scalar part1 = (M_k.cwiseProduct(H_part)).sum();
            Scalar dmu_k = 0.0;
            for (int i = 0; i < S; ++i) {
                dmu_k -= mg_part.pi(i) * dQ_k(i, i);
            }
            Scalar part2 = -(dmu_k / mu_part) * sum_part;
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

        auto build_dQ_nuc = [&](int u, int v, Scalar beta_val) -> Matrix {
            Matrix dQ = Matrix::Zero(S, S);
            for (int i = 0; i < S; ++i) {
                for (int j = 0; j < S; ++j) {
                    if (i == j) continue;
                    const auto& diff = gcode.diff_matrix[i][j];
                    if (diff.num_diffs == 1) {
                        int from = diff.diff_from[0];
                        int to = diff.diff_to[0];
                        if ((from == u && to == v) || (from == v && to == u)) {
                            Scalar rate_mod = diff.synonymous ? params.alpha : beta_val;
                            Scalar target_nuc_freq = aln.pos_nuc_frequencies(diff.diff_pos[0], to);
                            dQ(i, j) = rate_mod * 1.0 * target_nuc_freq;
                        }
                    }
                }
                dQ(i, i) = -dQ.row(i).sum();
            }
            return dQ;
        };

        auto compute_theta_grad = [&](int u, int v) -> Scalar {
            Scalar g = compute_part_grad(mg_test, mu_test, H_test, sum_tau_dL_dtau_test, build_dQ_nuc(u, v, params.beta));
            if (has_partition) {
                g += compute_part_grad(mg_bg, mu_bg, H_bg, sum_tau_dL_dtau_bg, build_dQ_nuc(u, v, omega_bg));
            }
            return g;
        };

        JointMG94Gradients res;
        res.log_likelihood = total_log_l;
        res.grad_branches.resize(num_nodes, 0.0);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                bool is_test = (!has_partition || (*is_test_branch)[node.id]);
                Scalar conv_b = is_test ? conv_test : conv_bg;
                res.grad_branches[node.id] = conv_b * dL_dtau[node.id];
            }
        }

        res.grad_omega = compute_part_grad(mg_test, mu_test, H_test, sum_tau_dL_dtau_test, dQ_omega);
        if (has_partition) {
            res.grad_omega_bg = compute_part_grad(mg_bg, mu_bg, H_bg, sum_tau_dL_dtau_bg, dQ_omega);
        }
        res.grad_theta_AC = compute_theta_grad(0, 1);
        res.grad_theta_AT = compute_theta_grad(0, 3);
        res.grad_theta_CG = compute_theta_grad(1, 2);
        res.grad_theta_CT = compute_theta_grad(1, 3);
        res.grad_theta_GT = compute_theta_grad(2, 3);

        return res;
    }

    /**
     * Compute stationary distribution for ExpCM at a single site.
     *
     * p_{r,x} = (π_{r,A(x)})^β · φ_{x1} φ_{x2} φ_{x3} / Z_r
     */
    static Vector61Codon compute_expcm_stationary(
        const GeneticCode& gcode,
        const Preferences& prefs,
        size_t site,
        Scalar beta,
        const Vector4& phi
    ) {
        Vector61Codon p = Vector61Codon::Zero();
        Scalar Z = 0.0;

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            char aa = gcode.sense_codon_aa[x];
            Scalar pref_beta = std::pow(prefs.at(site, aa), beta);

            const std::string& codon = gcode.sense_codons[x];
            int8_t n1 = char_to_nuc(codon[0]);
            int8_t n2 = char_to_nuc(codon[1]);
            int8_t n3 = char_to_nuc(codon[2]);
            Scalar nuc_freq = phi(n1) * phi(n2) * phi(n3);

            p(x) = pref_beta * nuc_freq;
            Z += p(x);
        }

        if (Z > 0.0) {
            p /= Z;
        }

        return p;
    }

    /**
     * Compute analytical derivative of site stationary distribution p_r w.r.t. β:
     *
     * ∂p_{r,x} / ∂β = p_{r,x} · [ ln(π_{r,A(x)}) - Σ_z p_{r,z} ln(π_{r,A(z)}) ]
     */
    static Vector61Codon compute_dp_dbeta(
        const GeneticCode& gcode,
        const Preferences& prefs,
        size_t site,
        const Vector61Codon& pi_site
    ) {
        Vector61Codon ln_pref = Vector61Codon::Zero();
        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            char aa = gcode.sense_codon_aa[x];
            ln_pref(x) = std::log(prefs.at(site, aa));
        }
        Scalar mean_ln_pref = pi_site.dot(ln_pref);
        return pi_site.array() * (ln_pref.array() - mean_ln_pref);
    }

    /**
     * Compute ExpCM log-likelihood at fixed parameters using per-site reversible
     * spectral decomposition.
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
        size_t num_sites = prefs.num_sites();
        std::vector<size_t> leaf_to_taxon = build_leaf_to_taxon_map(tree, aln);

        ExpCMQMatrix Q_matrix;
        Q_matrix.build(gcode, params);
        const Matrix61Codon& Q = Q_matrix.Q;

        std::vector<Scalar> site_lls(num_sites, 0.0);

        #pragma omp parallel for schedule(dynamic)
        for (size_t site = 0; site < num_sites; ++site) {
            size_t pattern_idx = aln.site_to_pattern[site];
            const auto& pattern = aln.patterns[pattern_idx];

            ExpCMFMatrix F_matrix;
            F_matrix.build(gcode, prefs, site, params);

            Matrix61Codon P = F_matrix.F.cwiseProduct(Q);
            for (int x = 0; x < 61; ++x) {
                P(x, x) = -P.row(x).sum();
            }

            Vector61Codon pi_site = compute_expcm_stationary(
                gcode, prefs, site, params.beta, params.phi
            );

            ExpCMSiteEigen site_eigen;
            site_eigen.decompose(P, pi_site);

            std::vector<Matrix61Codon> P_branches(num_nodes);
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar scaled_bl = node.branch_length * branch_scale;
                    P_branches[node.id] = site_eigen.transition_matrix(scaled_bl);
                }
            }

            std::vector<Vector61Codon> node_L(num_nodes);
            std::vector<double> node_log_scale(num_nodes, 0.0);

            for (int32_t node_id : tree.post_order) {
                const auto& node = tree.nodes[node_id];

                if (node.is_leaf) {
                    size_t t_idx = leaf_to_taxon[node_id];
                    init_leaf_likelihood(pattern, t_idx, 61, node_L[node_id]);
                    node_log_scale[node_id] = 0.0;
                } else {
                    double child_scale_sum = 0.0;
                    tree.eval_internal_node_with_model(node_id, P_branches, site_eigen, node_L, node_log_scale, node_L[node_id], child_scale_sum, branch_scale);
                    node_log_scale[node_id] = child_scale_sum;
                }
            }

            Scalar site_likelihood = pi_site.dot(node_L[tree.root_id]);
            if (site_likelihood > 0.0) {
                site_lls[site] = std::log(site_likelihood) + node_log_scale[tree.root_id];
            } else {
                site_lls[site] = -1e20;
            }
        }

        return neumaier_sum(site_lls);
    }

    /**
     * Compute log-likelihood for a single site using reversible spectral decomposition.
     */
    static Scalar compute_expcm_site_log_likelihood(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        size_t site,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0,
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr
    ) {
        params.validate();
        if (site >= prefs.num_sites()) {
            throw std::out_of_range("Site index " + std::to_string(site) +
                                   " exceeds num_sites " + std::to_string(prefs.num_sites()));
        }

        size_t num_nodes = tree.num_nodes();
        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = build_leaf_to_taxon_map(tree, aln));

        ExpCMQMatrix Q_matrix;
        Q_matrix.build(gcode, params);
        const Matrix61Codon& Q = Q_matrix.Q;

        size_t pattern_idx = aln.site_to_pattern[site];
        const auto& pattern = aln.patterns[pattern_idx];

        ExpCMFMatrix F_matrix;
        F_matrix.build(gcode, prefs, site, params);

        Matrix61Codon P = F_matrix.F.cwiseProduct(Q);
        for (int x = 0; x < 61; ++x) {
            P(x, x) = -P.row(x).sum();
        }

        Vector61Codon pi_site = compute_expcm_stationary(
            gcode, prefs, site, params.beta, params.phi
        );

        ExpCMSiteEigen site_eigen;
        site_eigen.decompose(P, pi_site);

        std::vector<Matrix61Codon> P_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar scaled_bl = node.branch_length * branch_scale;
                P_branches[node.id] = site_eigen.transition_matrix(scaled_bl);
            }
        }

        std::vector<Vector61Codon> node_L(num_nodes);
        std::vector<double> node_log_scale(num_nodes, 0.0);

        for (int32_t node_id : tree.post_order) {
            const auto& node = tree.nodes[node_id];

            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[node_id];
                init_leaf_likelihood(pattern, t_idx, 61, node_L[node_id]);
                node_log_scale[node_id] = 0.0;
            } else {
                double child_scale_sum = 0.0;
                tree.eval_internal_node_with_model(node_id, P_branches, site_eigen, node_L, node_log_scale, node_L[node_id], child_scale_sum, branch_scale);
                node_log_scale[node_id] = child_scale_sum;
            }
        }

        Scalar site_likelihood = pi_site.dot(node_L[tree.root_id]);
        if (site_likelihood > 0.0) {
            return std::log(site_likelihood) + node_log_scale[tree.root_id];
        } else {
            return -1e20;
        }
    }

    /**
     * Compute log-likelihood for all sites individually.
     */
    static std::vector<Scalar> compute_expcm_all_sites_log_likelihood(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0
    ) {
        size_t num_sites = prefs.num_sites();
        std::vector<Scalar> site_likelihoods(num_sites);
        std::vector<size_t> leaf_to_taxon = build_leaf_to_taxon_map(tree, aln);

        #pragma omp parallel for schedule(dynamic)
        for (size_t site = 0; site < num_sites; ++site) {
            site_likelihoods[site] = compute_expcm_site_log_likelihood(
                tree, aln, gcode, prefs, site, params, branch_scale, &leaf_to_taxon
            );
        }

        return site_likelihoods;
    }

    /**
     * Compute central finite-difference gradient of log-likelihood w.r.t. β.
     */
    static Scalar compute_expcm_beta_gradient_finite_diff(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0,
        Scalar epsilon = 1e-6
    ) {
        ExpCMParameters params_plus = params;
        params_plus.beta += epsilon;
        Scalar L_plus = compute_expcm_log_likelihood(tree, aln, gcode, prefs, params_plus, branch_scale);

        ExpCMParameters params_minus = params;
        params_minus.beta -= epsilon;
        Scalar L_minus = compute_expcm_log_likelihood(tree, aln, gcode, prefs, params_minus, branch_scale);

        return (L_plus - L_minus) / (2.0 * epsilon);
    }

    /**
     * Compute dF/dβ matrix (derivative of fixation bias w.r.t. β) for a single site.
     */
    static Matrix61Codon compute_dF_dbeta(
        const GeneticCode& gcode,
        const Preferences& prefs,
        size_t site,
        Scalar beta,
        Scalar omega
    ) {
        Matrix61Codon dF_dbeta = Matrix61Codon::Zero();
        const double EXPCM_ALMOST_ZERO_VAL = 1e-10;

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            for (int y = 0; y < gcode.num_sense_codons; ++y) {
                if (x == y) continue;

                char aa_x = gcode.sense_codon_aa[x];
                char aa_y = gcode.sense_codon_aa[y];

                if (aa_x == aa_y) {
                    dF_dbeta(x, y) = 0.0;
                } else {
                    Scalar pi_x = prefs.at(site, aa_x);
                    Scalar pi_y = prefs.at(site, aa_y);

                    Scalar ratio = pi_x / pi_y;
                    Scalar ratio_beta = std::pow(ratio, beta);
                    Scalar ln_ratio = std::log(ratio);

                    if (std::fabs(1.0 - ratio_beta) < EXPCM_ALMOST_ZERO_VAL) {
                        dF_dbeta(x, y) = omega * (-ln_ratio / 2.0);
                    } else {
                        Scalar denominator = (1.0 - ratio_beta);
                        Scalar denom_squared = denominator * denominator;
                        Scalar numerator = -ln_ratio * (1.0 - ratio_beta + beta * ratio_beta * ln_ratio);
                        dF_dbeta(x, y) = omega * (numerator / denom_squared);
                    }
                }
            }
        }

        dF_dbeta.diagonal().setZero();
        return dF_dbeta;
    }

    /**
     * Compute dF/dω matrix (derivative of fixation bias w.r.t. ω) for a single site.
     */
    static Matrix61Codon compute_dF_domega(
        const GeneticCode& gcode,
        const Preferences& prefs,
        size_t site,
        Scalar beta
    ) {
        Matrix61Codon dF_domega = Matrix61Codon::Zero();
        const double EXPCM_ALMOST_ZERO_VAL = 1e-10;

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            for (int y = 0; y < gcode.num_sense_codons; ++y) {
                if (x == y) continue;

                char aa_x = gcode.sense_codon_aa[x];
                char aa_y = gcode.sense_codon_aa[y];

                if (aa_x == aa_y) {
                    dF_domega(x, y) = 0.0;
                } else {
                    Scalar pi_x = prefs.at(site, aa_x);
                    Scalar pi_y = prefs.at(site, aa_y);

                    Scalar ratio = pi_x / pi_y;
                    Scalar ratio_beta = std::pow(ratio, beta);
                    Scalar ln_ratio_beta = std::log(ratio_beta);

                    if (std::fabs(1.0 - ratio_beta) < EXPCM_ALMOST_ZERO_VAL) {
                        dF_domega(x, y) = 1.0;
                    } else {
                        dF_domega(x, y) = -ln_ratio_beta / (1.0 - ratio_beta);
                    }
                }
            }
        }

        dF_domega.diagonal().setZero();
        return dF_domega;
    }

    /**
     * Check if nucleotide substitution is a transition (A↔G or C↔T).
     */
    static bool is_transition(int8_t nuc_from, int8_t nuc_to) {
        return (nuc_from == 0 && nuc_to == 2) ||
               (nuc_from == 2 && nuc_to == 0) ||
               (nuc_from == 1 && nuc_to == 3) ||
               (nuc_from == 3 && nuc_to == 1);
    }

    /**
     * Compute dQ/dκ matrix (derivative of mutation rates w.r.t. κ).
     */
    static Matrix61Codon compute_dQ_dkappa(
        const GeneticCode& gcode,
        const ExpCMParameters& params
    ) {
        Matrix61Codon dQ_dkappa = Matrix61Codon::Zero();

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            for (int y = 0; y < gcode.num_sense_codons; ++y) {
                if (x == y) continue;

                const auto& diff = gcode.diff_matrix[x][y];
                if (diff.num_diffs != 1) continue;

                int8_t nuc_to = diff.nuc_to;
                int8_t nuc_from = diff.nuc_from;

                if (is_transition(nuc_from, nuc_to)) {
                    dQ_dkappa(x, y) = params.phi(nuc_to);
                }
            }
        }

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            dQ_dkappa(x, x) = -dQ_dkappa.row(x).sum();
        }

        return dQ_dkappa;
    }

    struct ExpCMGradients {
        Scalar log_likelihood = 0.0;
        Scalar grad_beta = 0.0;
        Scalar grad_omega = 0.0;
        Scalar grad_kappa = 0.0;
    };

    /**
     * Jointly evaluate ExpCM log-likelihood and exact analytical gradients
     * (∂lnL/∂β, ∂lnL/∂ω, ∂lnL/∂κ) in a single pass over sites using
     * reversible spectral decomposition + 2-pass Inside-Outside + spectral Fréchet derivatives.
     */
    static ExpCMGradients compute_expcm_joint_frechet_gradients(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0
    ) {
        params.validate();

        size_t num_nodes = tree.num_nodes();
        size_t num_sites = prefs.num_sites();
        std::vector<size_t> leaf_to_taxon = build_leaf_to_taxon_map(tree, aln);

        ExpCMQMatrix Q_matrix;
        Q_matrix.build(gcode, params);
        const Matrix61Codon& Q = Q_matrix.Q;
        Matrix61Codon dQ_dkappa = compute_dQ_dkappa(gcode, params);

        std::vector<Scalar> site_ll(num_sites, 0.0);
        std::vector<Scalar> site_g_beta(num_sites, 0.0);
        std::vector<Scalar> site_g_omega(num_sites, 0.0);
        std::vector<Scalar> site_g_kappa(num_sites, 0.0);

        #pragma omp parallel
        {
            InsideOutsideWorkspace ws;
            ws.init(num_nodes, 61);
            std::vector<Matrix> P_branches(num_nodes, Matrix::Zero(61, 61));
            std::vector<std::array<Matrix, 3>> P_bridges(tree.bridges.size());

            #pragma omp for schedule(dynamic)
            for (size_t site = 0; site < num_sites; ++site) {
                size_t pattern_idx = aln.site_to_pattern[site];
                const auto& pattern = aln.patterns[pattern_idx];

                ExpCMFMatrix F_matrix;
                F_matrix.build(gcode, prefs, site, params);
                const Matrix61Codon& F = F_matrix.F;

                Matrix61Codon P = F.cwiseProduct(Q);
                for (int x = 0; x < 61; ++x) {
                    P(x, x) = -P.row(x).sum();
                }

                Vector61Codon pi_site = compute_expcm_stationary(
                    gcode, prefs, site, params.beta, params.phi
                );

                ExpCMSiteEigen site_eigen;
                site_eigen.decompose(P, pi_site);

                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) {
                        Scalar tb = node.branch_length * branch_scale;
                        P_branches[node.id] = site_eigen.transition_matrix(tb);
                    }
                }
                for (size_t b = 0; b < tree.bridges.size(); ++b) {
                    for (int r = 0; r < 3; ++r) {
                        Scalar tb = tree.bridges[b].branch_lengths[r] * branch_scale;
                        P_bridges[b][r] = site_eigen.transition_matrix(tb);
                    }
                }

                Vector pi_dyn = pi_site;
                compute_inside_outside_in_place(
                    tree, pattern, leaf_to_taxon, P_branches, pi_dyn, 61, ws,
                    tree.bridges.empty() ? nullptr : &P_bridges
                );

                Scalar L_r = ws.likelihood;
                site_ll[site] = ws.log_likelihood;
                if (L_r <= 0.0) continue;

                Scalar inv_L = 1.0 / L_r;

                // Accumulate spectral Fréchet sensitivity matrix H_r across all branches and bridges
                Matrix61Codon H_r = Matrix61Codon::Zero();
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                        int32_t bid = node.id;
                        Scalar tb = node.branch_length * branch_scale;

                        Vector61Codon u_v = site_eigen.U_pi.transpose() * ws.V[bid];
                        Vector61Codon v_d = site_eigen.V_pi * ws.D[bid];
                        Matrix61Codon A_tilde = inv_L * (u_v * v_d.transpose());

                        Matrix61Codon E;
                        for (int i = 0; i < 61; ++i) {
                            Scalar li = site_eigen.eigenvalues(i);
                            Scalar exp_li = std::exp(li * tb);
                            for (int j = 0; j < 61; ++j) {
                                Scalar lj = site_eigen.eigenvalues(j);
                                Scalar diff = li - lj;
                                if (std::fabs(diff) < 1e-9) {
                                    E(i, j) = tb * exp_li;
                                } else {
                                    E(i, j) = (exp_li - std::exp(lj * tb)) / diff;
                                }
                            }
                        }

                        H_r.noalias() += A_tilde.cwiseProduct(E);
                    }
                }
                for (size_t b = 0; b < tree.bridges.size(); ++b) {
                    const auto& wb = tree.bridges[b];
                    int32_t uid = wb.hub_node_id;
                    const Vector& M_A = ws.P_D[wb.children[0]];
                    const Vector& M_B = ws.P_D[wb.children[1]];
                    const Vector& M_C = ws.P_D[wb.children[2]];
                    const Vector& U_uid = ws.U[uid];
                    int total_boost = ws.boost_counts[uid];

                    for (int r = 0; r < 3; ++r) {
                        if (wb.weights[r] <= 0.0) continue;
                        Scalar tb = wb.branch_lengths[r] * branch_scale;

                        Vector D_br = (r == 0) ? M_A.cwiseProduct(M_B)
                                    : ((r == 1) ? M_A.cwiseProduct(M_C) : M_B.cwiseProduct(M_C));
                        Vector V_br = wb.weights[r] * U_uid.cwiseProduct(
                            (r == 0) ? M_C : ((r == 1) ? M_B : M_A));
                        int net_boost = total_boost;
                        while (net_boost > 0) {
                            V_br *= 1e150;
                            net_boost--;
                        }
                        while (net_boost < 0) {
                            V_br *= 1e-150;
                            net_boost++;
                        }

                        Vector61Codon u_v = site_eigen.U_pi.transpose() * V_br;
                        Vector61Codon v_d = site_eigen.V_pi * D_br;
                        Matrix61Codon A_tilde = inv_L * (u_v * v_d.transpose());

                        Matrix61Codon E;
                        for (int i = 0; i < 61; ++i) {
                            Scalar li = site_eigen.eigenvalues(i);
                            Scalar exp_li = std::exp(li * tb);
                            for (int j = 0; j < 61; ++j) {
                                Scalar lj = site_eigen.eigenvalues(j);
                                Scalar diff = li - lj;
                                if (std::fabs(diff) < 1e-9) {
                                    E(i, j) = tb * exp_li;
                                } else {
                                    E(i, j) = (exp_li - std::exp(lj * tb)) / diff;
                                }
                            }
                        }

                        H_r.noalias() += A_tilde.cwiseProduct(E);
                    }
                }

                // Back-project H_r to codon space: G_r = U_pi * H_r * V_pi
                // so that for any dP_r, sum( (V_pi * dP_r * U_pi) ⊙ H_r ) = sum( dP_r ⊙ (V_pi^T * H_r * U_pi^T) )
                Matrix61Codon W_r = site_eigen.V_pi.transpose() * H_r * site_eigen.U_pi.transpose();

                // 1. Beta gradient: includes branch Fréchet term + root stationary derivative term
                Matrix61Codon dF_dbeta = compute_dF_dbeta(gcode, prefs, site, params.beta, params.omega);
                Matrix61Codon dP_dbeta = dF_dbeta.cwiseProduct(Q);
                for (int x = 0; x < 61; ++x) {
                    dP_dbeta(x, x) = -dP_dbeta.row(x).sum();
                }
                Vector61Codon dp_dbeta = compute_dp_dbeta(gcode, prefs, site, pi_site);
                Scalar root_term_beta = dp_dbeta.dot(ws.D[tree.root_id]) * inv_L;
                site_g_beta[site] = root_term_beta + (dP_dbeta.cwiseProduct(W_r)).sum();

                // 2. Omega gradient
                Matrix61Codon dF_domega = compute_dF_domega(gcode, prefs, site, params.beta);
                Matrix61Codon dP_domega = dF_domega.cwiseProduct(Q);
                for (int x = 0; x < 61; ++x) {
                    dP_domega(x, x) = -dP_domega.row(x).sum();
                }
                site_g_omega[site] = (dP_domega.cwiseProduct(W_r)).sum();

                // 3. Kappa gradient
                Matrix61Codon dP_dkappa = F.cwiseProduct(dQ_dkappa);
                for (int x = 0; x < 61; ++x) {
                    dP_dkappa(x, x) = -dP_dkappa.row(x).sum();
                }
                site_g_kappa[site] = (dP_dkappa.cwiseProduct(W_r)).sum();
            }
        }

        ExpCMGradients res;
        res.log_likelihood = neumaier_sum(site_ll);
        res.grad_beta = neumaier_sum(site_g_beta);
        res.grad_omega = neumaier_sum(site_g_omega);
        res.grad_kappa = neumaier_sum(site_g_kappa);
        return res;
    }

    /**
     * Compute analytical beta gradient using Inside-Outside and spectral Fréchet derivatives.
     */
    static Scalar compute_expcm_beta_gradient_frechet(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0
    ) {
        return compute_expcm_joint_frechet_gradients(tree, aln, gcode, prefs, params, branch_scale).grad_beta;
    }

    /**
     * Compute analytical beta gradient (alias to exact spectral Fréchet derivative).
     */
    static Scalar compute_expcm_beta_gradient_analytical(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0
    ) {
        return compute_expcm_beta_gradient_frechet(tree, aln, gcode, prefs, params, branch_scale);
    }

    /**
     * Compute central finite-difference gradient of log-likelihood w.r.t. ω.
     */
    static Scalar compute_expcm_omega_gradient_finite_diff(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0,
        Scalar epsilon = 1e-6
    ) {
        ExpCMParameters params_plus = params;
        params_plus.omega += epsilon;
        Scalar L_plus = compute_expcm_log_likelihood(tree, aln, gcode, prefs, params_plus, branch_scale);

        ExpCMParameters params_minus = params;
        params_minus.omega -= epsilon;
        Scalar L_minus = compute_expcm_log_likelihood(tree, aln, gcode, prefs, params_minus, branch_scale);

        return (L_plus - L_minus) / (2.0 * epsilon);
    }

    /**
     * Compute analytical omega gradient using Inside-Outside and spectral Fréchet derivatives.
     */
    static Scalar compute_expcm_omega_gradient_frechet(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0
    ) {
        return compute_expcm_joint_frechet_gradients(tree, aln, gcode, prefs, params, branch_scale).grad_omega;
    }

    /**
     * Compute central finite-difference gradient of log-likelihood w.r.t. κ.
     */
    static Scalar compute_expcm_kappa_gradient_finite_diff(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0,
        Scalar epsilon = 1e-6
    ) {
        ExpCMParameters params_plus = params;
        params_plus.kappa += epsilon;
        Scalar L_plus = compute_expcm_log_likelihood(tree, aln, gcode, prefs, params_plus, branch_scale);

        ExpCMParameters params_minus = params;
        params_minus.kappa -= epsilon;
        Scalar L_minus = compute_expcm_log_likelihood(tree, aln, gcode, prefs, params_minus, branch_scale);

        return (L_plus - L_minus) / (2.0 * epsilon);
    }

    /**
     * Compute analytical kappa gradient using Inside-Outside and spectral Fréchet derivatives.
     */
    static Scalar compute_expcm_kappa_gradient_frechet(
        const Tree& tree,
        const Alignment& aln,
        const GeneticCode& gcode,
        const Preferences& prefs,
        const ExpCMParameters& params,
        Scalar branch_scale = 1.0
    ) {
        return compute_expcm_joint_frechet_gradients(tree, aln, gcode, prefs, params, branch_scale).grad_kappa;
    }

    /**
     * Compute matrix exponential using scaling and squaring method with Padé (5,5) approximant.
     * Retained for direct matrix-exponential validation tests.
     */
    static Matrix61Codon matrix_exp_scaling_squaring(const Matrix61Codon& M) {
        const int max_scale = 30;
        Scalar norm_inf = M.rowwise().lpNorm<1>().maxCoeff();

        if (norm_inf < 1e-10) {
            return Matrix61Codon::Identity() + M;
        }

        int s = 0;
        Scalar one_half_power = 1.0;
        while (norm_inf > 1.0 && s < max_scale) {
            norm_inf *= 0.5;
            one_half_power *= 0.5;
            s++;
        }

        Matrix61Codon M_scaled = M * one_half_power;
        Matrix61Codon result = pade_exp(M_scaled);

        for (int i = 0; i < s; ++i) {
            result = result * result;
        }

        return result;
    }

    /**
     * Padé (5,5) approximation for matrix exponential.
     */
    static Matrix61Codon pade_exp(const Matrix61Codon& X) {
        Matrix61Codon I = Matrix61Codon::Identity();
        Matrix61Codon X2 = X * X;
        Matrix61Codon X3 = X2 * X;
        Matrix61Codon X4 = X3 * X;
        Matrix61Codon X5 = X4 * X;

        Matrix61Codon num = 30240.0 * I + 15120.0 * X + 3360.0 * X2 + 420.0 * X3 + 30.0 * X4 + X5;
        Matrix61Codon den = 30240.0 * I - 15120.0 * X + 3360.0 * X2 - 420.0 * X3 + 30.0 * X4 - X5;

        return den.inverse() * num;
    }
};

} // namespace hyphy::core
