#pragma once

#include <hyphy/core/types.hpp>
#include <hyphy/core/alignment.hpp>
#include <hyphy/core/tree.hpp>
#include <algorithm>
#include <cstdint>
#include <vector>

namespace hyphy::core::ancestral {

/**
 * Result of a joint Viterbi (max-product) ancestral state reconstruction for a single site pattern.
 */
struct ViterbiReconstructionResult {
    bool is_constant = false;
    int8_t constant_state = -1;
    int8_t root_state = 0;
    std::vector<int8_t> node_states; // Indexed by node_id (0 .. num_nodes - 1)
};

/**
 * Joint maximum a posteriori (Viterbi / max-product) ancestral state reconstruction on a tree.
 *
 * - First checks if all observed (non-ambiguous) leaf states in [0, num_states) are identical;
 *   if so, returns `is_constant = true` and `constant_state = st` immediately without matrix products.
 * - Otherwise executes an upward post-order max-product pass with per-node underflow scaling,
 *   selects the root state maximizing `root_frequencies(c) * M_root(c)`, and traces back
 *   optimal states downward from the root to all nodes.
 */
inline ViterbiReconstructionResult reconstruct_viterbi_states(
    const Tree& tree,
    const SitePattern& pattern,
    const std::vector<size_t>& leaf_to_taxon,
    const std::vector<Matrix>& branch_transition_matrices,
    const Vector& root_frequencies,
    int num_states
) {
    ViterbiReconstructionResult res;
    size_t num_nodes = tree.num_nodes();
    res.node_states.assign(num_nodes, -1);

    // 1. Check if pattern is constant across all observed leaf taxa
    int8_t const_state = -1;
    bool is_constant = true;
    for (size_t nid = 0; nid < num_nodes; ++nid) {
        if (tree.nodes[nid].is_leaf) {
            size_t t_idx = leaf_to_taxon[nid];
            if (t_idx != static_cast<size_t>(-1)) {
                int8_t st = pattern.states[t_idx];
                if (st >= 0 && st < num_states) {
                    if (const_state == -1) {
                        const_state = st;
                    } else if (const_state != st) {
                        is_constant = false;
                        break;
                    }
                }
            }
        }
    }

    if (is_constant && const_state >= 0) {
        res.is_constant = true;
        res.constant_state = const_state;
        res.root_state = const_state;
        std::fill(res.node_states.begin(), res.node_states.end(), const_state);
        return res;
    }

    // 2. Upward max-product pass
    std::vector<std::vector<Scalar>> M_node(num_nodes, std::vector<Scalar>(num_states, 1.0));
    std::vector<std::vector<int8_t>> backpointer(num_nodes, std::vector<int8_t>(num_states, 0));
    std::vector<std::vector<Scalar>> val_node(num_nodes, std::vector<Scalar>(num_states, 1.0));

    for (int32_t nid : tree.post_order) {
        const auto& node = tree.nodes[nid];
        if (node.is_leaf) {
            size_t t_idx = leaf_to_taxon[nid];
            int8_t st = (t_idx != static_cast<size_t>(-1)) ? pattern.states[t_idx] : -1;
            if (st >= 0 && st < num_states) {
                for (int c = 0; c < num_states; ++c) {
                    M_node[nid][c] = (c == st) ? 1.0 : 0.0;
                }
            } else {
                for (int c = 0; c < num_states; ++c) {
                    M_node[nid][c] = 1.0;
                }
            }
        }

        if (nid != tree.root_id) {
            const Matrix& P = branch_transition_matrices[nid];
            Scalar max_val = 0.0;
            for (int p_st = 0; p_st < num_states; ++p_st) {
                Scalar best_score = -1.0;
                int8_t best_c = 0;
                for (int c = 0; c < num_states; ++c) {
                    Scalar score = P(p_st, c) * M_node[nid][c];
                    if (score > best_score) {
                        best_score = score;
                        best_c = static_cast<int8_t>(c);
                    }
                }
                backpointer[nid][p_st] = best_c;
                val_node[nid][p_st] = std::max<Scalar>(best_score, 0.0);
                if (val_node[nid][p_st] > max_val) {
                    max_val = val_node[nid][p_st];
                }
            }

            if (max_val > 0.0) {
                Scalar scale = 1.0 / max_val;
                for (int p_st = 0; p_st < num_states; ++p_st) {
                    val_node[nid][p_st] *= scale;
                }
            }

            int32_t parent_id = node.parent_id;
            Scalar p_max = 0.0;
            for (int p_st = 0; p_st < num_states; ++p_st) {
                M_node[parent_id][p_st] *= val_node[nid][p_st];
                if (M_node[parent_id][p_st] > p_max) {
                    p_max = M_node[parent_id][p_st];
                }
            }
            if (p_max > 0.0 && p_max < 1e-100) {
                Scalar p_scale = 1.0 / p_max;
                for (int p_st = 0; p_st < num_states; ++p_st) {
                    M_node[parent_id][p_st] *= p_scale;
                }
            }
        }
    }

    // 3. Optimal root state
    int8_t root_state = 0;
    Scalar best_root_score = -1.0;
    for (int c = 0; c < num_states; ++c) {
        Scalar score = root_frequencies(c) * M_node[tree.root_id][c];
        if (score > best_root_score) {
            best_root_score = score;
            root_state = static_cast<int8_t>(c);
        }
    }
    res.root_state = root_state;
    res.node_states[tree.root_id] = root_state;

    // 4. Downward traceback
    std::vector<int32_t> stack = {tree.root_id};
    while (!stack.empty()) {
        int32_t curr = stack.back();
        stack.pop_back();
        for (int32_t child_id : tree.nodes[curr].children) {
            int8_t parent_st = res.node_states[curr];
            res.node_states[child_id] = backpointer[child_id][parent_st];
            stack.push_back(child_id);
        }
    }

    return res;
}

} // namespace hyphy::core::ancestral
