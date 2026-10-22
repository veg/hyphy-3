#pragma once

#include "hyphy/core/types.hpp"
#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <regex>
#include <algorithm>
#include <cctype>

namespace hyphy::core {

struct TreeNode {
    int32_t id = INVALID_INDEX;
    std::string name;
    int32_t parent_id = INVALID_INDEX;
    std::vector<int32_t> children;
    Scalar branch_length = 0.0;
    bool is_leaf = false;
    std::string model_tag;
};

struct WheatstoneBridge {
    int32_t hub_node_id = INVALID_INDEX;       // Node u (collapsed 3-child hub: A, B, C)
    int32_t bridge_node_id = INVALID_INDEX;    // Original node v in backbone tree
    std::array<int32_t, 3> children{INVALID_INDEX, INVALID_INDEX, INVALID_INDEX}; // {id_A, id_B, id_C}
    std::array<Scalar, 3> weights{1.0, 0.0, 0.0};        // {w_1, w_2, w_3} summing to 1.0 (Ohmic currents)
    std::array<Scalar, 3> branch_lengths{0.01, 0.01, 0.01}; // {t_1, t_2, t_3} internal quartet bridge lengths
    Scalar ohmic_support = 1.0;                // I_Ohm = G_1 / (G_1 + G_2 + G_3)
    Scalar shunt_stress = 0.0;                 // tau_Ohm = G_2 / G_1
};

class Tree {
public:
    std::vector<TreeNode> nodes;
    int32_t root_id = INVALID_INDEX;
    std::vector<int32_t> post_order;
    std::unordered_map<std::string, int32_t> leaf_name_to_id;
    std::vector<WheatstoneBridge> bridges;
    std::vector<int32_t> node_to_bridge_idx;

    bool has_wheatstone_bridges() const {
        return !bridges.empty();
    }

    int32_t bridge_index(int32_t nid) const {
        return (nid >= 0 && static_cast<size_t>(nid) < node_to_bridge_idx.size())
            ? node_to_bridge_idx[nid]
            : -1;
    }

    void rebuild_bridge_index() {
        node_to_bridge_idx.assign(nodes.size(), -1);
        for (size_t k = 0; k < bridges.size(); ++k) {
            int32_t hid = bridges[k].hub_node_id;
            if (hid >= 0 && static_cast<size_t>(hid) < node_to_bridge_idx.size()) {
                node_to_bridge_idx[hid] = static_cast<int32_t>(k);
            }
        }
    }

    void scale_bridge_lengths(Scalar factor) {
        for (auto& wb : bridges) {
            for (int r = 0; r < 3; ++r) {
                wb.branch_lengths[r] *= factor;
            }
        }
    }

    template <typename VecType, typename MatType>
    void eval_internal_node(
        int32_t nid,
        const std::vector<MatType>& P_branches,
        const std::vector<VecType>& node_L,
        const std::vector<double>& node_log_scale,
        VecType& out_L,
        double& out_child_scale_sum,
        const std::vector<std::array<MatType, 3>>* P_bridges = nullptr,
        double* out_node_boost = nullptr,
        int* out_node_boost_count = nullptr
    ) const {
        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double LOG_SCALER = 345.38776394910684;

        int boost_count = 0;
        double boost_mult = 1.0;

        int32_t b_idx = bridge_index(nid);
        if (b_idx >= 0) {
            const auto& wb = bridges[static_cast<size_t>(b_idx)];
            int32_t id_A = wb.children[0];
            int32_t id_B = wb.children[1];
            int32_t id_C = wb.children[2];

            VecType M_A = P_branches[id_A] * node_L[id_A];
            VecType M_B = P_branches[id_B] * node_L[id_B];
            VecType M_C = P_branches[id_C] * node_L[id_C];
            out_child_scale_sum = node_log_scale[id_A] + node_log_scale[id_B] + node_log_scale[id_C];

            auto get_P_br = [&](int r) -> const MatType& {
                if (P_bridges != nullptr) {
                    return (*P_bridges)[static_cast<size_t>(b_idx)][r];
                }
                int32_t bid = wb.bridge_node_id;
                if (bid >= 0 && bid != nid && static_cast<size_t>(bid) < P_branches.size() &&
                    P_branches[bid].rows() == M_A.size()) {
                    return P_branches[bid];
                }
                if (nid >= 0 && static_cast<size_t>(nid) < P_branches.size() &&
                    P_branches[nid].rows() == M_A.size()) {
                    return P_branches[nid];
                }
                return P_branches[id_A];
            };

            out_L = VecType::Zero(M_A.size());
            if (wb.weights[0] > 0.0) {
                VecType Y1 = M_A.cwiseProduct(M_B);
                VecType D1 = (get_P_br(0) * Y1).cwiseProduct(M_C);
                out_L.noalias() += wb.weights[0] * D1;
            }
            if (wb.weights[1] > 0.0) {
                VecType Y2 = M_A.cwiseProduct(M_C);
                VecType D2 = (get_P_br(1) * Y2).cwiseProduct(M_B);
                out_L.noalias() += wb.weights[1] * D2;
            }
            if (wb.weights[2] > 0.0) {
                VecType Y3 = M_B.cwiseProduct(M_C);
                VecType D3 = (get_P_br(2) * Y3).cwiseProduct(M_A);
                out_L.noalias() += wb.weights[2] * D3;
            }
        } else {
            out_L.setOnes();
            out_child_scale_sum = 0.0;
            for (int32_t cid : nodes[nid].children) {
                out_L.array() *= (P_branches[cid] * node_L[cid]).array();
                out_child_scale_sum += node_log_scale[cid];
                double cur_max = out_L.maxCoeff();
                while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                    out_L *= SCALER_UP;
                    cur_max *= SCALER_UP;
                    boost_mult *= SCALER_UP;
                    boost_count++;
                }
            }
        }

        double max_val = out_L.maxCoeff();
        while (max_val < SCALER_THRESH && max_val > 0.0) {
            out_L *= SCALER_UP;
            max_val *= SCALER_UP;
            boost_mult *= SCALER_UP;
            boost_count++;
        }
        out_child_scale_sum -= boost_count * LOG_SCALER;
        if (out_node_boost != nullptr) {
            *out_node_boost = boost_mult;
        }
        if (out_node_boost_count != nullptr) {
            *out_node_boost_count = boost_count;
        }
    }

    template <typename VecType, typename MatType, typename ModelType>
    void eval_internal_node_with_model(
        int32_t nid,
        const std::vector<MatType>& P_branches,
        const ModelType& model,
        const std::vector<VecType>& node_L,
        const std::vector<double>& node_log_scale,
        VecType& out_L,
        double& out_child_scale_sum,
        Scalar branch_scale = 1.0,
        double* out_node_boost = nullptr,
        const std::vector<std::array<MatType, 3>>* P_bridges = nullptr
    ) const {
        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double LOG_SCALER = 345.38776394910684;

        int boost_count = 0;
        double boost_mult = 1.0;

        int32_t b_idx = bridge_index(nid);
        if (b_idx >= 0) {
            const auto& wb = bridges[static_cast<size_t>(b_idx)];
            int32_t id_A = wb.children[0];
            int32_t id_B = wb.children[1];
            int32_t id_C = wb.children[2];

            VecType M_A = P_branches[id_A] * node_L[id_A];
            VecType M_B = P_branches[id_B] * node_L[id_B];
            VecType M_C = P_branches[id_C] * node_L[id_C];
            out_child_scale_sum = node_log_scale[id_A] + node_log_scale[id_B] + node_log_scale[id_C];

            out_L = VecType::Zero(M_A.size());
            if (wb.weights[0] > 0.0) {
                MatType P1 = (P_bridges != nullptr)
                    ? (*P_bridges)[static_cast<size_t>(b_idx)][0]
                    : model.transition_matrix(wb.branch_lengths[0] * branch_scale);
                VecType Y1 = M_A.cwiseProduct(M_B);
                VecType D1 = (P1 * Y1).cwiseProduct(M_C);
                out_L.noalias() += wb.weights[0] * D1;
            }
            if (wb.weights[1] > 0.0) {
                MatType P2 = (P_bridges != nullptr)
                    ? (*P_bridges)[static_cast<size_t>(b_idx)][1]
                    : model.transition_matrix(wb.branch_lengths[1] * branch_scale);
                VecType Y2 = M_A.cwiseProduct(M_C);
                VecType D2 = (P2 * Y2).cwiseProduct(M_B);
                out_L.noalias() += wb.weights[1] * D2;
            }
            if (wb.weights[2] > 0.0) {
                MatType P3 = (P_bridges != nullptr)
                    ? (*P_bridges)[static_cast<size_t>(b_idx)][2]
                    : model.transition_matrix(wb.branch_lengths[2] * branch_scale);
                VecType Y3 = M_B.cwiseProduct(M_C);
                VecType D3 = (P3 * Y3).cwiseProduct(M_A);
                out_L.noalias() += wb.weights[2] * D3;
            }
        } else {
            out_L.setOnes();
            out_child_scale_sum = 0.0;
            for (int32_t cid : nodes[nid].children) {
                out_L.array() *= (P_branches[cid] * node_L[cid]).array();
                out_child_scale_sum += node_log_scale[cid];
                double cur_max = out_L.maxCoeff();
                while (cur_max < SCALER_THRESH && cur_max > 0.0) {
                    out_L *= SCALER_UP;
                    cur_max *= SCALER_UP;
                    boost_mult *= SCALER_UP;
                    boost_count++;
                }
            }
        }

        double max_val = out_L.maxCoeff();
        while (max_val < SCALER_THRESH && max_val > 0.0) {
            out_L *= SCALER_UP;
            max_val *= SCALER_UP;
            boost_mult *= SCALER_UP;
            boost_count++;
        }
        out_child_scale_sum -= boost_count * LOG_SCALER;
        if (out_node_boost != nullptr) {
            *out_node_boost = boost_mult;
        }
    }

    static Tree from_newick_file(const std::string& filepath, bool keep_rooted = true) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open Newick file: " + filepath);
        }
        std::stringstream buffer;
        buffer << file.rdbuf();
        return from_newick(buffer.str(), keep_rooted);
    }

    static Tree from_newick(std::string_view newick_sv, bool keep_rooted = true) {
        Tree tree;
        std::string newick(newick_sv);

        // Strip trailing semicolon and whitespace
        while (!newick.empty() && (newick.back() == ';' || newick.back() == ' ' || newick.back() == '\n' || newick.back() == '\r')) {
            newick.pop_back();
        }

        // Clean comments [comment], keep model tags {Tag}
        std::string clean;
        bool in_comment = false;
        for (char c : newick) {
            if (c == '[') { in_comment = true; continue; }
            if (c == ']') { in_comment = false; continue; }
            if (!in_comment && c != ' ' && c != '\t' && c != '\n' && c != '\r') {
                clean += c;
            }
        }

        // Parse Newick
        std::stack<int32_t> parent_stack;
        std::string token;
        int32_t current_id = INVALID_INDEX;

        auto create_node = [&]() -> int32_t {
            int32_t id = static_cast<int32_t>(tree.nodes.size());
            TreeNode node;
            node.id = id;
            tree.nodes.push_back(std::move(node));
            return id;
        };

        auto process_token = [&](int32_t target_node, std::string str) {
            if (str.empty()) return;
            std::string tag;
            auto open_brace = str.find('{');
            auto close_brace = str.find('}');
            if (open_brace != std::string::npos && close_brace != std::string::npos && close_brace > open_brace) {
                tag = str.substr(open_brace + 1, close_brace - open_brace - 1);
                str.erase(open_brace, close_brace - open_brace + 1);
            }
            auto colon = str.find(':');
            std::string name;
            Scalar length = 0.0;
            if (colon != std::string::npos) {
                name = str.substr(0, colon);
                try {
                    length = std::stod(str.substr(colon + 1));
                } catch (...) {
                    length = 0.0;
                }
            } else {
                name = str;
            }
            if (name.size() >= 2 && ((name.front() == '\'' && name.back() == '\'') || (name.front() == '"' && name.back() == '"'))) {
                name = name.substr(1, name.size() - 2);
            }
            tree.nodes[target_node].name = name;
            tree.nodes[target_node].branch_length = length;
            tree.nodes[target_node].model_tag = tag;
        };

        for (size_t i = 0; i < clean.size(); ++i) {
            char c = clean[i];
            if (c == '(') {
                int32_t internal_id = create_node();
                if (!parent_stack.empty()) {
                    int32_t p = parent_stack.top();
                    tree.nodes[p].children.push_back(internal_id);
                    tree.nodes[internal_id].parent_id = p;
                }
                parent_stack.push(internal_id);
            } else if (c == ',' || c == ')') {
                if (!token.empty()) {
                    int32_t leaf_id = create_node();
                    tree.nodes[leaf_id].is_leaf = true;
                    process_token(leaf_id, token);
                    token.clear();
                    if (!parent_stack.empty()) {
                        int32_t p = parent_stack.top();
                        tree.nodes[p].children.push_back(leaf_id);
                        tree.nodes[leaf_id].parent_id = p;
                    }
                }
                if (c == ')') {
                    if (!parent_stack.empty()) {
                        current_id = parent_stack.top();
                        parent_stack.pop();
                        // Peek forward for internal node name or branch length
                        size_t j = i + 1;
                        std::string int_token;
                        while (j < clean.size() && clean[j] != ',' && clean[j] != ')' && clean[j] != ';') {
                            int_token += clean[j];
                            ++j;
                        }
                        if (!int_token.empty()) {
                            process_token(current_id, int_token);
                            i = j - 1; // Advance loop
                        }
                    }
                }
            } else {
                token += c;
            }
        }

        if (tree.nodes.empty()) {
            throw std::runtime_error("Failed to parse Newick tree: empty tree");
        }

        // The root is the node with parent_id == INVALID_INDEX
        // Auto-names for unnamed nodes (and any colliding internal labels from
        // stitched subtrees) must be strictly unique so the JSON phylogeny
        // block never encounters duplicate node keys.
        std::unordered_set<std::string> used_names;
        used_names.reserve(tree.nodes.size() * 2);
        for (const auto& node : tree.nodes) {
            if (node.is_leaf && !node.name.empty()) {
                used_names.insert(node.name);
            }
        }
        int node_counter = 1;
        for (auto& node : tree.nodes) {
            if (node.parent_id == INVALID_INDEX) {
                tree.root_id = node.id;
            }
            if (!node.is_leaf) {
                if (node.name.empty() || used_names.count(node.name)) {
                    std::string candidate;
                    do {
                        candidate = "Node" + std::to_string(node_counter++);
                    } while (used_names.count(candidate));
                    node.name = candidate;
                }
                used_names.insert(node.name);
            } else if (node.name.empty()) {
                std::string candidate;
                do {
                    candidate = "Leaf" + std::to_string(node_counter++);
                } while (used_names.count(candidate));
                node.name = candidate;
                used_names.insert(node.name);
            }
        }

        if (!keep_rooted) {
            tree.unroot_if_bifurcating();
        }

        // Build leaf map and compute post-order traversal
        tree.leaf_name_to_id.clear();
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                tree.leaf_name_to_id[node.name] = node.id;
            }
        }

        tree.compute_post_order();
        return tree;
    }

    /**
     * @brief Collapse a degree-2 bifurcating root into an unrooted (trifurcating+) root,
     *        matching HyPhy's _TreeTopology behavior under reversible models.
     */
    void unroot_if_bifurcating() {
        if (root_id == INVALID_INDEX || nodes.size() <= 3) return;
        if (nodes[root_id].children.size() != 2) return;

        int32_t c0 = nodes[root_id].children[0];
        int32_t c1 = nodes[root_id].children[1];

        // Choose an internal child to absorb into the root; if both are internal,
        // prefer the one without an explicit branch model tag so tagged clades stay intact.
        int32_t absorb_id = INVALID_INDEX;
        int32_t other_id = INVALID_INDEX;
        if (!nodes[c0].is_leaf && !nodes[c1].is_leaf) {
            if (nodes[c1].model_tag.empty() && !nodes[c0].model_tag.empty()) {
                absorb_id = c1;
                other_id = c0;
            } else {
                absorb_id = c0;
                other_id = c1;
            }
        } else if (!nodes[c0].is_leaf) {
            absorb_id = c0;
            other_id = c1;
        } else if (!nodes[c1].is_leaf) {
            absorb_id = c1;
            other_id = c0;
        } else {
            return;
        }

        // Merge branch lengths across the former root split into other_id
        nodes[other_id].branch_length += nodes[absorb_id].branch_length;
        if (nodes[other_id].model_tag.empty() && !nodes[absorb_id].model_tag.empty()) {
            nodes[other_id].model_tag = nodes[absorb_id].model_tag;
        }

        // Replace absorb_id in root's children with absorb_id's children
        std::vector<int32_t> new_root_children;
        new_root_children.reserve(nodes[absorb_id].children.size() + 1);
        if (absorb_id == c0) {
            for (int32_t gc : nodes[absorb_id].children) {
                new_root_children.push_back(gc);
            }
            new_root_children.push_back(other_id);
        } else {
            new_root_children.push_back(other_id);
            for (int32_t gc : nodes[absorb_id].children) {
                new_root_children.push_back(gc);
            }
        }
        nodes[root_id].children = std::move(new_root_children);

        // Compact nodes array by removing absorb_id and remapping all node IDs
        int32_t old_n = static_cast<int32_t>(nodes.size());
        std::vector<int32_t> old_to_new(old_n, INVALID_INDEX);
        std::vector<TreeNode> new_nodes;
        new_nodes.reserve(old_n - 1);
        for (int32_t i = 0; i < old_n; ++i) {
            if (i == absorb_id) continue;
            int32_t new_id = static_cast<int32_t>(new_nodes.size());
            old_to_new[i] = new_id;
            TreeNode nd = std::move(nodes[i]);
            nd.id = new_id;
            new_nodes.push_back(std::move(nd));
        }

        root_id = old_to_new[root_id];
        for (auto& nd : new_nodes) {
            if (nd.id == root_id) {
                nd.parent_id = INVALID_INDEX;
                nd.branch_length = 0.0;
            }
            for (auto& cid : nd.children) {
                cid = old_to_new[cid];
                if (cid != INVALID_INDEX) {
                    new_nodes[cid].parent_id = nd.id;
                }
            }
        }
        nodes = std::move(new_nodes);

        // Update any Wheatstone bridge indices if present
        for (auto& wb : bridges) {
            if (wb.hub_node_id >= 0 && wb.hub_node_id < old_n) wb.hub_node_id = old_to_new[wb.hub_node_id];
            if (wb.bridge_node_id >= 0 && wb.bridge_node_id < old_n) wb.bridge_node_id = old_to_new[wb.bridge_node_id];
            for (int r = 0; r < 3; ++r) {
                if (wb.children[r] >= 0 && wb.children[r] < old_n) wb.children[r] = old_to_new[wb.children[r]];
            }
        }
        if (!bridges.empty()) {
            rebuild_bridge_index();
        }

        leaf_name_to_id.clear();
        for (const auto& nd : nodes) {
            if (nd.is_leaf) {
                leaf_name_to_id[nd.name] = nd.id;
            }
        }
        compute_post_order();
    }

    size_t num_leaves() const {
        return leaf_name_to_id.size();
    }

    size_t num_nodes() const {
        return nodes.size();
    }

    std::string to_newick(bool include_branch_lengths = true) const {
        if (root_id == INVALID_INDEX || nodes.empty()) return "";
        auto write_node = [&](auto& self, int32_t nid) -> std::string {
            const auto& node = nodes[nid];
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
            if (include_branch_lengths && nid != root_id) {
                std::ostringstream ss;
                ss << ":" << std::setprecision(8) << node.branch_length;
                res += ss.str();
            }
            return res;
        };
        return write_node(write_node, root_id) + ";";
    }

    /**
     * @brief Serialize tree/network to Cardona et al. (2008) Extended Newick (eNewick) with NHX annotations.
     *
     * For every active Wheatstone bridge at a 3-child hub u with children {A, B, C} and Ohmic weights
     * {w_1, w_2, w_3} (where w_1 pairs (A,B), w_2 shunts (A,C) with inheritance probability gamma = w_2/(w_1+w_2)):
     *   - Clade A is wrapped in a hybrid reticulation node `#H<k>` whose primary parent is the (A,B) bridge
     *     with weight `w_1` and whose secondary shunt parent is the (C, #H<k>) bridge with weight `w_2`.
     *   - Polytomy busbar nodes (>2 children without a bridge) are annotated with `[&&NHX:ohmic_busbar=true:degree=K]`.
     */
    std::string to_extended_newick(bool include_branch_lengths = true, bool include_nhx = true) const {
        if (root_id == INVALID_INDEX || nodes.empty()) return "";
        auto fmt_num = [](Scalar v, int prec = 6) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(prec) << v;
            return ss.str();
        };

        auto write_node = [&](auto& self, int32_t nid) -> std::string {
            const auto& node = nodes[nid];
            int32_t b_idx = bridge_index(nid);

            // Case 1: Active 4-Way Wheatstone Bridge Hub
            if (b_idx >= 0 && static_cast<size_t>(b_idx) < bridges.size()) {
                const auto& wb = bridges[static_cast<size_t>(b_idx)];
                int32_t id_A = wb.children[0];
                int32_t id_B = wb.children[1];
                int32_t id_C = wb.children[2];
                if (id_A >= 0 && id_B >= 0 && id_C >= 0 &&
                    static_cast<size_t>(id_A) < nodes.size() &&
                    static_cast<size_t>(id_B) < nodes.size() &&
                    static_cast<size_t>(id_C) < nodes.size()) {

                    std::string sub_A = self(self, id_A);
                    std::string sub_B = self(self, id_B);
                    std::string sub_C = self(self, id_C);

                    Scalar w1 = wb.weights[0];
                    Scalar w2 = wb.weights[1];
                    Scalar w3 = wb.weights[2];
                    Scalar w_tot = std::max<Scalar>(1e-12, w1 + w2);
                    Scalar gamma1 = w1 / w_tot;
                    Scalar gamma2 = w2 / w_tot;
                    Scalar t1 = std::max<Scalar>(1e-5, wb.branch_lengths[0]);
                    Scalar t2 = std::max<Scalar>(1e-5, wb.branch_lengths[1]);

                    std::string h_tag = "#H" + std::to_string(b_idx + 1);
                    std::string res;
                    if (w2 >= 0.01) {
                        // Explicit eNewick hybrid reticulation: clade A is reached via primary (A,B) bridge (gamma1)
                        // and secondary Wheatstone shunt (C, #H) bridge (gamma2)
                        res = "((" + sub_A + ")" + h_tag + ":0.0::" + fmt_num(gamma1, 4) + "," + sub_B + ")Bridge_" +
                              std::to_string(b_idx + 1) + "_AB:" + fmt_num(t1) + ",(" + sub_C + "," + h_tag +
                              ":0.0::" + fmt_num(gamma2, 4) + ")Bridge_" + std::to_string(b_idx + 1) + "_Shunt:" + fmt_num(t2) + ")";
                    } else {
                        res = "((" + sub_A + "," + sub_B + "):" + fmt_num(t1) + "," + sub_C + ")";
                    }
                    res += node.name;
                    if (include_nhx) {
                        res += "[&&NHX:wheatstone_hub=true:ohmic_support=" + fmt_num(wb.ohmic_support, 4) +
                               ":shunt_stress=" + fmt_num(wb.shunt_stress, 4) +
                               ":w1=" + fmt_num(w1, 4) + ":w2=" + fmt_num(w2, 4) + ":w3=" + fmt_num(w3, 4) + "]";
                    }
                    if (include_branch_lengths && nid != root_id) {
                        res += ":" + fmt_num(node.branch_length, 8);
                    }
                    return res;
                }
            }

            // Case 2: Standard Tree Node or Polytomy Busbar
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
            if (include_branch_lengths && nid != root_id) {
                res += ":" + fmt_num(node.branch_length, 8);
            }
            if (include_nhx && !node.is_leaf && node.children.size() > 2 && nid != root_id) {
                res += "[&&NHX:ohmic_busbar=true:degree=" + std::to_string(node.children.size()) + "]";
            }
            return res;
        };
        return write_node(write_node, root_id) + ";";
    }

    void compute_post_order() {
        post_order.clear();
        if (root_id == INVALID_INDEX) return;
        traverse_post_order(root_id);
    }

    void traverse_post_order(int32_t node_id) {
        for (int32_t child_id : nodes[node_id].children) {
            traverse_post_order(child_id);
        }
        post_order.push_back(node_id);
    }

    struct BranchSelection {
        std::vector<bool> is_selected;
        std::vector<std::string> selected_names;
        std::vector<std::string> unselected_names;
        size_t num_selected = 0;
        size_t num_unselected = 0;
    };

    BranchSelection select_branches(
        const std::string& spec = "All",
        const std::string& regex_str = "",
        const std::vector<std::string>& explicit_names = {},
        bool auto_detect_test_tags_when_all = false
    ) const {
        BranchSelection sel;
        size_t n = nodes.size();
        sel.is_selected.assign(n, false);

        auto to_lower = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return s;
        };
        auto trim = [](std::string s) {
            size_t start = s.find_first_not_of(" \t\r\n");
            if (start == std::string::npos) return std::string{};
            size_t end = s.find_last_not_of(" \t\r\n");
            return s.substr(start, end - start + 1);
        };
        auto strip_braces = [&](std::string s) {
            s = trim(std::move(s));
            if (s.size() >= 2 && s.front() == '{' && s.back() == '}') {
                return trim(s.substr(1, s.size() - 2));
            }
            return s;
        };
        auto is_fg_synonym = [](const std::string& s_low) {
            return s_low == "t" || s_low == "test" || s_low == "foreground" || s_low == "fg";
        };
        auto is_bg_synonym = [](const std::string& s_low) {
            return s_low == "r" || s_low == "ref" || s_low == "reference" || s_low == "background" || s_low == "bg";
        };

        bool has_tagged_test = false;
        if (auto_detect_test_tags_when_all) {
            for (const auto& node : nodes) {
                if (node.id != root_id && is_fg_synonym(to_lower(node.model_tag))) {
                    has_tagged_test = true;
                    break;
                }
            }
        }

        if (!regex_str.empty()) {
            std::regex re(regex_str, std::regex::icase);
            for (const auto& node : nodes) {
                if (node.id == root_id) continue;
                if (std::regex_search(node.name, re) || (!node.model_tag.empty() && std::regex_search(node.model_tag, re))) {
                    sel.is_selected[node.id] = true;
                }
            }
        } else if (!explicit_names.empty()) {
            std::unordered_set<std::string> lower_targets;
            for (const auto& item : explicit_names) {
                std::string cleaned = to_lower(strip_braces(item));
                if (!cleaned.empty()) lower_targets.insert(cleaned);
            }
            for (const auto& node : nodes) {
                if (node.id == root_id) continue;
                std::string name_low = to_lower(node.name);
                std::string tag_low = to_lower(node.model_tag);
                if (lower_targets.count(name_low) > 0 || (!tag_low.empty() && lower_targets.count(tag_low) > 0)) {
                    sel.is_selected[node.id] = true;
                }
            }
        } else {
            std::string mode = strip_braces(spec);
            std::string mode_low = to_lower(mode);

            if (mode.empty() || mode_low == "all") {
                for (const auto& node : nodes) {
                    if (node.id == root_id) continue;
                    if (auto_detect_test_tags_when_all && has_tagged_test) {
                        sel.is_selected[node.id] = is_fg_synonym(to_lower(node.model_tag));
                    } else {
                        sel.is_selected[node.id] = true;
                    }
                }
            } else if (mode_low == "all-branches" || mode_low == "all_branches") {
                for (const auto& node : nodes) {
                    if (node.id == root_id) continue;
                    sel.is_selected[node.id] = true;
                }
            } else if (mode_low == "internal") {
                for (const auto& node : nodes) {
                    if (node.id == root_id) continue;
                    sel.is_selected[node.id] = !node.is_leaf;
                }
            } else if (mode_low == "leaves" || mode_low == "terminal" || mode_low == "tips") {
                for (const auto& node : nodes) {
                    if (node.id == root_id) continue;
                    sel.is_selected[node.id] = node.is_leaf;
                }
            } else if (mode_low == "unlabeled" || mode_low == "untagged") {
                for (const auto& node : nodes) {
                    if (node.id == root_id) continue;
                    sel.is_selected[node.id] = node.model_tag.empty();
                }
            } else if (mode.find(',') != std::string::npos) {
                std::unordered_set<std::string> tokens;
                std::stringstream ss(mode);
                std::string item;
                while (std::getline(ss, item, ',')) {
                    std::string cleaned = to_lower(strip_braces(item));
                    if (!cleaned.empty()) tokens.insert(cleaned);
                }
                for (const auto& node : nodes) {
                    if (node.id == root_id) continue;
                    std::string name_low = to_lower(node.name);
                    std::string tag_low = to_lower(node.model_tag);
                    bool match = (tokens.count(name_low) > 0) || (!tag_low.empty() && tokens.count(tag_low) > 0);
                    if (!match && tokens.count("internal") > 0 && !node.is_leaf) match = true;
                    if (!match && (tokens.count("leaves") > 0 || tokens.count("terminal") > 0) && node.is_leaf) match = true;
                    if (!match && (tokens.count("unlabeled") > 0 || tokens.count("untagged") > 0) && node.model_tag.empty()) match = true;
                    if (!match && is_fg_synonym(tag_low)) {
                        for (const auto& tok : tokens) {
                            if (is_fg_synonym(tok)) { match = true; break; }
                        }
                    }
                    if (!match && is_bg_synonym(tag_low)) {
                        for (const auto& tok : tokens) {
                            if (is_bg_synonym(tok)) { match = true; break; }
                        }
                    }
                    sel.is_selected[node.id] = match;
                }
            } else {
                // Single token: check tag (including FG/BG synonyms) or branch name first
                bool any_exact = false;
                for (const auto& node : nodes) {
                    if (node.id == root_id) continue;
                    std::string name_low = to_lower(node.name);
                    std::string tag_low = to_lower(node.model_tag);
                    bool match = (!tag_low.empty() && tag_low == mode_low) || (name_low == mode_low);
                    if (!match && is_fg_synonym(mode_low) && is_fg_synonym(tag_low)) match = true;
                    if (!match && is_bg_synonym(mode_low) && is_bg_synonym(tag_low)) match = true;
                    if (match) {
                        sel.is_selected[node.id] = true;
                        any_exact = true;
                    }
                }
                // Fallback to regex only if mode contains regex metacharacters
                if (!any_exact && spec.find_first_of(".*+?^$[]()|\\") != std::string::npos) {
                    try {
                        std::regex re(spec, std::regex::icase);
                        for (const auto& node : nodes) {
                            if (node.id == root_id) continue;
                            if (std::regex_search(node.name, re) || (!node.model_tag.empty() && std::regex_search(node.model_tag, re))) {
                                sel.is_selected[node.id] = true;
                            }
                        }
                    } catch (...) {
                    }
                }
            }
        }

        for (const auto& node : nodes) {
            if (node.id == root_id) continue;
            if (sel.is_selected[node.id]) {
                sel.selected_names.push_back(node.name);
                sel.num_selected++;
            } else {
                sel.unselected_names.push_back(node.name);
                sel.num_unselected++;
            }
        }

        if (sel.num_selected == 0 && n > 1) {
            std::string requested = !regex_str.empty() ? regex_str : spec;
            throw std::runtime_error("No branches matched '--branches " + requested + "'.");
        }

        return sel;
    }
};

} // namespace hyphy::core
