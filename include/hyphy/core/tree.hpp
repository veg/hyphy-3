#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/nexus_utils.hpp"
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
                if (P_bridges == nullptr || static_cast<size_t>(b_idx) >= P_bridges->size()) {
                    throw std::logic_error("Wheatstone bridge transition matrices P_bridges required when tree has bridges");
                }
                return (*P_bridges)[static_cast<size_t>(b_idx)][r];
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
            if (!nodes[nid].children.empty()) {
                out_L.setOnes(P_branches[nodes[nid].children[0]].rows());
            } else {
                out_L.setOnes();
            }
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
        const std::vector<std::array<MatType, 3>>* P_bridges = nullptr,
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
            if (!nodes[nid].children.empty()) {
                out_L.setOnes(P_branches[nodes[nid].children[0]].rows());
            } else {
                out_L.setOnes();
            }
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
        std::string newick;
        if (looks_like_nexus_tree_input(newick_sv)) {
            std::string extracted = extract_nexus_tree_newick(strip_nexus_comments(newick_sv), 0);
            newick = extracted.empty() ? std::string(newick_sv) : std::move(extracted);
        } else {
            newick.assign(newick_sv);
        }

        // Strip trailing semicolon and whitespace
        while (!newick.empty() && (newick.back() == ';' || newick.back() == ' ' || newick.back() == '\t' || newick.back() == '\n' || newick.back() == '\r')) {
            newick.pop_back();
        }

        // Clean comments [comment] outside quotes, keep model tags {Tag}, preserve characters inside quotes
        std::string clean;
        clean.reserve(newick.size());
        bool in_comment = false;
        char quote_char = '\0';
        for (size_t k = 0; k < newick.size(); ++k) {
            char c = newick[k];
            if (quote_char != '\0') {
                clean += c;
                if (c == quote_char) {
                    if (quote_char == '\'' && k + 1 < newick.size() && newick[k + 1] == '\'') {
                        clean += newick[++k];
                    } else {
                        quote_char = '\0';
                    }
                }
                continue;
            }
            if (in_comment) {
                if (c == ']') in_comment = false;
                continue;
            }
            if (c == '[') {
                in_comment = true;
                continue;
            }
            if (c == '\'' || c == '"') {
                quote_char = c;
                clean += c;
                continue;
            }
            if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
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
            size_t open_brace = std::string::npos;
            size_t close_brace = std::string::npos;
            {
                char q = '\0';
                for (size_t k = 0; k < str.size(); ++k) {
                    char ch = str[k];
                    if (q != '\0') {
                        if (ch == q) {
                            if (q == '\'' && k + 1 < str.size() && str[k + 1] == '\'') {
                                ++k;
                            } else {
                                q = '\0';
                            }
                        }
                        continue;
                    }
                    if (ch == '\'' || ch == '"') {
                        q = ch;
                        continue;
                    }
                    if (ch == '{' && open_brace == std::string::npos) {
                        open_brace = k;
                    } else if (ch == '}' && open_brace != std::string::npos && close_brace == std::string::npos) {
                        close_brace = k;
                        break;
                    }
                }
            }
            if (open_brace != std::string::npos && close_brace != std::string::npos && close_brace > open_brace) {
                tag = str.substr(open_brace + 1, close_brace - open_brace - 1);
                str.erase(open_brace, close_brace - open_brace + 1);
            }

            size_t colon = std::string::npos;
            {
                char q = '\0';
                for (size_t k = 0; k < str.size(); ++k) {
                    char ch = str[k];
                    if (q != '\0') {
                        if (ch == q) {
                            if (q == '\'' && k + 1 < str.size() && str[k + 1] == '\'') {
                                ++k;
                            } else {
                                q = '\0';
                            }
                        }
                        continue;
                    }
                    if (ch == '\'' || ch == '"') {
                        q = ch;
                        continue;
                    }
                    if (ch == ':') {
                        colon = k;
                        break;
                    }
                }
            }

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
                char q = name.front();
                std::string inner = name.substr(1, name.size() - 2);
                if (q == '\'') {
                    std::string unescaped;
                    unescaped.reserve(inner.size());
                    for (size_t k = 0; k < inner.size(); ++k) {
                        if (inner[k] == '\'' && k + 1 < inner.size() && inner[k + 1] == '\'') {
                            unescaped += '\'';
                            ++k;
                        } else {
                            unescaped += inner[k];
                        }
                    }
                    name = std::move(unescaped);
                } else {
                    name = std::move(inner);
                }
            }
            tree.nodes[target_node].name = name;
            tree.nodes[target_node].branch_length = length;
            tree.nodes[target_node].model_tag = tag;
        };

        char quote = '\0';
        for (size_t i = 0; i < clean.size(); ++i) {
            char c = clean[i];
            if (quote != '\0') {
                token += c;
                if (c == quote) {
                    if (quote == '\'' && i + 1 < clean.size() && clean[i + 1] == '\'') {
                        token += clean[++i];
                    } else {
                        quote = '\0';
                    }
                }
                continue;
            }
            if (c == '\'' || c == '"') {
                quote = c;
                token += c;
                continue;
            }
            if (c == ';') {
                break;
            }
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
                    if (parent_stack.empty()) {
                        throw std::runtime_error("Invalid Newick tree: unbalanced parentheses");
                    }
                    current_id = parent_stack.top();
                    parent_stack.pop();
                    // Peek forward for internal node name or branch length (respecting quotes)
                    size_t j = i + 1;
                    std::string int_token;
                    char int_quote = '\0';
                    while (j < clean.size()) {
                        char cj = clean[j];
                        if (int_quote != '\0') {
                            int_token += cj;
                            if (cj == int_quote) {
                                if (int_quote == '\'' && j + 1 < clean.size() && clean[j + 1] == '\'') {
                                    int_token += clean[++j];
                                } else {
                                    int_quote = '\0';
                                }
                            }
                            ++j;
                            continue;
                        }
                        if (cj == '\'' || cj == '"') {
                            int_quote = cj;
                            int_token += cj;
                            ++j;
                            continue;
                        }
                        if (cj == ',' || cj == ')' || cj == '(' || cj == ';') {
                            break;
                        }
                        int_token += cj;
                        ++j;
                    }
                    if (int_quote != '\0') {
                        throw std::runtime_error("Invalid Newick tree: unterminated quoted label");
                    }
                    if (!int_token.empty()) {
                        process_token(current_id, int_token);
                        i = j - 1; // Advance loop
                    }
                }
            } else {
                token += c;
            }
        }

        if (quote != '\0') {
            throw std::runtime_error("Invalid Newick tree: unterminated quoted label");
        }
        if (!parent_stack.empty()) {
            throw std::runtime_error("Invalid Newick tree: unbalanced parentheses");
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
     * @brief Prune tree down to the specified leaf names (case-insensitive matching),
     * collapsing degree-1 internal nodes and accumulating branch lengths.
     */
    Tree prune_to_leaves(const std::vector<std::string>& keep_leaf_names, bool unroot = true) const {
        if (root_id == INVALID_INDEX || nodes.empty()) return *this;

        auto to_lower_str = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return s;
        };
        std::unordered_set<std::string> keep_lower;
        keep_lower.reserve(keep_leaf_names.size());
        for (const auto& name : keep_leaf_names) {
            keep_lower.insert(to_lower_str(name));
        }

        auto serialize_pruned = [&](auto& self, int32_t nid) -> std::pair<std::string, Scalar> {
            const auto& nd = nodes[nid];
            if (nd.is_leaf) {
                if (keep_lower.count(to_lower_str(nd.name)) > 0) {
                    std::string token = nd.name;
                    if (!nd.model_tag.empty()) token += "{" + nd.model_tag + "}";
                    return {token, nd.branch_length};
                }
                return {"", 0.0};
            }
            std::vector<std::pair<std::string, Scalar>> kept_children;
            for (int32_t cid : nd.children) {
                auto sub = self(self, cid);
                if (!sub.first.empty()) {
                    kept_children.push_back(std::move(sub));
                }
            }
            if (kept_children.empty()) {
                return {"", 0.0};
            }
            if (kept_children.size() == 1) {
                Scalar merged_bl = kept_children[0].second + ((nid == root_id) ? 0.0 : nd.branch_length);
                return {kept_children[0].first, merged_bl};
            }
            std::ostringstream oss;
            oss << "(";
            for (size_t i = 0; i < kept_children.size(); ++i) {
                if (i > 0) oss << ",";
                oss << kept_children[i].first << ":" << std::setprecision(10) << kept_children[i].second;
            }
            oss << ")" << nd.name;
            if (!nd.model_tag.empty()) oss << "{" + nd.model_tag + "}";
            return {oss.str(), (nid == root_id) ? 0.0 : nd.branch_length};
        };

        auto root_pair = serialize_pruned(serialize_pruned, root_id);
        if (root_pair.first.empty() || root_pair.first.find('(') == std::string::npos) {
            throw std::runtime_error("Pruning left fewer than 2 tips in the tree.");
        }
        return Tree::from_newick(root_pair.first + ";", unroot);
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
        bool auto_detect_test_tags_when_all = false,
        bool throw_on_empty = true
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

        if (throw_on_empty && sel.num_selected == 0 && n > 1) {
            std::string requested = !regex_str.empty() ? regex_str : spec;
            throw std::runtime_error("No branches matched '--branches " + requested + "'.");
        }

        return sel;
    }

    static std::vector<std::string> parse_taxon_list_or_file(const std::string& spec) {
        auto trim_ws = [](const std::string& s) -> std::string {
            size_t a = s.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) return {};
            size_t b = s.find_last_not_of(" \t\r\n");
            return s.substr(a, b - a + 1);
        };
        auto strip_quotes = [](const std::string& s) -> std::string {
            if (s.size() >= 2 &&
                ((s.front() == '\'' && s.back() == '\'') ||
                 (s.front() == '"' && s.back() == '"'))) {
                return s.substr(1, s.size() - 2);
            }
            return s;
        };

        std::vector<std::string> items;
        if (spec.empty()) return items;

        if (spec.find('\n') == std::string::npos && spec.find('\r') == std::string::npos) {
            std::ifstream ifs(spec);
            if (ifs.is_open()) {
                std::string line;
                bool is_fasta = false;
                while (std::getline(ifs, line)) {
                    std::string t = trim_ws(line);
                    if (t.empty() || t[0] == '#') continue;
                    if (t[0] == '>') {
                        is_fasta = true;
                        std::string hdr = trim_ws(t.substr(1));
                        size_t sp = hdr.find_first_of(" \t");
                        if (sp != std::string::npos) hdr = hdr.substr(0, sp);
                        hdr = strip_quotes(trim_ws(hdr));
                        if (!hdr.empty()) items.push_back(std::move(hdr));
                    } else if (!is_fasta) {
                        std::stringstream ss(t);
                        std::string tok;
                        while (std::getline(ss, tok, ',')) {
                            tok = strip_quotes(trim_ws(tok));
                            if (!tok.empty()) items.push_back(std::move(tok));
                        }
                    }
                }
                return items;
            }
        }

        std::stringstream ss(spec);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            tok = strip_quotes(trim_ws(tok));
            if (!tok.empty()) items.push_back(std::move(tok));
        }
        return items;
    }

    /**
     * @brief Parse a 2-column CSV/TSV/semicolon branch label map file
     *        (`Taxon,Group` / `Taxon\tGroup`) preserving first-seen group order.
     */
    static std::vector<std::pair<std::string, std::vector<std::string>>> parse_label_map_file(
        const std::string& filepath
    ) {
        auto trim_ws = [](const std::string& s) -> std::string {
            size_t a = s.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) return {};
            size_t b = s.find_last_not_of(" \t\r\n");
            return s.substr(a, b - a + 1);
        };
        auto to_low = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return s;
        };
        auto strip_quotes = [](const std::string& s) -> std::string {
            if (s.size() >= 2 &&
                ((s.front() == '\'' && s.back() == '\'') ||
                 (s.front() == '"' && s.back() == '"'))) {
                return s.substr(1, s.size() - 2);
            }
            return s;
        };
        auto strip_braces = [&](const std::string& s) -> std::string {
            if (s.size() >= 2 && s.front() == '{' && s.back() == '}') {
                return trim_ws(s.substr(1, s.size() - 2));
            }
            return s;
        };

        std::ifstream ifs(filepath);
        if (!ifs.is_open()) {
            throw std::runtime_error("Could not open branch label file: " + filepath);
        }

        const std::unordered_set<std::string> col0_headers = {
            "taxon", "taxa", "species", "tip", "tips", "branch", "branches",
            "sequence", "name", "sample", "id"
        };
        const std::unordered_set<std::string> col1_headers = {
            "label", "labels", "group", "groups", "regime", "regimes", "tag",
            "tags", "clade", "set", "partition", "category", "phenotype", "trait"
        };

        std::vector<std::pair<std::string, std::vector<std::string>>> groups;
        std::unordered_map<std::string, size_t> group_idx_by_low;
        std::unordered_map<std::string, std::string> taxon_to_group_low;

        std::string raw_line;
        bool seen_first_data_row = false;
        while (std::getline(ifs, raw_line)) {
            std::string line = trim_ws(raw_line);
            if (line.empty() || line[0] == '#') continue;

            char delim = '\0';
            size_t dpos = line.find('\t');
            if (dpos != std::string::npos) {
                delim = '\t';
            } else if ((dpos = line.find(',')) != std::string::npos) {
                delim = ',';
            } else if ((dpos = line.find(';')) != std::string::npos) {
                delim = ';';
            }
            if (delim == '\0') {
                throw std::runtime_error(
                    "Invalid row in branch label file '" + filepath +
                    "' (expected 2 columns separated by tab, comma, or semicolon): " + line
                );
            }

            std::string col0 = strip_quotes(trim_ws(line.substr(0, dpos)));
            size_t next_dpos = line.find(delim, dpos + 1);
            std::string col1_raw = (next_dpos == std::string::npos)
                ? line.substr(dpos + 1)
                : line.substr(dpos + 1, next_dpos - (dpos + 1));
            std::string col1 = trim_ws(strip_braces(strip_quotes(trim_ws(col1_raw))));

            if (!seen_first_data_row) {
                seen_first_data_row = true;
                if (col0_headers.count(to_low(col0)) > 0 && col1_headers.count(to_low(col1)) > 0) {
                    continue;
                }
            }

            if (col0.empty()) {
                throw std::runtime_error(
                    "Empty taxon/branch name in branch label file '" + filepath + "': " + line
                );
            }
            if (col1.empty()) {
                continue; // Unlabeled row
            }

            std::string clean_label;
            clean_label.reserve(col1.size());
            bool prev_ws = false;
            for (unsigned char ch : col1) {
                if (std::isspace(ch)) {
                    if (!prev_ws) {
                        clean_label.push_back('_');
                        prev_ws = true;
                    }
                } else {
                    clean_label.push_back(static_cast<char>(ch));
                    prev_ws = false;
                }
            }

            if (clean_label.find_first_of(",():;{}[]") != std::string::npos) {
                throw std::runtime_error(
                    "Invalid Newick character in group label '" + clean_label + "' from '" + filepath + "'."
                );
            }

            std::string tax_low = to_low(col0);
            std::string grp_low = to_low(clean_label);
            auto it_existing = taxon_to_group_low.find(tax_low);
            if (it_existing != taxon_to_group_low.end()) {
                if (it_existing->second != grp_low) {
                    throw std::runtime_error(
                        "Conflicting group labels for taxon '" + col0 + "' in '" + filepath + "'."
                    );
                }
                continue;
            }
            taxon_to_group_low.emplace(tax_low, grp_low);

            auto it_g = group_idx_by_low.find(grp_low);
            if (it_g == group_idx_by_low.end()) {
                size_t new_idx = groups.size();
                group_idx_by_low.emplace(grp_low, new_idx);
                groups.emplace_back(clean_label, std::vector<std::string>{col0});
            } else {
                groups[it_g->second].second.push_back(col0);
            }
        }

        if (groups.empty()) {
            throw std::runtime_error("Branch label file '" + filepath + "' contained no taxon-to-group mappings.");
        }
        return groups;
    }

    static std::string normalize_internal_label_mode(const std::string& internal_mode_raw) {
        std::string m;
        size_t a = internal_mode_raw.find_first_not_of(" \t\r\n");
        if (a != std::string::npos) {
            size_t b = internal_mode_raw.find_last_not_of(" \t\r\n");
            m = internal_mode_raw.substr(a, b - a + 1);
        }
        std::transform(m.begin(), m.end(), m.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (m.empty() || m == "all" || m == "clade" || m == "monophyletic") return "all";
        if (m == "parsimony" || m == "fitch") return "parsimony";
        if (m == "mrca") return "mrca";
        if (m == "paraphyletic") return "paraphyletic";
        if (m == "only" || m == "internal" || m == "internal-only") return "only";
        if (m == "none" || m == "leaves" || m == "tips" || m == "terminal") return "none";
        throw std::runtime_error(
            "Unknown --internal-nodes strategy '" + internal_mode_raw +
            "'. Expected All, Parsimony, MRCA, Paraphyletic, Only, or None."
        );
    }

    /**
     * @brief Assign `model_tag = label_tag` to selected branches and propagate
     *        to internal branches according to `internal_mode_raw`.
     */
    size_t apply_branch_labels(
        const std::string& label_tag,
        const std::string& list_spec,
        const std::string& internal_mode_raw = "All",
        bool allow_empty_match = false
    ) {
        auto trim_ws = [](const std::string& s) -> std::string {
            size_t a = s.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) return {};
            size_t b = s.find_last_not_of(" \t\r\n");
            return s.substr(a, b - a + 1);
        };
        std::string clean_tag = trim_ws(label_tag);
        if (clean_tag.size() >= 2 && clean_tag.front() == '{' && clean_tag.back() == '}') {
            clean_tag = trim_ws(clean_tag.substr(1, clean_tag.size() - 2));
        }
        if (clean_tag.empty()) {
            throw std::runtime_error("Branch label tag cannot be empty.");
        }
        std::string imode = normalize_internal_label_mode(internal_mode_raw);
        std::vector<std::string> explicit_items = parse_taxon_list_or_file(list_spec);
        bool use_explicit = (explicit_items.size() > 1) ||
                            (explicit_items.size() == 1 && list_spec != explicit_items[0]);
        return apply_branch_labels_resolved(
            clean_tag,
            use_explicit ? explicit_items : std::vector<std::string>{},
            use_explicit ? std::string{} : (list_spec.empty() ? std::string{"All"} : list_spec),
            imode,
            allow_empty_match
        );
    }

    /**
     * @brief Apply multiple named branch label groups (e.g. parsed from CSV/TSV)
     *        onto the tree with clade/parsimony internal-branch propagation.
     */
    size_t apply_label_map(
        const std::vector<std::pair<std::string, std::vector<std::string>>>& groups,
        const std::string& internal_mode_raw = "All",
        const std::string& source_name = "--labels"
    ) {
        std::string imode = normalize_internal_label_mode(internal_mode_raw);
        if (groups.empty() || nodes.size() <= 1 || root_id == INVALID_INDEX) return 0;

        auto to_low = [](std::string s) {
            size_t a = s.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) return std::string{};
            size_t b = s.find_last_not_of(" \t\r\n");
            s = s.substr(a, b - a + 1);
            if (s.size() >= 2 &&
                ((s.front() == '\'' && s.back() == '\'') ||
                 (s.front() == '"' && s.back() == '"'))) {
                s = s.substr(1, s.size() - 2);
            }
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return s;
        };

        const size_t G = groups.size();
        const size_t N = nodes.size();
        std::unordered_map<std::string, size_t> taxon_to_g;
        for (size_t g = 0; g < G; ++g) {
            for (const auto& t : groups[g].second) {
                std::string tl = to_low(t);
                if (!tl.empty()) taxon_to_g[tl] = g;
            }
        }

        std::vector<int> node_group(N, -1);
        std::vector<std::vector<std::string>> present_items(G);
        size_t matched_nodes = 0;
        for (const auto& nd : nodes) {
            if (nd.id == root_id) continue;
            std::string nl = to_low(nd.name);
            auto it = taxon_to_g.find(nl);
            if (it != taxon_to_g.end()) {
                node_group[nd.id] = static_cast<int>(it->second);
                present_items[it->second].push_back(nd.name);
                matched_nodes++;
            }
        }

        if (matched_nodes == 0) {
            throw std::runtime_error(
                "No taxa in '" + source_name + "' matched any branch or tip in the tree."
            );
        }

        compute_post_order();

        if (imode == "all") {
            std::vector<int> total_leaves(N, 0);
            std::vector<std::vector<int>> group_leaves(G, std::vector<int>(N, 0));
            for (int32_t nid : post_order) {
                const auto& nd = nodes[nid];
                if (nd.is_leaf) {
                    total_leaves[nid] = 1;
                    if (node_group[nid] >= 0) {
                        group_leaves[node_group[nid]][nid] = 1;
                    }
                } else {
                    for (int32_t cid : nd.children) {
                        total_leaves[nid] += total_leaves[cid];
                        for (size_t g = 0; g < G; ++g) {
                            group_leaves[g][nid] += group_leaves[g][cid];
                        }
                    }
                }
            }
            for (const auto& nd : nodes) {
                if (nd.id == root_id || nd.is_leaf || total_leaves[nd.id] <= 0) continue;
                for (size_t g = 0; g < G; ++g) {
                    if (group_leaves[g][nd.id] == total_leaves[nd.id]) {
                        node_group[nd.id] = static_cast<int>(g);
                        break;
                    }
                }
            }
        } else if (imode == "parsimony") {
            const size_t num_states = G + 1; // State G == unlabeled
            std::vector<std::vector<uint8_t>> fitch(N, std::vector<uint8_t>(num_states, 0));
            for (int32_t nid : post_order) {
                const auto& nd = nodes[nid];
                if (nd.is_leaf) {
                    int s = (node_group[nid] >= 0 && static_cast<size_t>(node_group[nid]) < G)
                          ? node_group[nid]
                          : static_cast<int>(G);
                    fitch[nid][s] = 1;
                } else if (!nd.children.empty()) {
                    std::vector<size_t> counts(num_states, 0);
                    size_t max_c = 0;
                    for (int32_t cid : nd.children) {
                        for (size_t s = 0; s < num_states; ++s) {
                            if (fitch[cid][s]) {
                                counts[s]++;
                                if (counts[s] > max_c) max_c = counts[s];
                            }
                        }
                    }
                    for (size_t s = 0; s < num_states; ++s) {
                        if (counts[s] == max_c && max_c > 0) fitch[nid][s] = 1;
                    }
                } else {
                    fitch[nid][G] = 1;
                }
            }

            auto is_ref_name = [&](const std::string& raw) -> bool {
                std::string rl = to_low(raw);
                return rl == "reference" || rl == "ref" || rl == "background" || rl == "bg" || rl == "r";
            };
            std::vector<int> state(N, -1);
            if (root_id >= 0 && static_cast<size_t>(root_id) < N) {
                size_t n_opt = 0;
                int single_s = -1;
                int ref_s = -1;
                for (size_t s = 0; s < num_states; ++s) {
                    if (fitch[root_id][s]) {
                        n_opt++;
                        single_s = static_cast<int>(s);
                        if (s < G && ref_s < 0 && is_ref_name(groups[s].first)) {
                            ref_s = static_cast<int>(s);
                        }
                    }
                }
                if (n_opt == 1) {
                    state[root_id] = single_s;
                } else if (fitch[root_id][G]) {
                    state[root_id] = static_cast<int>(G);
                } else if (ref_s >= 0) {
                    state[root_id] = ref_s;
                } else {
                    state[root_id] = -1;
                }
            }

            for (auto it = post_order.rbegin(); it != post_order.rend(); ++it) {
                int32_t nid = *it;
                if (nid == root_id) continue;
                int32_t pid = nodes[nid].parent_id;
                int p_state = (pid >= 0 && static_cast<size_t>(pid) < N) ? state[pid] : -1;
                if (p_state >= 0 && fitch[nid][static_cast<size_t>(p_state)]) {
                    state[nid] = p_state;
                } else {
                    size_t n_opt = 0;
                    int single_s = -1;
                    for (size_t s = 0; s < num_states; ++s) {
                        if (fitch[nid][s]) {
                            n_opt++;
                            single_s = static_cast<int>(s);
                        }
                    }
                    state[nid] = (n_opt == 1) ? single_s : -1;
                }
                if (!nodes[nid].is_leaf && node_group[nid] < 0 &&
                    state[nid] >= 0 && static_cast<size_t>(state[nid]) < G) {
                    node_group[nid] = state[nid];
                }
            }
        } else {
            size_t total_assigned = 0;
            for (size_t g = 0; g < G; ++g) {
                if (present_items[g].empty()) continue;
                total_assigned += apply_branch_labels_resolved(
                    groups[g].first, present_items[g], "", imode, /*allow_empty_match=*/true
                );
            }
            return total_assigned;
        }

        size_t assigned = 0;
        for (auto& nd : nodes) {
            if (nd.id == root_id) continue;
            if (node_group[nd.id] >= 0 && static_cast<size_t>(node_group[nd.id]) < G) {
                nd.model_tag = groups[static_cast<size_t>(node_group[nd.id])].first;
                assigned++;
            }
        }
        return assigned;
    }

private:
    size_t apply_branch_labels_resolved(
        const std::string& clean_tag,
        const std::vector<std::string>& explicit_items,
        const std::string& fallback_spec,
        const std::string& imode,
        bool allow_empty_match
    ) {
        BranchSelection sel = select_branches(
            explicit_items.empty() ? fallback_spec : std::string{},
            "",
            explicit_items,
            false,
            !allow_empty_match
        );
        if (sel.num_selected == 0) return 0;

        size_t N = nodes.size();
        std::vector<bool> target = sel.is_selected;

        compute_post_order();
        std::vector<int> total_leaves(N, 0);
        std::vector<int> selected_leaves(N, 0);
        int global_selected_leaves = 0;
        for (int32_t nid : post_order) {
            const auto& nd = nodes[nid];
            if (nd.is_leaf) {
                total_leaves[nid] = 1;
                selected_leaves[nid] = target[nid] ? 1 : 0;
                global_selected_leaves += selected_leaves[nid];
            } else {
                for (int32_t cid : nd.children) {
                    total_leaves[nid] += total_leaves[cid];
                    selected_leaves[nid] += selected_leaves[cid];
                }
            }
        }

        if (global_selected_leaves > 0) {
            if (imode == "all") {
                for (const auto& nd : nodes) {
                    if (nd.id == root_id || nd.is_leaf) continue;
                    if (total_leaves[nd.id] > 0 && selected_leaves[nd.id] == total_leaves[nd.id]) {
                        target[nd.id] = true;
                    }
                }
            } else if (imode == "parsimony") {
                std::vector<uint8_t> fitch(N, 0); // bit 0: unselected, bit 1: selected
                for (int32_t nid : post_order) {
                    const auto& nd = nodes[nid];
                    if (nd.is_leaf) {
                        fitch[nid] = target[nid] ? 2 : 1;
                    } else if (!nd.children.empty()) {
                        size_t c0 = 0, c1 = 0;
                        for (int32_t cid : nd.children) {
                            if (fitch[cid] & 1) c0++;
                            if (fitch[cid] & 2) c1++;
                        }
                        if (c0 > c1) fitch[nid] = 1;
                        else if (c1 > c0) fitch[nid] = 2;
                        else fitch[nid] = 3;
                    } else {
                        fitch[nid] = 1;
                    }
                }
                std::vector<int> state(N, 0);
                if (root_id >= 0 && static_cast<size_t>(root_id) < N) {
                    // Break root tie toward unselected (state 0)
                    state[root_id] = (fitch[root_id] == 2) ? 1 : 0;
                }
                for (auto it = post_order.rbegin(); it != post_order.rend(); ++it) {
                    int32_t nid = *it;
                    if (nid == root_id) continue;
                    int32_t pid = nodes[nid].parent_id;
                    int p_state = (pid >= 0 && static_cast<size_t>(pid) < N) ? state[pid] : 0;
                    uint8_t p_bit = static_cast<uint8_t>(1u << static_cast<unsigned>(p_state));
                    if (fitch[nid] & p_bit) {
                        state[nid] = p_state;
                    } else {
                        state[nid] = (fitch[nid] == 2) ? 1 : 0;
                    }
                    if (!nodes[nid].is_leaf) {
                        target[nid] = (state[nid] == 1);
                    }
                }
            } else if (imode == "only") {
                for (const auto& nd : nodes) {
                    if (nd.id == root_id) continue;
                    if (nd.is_leaf) {
                        target[nd.id] = false;
                    } else if (total_leaves[nd.id] > 0 && selected_leaves[nd.id] == total_leaves[nd.id]) {
                        target[nd.id] = true;
                    }
                }
            } else if (imode == "mrca") {
                int32_t mrca_id = root_id;
                for (int32_t nid : post_order) {
                    if (selected_leaves[nid] == global_selected_leaves) {
                        mrca_id = nid;
                        break;
                    }
                }
                auto mark_clade = [&](auto& self, int32_t u) -> void {
                    if (u != root_id) target[u] = true;
                    for (int32_t c : nodes[u].children) self(self, c);
                };
                mark_clade(mark_clade, mrca_id);
            } else if (imode == "paraphyletic") {
                for (const auto& nd : nodes) {
                    if (nd.id == root_id || nd.is_leaf) continue;
                    int children_with_sel = 0;
                    for (int32_t cid : nd.children) {
                        if (selected_leaves[cid] > 0) children_with_sel++;
                    }
                    bool outside_has_sel = (global_selected_leaves - selected_leaves[nd.id]) > 0;
                    if ((children_with_sel >= 2) ||
                        (children_with_sel >= 1 && outside_has_sel && selected_leaves[nd.id] == total_leaves[nd.id])) {
                        target[nd.id] = true;
                    }
                }
            } else if (imode == "none") {
                for (const auto& nd : nodes) {
                    if (!nd.is_leaf && !sel.is_selected[nd.id]) {
                        target[nd.id] = false;
                    }
                }
            }
        }

        size_t assigned = 0;
        for (auto& nd : nodes) {
            if (nd.id == root_id) continue;
            if (target[nd.id]) {
                nd.model_tag = clean_tag;
                assigned++;
            }
        }
        return assigned;
    }
};

} // namespace hyphy::core
