#pragma once

#include "hyphy/core/types.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <stack>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>

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

class Tree {
public:
    std::vector<TreeNode> nodes;
    int32_t root_id = INVALID_INDEX;
    std::vector<int32_t> post_order;
    std::unordered_map<std::string, int32_t> leaf_name_to_id;

    static Tree from_newick_file(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open Newick file: " + filepath);
        }
        std::stringstream buffer;
        buffer << file.rdbuf();
        return from_newick(buffer.str());
    }

    static Tree from_newick(std::string_view newick_sv) {
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
        int node_counter = 1;
        for (auto& node : tree.nodes) {
            if (node.parent_id == INVALID_INDEX) {
                tree.root_id = node.id;
            }
            if (node.name.empty()) {
                node.name = "Node" + std::to_string(node_counter++);
            }
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

private:
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
};

} // namespace hyphy::core
