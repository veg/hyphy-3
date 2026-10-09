#pragma once

#include "hyphy/core/types.hpp"
#include <memory>
#include <vector>
#include <functional>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace hyphy::autograd {

using hyphy::core::Scalar;

class Node {
public:
    Scalar val = 0.0;
    Scalar grad = 0.0;
    bool requires_grad = false;
    std::vector<std::shared_ptr<Node>> parents;
    std::function<void(Scalar)> backward_fn;

    Node(Scalar v = 0.0, bool req_grad = false)
        : val(v), grad(0.0), requires_grad(req_grad) {}
};

class Var {
public:
    std::shared_ptr<Node> node;

    Var() : node(std::make_shared<Node>(0.0, false)) {}
    Var(Scalar v, bool req_grad = false) : node(std::make_shared<Node>(v, req_grad)) {}
    Var(std::shared_ptr<Node> n) : node(std::move(n)) {}

    Scalar val() const { return node ? node->val : 0.0; }
    Scalar grad() const { return node ? node->grad : 0.0; }
    void set_val(Scalar v) { if (node) node->val = v; }
    void set_grad(Scalar g) { if (node) node->grad = g; }
    bool requires_grad() const { return node ? node->requires_grad : false; }
    void zero_grad() { if (node) node->grad = 0.0; }

    void backward(Scalar grad_output = 1.0) {
        if (!node) return;

        // 1. Topological sort via post-order DFS
        std::vector<std::shared_ptr<Node>> topo;
        std::unordered_set<Node*> visited;

        std::function<void(const std::shared_ptr<Node>&)> build_topo =
            [&](const std::shared_ptr<Node>& n) {
                if (!n || visited.count(n.get())) return;
                visited.insert(n.get());
                for (const auto& parent : n->parents) {
                    build_topo(parent);
                }
                topo.push_back(n);
            };

        build_topo(node);

        // 2. Initialize root gradient
        node->grad = grad_output;

        // 3. Reverse topological traversal calling backward_fn
        for (auto it = topo.rbegin(); it != topo.rend(); ++it) {
            auto& curr = *it;
            if (curr->backward_fn && std::abs(curr->grad) > 0.0) {
                curr->backward_fn(curr->grad);
            }
        }
    }
};

// Operator Overloads

// Addition: a + b
inline Var operator+(const Var& a, const Var& b) {
    bool req = a.requires_grad() || b.requires_grad();
    auto res_node = std::make_shared<Node>(a.val() + b.val(), req);
    if (req) {
        res_node->parents = {a.node, b.node};
        res_node->backward_fn = [a_node = a.node, b_node = b.node](Scalar g) {
            if (a_node->requires_grad) a_node->grad += g;
            if (b_node->requires_grad) b_node->grad += g;
        };
    }
    return Var(res_node);
}

inline Var operator+(const Var& a, Scalar b) { return a + Var(b); }
inline Var operator+(Scalar a, const Var& b) { return Var(a) + b; }

// Subtraction: a - b
inline Var operator-(const Var& a, const Var& b) {
    bool req = a.requires_grad() || b.requires_grad();
    auto res_node = std::make_shared<Node>(a.val() - b.val(), req);
    if (req) {
        res_node->parents = {a.node, b.node};
        res_node->backward_fn = [a_node = a.node, b_node = b.node](Scalar g) {
            if (a_node->requires_grad) a_node->grad += g;
            if (b_node->requires_grad) b_node->grad -= g;
        };
    }
    return Var(res_node);
}

inline Var operator-(const Var& a, Scalar b) { return a - Var(b); }
inline Var operator-(Scalar a, const Var& b) { return Var(a) - b; }

// Unary negation: -a
inline Var operator-(const Var& a) {
    bool req = a.requires_grad();
    auto res_node = std::make_shared<Node>(-a.val(), req);
    if (req) {
        res_node->parents = {a.node};
        res_node->backward_fn = [a_node = a.node](Scalar g) {
            if (a_node->requires_grad) a_node->grad -= g;
        };
    }
    return Var(res_node);
}

// Multiplication: a * b
inline Var operator*(const Var& a, const Var& b) {
    bool req = a.requires_grad() || b.requires_grad();
    auto res_node = std::make_shared<Node>(a.val() * b.val(), req);
    if (req) {
        res_node->parents = {a.node, b.node};
        Scalar a_val = a.val();
        Scalar b_val = b.val();
        res_node->backward_fn = [a_node = a.node, b_node = b.node, a_val, b_val](Scalar g) {
            if (a_node->requires_grad) a_node->grad += g * b_val;
            if (b_node->requires_grad) b_node->grad += g * a_val;
        };
    }
    return Var(res_node);
}

inline Var operator*(const Var& a, Scalar b) { return a * Var(b); }
inline Var operator*(Scalar a, const Var& b) { return Var(a) * b; }

// Division: a / b
inline Var operator/(const Var& a, const Var& b) {
    bool req = a.requires_grad() || b.requires_grad();
    Scalar b_val = b.val();
    if (std::abs(b_val) < 1e-15) b_val = (b_val >= 0 ? 1e-15 : -1e-15);
    auto res_node = std::make_shared<Node>(a.val() / b_val, req);
    if (req) {
        res_node->parents = {a.node, b.node};
        Scalar a_val = a.val();
        res_node->backward_fn = [a_node = a.node, b_node = b.node, a_val, b_val](Scalar g) {
            if (a_node->requires_grad) a_node->grad += g / b_val;
            if (b_node->requires_grad) b_node->grad -= g * a_val / (b_val * b_val);
        };
    }
    return Var(res_node);
}

inline Var operator/(const Var& a, Scalar b) { return a / Var(b); }
inline Var operator/(Scalar a, const Var& b) { return Var(a) / b; }

// Exponential: exp(a)
inline Var exp(const Var& a) {
    Scalar val = std::exp(a.val());
    bool req = a.requires_grad();
    auto res_node = std::make_shared<Node>(val, req);
    if (req) {
        res_node->parents = {a.node};
        res_node->backward_fn = [a_node = a.node, val](Scalar g) {
            if (a_node->requires_grad) a_node->grad += g * val;
        };
    }
    return Var(res_node);
}

// Natural logarithm: log(a)
inline Var log(const Var& a) {
    Scalar a_val = std::max(a.val(), 1e-15);
    Scalar val = std::log(a_val);
    bool req = a.requires_grad();
    auto res_node = std::make_shared<Node>(val, req);
    if (req) {
        res_node->parents = {a.node};
        res_node->backward_fn = [a_node = a.node, a_val](Scalar g) {
            if (a_node->requires_grad) a_node->grad += g / a_val;
        };
    }
    return Var(res_node);
}

// Power: pow(a, p)
inline Var pow(const Var& a, Scalar p) {
    Scalar a_val = a.val();
    Scalar val = std::pow(a_val, p);
    bool req = a.requires_grad();
    auto res_node = std::make_shared<Node>(val, req);
    if (req) {
        res_node->parents = {a.node};
        res_node->backward_fn = [a_node = a.node, a_val, p](Scalar g) {
            if (a_node->requires_grad) a_node->grad += g * p * std::pow(a_val, p - 1.0);
        };
    }
    return Var(res_node);
}

// Square root: sqrt(a)
inline Var sqrt(const Var& a) {
    return pow(a, 0.5);
}

} // namespace hyphy::autograd
