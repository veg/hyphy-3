#pragma once

#include "hyphy/autograd/var.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"

namespace hyphy::autograd {

using hyphy::core::Tree;
using hyphy::core::Alignment;
using hyphy::core::MG94Parameters;
using hyphy::core::LikelihoodEngine;

// Differentiable Phylogenetic Tree Log-Likelihood Function
// Forward: evaluates ln L(t) under MG94 codon model
// Backward: evaluates exact analytical Inside-Outside adjoint gradients d ln L / d t_b
inline Var tree_log_likelihood(
    const Tree& tree,
    const Alignment& aln,
    const MG94Parameters& params,
    const std::vector<Var>& branch_lengths
) {
    size_t num_nodes = tree.num_nodes();
    std::vector<Scalar> bl(num_nodes, 0.0);
    bool any_req = false;
    for (size_t i = 0; i < num_nodes && i < branch_lengths.size(); ++i) {
        bl[i] = branch_lengths[i].val();
        if (branch_lengths[i].requires_grad()) any_req = true;
    }

    // Call high-performance C++ Inside-Outside adjoint solver
    auto [total_log_l, grads] = LikelihoodEngine::compute_branch_length_gradients(
        tree, aln, params, bl
    );

    auto res_node = std::make_shared<Node>(total_log_l, any_req);

    if (any_req) {
        for (size_t i = 0; i < num_nodes && i < branch_lengths.size(); ++i) {
            if (branch_lengths[i].requires_grad()) {
                res_node->parents.push_back(branch_lengths[i].node);
            }
        }

        res_node->backward_fn = [branch_lengths, grads, root_id = tree.root_id](Scalar g) {
            for (size_t i = 0; i < branch_lengths.size() && i < grads.size(); ++i) {
                if (i != static_cast<size_t>(root_id) && branch_lengths[i].requires_grad()) {
                    branch_lengths[i].node->grad += g * grads[i];
                }
            }
        };
    }

    return Var(res_node);
}

} // namespace hyphy::autograd
