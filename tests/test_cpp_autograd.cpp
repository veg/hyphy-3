#include "hyphy/autograd/var.hpp"
#include "hyphy/autograd/tree_likelihood.hpp"
#include "hyphy/autograd/optimizer.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>

using namespace hyphy::core;
using namespace hyphy::autograd;

void test_basic_graph_derivatives() {
    std::cout << "[Test 1] Testing basic C++ autograd computational graph..." << std::endl;

    // f(x, y) = x * x * y + exp(x) - log(y)
    // df/dx = 2 * x * y + exp(x)
    // df/dy = x * x - 1 / y
    Scalar x_val = 1.5;
    Scalar y_val = 2.0;

    Var x(x_val, true);
    Var y(y_val, true);

    Var f = x * x * y + exp(x) - log(y);

    Scalar expected_f = x_val * x_val * y_val + std::exp(x_val) - std::log(y_val);
    assert(std::abs(f.val() - expected_f) < 1e-12);

    f.backward();

    Scalar expected_df_dx = 2.0 * x_val * y_val + std::exp(x_val);
    Scalar expected_df_dy = x_val * x_val - 1.0 / y_val;

    Scalar err_x = std::abs(x.grad() - expected_df_dx);
    Scalar err_y = std::abs(y.grad() - expected_df_dy);

    std::cout << "  df/dx: computed = " << x.grad() << ", expected = " << expected_df_dx << ", err = " << err_x << std::endl;
    std::cout << "  df/dy: computed = " << y.grad() << ", expected = " << expected_df_dy << ", err = " << err_y << std::endl;

    assert(err_x < 1e-12);
    assert(err_y < 1e-12);
    std::cout << "  -> PASSED!" << std::endl;
}

void test_tree_likelihood_adjoint_gradients(const std::string& filepath) {
    std::cout << "[Test 2] Testing tree log-likelihood analytical adjoint gradients vs finite differences..." << std::endl;

    Alignment aln = Alignment::load(filepath);
    Tree tree = Tree::from_newick(aln.embedded_tree_newick);
    MG94Parameters params;
    params.alpha = 1.0;
    params.beta = 0.4;

    size_t num_nodes = tree.num_nodes();
    std::vector<Scalar> init_lengths(num_nodes, 0.05);
    init_lengths[tree.root_id] = 0.0;

    std::vector<Var> branch_vars;
    for (size_t i = 0; i < num_nodes; ++i) {
        branch_vars.emplace_back(init_lengths[i], i != static_cast<size_t>(tree.root_id));
    }

    Var ll = tree_log_likelihood(tree, aln, params, branch_vars);
    std::cout << "  Tree log-likelihood at t=0.05: " << ll.val() << std::endl;

    ll.backward();

    // Check against central finite difference for each non-root branch
    Scalar eps = 1e-6;
    Scalar max_rel_err = 0.0;

    for (size_t nid = 0; nid < num_nodes; ++nid) {
        if (nid == static_cast<size_t>(tree.root_id)) continue;

        auto bl_plus = init_lengths;
        bl_plus[nid] += eps;
        auto [ll_plus, _] = LikelihoodEngine::compute_branch_length_gradients(tree, aln, params, bl_plus);

        auto bl_minus = init_lengths;
        bl_minus[nid] -= eps;
        auto [ll_minus, __] = LikelihoodEngine::compute_branch_length_gradients(tree, aln, params, bl_minus);

        Scalar num_grad = (ll_plus - ll_minus) / (2.0 * eps);
        Scalar ana_grad = branch_vars[nid].grad();

        Scalar denom = std::max({std::abs(ana_grad), std::abs(num_grad), 1e-4});
        Scalar rel_err = std::abs(ana_grad - num_grad) / denom;

        if (rel_err > max_rel_err) max_rel_err = rel_err;

        if (nid < 8) {
            std::cout << "  Branch " << nid << " (" << tree.nodes[nid].name << "): ana = "
                      << ana_grad << ", num = " << num_grad << ", rel_err = " << rel_err << std::endl;
        }

        assert(rel_err < 1e-5);
    }

    std::cout << "  Max relative error across all branches: " << max_rel_err << std::endl;
    assert(max_rel_err < 1e-5);
    std::cout << "  -> PASSED!" << std::endl;
}

void test_composite_regularized_loss(const std::string& filepath) {
    std::cout << "[Test 3] Testing composite regularized loss with prior in C++ autograd..." << std::endl;

    Alignment aln = Alignment::load(filepath);
    Tree tree = Tree::from_newick(aln.embedded_tree_newick);
    MG94Parameters params;

    size_t num_nodes = tree.num_nodes();
    std::vector<Var> branch_vars;
    for (size_t i = 0; i < num_nodes; ++i) {
        branch_vars.emplace_back(0.05, i != static_cast<size_t>(tree.root_id));
    }

    Var ll = tree_log_likelihood(tree, aln, params, branch_vars);

    // Add L2 penalty: loss = -ll + lambda * sum(t_i^2)
    Scalar lambda = 10.0;
    Var penalty = 0.0;
    for (size_t i = 0; i < num_nodes; ++i) {
        if (i != static_cast<size_t>(tree.root_id)) {
            penalty = penalty + branch_vars[i] * branch_vars[i];
        }
    }

    Var loss = -ll + lambda * penalty;
    loss.backward();

    // Check that grad = -dL_dt + 2 * lambda * t
    auto [_, ana_grads] = LikelihoodEngine::compute_branch_length_gradients(
        tree, aln, params, std::vector<Scalar>(num_nodes, 0.05)
    );

    for (size_t i = 0; i < num_nodes; ++i) {
        if (i == static_cast<size_t>(tree.root_id)) continue;
        Scalar expected_grad = -ana_grads[i] + 2.0 * lambda * branch_vars[i].val();
        Scalar err = std::abs(branch_vars[i].grad() - expected_grad);
        assert(err < 1e-10);
    }

    std::cout << "  -> PASSED!" << std::endl;
}

void test_cpp_adam_optimization(const std::string& filepath) {
    std::cout << "[Test 4] Testing native C++ Adam branch length optimizer..." << std::endl;

    Alignment aln = Alignment::load(filepath);
    Tree tree = Tree::from_newick(aln.embedded_tree_newick);
    MG94Parameters params;
    params.alpha = 1.0;
    params.beta = 0.4;

    size_t num_nodes = tree.num_nodes();
    // Parameterize in log space: t = exp(log_t)
    std::vector<Var> log_t;
    std::vector<Var*> opt_params;

    for (size_t i = 0; i < num_nodes; ++i) {
        bool req = (i != static_cast<size_t>(tree.root_id));
        log_t.emplace_back(std::log(0.05), req);
    }
    for (size_t i = 0; i < num_nodes; ++i) {
        if (i != static_cast<size_t>(tree.root_id)) {
            opt_params.push_back(&log_t[i]);
        }
    }

    Adam optimizer(opt_params, /*lr=*/0.05);

    // Initial forward pass
    std::vector<Var> init_branches;
    for (size_t i = 0; i < num_nodes; ++i) {
        init_branches.push_back(i != static_cast<size_t>(tree.root_id) ? exp(log_t[i]) : Var(0.0));
    }
    Var init_ll = tree_log_likelihood(tree, aln, params, init_branches);
    Scalar initial_val = init_ll.val();
    std::cout << "  Initial log-likelihood: " << initial_val << std::endl;

    Scalar final_val = initial_val;
    for (int step = 1; step <= 10; ++step) {
        optimizer.zero_grad();

        std::vector<Var> branches;
        branches.reserve(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i) {
            branches.push_back(i != static_cast<size_t>(tree.root_id) ? exp(log_t[i]) : Var(0.0));
        }

        Var loss = -tree_log_likelihood(tree, aln, params, branches);
        final_val = -loss.val();
        loss.backward();

        optimizer.step();

        if (step % 2 == 0) {
            std::cout << "  Step " << step << ": log-likelihood = " << final_val << std::endl;
        }
    }

    std::cout << "  Final log-likelihood: " << final_val << " (delta: " << (final_val - initial_val) << ")" << std::endl;
    assert(final_val > initial_val + 50.0);
    std::cout << "  -> PASSED!" << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "   HYPHY 3 C++ AUTOGRAD TEST SUITE      " << std::endl;
    std::cout << "========================================" << std::endl;

    test_basic_graph_derivatives();

    std::string test_file = "tests/data/adh.nex";
    if (!std::filesystem::exists(test_file)) {
        test_file = "../tests/data/adh.nex";
    }
    if (!std::filesystem::exists(test_file)) {
        test_file = "../../tests/data/adh.nex";
    }

    if (std::filesystem::exists(test_file)) {
        test_tree_likelihood_adjoint_gradients(test_file);
        test_composite_regularized_loss(test_file);
        test_cpp_adam_optimization(test_file);
    } else {
        std::cerr << "Warning: Could not find adh.nex test file at " << test_file << std::endl;
    }

    std::cout << "\nALL C++ AUTOGRAD TESTS PASSED SUCCESSFULLY!" << std::endl;
    return 0;
}
