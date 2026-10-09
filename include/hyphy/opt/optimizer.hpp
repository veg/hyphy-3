#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/param_map.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "LBFGSB.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace hyphy::opt {

using namespace hyphy::core;

struct GTRFitResult {
    Scalar log_likelihood = 0.0;
    Scalar aicc = 0.0;
    GTRParameters params;
    Tree tree;
    int iterations = 0;
};

// Fits nucleotide GTR model (substitution rates + branch lengths) on an alignment
class GTRFitter {
public:
    Tree tree;
    Alignment aln;
    Alignment nuc_aln;

    GTRFitter(Tree t, const Alignment& a)
        : tree(std::move(t)), aln(a), nuc_aln(a.to_nucleotide_alignment()) {}

    GTRFitResult fit() {
        // Collect non-root branch node indices
        std::vector<int32_t> branch_node_ids;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                branch_node_ids.push_back(node.id);
            }
        }
        size_t num_branches = branch_node_ids.size();

        GTRParameters params;
        params.theta_AC = 0.5;
        params.theta_AT = 0.3;
        params.theta_CG = 0.5;
        params.theta_CT = 1.0;
        params.theta_GT = 0.3;

        // Initialize branch lengths to 0.05 if not set or zero
        Tree curr_tree = tree;
        for (size_t b = 0; b < num_branches; ++b) {
            int32_t nid = branch_node_ids[b];
            if (curr_tree.nodes[nid].branch_length < 1e-4) {
                curr_tree.nodes[nid].branch_length = 0.05;
            }
        }

        auto eval_tree_lnl = [&](const Tree& t, const GTRParameters& p) -> Scalar {
            GTRMatrix gtr;
            gtr.update(p, nuc_aln.nuc_frequencies);
            return LikelihoodEngine::compute_gtr_log_likelihood(t, nuc_aln, gtr);
        };

        // Stage 1: Optimize global tree scale factor s
        auto scale_obj = [&](Scalar s) -> Scalar {
            Tree t_scale = curr_tree;
            for (size_t b = 0; b < num_branches; ++b) {
                int32_t nid = branch_node_ids[b];
                t_scale.nodes[nid].branch_length *= s;
            }
            return -eval_tree_lnl(t_scale, params);
        };

        auto [opt_s, opt_s_nll] = Brent1D::minimize(scale_obj, 0.01, 1.0, 50.0, 1e-3, 30);
        for (size_t b = 0; b < num_branches; ++b) {
            int32_t nid = branch_node_ids[b];
            curr_tree.nodes[nid].branch_length *= opt_s;
        }

        // Stage 2 & 3: Alternating coordinate-wise optimization of GTR rates and individual branch lengths
        int total_evals = 0;
        for (int cycle = 0; cycle < 3; ++cycle) {
            // Optimize 5 GTR rate parameters coordinate-wise
            auto opt_rate = [&](Scalar& r_ref) {
                auto rate_obj = [&](Scalar val) -> Scalar {
                    Scalar old_val = r_ref;
                    r_ref = val;
                    Scalar lnl = eval_tree_lnl(curr_tree, params);
                    r_ref = old_val;
                    return -lnl;
                };
                auto [best_r, _] = Brent1D::minimize(rate_obj, 1e-4, r_ref, 20.0, 1e-3, 20);
                r_ref = best_r;
                total_evals += 20;
            };

            opt_rate(params.theta_AC);
            opt_rate(params.theta_AT);
            opt_rate(params.theta_CG);
            opt_rate(params.theta_CT);
            opt_rate(params.theta_GT);

            // Optimize individual branch lengths coordinate-wise
            for (size_t b = 0; b < num_branches; ++b) {
                int32_t nid = branch_node_ids[b];
                Scalar init_bl = curr_tree.nodes[nid].branch_length;
                auto bl_obj = [&](Scalar bl) -> Scalar {
                    Scalar old_bl = curr_tree.nodes[nid].branch_length;
                    curr_tree.nodes[nid].branch_length = bl;
                    Scalar lnl = eval_tree_lnl(curr_tree, params);
                    curr_tree.nodes[nid].branch_length = old_bl;
                    return -lnl;
                };
                auto [best_bl, _] = Brent1D::minimize(bl_obj, 1e-6, std::max(init_bl, 1e-4), 5.0, 1e-3, 20);
                curr_tree.nodes[nid].branch_length = best_bl;
                total_evals += 20;
            }
        }

        Scalar final_lnl = eval_tree_lnl(curr_tree, params);

        GTRFitResult res;
        res.log_likelihood = final_lnl;
        res.params = params;
        res.tree = curr_tree;
        res.iterations = 3;

        size_t K = 5 + num_branches;
        size_t N = nuc_aln.num_nucleotides;
        res.aicc = -2.0 * res.log_likelihood + 2.0 * K * N / (N - K - 1.0);
        return res;
    }
};

struct FitResult {
    Scalar log_likelihood = 0.0;
    Vector x_opt;
    std::vector<std::string> param_names;
    int iterations = 0;
    bool converged = false;
    Tree tree;
};

// Fits global parameters (e.g. omega, transition rates) on a fixed tree
class MG94Fitter {
public:
    Tree tree;
    Alignment aln;

    MG94Fitter(Tree t, Alignment a) : tree(std::move(t)), aln(std::move(a)) {}

    // Functor for L-BFGS-B: minimizes negative log-likelihood
    struct Objective {
        const Tree& tree;
        const Alignment& aln;
        MG94Parameters base_params;
        bool fit_rates = false;

        Scalar operator()(const Vector& x, Vector& grad) {
            Scalar omega = x(0);

            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omega;
            if (fit_rates && x.size() >= 6) {
                p.theta_AC = x(1);
                p.theta_AT = x(2);
                p.theta_CG = x(3);
                p.theta_CT = x(4);
                p.theta_GT = x(5);
            }

            auto eval_lnl = [&](const MG94Parameters& params) -> Scalar {
                MG94Matrix mg;
                mg.update(params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *(aln.code ? aln.code : GeneticCode::universal()));
                return LikelihoodEngine::compute_mg94_log_likelihood(tree, aln, mg);
            };

            Scalar f_val = -eval_lnl(p);

            grad.resize(x.size());
            for (int i = 0; i < x.size(); ++i) {
                Scalar h = std::max(1e-6, 1e-5 * std::abs(x(i)));
                Vector x_plus = x;
                Vector x_minus = x;
                x_plus(i) += h;
                x_minus(i) -= h;

                MG94Parameters p_plus = p;
                MG94Parameters p_minus = p;

                if (i == 0) {
                    p_plus.beta = x_plus(0);
                    p_minus.beta = x_minus(0);
                } else if (fit_rates) {
                    if (i == 1) { p_plus.theta_AC = x_plus(1); p_minus.theta_AC = x_minus(1); }
                    if (i == 2) { p_plus.theta_AT = x_plus(2); p_minus.theta_AT = x_minus(2); }
                    if (i == 3) { p_plus.theta_CG = x_plus(3); p_minus.theta_CG = x_minus(3); }
                    if (i == 4) { p_plus.theta_CT = x_plus(4); p_minus.theta_CT = x_minus(4); }
                    if (i == 5) { p_plus.theta_GT = x_plus(5); p_minus.theta_GT = x_minus(5); }
                }

                Scalar f_plus = -eval_lnl(p_plus);
                Scalar f_minus = -eval_lnl(p_minus);
                grad(i) = (f_plus - f_minus) / (2.0 * h);
            }

            return f_val;
        }
    };

    FitResult fit_omega(Scalar init_omega = 1.0, MG94Parameters base_params = MG94Parameters{}) {
        LBFGSpp::LBFGSBParam<Scalar> param;
        param.m = 6;
        param.epsilon = 1e-5;
        param.max_iterations = 50;

        LBFGSpp::LBFGSBSolver<Scalar> solver(param);

        Objective obj{tree, aln, base_params, false};

        Vector x(1);
        x(0) = init_omega;

        Vector lb(1), ub(1);
        lb(0) = 1e-4;
        ub(0) = 100.0;

        Scalar fx = 0.0;
        int niter = 0;
        try {
            niter = solver.minimize(obj, x, fx, lb, ub);
        } catch (const std::exception& e) {
            std::cerr << "L-BFGS-B warning: " << e.what() << "\n";
        }

        FitResult res;
        res.log_likelihood = -fx;
        res.x_opt = x;
        res.param_names = {"omega"};
        res.iterations = niter;
        res.converged = true;
        res.tree = tree;
        return res;
    }

    // Jointly fits omega and global tree branch length scaling factor s
    FitResult fit_omega_and_scale(Scalar init_omega = 1.0, MG94Parameters base_params = MG94Parameters{}) {
        auto eval_at = [&](Scalar omega, Scalar s) -> Scalar {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omega;

            MG94Matrix mg;
            mg.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *(aln.code ? aln.code : GeneticCode::universal()));
            Scalar conv = (mg.scale_factor > 1e-12) ? (3.0 / mg.scale_factor) : 1.0;

            Tree t_eval = tree;
            for (auto& node : t_eval.nodes) {
                if (node.id != t_eval.root_id) {
                    node.branch_length *= (s * conv);
                }
            }

            return -LikelihoodEngine::compute_mg94_log_likelihood(t_eval, aln, mg);
        };

        auto opt_pt = NelderMead2D::minimize(eval_at, init_omega, 1.0, 1e-4, 50.0, 1e-4, 80);

        FitResult res;
        res.log_likelihood = -opt_pt.f;
        res.x_opt.resize(2);
        res.x_opt(0) = opt_pt.a; // omega
        res.x_opt(1) = opt_pt.b; // scale s
        res.param_names = {"omega", "tree_scale"};
        res.iterations = 50;
        res.converged = true;

        res.tree = tree;
        for (auto& node : res.tree.nodes) {
            if (node.id != res.tree.root_id) {
                node.branch_length *= opt_pt.b;
            }
        }
        return res;
    }
};

} // namespace hyphy::opt
