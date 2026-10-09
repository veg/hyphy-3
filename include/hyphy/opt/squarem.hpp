#pragma once

#include "hyphy/core/types.hpp"
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace hyphy::opt {

using namespace hyphy::core;

/**
 * SQUAREM: Squared Polynomial Extrapolation of Expectation-Maximization
 * (Varadhan and Roland, 2008, Scandinavian Journal of Statistics)
 * 
 * Accelerates any contractive fixed-point map F(x) (such as an EM step)
 * using Cauchy / Barzilai-Borwein step lengths with backtracking stabilization.
 */
class SQUAREM {
public:
    struct Result {
        Vector x_opt;
        Scalar log_likelihood = -1e300;
        int iterations = 0;
        bool converged = false;
    };

    /**
     * Minimize negative log-likelihood (maximize log-likelihood) using SQUAREM.
     * @param F EM fixed-point operator x -> F(x)
     * @param eval_lnL Objective function evaluating log-likelihood ln L(x)
     * @param x0 Initial parameter vector
     * @param lb Lower bound vector
     * @param ub Upper bound vector
     * @param tol Convergence tolerance on |ln L_new - ln L_old| or |x_new - x_old|
     * @param max_iter Maximum number of SQUAREM macro-iterations
     */
    static Result optimize(
        std::function<Vector(const Vector&)> F,
        std::function<Scalar(const Vector&)> eval_lnL,
        const Vector& x0,
        const Vector& lb,
        const Vector& ub,
        Scalar tol = 1e-5,
        int max_iter = 100
    ) {
        Result res;
        Vector x = x0;
        // Clamp to bounds
        for (Eigen::Index i = 0; i < x.size(); ++i) {
            x(i) = std::clamp(x(i), lb(i), ub(i));
        }

        Scalar cur_lnL = eval_lnL(x);
        res.x_opt = x;
        res.log_likelihood = cur_lnL;

        for (int iter = 0; iter < max_iter; ++iter) {
            res.iterations = iter + 1;

            // Two EM steps
            Vector x1 = F(x);
            for (Eigen::Index i = 0; i < x1.size(); ++i) {
                x1(i) = std::clamp(x1(i), lb(i), ub(i));
            }

            Vector x2 = F(x1);
            for (Eigen::Index i = 0; i < x2.size(); ++i) {
                x2(i) = std::clamp(x2(i), lb(i), ub(i));
            }

            Vector r = x1 - x;
            Vector v = (x2 - x1) - r;

            Scalar r_norm2 = r.squaredNorm();
            Scalar v_norm2 = v.squaredNorm();

            if (r_norm2 < 1e-16) {
                // Already at fixed point
                res.converged = true;
                res.x_opt = x2;
                res.log_likelihood = eval_lnL(x2);
                break;
            }

            Vector x_acc = x2;
            if (v_norm2 > 1e-16) {
                // Barzilai-Borwein / Cauchy step size (SQUAREM-1 / S1)
                Scalar alpha = -std::sqrt(r_norm2 / v_norm2);

                // Quadratic extrapolation: x_acc = x - 2 * alpha * r + alpha^2 * v
                x_acc = x - 2.0 * alpha * r + alpha * alpha * v;
                for (Eigen::Index i = 0; i < x_acc.size(); ++i) {
                    x_acc(i) = std::clamp(x_acc(i), lb(i), ub(i));
                }

                // Stabilize with one EM step from the extrapolated point
                Vector x_stab = F(x_acc);
                for (Eigen::Index i = 0; i < x_stab.size(); ++i) {
                    x_stab(i) = std::clamp(x_stab(i), lb(i), ub(i));
                }

                Scalar stab_lnL = eval_lnL(x_stab);
                if (stab_lnL >= cur_lnL) {
                    x = x_stab;
                    Scalar delta_L = stab_lnL - cur_lnL;
                    cur_lnL = stab_lnL;
                    res.x_opt = x;
                    res.log_likelihood = cur_lnL;
                    if (delta_L < tol) {
                        res.converged = true;
                        break;
                    }
                    continue;
                }
            }

            // Fallback: standard EM update x2
            Scalar x2_lnL = eval_lnL(x2);
            Scalar delta_L = x2_lnL - cur_lnL;
            cur_lnL = x2_lnL;
            x = x2;
            res.x_opt = x;
            res.log_likelihood = cur_lnL;

            if (delta_L < tol && delta_L >= 0.0) {
                res.converged = true;
                break;
            }
        }

        return res;
    }
};

} // namespace hyphy::opt
