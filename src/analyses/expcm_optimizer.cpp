#include "hyphy/analyses/expcm_optimizer.hpp"
#include "hyphy/core/likelihood.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace hyphy::analyses {

using hyphy::core::Scalar;
using hyphy::core::GeneticCode;
using hyphy::core::Preferences;
using hyphy::core::ExpCMParameters;
using hyphy::core::Tree;
using hyphy::core::Alignment;
using hyphy::core::LikelihoodEngine;

ExpCMOptimizer::ExpCMOptimizer(
    const Tree& tree,
    const Alignment& aln,
    const GeneticCode& gcode,
    const Preferences& prefs
) : tree_(tree), aln_(aln), gcode_(gcode), prefs_(prefs) {}

ExpCMParameters ExpCMOptimizer::log_to_params(const Eigen::Vector3d& log_params) {
    ExpCMParameters params;
    params.beta = std::exp(log_params(0));
    params.omega = std::exp(log_params(1));
    params.kappa = std::exp(log_params(2));
    params.phi = {0.35, 0.20, 0.25, 0.20};
    return params;
}

Eigen::Vector3d ExpCMOptimizer::params_to_log(const ExpCMParameters& params) {
    Eigen::Vector3d log_params;
    log_params(0) = std::log(params.beta);
    log_params(1) = std::log(params.omega);
    log_params(2) = std::log(params.kappa);
    return log_params;
}

void ExpCMOptimizer::enforce_bounds(Eigen::Vector3d& log_params) {
    log_params(0) = std::clamp(log_params(0), std::log(bounds_.beta_min), std::log(bounds_.beta_max));
    log_params(1) = std::clamp(log_params(1), std::log(bounds_.omega_min), std::log(bounds_.omega_max));
    log_params(2) = std::clamp(log_params(2), std::log(bounds_.kappa_min), std::log(bounds_.kappa_max));
}

Scalar ExpCMOptimizer::evaluate_objective(const Eigen::Vector3d& log_params) {
    ExpCMParameters params = log_to_params(log_params);
    Scalar lnL = LikelihoodEngine::compute_expcm_log_likelihood(tree_, aln_, gcode_, prefs_, params);
    return -lnL;
}

Eigen::Vector3d ExpCMOptimizer::evaluate_gradient(const Eigen::Vector3d& log_params) {
    ExpCMParameters params = log_to_params(log_params);

    Scalar grad_beta = LikelihoodEngine::compute_expcm_beta_gradient_finite_diff(
        tree_, aln_, gcode_, prefs_, params, 1.0, 1e-5
    );
    Scalar grad_omega = LikelihoodEngine::compute_expcm_omega_gradient_finite_diff(
        tree_, aln_, gcode_, prefs_, params, 1.0, 1e-5
    );
    Scalar grad_kappa = LikelihoodEngine::compute_expcm_kappa_gradient_finite_diff(
        tree_, aln_, gcode_, prefs_, params, 1.0, 1e-5
    );

    Eigen::Vector3d gradient;
    gradient(0) = -grad_beta * params.beta;
    gradient(1) = -grad_omega * params.omega;
    gradient(2) = -grad_kappa * params.kappa;

    return gradient;
}

Eigen::Vector3d ExpCMOptimizer::compute_lbfgs_direction(const Eigen::Vector3d& gradient) {
    int m = std::min(static_cast<int>(past_steps_.size()), static_cast<int>(criteria_.memory_size));

    Eigen::Vector3d q = gradient;
    std::vector<double> alphas(m);
    for (int i = m - 1; i >= 0; --i) {
        alphas[i] = rho_[i] * past_steps_[i].dot(q);
        q = q - alphas[i] * past_gradients_[i];
    }

    double gamma = 1.0;
    if (m > 0) {
        double ys = past_gradients_[m - 1].dot(past_steps_[m - 1]);
        double yy = past_gradients_[m - 1].dot(past_gradients_[m - 1]);
        if (yy > 1e-16) gamma = ys / yy;
    }

    q = gamma * q;

    for (int i = 0; i < m; ++i) {
        double beta = rho_[i] * past_gradients_[i].dot(q);
        q = q + (alphas[i] - beta) * past_steps_[i];
    }

    return -q;
}

double ExpCMOptimizer::line_search_strong_wolfe(
    const Eigen::Vector3d& direction,
    Scalar f0,
    const Eigen::Vector3d& grad0
) {
    const double c1 = 1e-4;
    const double c2 = 0.9;
    const int max_evals = 10;

    double step_size = 1.0;

    for (int i = 0; i < max_evals; ++i) {
        Eigen::Vector3d new_params = current_log_params_ + step_size * direction;
        enforce_bounds(new_params);

        Scalar f_new = evaluate_objective(new_params);

        if (f_new <= f0 + c1 * step_size * grad0.dot(direction)) {
            Eigen::Vector3d grad_new = evaluate_gradient(new_params);

            if (std::fabs(grad_new.dot(direction)) <= c2 * std::fabs(grad0.dot(direction))) {
                return step_size;
            }
        }

        step_size *= 0.5;
    }

    return step_size;
}

void ExpCMOptimizer::update_lbfgs_state(
    const Eigen::Vector3d& s,
    const Eigen::Vector3d& y
) {
    double sy = s.dot(y);

    if (sy <= 1e-16) {
        return;
    }

    if (past_steps_.size() >= criteria_.memory_size) {
        past_steps_.erase(past_steps_.begin());
        past_gradients_.erase(past_gradients_.begin());
        rho_.erase(rho_.begin());
    }

    past_steps_.push_back(s);
    past_gradients_.push_back(y);
    rho_.push_back(1.0 / sy);
}

ExpCMOptimizer::Result ExpCMOptimizer::optimize(
    const ExpCMParameters& initial_params,
    const Bounds& bounds,
    const ConvergenceCriteria& criteria
) {
    bounds_ = bounds;
    criteria_ = criteria;

    Result result;
    result.likelihood_history.clear();
    result.gradient_magnitude_history.clear();

    current_log_params_ = params_to_log(initial_params);
    enforce_bounds(current_log_params_);

    current_likelihood_ = evaluate_objective(current_log_params_);
    current_gradient_ = evaluate_gradient(current_log_params_);

    std::cout << "L-BFGS Optimization Starting\n";
    std::cout << "Initial parameters: β=" << std::exp(current_log_params_(0))
              << ", ω=" << std::exp(current_log_params_(1))
              << ", κ=" << std::exp(current_log_params_(2)) << "\n";
    std::cout << "Initial objective (negative ln L): " << current_likelihood_ << "\n";
    std::cout << std::string(80, '-') << "\n";

    for (size_t iteration = 0; iteration < criteria_.max_iterations; ++iteration) {
        Eigen::Vector3d direction = compute_lbfgs_direction(current_gradient_);

        double step_size = line_search_strong_wolfe(direction, current_likelihood_, current_gradient_);

        Eigen::Vector3d new_log_params = current_log_params_ + step_size * direction;
        enforce_bounds(new_log_params);

        Scalar new_likelihood = evaluate_objective(new_log_params);
        Eigen::Vector3d new_gradient = evaluate_gradient(new_log_params);

        Eigen::Vector3d s = new_log_params - current_log_params_;
        Eigen::Vector3d y = new_gradient - current_gradient_;
        update_lbfgs_state(s, y);

        result.likelihood_history.push_back(-new_likelihood);
        double grad_magnitude = new_gradient.norm();
        result.gradient_magnitude_history.push_back(grad_magnitude);

        if (iteration % 10 == 0 || iteration == 0) {
            std::cout << "Iter " << std::setw(4) << iteration
                      << ": ln L = " << std::setw(10) << -new_likelihood
                      << ", ||grad|| = " << std::setw(12) << grad_magnitude
                      << ", β=" << std::setw(8) << std::exp(new_log_params(0))
                      << ", ω=" << std::setw(8) << std::exp(new_log_params(1))
                      << ", κ=" << std::setw(8) << std::exp(new_log_params(2)) << "\n";
        }

        double grad_norm = new_gradient.norm();
        if (grad_norm < criteria_.gradient_tolerance) {
            current_log_params_ = new_log_params;
            current_likelihood_ = new_likelihood;
            current_gradient_ = new_gradient;
            result.convergence_reason = "Gradient norm below tolerance";
            result.iterations = iteration + 1;
            std::cout << "\n✓ Converged: " << result.convergence_reason << "\n";
            break;
        }

        double relative_change = std::fabs(current_likelihood_ - new_likelihood) / (std::fabs(current_likelihood_) + 1e-10);
        if (relative_change < criteria_.relative_tol) {
            current_log_params_ = new_log_params;
            current_likelihood_ = new_likelihood;
            current_gradient_ = new_gradient;
            result.convergence_reason = "Relative change below tolerance";
            result.iterations = iteration + 1;
            std::cout << "\n✓ Converged: " << result.convergence_reason << "\n";
            break;
        }

        current_log_params_ = new_log_params;
        current_likelihood_ = new_likelihood;
        current_gradient_ = new_gradient;
    }

    result.optimal_params = current_log_params_;
    result.optimal_likelihood = -current_likelihood_;

    if (result.iterations == 0) {
        result.convergence_reason = "Max iterations reached";
        result.iterations = criteria_.max_iterations;
    }

    std::cout << std::string(80, '-') << "\n";
    std::cout << "Final parameters: β=" << std::setw(10) << std::exp(current_log_params_(0))
              << ", ω=" << std::setw(10) << std::exp(current_log_params_(1))
              << ", κ=" << std::setw(10) << std::exp(current_log_params_(2)) << "\n";
    std::cout << "Final objective (ln L): " << std::setw(12) << -current_likelihood_ << "\n";
    std::cout << "Final gradient norm: " << std::setw(12) << current_gradient_.norm() << "\n";
    std::cout << "Converged in " << result.iterations << " iterations\n";

    return result;
}

} // namespace hyphy::analyses
