#pragma once

#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/types.hpp"
#include <Eigen/Core>
#include <vector>
#include <string>
#include <cmath>

namespace hyphy::analyses {

class ExpCMOptimizer {
public:
    // Result structure for optimization
    struct Result {
        Eigen::Vector3d optimal_params;        // [log(β), log(ω), log(κ)] in log-space
        hyphy::core::Scalar optimal_likelihood;
        size_t iterations;
        std::string convergence_reason;
        std::vector<hyphy::core::Scalar> likelihood_history;
        std::vector<hyphy::core::Scalar> gradient_magnitude_history;

        // Convert back to original parameter space
        Eigen::Vector3d get_real_params() const {
            return optimal_params.array().exp();
        }
    };

    // Parameter bounds
    struct Bounds {
        double beta_min;
        double beta_max;
        double omega_min;
        double omega_max;
        double kappa_min;
        double kappa_max;

        Bounds()
            : beta_min(1e-5), beta_max(10.0),
              omega_min(1e-5), omega_max(100.0),
              kappa_min(1e-5), kappa_max(100.0) {}
    };

    // Convergence criteria
    struct ConvergenceCriteria {
        double gradient_tolerance;
        double relative_tol;
        size_t max_iterations;
        double line_search_tolerance;
        size_t memory_size;

        ConvergenceCriteria()
            : gradient_tolerance(1e-6),
              relative_tol(1e-10),
              max_iterations(1000),
              line_search_tolerance(1e-4),
              memory_size(10) {}
    };

    // Constructor
    ExpCMOptimizer(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs
    );

    // Main optimization method
    Result optimize(
        const hyphy::core::ExpCMParameters& initial_params,
        const Bounds& bounds = Bounds(),
        const ConvergenceCriteria& criteria = ConvergenceCriteria()
    );

    // For testing/debugging: evaluate objective and gradient at given parameters
    hyphy::core::Scalar evaluate_objective(const Eigen::Vector3d& log_params);
    Eigen::Vector3d evaluate_gradient(const Eigen::Vector3d& log_params);

private:
    // References to input data
    const hyphy::core::Tree& tree_;
    const hyphy::core::Alignment& aln_;
    const hyphy::core::GeneticCode& gcode_;
    const hyphy::core::Preferences& prefs_;

    // Configuration
    Bounds bounds_;
    ConvergenceCriteria criteria_;

    // Working storage
    Eigen::Vector3d current_log_params_;
    hyphy::core::Scalar current_likelihood_;
    Eigen::Vector3d current_gradient_;

    // L-BFGS state
    std::vector<Eigen::Vector3d> past_steps_;    // {s_k} vectors
    std::vector<Eigen::Vector3d> past_gradients_; // {y_k} vectors
    std::vector<double> rho_;                     // {ρ_k} = 1/(y_k^T s_k)

    // Helper methods
    hyphy::core::ExpCMParameters log_to_params(const Eigen::Vector3d& log_params);
    Eigen::Vector3d params_to_log(const hyphy::core::ExpCMParameters& params);

    void enforce_bounds(Eigen::Vector3d& log_params);

    Eigen::Vector3d compute_lbfgs_direction(const Eigen::Vector3d& gradient);
    double line_search_strong_wolfe(
        const Eigen::Vector3d& direction,
        hyphy::core::Scalar f0,
        const Eigen::Vector3d& grad0
    );

    void update_lbfgs_state(
        const Eigen::Vector3d& s,
        const Eigen::Vector3d& y
    );
};

} // namespace hyphy::analyses
