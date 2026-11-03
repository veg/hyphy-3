#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "hyphy/analyses/expcm_optimizer.hpp"
#include "hyphy/core/likelihood.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace hyphy::core;
using namespace hyphy::analyses;

const std::string ALIGNMENT_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/alignment.fasta";
const std::string TREE_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/tree.newick";
const std::string PREFS_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/prefs.csv";

const double BETA_REF = 3.38189;
const double OMEGA_REF = 0.977468;
const double KAPPA_REF = 8.75598;

TEST_SUITE("Gradient Validation") {

    TEST_CASE("Gradient direction check at reference point") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        // Test at reference point
        ExpCMParameters ref_params;
        ref_params.beta = BETA_REF;
        ref_params.omega = OMEGA_REF;
        ref_params.kappa = KAPPA_REF;
        ref_params.phi = {0.355409, 0.201765, 0.249379, 1.0 - 0.355409 - 0.201765 - 0.249379};

        Eigen::Vector3d log_ref = optimizer.params_to_log(ref_params);
        Eigen::Vector3d grad_ref = optimizer.evaluate_gradient(log_ref);

        std::cout << "\n=== GRADIENT AT REFERENCE POINT ===\n";
        std::cout << "Reference parameters:\n";
        std::cout << "  β    = " << BETA_REF << "\n";
        std::cout << "  ω    = " << OMEGA_REF << "\n";
        std::cout << "  κ    = " << KAPPA_REF << "\n";
        std::cout << "  ln L = " << -optimizer.evaluate_objective(log_ref) << "\n\n";
        
        std::cout << "Gradient at reference:\n";
        std::cout << "  ∂(-lnL)/∂(log β) = " << std::setw(12) << grad_ref(0) << "\n";
        std::cout << "  ∂(-lnL)/∂(log ω) = " << std::setw(12) << grad_ref(1) << "\n";
        std::cout << "  ∂(-lnL)/∂(log κ) = " << std::setw(12) << grad_ref(2) << "\n";
        std::cout << "  ||grad|| = " << grad_ref.norm() << "\n\n";

        // Check gradient magnitude - at optimum should be small
        CHECK(grad_ref.norm() < 1.5);  // Reasonable tolerance for near-optimum
        
        // Test directional derivatives by moving in parameter space
        std::cout << "Testing directional derivatives:\n";
        
        // Increase beta slightly
        Eigen::Vector3d log_beta_up = log_ref;
        log_beta_up(0) += 0.01;
        double f_ref = optimizer.evaluate_objective(log_ref);
        double f_beta_up = optimizer.evaluate_objective(log_beta_up);
        double deriv_beta_fd = (f_beta_up - f_ref) / 0.01;  // Finite difference
        double deriv_beta_grad = grad_ref(0);  // From gradient
        
        std::cout << "  β direction: FD = " << std::setw(10) << deriv_beta_fd 
                  << ", Grad = " << std::setw(10) << deriv_beta_grad
                  << ", Error = " << std::fabs(deriv_beta_fd - deriv_beta_grad) << "\n";

        // Increase omega slightly
        Eigen::Vector3d log_omega_up = log_ref;
        log_omega_up(1) += 0.01;
        double f_omega_up = optimizer.evaluate_objective(log_omega_up);
        double deriv_omega_fd = (f_omega_up - f_ref) / 0.01;
        double deriv_omega_grad = grad_ref(1);
        
        std::cout << "  ω direction: FD = " << std::setw(10) << deriv_omega_fd 
                  << ", Grad = " << std::setw(10) << deriv_omega_grad
                  << ", Error = " << std::fabs(deriv_omega_fd - deriv_omega_grad) << "\n";

        // Increase kappa slightly
        Eigen::Vector3d log_kappa_up = log_ref;
        log_kappa_up(2) += 0.01;
        double f_kappa_up = optimizer.evaluate_objective(log_kappa_up);
        double deriv_kappa_fd = (f_kappa_up - f_ref) / 0.01;
        double deriv_kappa_grad = grad_ref(2);
        
        std::cout << "  κ direction: FD = " << std::setw(10) << deriv_kappa_fd 
                  << ", Grad = " << std::setw(10) << deriv_kappa_grad
                  << ", Error = " << std::fabs(deriv_kappa_fd - deriv_kappa_grad) << "\n";

        // All finite differences should match gradient reasonably
        CHECK(std::fabs(deriv_beta_fd - deriv_beta_grad) < 1.0);
        CHECK(std::fabs(deriv_omega_fd - deriv_omega_grad) < 1.0);
        CHECK(std::fabs(deriv_kappa_fd - deriv_kappa_grad) < 1.0);
    }

    TEST_CASE("Gradient direction check at offset point") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        // Test at offset point (from failed Phase 3C)
        ExpCMParameters offset_params;
        offset_params.beta = 3.0;
        offset_params.omega = 0.8;
        offset_params.kappa = 7.0;
        offset_params.phi = {0.35, 0.20, 0.25, 0.20};

        Eigen::Vector3d log_offset = optimizer.params_to_log(offset_params);
        Eigen::Vector3d grad_offset = optimizer.evaluate_gradient(log_offset);
        double f_offset = optimizer.evaluate_objective(log_offset);

        std::cout << "\n=== GRADIENT AT OFFSET POINT ===\n";
        std::cout << "Offset parameters:\n";
        std::cout << "  β    = " << offset_params.beta << "\n";
        std::cout << "  ω    = " << offset_params.omega << "\n";
        std::cout << "  κ    = " << offset_params.kappa << "\n";
        std::cout << "  ln L = " << -f_offset << "\n\n";
        
        std::cout << "Gradient at offset:\n";
        std::cout << "  ∂(-lnL)/∂(log β) = " << std::setw(12) << grad_offset(0) << "\n";
        std::cout << "  ∂(-lnL)/∂(log ω) = " << std::setw(12) << grad_offset(1) << "\n";
        std::cout << "  ∂(-lnL)/∂(log κ) = " << std::setw(12) << grad_offset(2) << "\n";
        std::cout << "  ||grad|| = " << grad_offset.norm() << "\n\n";

        // Check direction toward reference
        std::cout << "Direction from offset toward reference:\n";
        double delta_log_beta = std::log(BETA_REF) - log_offset(0);
        double delta_log_omega = std::log(OMEGA_REF) - log_offset(1);
        double delta_log_kappa = std::log(KAPPA_REF) - log_offset(2);
        
        std::cout << "  Δ log β = " << std::setw(10) << delta_log_beta << "\n";
        std::cout << "  Δ log ω = " << std::setw(10) << delta_log_omega << "\n";
        std::cout << "  Δ log κ = " << std::setw(10) << delta_log_kappa << "\n\n";
        
        // Dot product with gradient - should be negative (opposite direction) if gradient points to optimum
        double dot_product = grad_offset(0) * delta_log_beta + 
                           grad_offset(1) * delta_log_omega + 
                           grad_offset(2) * delta_log_kappa;
        
        std::cout << "Dot product (grad · direction_to_ref) = " << dot_product << "\n";
        if (dot_product < 0) {
            std::cout << "✓ CORRECT: Gradient opposes direction to reference (points to different optimum)\n";
        } else {
            std::cout << "✗ WRONG: Gradient points toward reference (but we converged away!)\n";
        }
    }

    TEST_CASE("Analytical Frechet gradients match central finite differences to < 1e-5") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMParameters params;
        params.beta = 3.0;
        params.omega = 0.8;
        params.kappa = 7.0;
        params.phi = {0.35, 0.20, 0.25, 0.20};

        Scalar eps = 1e-6;
        Scalar dL_dbeta_frechet = LikelihoodEngine::compute_expcm_beta_gradient_frechet(tree, aln, gcode, prefs, params);
        Scalar dL_dbeta_fd = LikelihoodEngine::compute_expcm_beta_gradient_finite_diff(tree, aln, gcode, prefs, params, 1.0, eps);

        Scalar dL_domega_frechet = LikelihoodEngine::compute_expcm_omega_gradient_frechet(tree, aln, gcode, prefs, params);
        Scalar dL_domega_fd = LikelihoodEngine::compute_expcm_omega_gradient_finite_diff(tree, aln, gcode, prefs, params, 1.0, eps);

        Scalar dL_dkappa_frechet = LikelihoodEngine::compute_expcm_kappa_gradient_frechet(tree, aln, gcode, prefs, params);
        Scalar dL_dkappa_fd = LikelihoodEngine::compute_expcm_kappa_gradient_finite_diff(tree, aln, gcode, prefs, params, 1.0, eps);

        std::cout << "\n=== ANALYTICAL FRECHET VS CENTRAL FINITE DIFFERENCE ===\n";
        std::cout << "  dL/dbeta:  Frechet = " << dL_dbeta_frechet << ", FD = " << dL_dbeta_fd
                  << ", |diff| = " << std::fabs(dL_dbeta_frechet - dL_dbeta_fd) << "\n";
        std::cout << "  dL/domega: Frechet = " << dL_domega_frechet << ", FD = " << dL_domega_fd
                  << ", |diff| = " << std::fabs(dL_domega_frechet - dL_domega_fd) << "\n";
        std::cout << "  dL/dkappa: Frechet = " << dL_dkappa_frechet << ", FD = " << dL_dkappa_fd
                  << ", |diff| = " << std::fabs(dL_dkappa_frechet - dL_dkappa_fd) << "\n";

        CHECK(std::fabs(dL_dbeta_frechet - dL_dbeta_fd) < 1e-5);
        CHECK(std::fabs(dL_domega_frechet - dL_domega_fd) < 1e-5);
        CHECK(std::fabs(dL_dkappa_frechet - dL_dkappa_fd) < 1e-5);
    }
}

