#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "hyphy/analyses/expcm_optimizer.hpp"
#include "hyphy/core/likelihood.hpp"
#include <iostream>
#include <cmath>

using namespace hyphy::core;
using namespace hyphy::analyses;

// Setup test data paths
const std::string ALIGNMENT_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/alignment.fasta";
const std::string TREE_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/tree.newick";
const std::string PREFS_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/prefs.csv";

// Reference parameters (phydms optimal)
const double BETA_REF = 3.38189;
const double OMEGA_REF = 0.977468;
const double KAPPA_REF = 8.75598;
const double LNL_REF = -775.21;

TEST_SUITE("ExpCM L-BFGS Optimizer") {

    TEST_CASE("Load test data") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        CHECK(aln.sequences.size() == 19);
        CHECK(tree.num_nodes() == 36);
        CHECK(prefs.num_sites() > 0);
    }

    TEST_CASE("Optimizer construction") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);
        CHECK(true);  // Just verify construction
    }

    TEST_CASE("Objective function evaluation") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        Eigen::Vector3d log_params;
        log_params(0) = std::log(BETA_REF);
        log_params(1) = std::log(OMEGA_REF);
        log_params(2) = std::log(KAPPA_REF);

        Scalar obj = optimizer.evaluate_objective(log_params);

        // Expected: -ln L = 775.21
        CHECK(std::fabs(obj - 775.21) < 0.2);
    }

    TEST_CASE("Gradient evaluation returns finite values") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        Eigen::Vector3d log_params;
        log_params(0) = std::log(BETA_REF);
        log_params(1) = std::log(OMEGA_REF);
        log_params(2) = std::log(KAPPA_REF);

        Eigen::Vector3d grad = optimizer.evaluate_gradient(log_params);

        CHECK(!std::isnan(grad(0)));
        CHECK(!std::isnan(grad(1)));
        CHECK(!std::isnan(grad(2)));
        CHECK(!std::isinf(grad(0)));
        CHECK(!std::isinf(grad(1)));
        CHECK(!std::isinf(grad(2)));
    }

    TEST_CASE("Single-parameter optimization (beta only)") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        ExpCMParameters initial;
        initial.beta = 3.2;  // Close to optimal
        initial.omega = OMEGA_REF;
        initial.kappa = KAPPA_REF;
        initial.phi = {0.35, 0.20, 0.25, 0.20};

        ExpCMOptimizer::Bounds bounds;
        bounds.beta_min = 0.5;
        bounds.beta_max = 5.0;
        bounds.omega_min = OMEGA_REF * 0.9;
        bounds.omega_max = OMEGA_REF * 1.1;
        bounds.kappa_min = KAPPA_REF * 0.9;
        bounds.kappa_max = KAPPA_REF * 1.1;

        ExpCMOptimizer::ConvergenceCriteria criteria;
        criteria.gradient_tolerance = 1e-3;
        criteria.max_iterations = 50;

        auto result = optimizer.optimize(initial, bounds, criteria);

        CHECK(result.iterations < criteria.max_iterations);
        double final_beta = result.get_real_params()(0);
        CHECK(std::fabs(final_beta - BETA_REF) < 0.3);
    }

    TEST_CASE("Joint optimization (all three parameters)") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        ExpCMParameters initial;
        initial.beta = 3.0;
        initial.omega = 0.9;
        initial.kappa = 8.0;
        initial.phi = {0.35, 0.20, 0.25, 0.20};

        ExpCMOptimizer::Bounds bounds;
        bounds.beta_min = 0.5;
        bounds.beta_max = 5.0;
        bounds.omega_min = 0.5;
        bounds.omega_max = 2.0;
        bounds.kappa_min = 1.0;
        bounds.kappa_max = 15.0;

        ExpCMOptimizer::ConvergenceCriteria criteria;
        criteria.gradient_tolerance = 1e-3;
        criteria.max_iterations = 100;

        auto result = optimizer.optimize(initial, bounds, criteria);

        CHECK(result.iterations < criteria.max_iterations);
        auto final_params = result.get_real_params();

        // Less strict tolerances for joint optimization from distant start
        CHECK(std::fabs(final_params(0) - BETA_REF) < 0.5);
        CHECK(std::fabs(final_params(1) - OMEGA_REF) < 0.2);
        CHECK(std::fabs(final_params(2) - KAPPA_REF) < 0.5);
    }

    TEST_CASE("Convergence from multiple starting points") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::vector<Eigen::Vector3d> final_params_list;

        std::vector<std::tuple<double, double, double>> starts = {
            {3.0, 0.9, 8.0},
            {3.5, 1.0, 9.0},
            {3.2, 0.95, 8.5}
        };

        for (const auto& start : starts) {
            ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

            ExpCMParameters initial;
            initial.beta = std::get<0>(start);
            initial.omega = std::get<1>(start);
            initial.kappa = std::get<2>(start);
            initial.phi = {0.35, 0.20, 0.25, 0.20};

            ExpCMOptimizer::Bounds bounds;
            bounds.beta_min = 0.5;
            bounds.beta_max = 5.0;
            bounds.omega_min = 0.5;
            bounds.omega_max = 2.0;
            bounds.kappa_min = 1.0;
            bounds.kappa_max = 15.0;

            ExpCMOptimizer::ConvergenceCriteria criteria;
            criteria.max_iterations = 100;

            auto result = optimizer.optimize(initial, bounds, criteria);
            final_params_list.push_back(result.get_real_params());
        }

        // All should be reasonably close to each other
        for (size_t i = 1; i < final_params_list.size(); ++i) {
            Eigen::Vector3d diff = final_params_list[i] - final_params_list[0];
            CHECK(diff.norm() < 0.5);
        }
    }

    TEST_CASE("Likelihood increases monotonically") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        ExpCMParameters initial;
        initial.beta = 3.0;
        initial.omega = 0.9;
        initial.kappa = 8.0;
        initial.phi = {0.35, 0.20, 0.25, 0.20};

        ExpCMOptimizer::Bounds bounds;
        bounds.beta_min = 0.5;
        bounds.beta_max = 5.0;
        bounds.omega_min = 0.5;
        bounds.omega_max = 2.0;
        bounds.kappa_min = 1.0;
        bounds.kappa_max = 15.0;

        ExpCMOptimizer::ConvergenceCriteria criteria;
        criteria.max_iterations = 100;

        auto result = optimizer.optimize(initial, bounds, criteria);

        // Check monotonicity (allowing for numerical rounding)
        for (size_t i = 1; i < result.likelihood_history.size(); ++i) {
            CHECK(result.likelihood_history[i] >= result.likelihood_history[i-1] - 1e-6);
        }
    }

    TEST_CASE("Optimization near optimum") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        // Start very close to the optimum
        ExpCMParameters initial;
        initial.beta = BETA_REF * 0.99;
        initial.omega = OMEGA_REF * 0.99;
        initial.kappa = KAPPA_REF * 0.99;
        initial.phi = {0.35, 0.20, 0.25, 0.20};

        ExpCMOptimizer::Bounds bounds;
        bounds.beta_min = 0.5;
        bounds.beta_max = 5.0;
        bounds.omega_min = 0.5;
        bounds.omega_max = 2.0;
        bounds.kappa_min = 1.0;
        bounds.kappa_max = 15.0;

        ExpCMOptimizer::ConvergenceCriteria criteria;
        criteria.gradient_tolerance = 1e-5;
        criteria.max_iterations = 100;

        auto result = optimizer.optimize(initial, bounds, criteria);

        auto final_params = result.get_real_params();

        // Should converge to optimum quickly from nearby point
        CHECK(std::fabs(final_params(0) - BETA_REF) < 0.05);
        CHECK(std::fabs(final_params(1) - OMEGA_REF) < 0.01);
        CHECK(std::fabs(final_params(2) - KAPPA_REF) < 0.05);
        CHECK(std::fabs(result.optimal_likelihood - LNL_REF) < 0.2);
    }

    TEST_CASE("Likelihood landscape check") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== LIKELIHOOD LANDSCAPE CHECK ===\n";

        // Test at reference point
        ExpCMParameters ref_params{BETA_REF, OMEGA_REF, KAPPA_REF, {0.35, 0.20, 0.25, 0.20}};
        auto lnl_ref = LikelihoodEngine::compute_expcm_log_likelihood(tree, aln, gcode, prefs, ref_params);

        // Test at starting point
        ExpCMParameters start_params{3.0, 0.8, 7.0, {0.35, 0.20, 0.25, 0.20}};
        auto lnl_start = LikelihoodEngine::compute_expcm_log_likelihood(tree, aln, gcode, prefs, start_params);

        // Test at convergence point
        ExpCMParameters conv_params{3.76814, 0.888331, 9.72467, {0.35, 0.20, 0.25, 0.20}};
        auto lnl_conv = LikelihoodEngine::compute_expcm_log_likelihood(tree, aln, gcode, prefs, conv_params);

        std::cout << "Reference (3.382, 0.977, 8.756): ln L = " << -lnl_ref << "\n";
        std::cout << "Starting (3.0, 0.8, 7.0):         ln L = " << -lnl_start << "\n";
        std::cout << "Convergence (3.768, 0.888, 9.725): ln L = " << -lnl_conv << "\n\n";

        std::cout << "Comparison:\n";
        std::cout << "  Start > Reference? " << (-lnl_start > -lnl_ref ? "YES" : "NO") << "\n";
        std::cout << "  Conv > Reference? " << (-lnl_conv > -lnl_ref ? "YES" : "NO") << "\n";
        std::cout << "  Conv > Start? " << (-lnl_conv > -lnl_start ? "YES" : "NO") << "\n\n";

        if (-lnl_ref < -lnl_conv) {
            std::cout << "✓ Reference IS at a better minimum than convergence point.\n";
            std::cout << "  This means the optimizer got stuck at a local minimum!\n";
            std::cout << "  Root cause: Likely gradient or line search issue.\n";
        } else {
            std::cout << "✗ Convergence point is BETTER than reference!\n";
            std::cout << "  This means the reference values may be wrong, or\n";
            std::cout << "  the likelihood function differs from phydms.\n";
        }
    }
}
