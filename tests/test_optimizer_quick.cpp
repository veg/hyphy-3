#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "hyphy/analyses/expcm_optimizer.hpp"
#include "hyphy/core/likelihood.hpp"
#include <iostream>

using namespace hyphy::core;
using namespace hyphy::analyses;

const std::string ALIGNMENT_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/alignment.fasta";
const std::string TREE_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/tree.newick";
const std::string PREFS_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/prefs.csv";

TEST_SUITE("Quick Optimizer Verification") {

    TEST_CASE("Data loads successfully") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        CHECK(aln.sequences.size() == 19);
        CHECK(tree.num_nodes() == 36);
        CHECK(prefs.num_sites() > 0);

        std::cout << "✓ Data loaded: " << aln.sequences.size() << " taxa, "
                  << tree.num_nodes() << " nodes, "
                  << prefs.num_sites() << " sites\n";
    }

    TEST_CASE("Optimizer constructs") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);
        std::cout << "✓ Optimizer constructed successfully\n";
        CHECK(true);
    }

    TEST_CASE("Objective function works") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        Eigen::Vector3d log_params;
        log_params(0) = std::log(3.38189);
        log_params(1) = std::log(0.977468);
        log_params(2) = std::log(8.75598);

        std::cout << "Computing objective (this may take 10-30 seconds)...\n";
        Scalar obj = optimizer.evaluate_objective(log_params);

        std::cout << "✓ Objective computed: " << obj << " (expected ~775.21)\n";
        CHECK(std::fabs(obj - 775.21) < 1.0);  // Loose tolerance
    }

    TEST_CASE("Gradient function works") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        Eigen::Vector3d log_params;
        log_params(0) = std::log(3.38189);
        log_params(1) = std::log(0.977468);
        log_params(2) = std::log(8.75598);

        std::cout << "Computing gradient (this may take 30-90 seconds)...\n";
        Eigen::Vector3d grad = optimizer.evaluate_gradient(log_params);

        std::cout << "✓ Gradient computed: [" << grad(0) << ", "
                  << grad(1) << ", " << grad(2) << "]\n";

        CHECK(!std::isnan(grad(0)));
        CHECK(!std::isnan(grad(1)));
        CHECK(!std::isnan(grad(2)));
        CHECK(!std::isinf(grad(0)));
        CHECK(!std::isinf(grad(1)));
        CHECK(!std::isinf(grad(2)));
    }

    TEST_CASE("Quick optimization (10 iterations max)") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        ExpCMOptimizer optimizer(tree, aln, gcode, prefs);

        ExpCMParameters initial;
        initial.beta = 3.3;
        initial.omega = 0.97;
        initial.kappa = 8.7;
        initial.phi = {0.35, 0.20, 0.25, 0.20};

        // Use conservative bounds to avoid line search issues
        ExpCMOptimizer::Bounds bounds;
        bounds.beta_min = 2.0;
        bounds.beta_max = 4.5;
        bounds.omega_min = 0.5;
        bounds.omega_max = 1.5;
        bounds.kappa_min = 7.0;
        bounds.kappa_max = 10.0;

        ExpCMOptimizer::ConvergenceCriteria criteria;
        criteria.max_iterations = 3;  // Just 3 iterations for quick check

        std::cout << "Running 3 optimization iterations (this may take 3-5 minutes)...\n";
        auto result = optimizer.optimize(initial, bounds, criteria);

        std::cout << "✓ Optimization completed in " << result.iterations << " iterations\n";
        std::cout << "  Final ln L: " << result.optimal_likelihood << "\n";
        std::cout << "  Final parameters: β=" << result.get_real_params()(0)
                  << ", ω=" << result.get_real_params()(1)
                  << ", κ=" << result.get_real_params()(2) << "\n";
        std::cout << "  Gradient norm: " << result.gradient_magnitude_history.back() << "\n";

        CHECK(result.iterations <= 3);
        CHECK(!std::isnan(result.optimal_likelihood));
        CHECK(!std::isinf(result.optimal_likelihood));
    }
}
