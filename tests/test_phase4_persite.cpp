#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "hyphy/core/likelihood.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace hyphy::core;

// Test data paths
const std::string ALIGNMENT_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/alignment.fasta";
const std::string TREE_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/tree.newick";
const std::string PREFS_FILE = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/prefs.csv";

// Reference parameters from Phase 3C
const double BETA_REF = 3.367;
const double OMEGA_REF = 0.979;
const double KAPPA_REF = 8.829;
const double TOTAL_LNL_REF = -775.282;

TEST_SUITE("Phase 4: Per-Site Analysis") {

    TEST_CASE("Load and validate data") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== DATA VALIDATION ===\n";
        std::cout << "Alignment: " << aln.sequences.size() << " taxa, "
                  << prefs.num_sites() << " sites\n";
        std::cout << "Tree: " << tree.num_nodes() << " nodes\n";
        std::cout << "Preferences: " << prefs.num_sites() << " sites loaded\n";

        CHECK(aln.sequences.size() == 19);
        CHECK(prefs.num_sites() == 130);
        CHECK(tree.num_nodes() == 36);
    }

    TEST_CASE("Compute single site likelihood") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== SINGLE SITE LIKELIHOOD TEST ===\n";

        // Test parameters from Phase 3C
        ExpCMParameters params;
        params.beta = BETA_REF;
        params.omega = OMEGA_REF;
        params.kappa = KAPPA_REF;
        params.phi = {0.35, 0.20, 0.25, 0.20};

        // Compute likelihood for first 5 sites
        std::cout << "Site     ln L\n";
        std::cout << "-----  ----------\n";
        for (size_t site = 0; site < 5; ++site) {
            Scalar site_lnl = LikelihoodEngine::compute_expcm_site_log_likelihood(
                tree, aln, gcode, prefs, site, params
            );
            std::cout << std::setw(5) << site << "  "
                      << std::setw(10) << std::fixed << std::setprecision(4) << site_lnl << "\n";
            CHECK(!std::isnan(site_lnl));
            CHECK(!std::isinf(site_lnl));
            CHECK(site_lnl < 0);  // Log likelihood should be negative
        }
    }

    TEST_CASE("Compute all sites log likelihood and verify total") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== ALL SITES LIKELIHOOD TEST ===\n";

        ExpCMParameters params;
        params.beta = BETA_REF;
        params.omega = OMEGA_REF;
        params.kappa = KAPPA_REF;
        params.phi = {0.35, 0.20, 0.25, 0.20};

        // Compute all site likelihoods
        auto site_lnls = LikelihoodEngine::compute_expcm_all_sites_log_likelihood(
            tree, aln, gcode, prefs, params
        );

        std::cout << "Number of sites: " << site_lnls.size() << "\n";
        CHECK(site_lnls.size() == 130);

        // Sum per-site likelihoods
        double sum_lnl = 0.0;
        double min_lnl = site_lnls[0];
        double max_lnl = site_lnls[0];

        for (size_t i = 0; i < site_lnls.size(); ++i) {
            sum_lnl += site_lnls[i];
            min_lnl = std::min(min_lnl, site_lnls[i]);
            max_lnl = std::max(max_lnl, site_lnls[i]);

            // Verify all are valid numbers
            CHECK(!std::isnan(site_lnls[i]));
            CHECK(!std::isinf(site_lnls[i]));
            CHECK(site_lnls[i] < 0);
        }

        // Compute full likelihood for comparison
        Scalar full_lnl = LikelihoodEngine::compute_expcm_log_likelihood(
            tree, aln, gcode, prefs, params
        );

        std::cout << "\n--- Summary Statistics ---\n";
        std::cout << "Full likelihood:     " << std::setw(12) << std::fixed << std::setprecision(6) << full_lnl << "\n";
        std::cout << "Sum of per-site:     " << std::setw(12) << std::fixed << std::setprecision(6) << sum_lnl << "\n";
        std::cout << "Difference:          " << std::setw(12) << std::fixed << std::setprecision(8)
                  << std::fabs(full_lnl - sum_lnl) << "\n";
        std::cout << "Min site ln L:       " << std::setw(12) << std::fixed << std::setprecision(6) << min_lnl << "\n";
        std::cout << "Max site ln L:       " << std::setw(12) << std::fixed << std::setprecision(6) << max_lnl << "\n";
        std::cout << "Mean site ln L:      " << std::setw(12) << std::fixed << std::setprecision(6)
                  << sum_lnl / site_lnls.size() << "\n";

        // The sum should match the full likelihood (within floating point error)
        double rel_error = std::fabs(full_lnl - sum_lnl) / std::fabs(full_lnl);
        std::cout << "\nRelative error: " << std::scientific << std::setprecision(2) << rel_error << "\n";

        CHECK(rel_error < 1e-10);  // Should match very closely (floating point precision)

        std::cout << "\n✓ Per-site likelihoods sum correctly to full likelihood\n";
    }

    TEST_CASE("Distribution of per-site likelihoods") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== PER-SITE LIKELIHOOD DISTRIBUTION ===\n";

        ExpCMParameters params;
        params.beta = BETA_REF;
        params.omega = OMEGA_REF;
        params.kappa = KAPPA_REF;
        params.phi = {0.35, 0.20, 0.25, 0.20};

        auto site_lnls = LikelihoodEngine::compute_expcm_all_sites_log_likelihood(
            tree, aln, gcode, prefs, params
        );

        // Compute percentiles
        std::vector<double> sorted_lnls = std::vector<double>(site_lnls.begin(), site_lnls.end());
        std::sort(sorted_lnls.begin(), sorted_lnls.end());

        size_t p10_idx = (10 * sorted_lnls.size()) / 100;
        size_t p25_idx = (25 * sorted_lnls.size()) / 100;
        size_t p50_idx = (50 * sorted_lnls.size()) / 100;
        size_t p75_idx = (75 * sorted_lnls.size()) / 100;
        size_t p90_idx = (90 * sorted_lnls.size()) / 100;

        std::cout << "Percentiles of per-site ln L:\n";
        std::cout << "10th:  " << std::setw(10) << std::fixed << std::setprecision(4) << sorted_lnls[p10_idx] << "\n";
        std::cout << "25th:  " << std::setw(10) << std::fixed << std::setprecision(4) << sorted_lnls[p25_idx] << "\n";
        std::cout << "50th:  " << std::setw(10) << std::fixed << std::setprecision(4) << sorted_lnls[p50_idx] << "\n";
        std::cout << "75th:  " << std::setw(10) << std::fixed << std::setprecision(4) << sorted_lnls[p75_idx] << "\n";
        std::cout << "90th:  " << std::setw(10) << std::fixed << std::setprecision(4) << sorted_lnls[p90_idx] << "\n";

        // Count sites by likelihood ranges
        int count_below_minus6 = 0, count_6_to_5 = 0, count_5_to_4 = 0, count_4_to_0 = 0;
        for (auto lnl : site_lnls) {
            if (lnl < -6) count_below_minus6++;
            else if (lnl < -5) count_6_to_5++;
            else if (lnl < -4) count_5_to_4++;
            else count_4_to_0++;
        }

        std::cout << "\nSites by ln L range:\n";
        std::cout << "< -6:        " << std::setw(4) << count_below_minus6 << " sites\n";
        std::cout << "-6 to -5:    " << std::setw(4) << count_6_to_5 << " sites\n";
        std::cout << "-5 to -4:    " << std::setw(4) << count_5_to_4 << " sites\n";
        std::cout << "-4 to 0:     " << std::setw(4) << count_4_to_0 << " sites\n";

        CHECK(site_lnls.size() == 130);
        CHECK(count_below_minus6 + count_6_to_5 + count_5_to_4 + count_4_to_0 == 130);
    }

    TEST_CASE("Per-site likelihood consistency across runs") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== CONSISTENCY TEST ===\n";

        ExpCMParameters params;
        params.beta = BETA_REF;
        params.omega = OMEGA_REF;
        params.kappa = KAPPA_REF;
        params.phi = {0.35, 0.20, 0.25, 0.20};

        // Compute twice and compare
        auto sites1 = LikelihoodEngine::compute_expcm_all_sites_log_likelihood(
            tree, aln, gcode, prefs, params
        );

        auto sites2 = LikelihoodEngine::compute_expcm_all_sites_log_likelihood(
            tree, aln, gcode, prefs, params
        );

        std::cout << "Checking consistency across runs...\n";
        for (size_t i = 0; i < sites1.size(); ++i) {
            double diff = std::fabs(sites1[i] - sites2[i]);
            CHECK(diff < 1e-14);  // Should be identical (same computation)
        }

        std::cout << "✓ All per-site likelihoods consistent\n";
    }
}
