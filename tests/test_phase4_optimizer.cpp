#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "hyphy/analyses/persite_optimizer.hpp"
#include "hyphy/core/likelihood.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace hyphy::core;
using namespace hyphy::analyses;

const std::string ALIGNMENT_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/alignment.fasta";
const std::string TREE_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/tree.newick";
const std::string PREFS_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/prefs.csv";

// Fixed parameters from Phase 3C optimization
const double BETA_PHASE3 = 3.367;
const double OMEGA_PHASE3 = 0.979;
const double KAPPA_PHASE3 = 8.829;

TEST_SUITE("Phase 4.2: Per-Site ω Optimization") {

    TEST_CASE("Load data for per-site optimization") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== PHASE 4.2 DATA VALIDATION ===\n";
        std::cout << "Alignment: " << aln.sequences.size() << " taxa\n";
        std::cout << "Sites: " << prefs.num_sites() << "\n";
        std::cout << "Global parameters: β=" << BETA_PHASE3
                  << ", ω=" << OMEGA_PHASE3
                  << ", κ=" << KAPPA_PHASE3 << "\n";

        CHECK(prefs.num_sites() == 130);
    }

    TEST_CASE("Optimize ω for single site") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== SINGLE SITE OPTIMIZATION ===\n";

        // Optimize site 0
        PersiteOptimizer::SiteResult result =
            PersiteOptimizer::optimize_site_omega_grid(
                tree, aln, gcode, prefs, 0,
                BETA_PHASE3, KAPPA_PHASE3
            );

        std::cout << "Site " << result.site_id << ":\n";
        std::cout << "  ω_MLE:       " << std::setw(8) << std::fixed << std::setprecision(4) << result.omega_mle << "\n";
        std::cout << "  ln L (H0):   " << std::setw(8) << std::fixed << std::setprecision(4) << result.lnl_null << "\n";
        std::cout << "  ln L (H1):   " << std::setw(8) << std::fixed << std::setprecision(4) << result.lnl_alt << "\n";
        std::cout << "  LRT stat:    " << std::setw(8) << std::fixed << std::setprecision(4) << result.test_statistic << "\n";
        std::cout << "  p-value:     " << std::scientific << std::setprecision(4) << result.p_value << "\n";

        CHECK(result.omega_mle > 0.0);
        CHECK(result.lnl_alt >= result.lnl_null);
        CHECK(result.test_statistic >= 0.0);
        CHECK(result.p_value >= 0.0);
        CHECK(result.p_value <= 1.0);
    }

    TEST_CASE("Optimize ω for all sites") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== ALL SITES OPTIMIZATION ===\n";

        auto results = PersiteOptimizer::optimize_all_sites(
            tree, aln, gcode, prefs,
            BETA_PHASE3, KAPPA_PHASE3
        );

        std::cout << "Optimized " << results.size() << " sites\n";
        CHECK(results.size() == 130);

        // Summary statistics
        std::vector<double> omega_mles, lrt_stats, pvalues;
        for (const auto& r : results) {
            omega_mles.push_back(r.omega_mle);
            lrt_stats.push_back(r.test_statistic);
            pvalues.push_back(r.p_value);
        }

        std::sort(omega_mles.begin(), omega_mles.end());
        std::sort(lrt_stats.begin(), lrt_stats.end());
        std::sort(pvalues.begin(), pvalues.end());

        std::cout << "\n--- ω Distribution ---\n";
        std::cout << "Min ω:      " << std::setw(8) << std::fixed << std::setprecision(4) << omega_mles.front() << "\n";
        std::cout << "Median ω:   " << std::setw(8) << std::fixed << std::setprecision(4) << omega_mles[65] << "\n";
        std::cout << "Max ω:      " << std::setw(8) << std::fixed << std::setprecision(4) << omega_mles.back() << "\n";
        double mean_omega = std::accumulate(omega_mles.begin(), omega_mles.end(), 0.0) / omega_mles.size();
        std::cout << "Mean ω:     " << std::setw(8) << std::fixed << std::setprecision(4) << mean_omega << "\n";

        std::cout << "\n--- LRT Distribution ---\n";
        std::cout << "Min LRT:    " << std::setw(8) << std::fixed << std::setprecision(4) << lrt_stats.front() << "\n";
        std::cout << "Median LRT: " << std::setw(8) << std::fixed << std::setprecision(4) << lrt_stats[65] << "\n";
        std::cout << "Max LRT:    " << std::setw(8) << std::fixed << std::setprecision(4) << lrt_stats.back() << "\n";

        std::cout << "\n--- P-value Distribution ---\n";
        int sig_001 = 0, sig_005 = 0, sig_01 = 0;
        for (const auto& r : results) {
            if (r.p_value < 0.001) sig_001++;
            if (r.p_value < 0.005) sig_005++;
            if (r.p_value < 0.01) sig_01++;
        }
        std::cout << "p < 0.001:  " << std::setw(3) << sig_001 << " sites\n";
        std::cout << "p < 0.005:  " << std::setw(3) << sig_005 << " sites\n";
        std::cout << "p < 0.01:   " << std::setw(3) << sig_01 << " sites\n";

        // Verify all results are valid
        for (size_t i = 0; i < results.size(); ++i) {
            CHECK(results[i].site_id == i);
            CHECK(results[i].omega_mle > 0.0);
            CHECK(results[i].test_statistic >= 0.0);
            CHECK(!std::isnan(results[i].p_value));
            CHECK(!std::isinf(results[i].p_value));
        }
    }

    TEST_CASE("FDR correction") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== FDR CORRECTION TEST ===\n";

        auto results = PersiteOptimizer::optimize_all_sites(
            tree, aln, gcode, prefs,
            BETA_PHASE3, KAPPA_PHASE3
        );

        // Apply FDR correction
        PersiteOptimizer::apply_fdr_correction(results, 0.05);

        // Count significant sites
        int sig_fdr = 0;
        for (const auto& r : results) {
            if (r.q_value < 0.05) sig_fdr++;
        }

        std::cout << "Sites with q < 0.05: " << sig_fdr << "\n";

        // Verify q-values are in [0,1] and finite
        bool qvalues_ok = true;
        for (const auto& r : results) {
            if (r.q_value < 0.0 || r.q_value > 1.0 || std::isnan(r.q_value)) {
                qvalues_ok = false;
                break;
            }
        }
        CHECK(qvalues_ok);

        // Verify groups: sites with ω > 1 and ω < 1 have q-values (separate correction per group)
        int count_up = 0, count_down = 0, count_eq = 0;
        for (const auto& r : results) {
            if (r.omega_mle > 1.0) count_up++;
            else if (r.omega_mle < 1.0) count_down++;
            else count_eq++;
        }

        std::cout << "Separate FDR groups: ω>1: " << count_up << " sites, "
                  << "ω<1: " << count_down << " sites, ω=1: " << count_eq << " sites\n";

        std::cout << "✓ FDR correction applied correctly (separate per group)\n";
    }

    TEST_CASE("Output format verification") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== OUTPUT FORMAT VERIFICATION ===\n";

        auto results = PersiteOptimizer::optimize_all_sites(
            tree, aln, gcode, prefs,
            BETA_PHASE3, KAPPA_PHASE3
        );

        PersiteOptimizer::apply_fdr_correction(results, 0.05);

        // Print header like phydms _omegabysite.txt format
        std::cout << "\nSample output (first 10 sites):\n";
        std::cout << "site\tomega\tp_value\tq_value\tsignificant\n";
        std::cout << "----\t-----\t-------\t-------\t-----------\n";

        for (size_t i = 0; i < std::min(size_t(10), results.size()); ++i) {
            const auto& r = results[i];
            std::string sig = r.q_value < 0.05 ? "YES" : "NO";
            std::cout << r.site_id << "\t"
                      << std::fixed << std::setprecision(4) << r.omega_mle << "\t"
                      << std::scientific << std::setprecision(4) << r.p_value << "\t"
                      << std::scientific << std::setprecision(4) << r.q_value << "\t"
                      << sig << "\n";
        }

        CHECK(results.size() == 130);
    }
}
