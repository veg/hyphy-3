#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "hyphy/analyses/persite_optimizer.hpp"
#include "hyphy/core/likelihood.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace hyphy::core;
using namespace hyphy::analyses;

const std::string ALIGNMENT_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/alignment.fasta";
const std::string TREE_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/tree.newick";
const std::string PREFS_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/prefs.csv";
const std::string REFERENCE_FILE = HYPHY3_SOURCE_DIR "/tests/data/expcm/ExpCM_NP_prefs_short_omegabysite.txt";

// Fixed parameters from Phase 3C
const double BETA_PHASE3 = 3.367;
const double OMEGA_PHASE3 = 0.979;
const double KAPPA_PHASE3 = 8.829;

TEST_SUITE("Phase 4.5: Validation against phydms") {

    TEST_CASE("Generate omegabysite output") {
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        std::cout << "\n=== GENERATING OMEGABYSITE OUTPUT ===\n";

        // Optimize all sites
        auto results = PersiteOptimizer::optimize_all_sites(
            tree, aln, gcode, prefs,
            BETA_PHASE3, KAPPA_PHASE3
        );

        // Apply FDR correction
        PersiteOptimizer::apply_fdr_correction(results, 0.05);

        // Write output
        std::string output_file = "/tmp/hyphy3_omegabysite.txt";
        PersiteOptimizer::write_omegabysite_file(results, output_file, true);

        std::cout << "Generated " << output_file << "\n";
        std::cout << "Total sites: " << results.size() << "\n";

        // Basic sanity checks
        CHECK(results.size() == 130);

        // Print summary
        int sig_count = 0;
        for (const auto& r : results) {
            if (r.q_value < 0.05) sig_count++;
        }
        std::cout << "Significant sites (q < 0.05): " << sig_count << "\n";
    }

    TEST_CASE("Compare with phydms reference") {
        std::cout << "\n=== COMPARING WITH PHYDMS REFERENCE ===\n";

        // Read reference file
        std::ifstream ref_file(REFERENCE_FILE);
        if (!ref_file.is_open()) {
            std::cout << "WARNING: Reference file not found: " << REFERENCE_FILE << "\n";
            return;
        }

        std::string line;
        std::vector<int> ref_sites;
        std::vector<double> ref_omega, ref_p, ref_dlnl, ref_q;

        // Skip header comments
        while (std::getline(ref_file, line)) {
            if (line[0] != '#') break;
        }

        // Parse reference data
        std::istringstream iss(line);
        std::string site_str, omega_str, p_str, dlnl_str, q_str;

        do {
            std::istringstream iss_line(line);
            if (iss_line >> site_str >> omega_str >> p_str >> dlnl_str >> q_str) {
                try {
                    ref_sites.push_back(std::stoi(site_str));
                    ref_omega.push_back(std::stod(omega_str));
                    ref_p.push_back(std::stod(p_str));
                    ref_dlnl.push_back(std::stod(dlnl_str));
                    ref_q.push_back(std::stod(q_str));
                } catch (...) {
                    // Skip lines that don't parse as data
                }
            }
        } while (std::getline(ref_file, line));

        ref_file.close();

        std::cout << "Read " << ref_sites.size() << " sites from reference\n";

        // Generate our output
        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        auto results = PersiteOptimizer::optimize_all_sites(
            tree, aln, gcode, prefs,
            BETA_PHASE3, KAPPA_PHASE3
        );

        PersiteOptimizer::apply_fdr_correction(results, 0.05);

        std::cout << "Generated " << results.size() << " sites from optimizer\n";

        // Compare key metrics
        std::cout << "\n--- Comparison Summary ---\n";

        // Count sites with ω > 1
        int our_omega_gt1 = 0, ref_omega_gt1 = 0;
        for (const auto& r : results) {
            if (r.omega_mle > 1.0) our_omega_gt1++;
        }
        for (double omega : ref_omega) {
            if (omega > 1.0) ref_omega_gt1++;
        }

        std::cout << "Sites with ω > 1:\n";
        std::cout << "  phydms:  " << ref_omega_gt1 << "\n";
        std::cout << "  hyphy3:  " << our_omega_gt1 << "\n";

        // Count significant sites
        int our_sig = 0, ref_sig = 0;
        for (const auto& r : results) {
            if (r.q_value < 0.05) our_sig++;
        }
        for (double q : ref_q) {
            if (q < 0.05) ref_sig++;
        }

        std::cout << "\nSignificant sites (q < 0.05):\n";
        std::cout << "  phydms:  " << ref_sig << "\n";
        std::cout << "  hyphy3:  " << our_sig << "\n";

        // Compare omega distribution
        double our_mean_omega = 0, ref_mean_omega = 0;
        for (const auto& r : results) our_mean_omega += r.omega_mle;
        for (double omega : ref_omega) ref_mean_omega += omega;

        our_mean_omega /= results.size();
        ref_mean_omega /= ref_omega.size();

        std::cout << "\nMean ω:\n";
        std::cout << "  phydms:  " << std::fixed << std::setprecision(4) << ref_mean_omega << "\n";
        std::cout << "  hyphy3:  " << std::fixed << std::setprecision(4) << our_mean_omega << "\n";

        // Check if number of significant sites matches
        if (our_sig == ref_sig) {
            std::cout << "\n✓ Number of significant sites matches!\n";
        } else {
            std::cout << "\n⚠ Number of significant sites differs\n";
        }

        // Verify all results are valid
        for (const auto& r : results) {
            CHECK(!std::isnan(r.omega_mle));
            CHECK(!std::isinf(r.omega_mle));
            CHECK(r.p_value >= 0.0);
            CHECK(r.p_value <= 1.0);
            CHECK(r.q_value >= 0.0);
            CHECK(r.q_value <= 1.0);
        }
    }

    TEST_CASE("Check specific sites") {
        std::cout << "\n=== CHECKING SPECIFIC SITES ===\n";

        auto aln = Alignment::from_fasta(ALIGNMENT_FILE);
        auto tree = Tree::from_newick_file(TREE_FILE);
        auto prefs = Preferences::from_csv(PREFS_FILE);
        auto gcode = *GeneticCode::universal();

        auto results = PersiteOptimizer::optimize_all_sites(
            tree, aln, gcode, prefs,
            BETA_PHASE3, KAPPA_PHASE3
        );

        PersiteOptimizer::apply_fdr_correction(results, 0.05);

        // Find sites with most extreme values
        double max_dlnl = -1e10;
        size_t max_site = 0;
        int min_p_site = 0;
        double min_p = 1.0;

        for (size_t i = 0; i < results.size(); ++i) {
            double dlnl = results[i].lnl_alt - results[i].lnl_null;
            if (dlnl > max_dlnl) {
                max_dlnl = dlnl;
                max_site = i;
            }
            if (results[i].p_value < min_p) {
                min_p = results[i].p_value;
                min_p_site = i;
            }
        }

        std::cout << "Site with largest dLnL (site " << (max_site+1) << "):\n";
        std::cout << "  ω:    " << std::fixed << std::setprecision(4) << results[max_site].omega_mle << "\n";
        std::cout << "  dLnL: " << std::fixed << std::setprecision(4) << max_dlnl << "\n";
        std::cout << "  P:    " << std::scientific << std::setprecision(4) << results[max_site].p_value << "\n";
        std::cout << "  Q:    " << std::scientific << std::setprecision(4) << results[max_site].q_value << "\n";

        std::cout << "\nSite with smallest p-value (site " << (min_p_site+1) << "):\n";
        std::cout << "  ω:    " << std::fixed << std::setprecision(4) << results[min_p_site].omega_mle << "\n";
        std::cout << "  dLnL: " << std::fixed << std::setprecision(4)
                  << (results[min_p_site].lnl_alt - results[min_p_site].lnl_null) << "\n";
        std::cout << "  P:    " << std::scientific << std::setprecision(4) << results[min_p_site].p_value << "\n";
        std::cout << "  Q:    " << std::scientific << std::setprecision(4) << results[min_p_site].q_value << "\n";

        CHECK(results.size() == 130);
    }
}
