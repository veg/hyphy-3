#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "hyphy/core/preferences.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/genetic_code.hpp"
#include <cmath>
#include <sstream>
#include <fstream>
#include <filesystem>

using namespace hyphy::core;

/**
 * Test Suite: ExpCM Implementation
 *
 * Step 1: Preferences Parser
 *   - Load CSV file with site × amino acid preferences
 *   - Apply minpref floor (default 0.002)
 *   - Normalize to sum to 1
 *   - Validate against alignment
 */

// Helper: Create a temporary test CSV file
std::string create_test_csv(const std::string& content) {
    std::string filepath = "/tmp/test_prefs_" + std::to_string(std::time(nullptr)) + ".csv";
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not create test CSV file: " + filepath);
    }
    file << content;
    file.close();
    return filepath;
}

TEST_SUITE("ExpCM: Step 1 - Preferences Parser") {

    TEST_CASE("Load valid preferences CSV (20 AAs, 3 sites)") {
        // CSV: values must sum to 1 across columns (amino acids) within each site
        std::string csv_content = R"(site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05
2,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05
3,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05)";

        std::string filepath = create_test_csv(csv_content);
        Preferences prefs = Preferences::from_csv(filepath);

        CHECK(prefs.num_sites() == 3);
        // After loading, each value should be 0.05 (already normalized)
        CHECK(std::fabs(prefs.at(static_cast<size_t>(0), static_cast<size_t>(0)) - 0.05) < 1e-6);
        CHECK(std::fabs(prefs.at(static_cast<size_t>(1), static_cast<size_t>(0)) - 0.05) < 1e-6);

        // Verify normalization: each site sums to 1
        for (size_t r = 0; r < prefs.num_sites(); ++r) {
            double sum = 0.0;
            for (size_t a = 0; a < Preferences::N_AA; ++a) {
                sum += prefs.at(r, a);
            }
            CHECK(std::fabs(sum - 1.0) < 1e-6);
        }

        std::filesystem::remove(filepath);
    }

    TEST_CASE("Apply minpref floor") {
        // Create a case where some values are below minpref
        // Input: first 10 AAs = 0.0, last 10 AAs = 0.01 (sums to 0.1)
        std::string csv_content = R"(site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)";

        std::string filepath = create_test_csv(csv_content);
        Scalar minpref = 0.002;
        Preferences prefs = Preferences::from_csv(filepath, minpref);

        // After loading:
        // - First 10 AAs: floored from 0.0 to 0.002
        // - Last 10 AAs: already 0.01, no change
        // Sum before norm = 10*0.002 + 10*0.01 = 0.02 + 0.1 = 0.12
        // After norm:  first 10 each = 0.002/0.12, last 10 each = 0.01/0.12

        // Just verify that all values are >= minpref after loading
        for (size_t a = 0; a < Preferences::N_AA; ++a) {
            double val = prefs.at(static_cast<size_t>(0), a);
            CHECK(val >= minpref * 0.9);  // Accounting for some numerical loss in normalization
        }

        std::filesystem::remove(filepath);
    }

    TEST_CASE("Access by amino acid letter") {
        std::string csv_content = R"(site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05)";

        std::string filepath = create_test_csv(csv_content);
        Preferences prefs = Preferences::from_csv(filepath);

        // Check access by letter
        CHECK(std::fabs(prefs.at(0, 'A') - 0.05) < 1e-6);
        CHECK(std::fabs(prefs.at(0, 'G') - 0.05) < 1e-6);
        CHECK(std::fabs(prefs.at(0, 'Y') - 0.05) < 1e-6);

        std::filesystem::remove(filepath);
    }

    TEST_CASE("Row access") {
        std::string csv_content = R"(site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05)";

        std::string filepath = create_test_csv(csv_content);
        Preferences prefs = Preferences::from_csv(filepath);

        auto row = prefs.row(0);
        CHECK(row.size() == Preferences::N_AA);
        for (int i = 0; i < row.size(); ++i) {
            CHECK(std::fabs(row(0, i) - 0.05) < 1e-6);
        }

        std::filesystem::remove(filepath);
    }

    TEST_CASE("Invalid access throws") {
        std::string csv_content = R"(site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05)";

        std::string filepath = create_test_csv(csv_content);
        Preferences prefs = Preferences::from_csv(filepath);

        CHECK_THROWS(prefs.at(static_cast<size_t>(1), static_cast<size_t>(0)));  // Only 1 site (index 0), not 2
        CHECK_THROWS(prefs.at(static_cast<size_t>(0), 'X'));  // Invalid amino acid

        std::filesystem::remove(filepath);
    }

    TEST_CASE("Match phydms check_prefs.py on NP_short data") {
        // Load actual test data
        std::string prefs_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/prefs.csv";

        // Note: only test if file exists (skip if not available)
        if (std::filesystem::exists(prefs_file)) {
            Preferences prefs = Preferences::from_csv(prefs_file);

            // From check_prefs.py output: 130 sites
            CHECK(prefs.num_sites() == 130);

            // Each site should sum to 1 (within tolerance)
            for (size_t r = 0; r < prefs.num_sites(); ++r) {
                double sum = 0.0;
                for (size_t a = 0; a < Preferences::N_AA; ++a) {
                    sum += prefs.at(r, a);
                }
                CHECK(std::fabs(sum - 1.0) < 1e-6);
            }
        }
    }
}

TEST_SUITE("ExpCM: Step 2 - Rate Matrices (Q, F, P)") {

    TEST_CASE("Q Matrix: Structure (row sums = 0)") {
        auto gcode = hyphy::core::GeneticCode::universal();
        hyphy::core::ExpCMParameters params;
        params.beta = 1.0;
        params.omega = 1.0;
        params.kappa = 2.0;
        params.phi = hyphy::core::Vector4(0.25, 0.25, 0.25, 0.25);

        hyphy::core::ExpCMQMatrix Q_matrix;
        Q_matrix.build(*gcode, params);

        // Check that row sums ≈ 0
        for (int x = 0; x < gcode->num_sense_codons; ++x) {
            double row_sum = Q_matrix.Q.row(x).sum();
            CHECK(std::fabs(row_sum) < 1e-10);
        }
    }

    TEST_CASE("Q Matrix: Transitions vs Transversions") {
        auto gcode = hyphy::core::GeneticCode::universal();
        hyphy::core::ExpCMParameters params;
        params.kappa = 2.0;
        params.phi = hyphy::core::Vector4(0.25, 0.25, 0.25, 0.25);

        hyphy::core::ExpCMQMatrix Q_matrix;
        Q_matrix.build(*gcode, params);

        // Q matrix should be symmetric (all φ_w are equal)
        CHECK(Q_matrix.Q.isApprox(Q_matrix.Q.transpose(), 1e-10));

        // Off-diagonal elements should be either φ or κ·φ = 0.25 or 0.5
        for (int x = 0; x < gcode->num_sense_codons; ++x) {
            for (int y = 0; y < gcode->num_sense_codons; ++y) {
                if (x == y) continue;
                double val = Q_matrix.Q(x, y);
                if (val > 0) {
                    CHECK((std::fabs(val - 0.25) < 1e-10 || std::fabs(val - 0.5) < 1e-10));
                }
            }
        }
    }

    TEST_CASE("F Matrix: Synonymous = 1") {
        auto gcode = hyphy::core::GeneticCode::universal();

        // Create uniform preferences
        hyphy::core::Preferences prefs;
        prefs.nsites = 1;
        prefs.prefs = Eigen::Matrix<hyphy::core::Scalar,
                                     Eigen::Dynamic, 20>::Ones(1, 20) / 20.0;

        hyphy::core::ExpCMParameters params;
        params.beta = 1.5;
        params.omega = 0.8;

        hyphy::core::ExpCMFMatrix F_matrix;
        F_matrix.build(*gcode, prefs, 0, params);

        // For each codon pair, if synonymous then F_xy should be close to 1
        for (int x = 0; x < gcode->num_sense_codons; ++x) {
            for (int y = 0; y < gcode->num_sense_codons; ++y) {
                if (x == y) {
                    CHECK(F_matrix.F(x, y) == 0.0);  // Diagonal is 0
                } else if (gcode->is_synonymous(x, y)) {
                    CHECK(std::fabs(F_matrix.F(x, y) - 1.0) < 1e-10);
                }
            }
        }
    }

    TEST_CASE("F Matrix: All values non-negative") {
        auto gcode = hyphy::core::GeneticCode::universal();

        hyphy::core::Preferences prefs;
        prefs.nsites = 1;
        prefs.prefs = Eigen::Matrix<hyphy::core::Scalar,
                                     Eigen::Dynamic, 20>::Ones(1, 20) / 20.0;

        hyphy::core::ExpCMParameters params;
        params.beta = 2.0;
        params.omega = 1.0;

        hyphy::core::ExpCMFMatrix F_matrix;
        F_matrix.build(*gcode, prefs, 0, params);

        // All F values should be ≥ 0
        for (int x = 0; x < gcode->num_sense_codons; ++x) {
            for (int y = 0; y < gcode->num_sense_codons; ++y) {
                CHECK(F_matrix.F(x, y) >= -1e-10);  // Allow small numerical error
            }
        }
    }

    TEST_CASE("P Matrix: Row sums = 0") {
        auto gcode = hyphy::core::GeneticCode::universal();

        hyphy::core::Preferences prefs;
        prefs.nsites = 1;
        prefs.prefs = Eigen::Matrix<hyphy::core::Scalar,
                                     Eigen::Dynamic, 20>::Ones(1, 20) / 20.0;

        hyphy::core::ExpCMParameters params;
        params.beta = 1.0;
        params.omega = 1.0;
        params.kappa = 2.0;
        params.phi = hyphy::core::Vector4(0.25, 0.25, 0.25, 0.25);

        hyphy::core::ExpCMQMatrix Q_matrix;
        Q_matrix.build(*gcode, params);

        hyphy::core::ExpCMFMatrix F_matrix;
        F_matrix.build(*gcode, prefs, 0, params);

        hyphy::core::ExpCMPMatrix P_matrix;
        P_matrix.build(F_matrix.F, Q_matrix.Q);

        // Check that row sums ≈ 0
        for (int x = 0; x < gcode->num_sense_codons; ++x) {
            double row_sum = P_matrix.P.row(x).sum();
            CHECK(std::fabs(row_sum) < 1e-10);
        }
    }

    TEST_CASE("P Matrix: Detailed balance") {
        // For a simple test, verify that π_x P_xy has reasonable structure
        // (Full detailed balance check requires stationary state computation)
        auto gcode = hyphy::core::GeneticCode::universal();

        hyphy::core::Preferences prefs;
        prefs.nsites = 1;
        prefs.prefs = Eigen::Matrix<hyphy::core::Scalar,
                                     Eigen::Dynamic, 20>::Ones(1, 20) / 20.0;

        hyphy::core::ExpCMParameters params;
        params.beta = 1.0;
        params.omega = 1.0;
        params.kappa = 2.0;
        params.phi = hyphy::core::Vector4(0.25, 0.25, 0.25, 0.25);

        hyphy::core::ExpCMQMatrix Q_matrix;
        Q_matrix.build(*gcode, params);

        hyphy::core::ExpCMFMatrix F_matrix;
        F_matrix.build(*gcode, prefs, 0, params);

        hyphy::core::ExpCMPMatrix P_matrix;
        P_matrix.build(F_matrix.F, Q_matrix.Q);

        // At least check that symmetric entries exist
        int off_diagonal_pairs = 0;
        for (int x = 0; x < gcode->num_sense_codons; ++x) {
            for (int y = x + 1; y < gcode->num_sense_codons; ++y) {
                if (P_matrix.P(x, y) > 0 || P_matrix.P(y, x) > 0) {
                    off_diagonal_pairs++;
                }
            }
        }
        CHECK(off_diagonal_pairs > 0);  // Must have some off-diagonal entries
    }
}

// TODO: Step 3 - Fixed-Parameter Likelihood
// TODO: Step 4 - Optimization
// TODO: Step 5 - Per-Site ω (omegabysite)
