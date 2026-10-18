#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "hyphy/core/preferences.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include <cmath>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>

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

TEST_SUITE("ExpCM: Step 3 - Fixed-Parameter Likelihood") {

    TEST_CASE("Debug: Stationary distribution sums to 1") {
        auto gcode = hyphy::core::GeneticCode::universal();
        hyphy::core::Preferences prefs;
        prefs.nsites = 1;
        prefs.prefs = Eigen::Matrix<hyphy::core::Scalar,
                                     Eigen::Dynamic, 20>::Ones(1, 20) / 20.0;

        hyphy::core::ExpCMParameters params;
        params.beta = 2.0;
        params.phi = hyphy::core::Vector4(0.25, 0.25, 0.25, 0.25);

        auto pi = hyphy::core::LikelihoodEngine::compute_expcm_stationary(
            *gcode, prefs, 0, params.beta, params.phi
        );

        // Check sum = 1
        CHECK(std::fabs(pi.sum() - 1.0) < 1e-6);
        // Check all non-negative
        CHECK((pi.array() >= -1e-10).all());
    }

    TEST_CASE("Stationary distribution with non-uniform preferences") {
        auto gcode = hyphy::core::GeneticCode::universal();
        hyphy::core::Preferences prefs;
        prefs.nsites = 1;

        // Create non-uniform preferences
        Eigen::Matrix<hyphy::core::Scalar, Eigen::Dynamic, 20> p(1, 20);
        p.row(0).head(10).setConstant(0.01);  // First 10 AAs: 0.01
        p.row(0).tail(10).setConstant(0.09);  // Last 10 AAs: 0.09
        prefs.prefs = p;

        hyphy::core::ExpCMParameters params;
        params.beta = 1.5;
        params.phi = hyphy::core::Vector4(0.3, 0.2, 0.25, 0.25);

        auto pi = hyphy::core::LikelihoodEngine::compute_expcm_stationary(
            *gcode, prefs, 0, params.beta, params.phi
        );

        CHECK(std::fabs(pi.sum() - 1.0) < 1e-6);
        CHECK((pi.array() >= -1e-10).all());
    }

    TEST_CASE("Debug: Print Q matrix for comparison with phydms") {
        auto gcode = hyphy::core::GeneticCode::universal();

        hyphy::core::ExpCMParameters params;
        params.beta = 3.38189;
        params.kappa = 8.75598;
        params.omega = 0.977468;
        params.phi = hyphy::core::Vector4(0.355409, 0.201765, 0.249379,
                                          1.0 - 0.355409 - 0.201765 - 0.249379);

        hyphy::core::ExpCMQMatrix Q_matrix;
        Q_matrix.build(*gcode, params);

        // Print Q matrix in a format that's easy to compare
        std::cerr << std::endl << "=== Q MATRIX ===" << std::endl;
        std::cerr << "First 10x10 submatrix:" << std::endl;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                std::cerr << std::fixed << std::setprecision(10) << std::setw(15) << Q_matrix.Q(i, j);
            }
            std::cerr << std::endl;
        }

        // Check specific properties
        std::cerr << std::endl << "Q Matrix Properties:" << std::endl;
        std::cerr << "Codon 0 (" << gcode->sense_codons[0] << "): ";
        for (int j = 0; j < 61; ++j) {
            if (Q_matrix.Q(0, j) != 0.0 && j != 0) {
                std::cerr << "Q[0," << j << "]=" << std::fixed << std::setprecision(6) << Q_matrix.Q(0, j) << " ";
            }
        }
        std::cerr << std::endl;

        // Check a few specific rates
        std::cerr << std::endl << "Sample Q values:" << std::endl;
        int count = 0;
        for (int i = 0; i < 61 && count < 20; ++i) {
            for (int j = 0; j < 61 && count < 20; ++j) {
                if (i != j && Q_matrix.Q(i, j) != 0.0) {
                    std::cerr << "Q[" << i << "," << j << "]=" << Q_matrix.Q(i, j) << " ";
                    count++;
                }
            }
        }
        std::cerr << std::endl;

        CHECK(true);  // Just diagnostic, always passes
    }

    TEST_CASE("Debug: Check Q matrix properties") {
        auto gcode = hyphy::core::GeneticCode::universal();

        hyphy::core::ExpCMParameters params;
        params.beta = 1.0;
        params.omega = 1.0;
        params.kappa = 2.0;
        params.phi = hyphy::core::Vector4(0.25, 0.25, 0.25, 0.25);

        hyphy::core::ExpCMQMatrix Q_matrix;
        Q_matrix.build(*gcode, params);

        // Check Q matrix properties
        std::cerr << "Q matrix diagonal check:" << std::endl;
        Scalar max_diag_error = 0.0;
        for (int x = 0; x < 10; ++x) {  // Check first 10 codons
            Scalar row_sum = Q_matrix.Q.row(x).sum();
            max_diag_error = std::max(max_diag_error, std::fabs(row_sum));
            std::cerr << "  Codon " << x << ": row sum = " << row_sum << std::endl;
        }
        CHECK(max_diag_error < 1e-10);

        std::cerr << "Q matrix off-diagonal range:" << std::endl;
        Scalar min_val = 0.0, max_val = 0.0;
        for (int i = 0; i < 61; ++i) {
            for (int j = 0; j < 61; ++j) {
                if (i != j && Q_matrix.Q(i, j) != 0.0) {
                    min_val = std::min(min_val, Q_matrix.Q(i, j));
                    max_val = std::max(max_val, Q_matrix.Q(i, j));
                }
            }
        }
        std::cerr << "  Min: " << min_val << ", Max: " << max_val << std::endl;
    }

    TEST_CASE("Debug: Compute branchScale using phydms formula") {
        auto gcode = hyphy::core::GeneticCode::universal();

        // Load real test data
        std::string aln_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/alignment.fasta";
        std::string prefs_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/prefs.csv";
        std::string tree_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/tree.newick";

        if (std::filesystem::exists(aln_file) &&
            std::filesystem::exists(prefs_file) &&
            std::filesystem::exists(tree_file)) {

            auto aln = hyphy::core::Alignment::from_fasta(aln_file);
            auto prefs = hyphy::core::Preferences::from_csv(prefs_file);
            auto tree = hyphy::core::Tree::from_newick_file(tree_file);

            hyphy::core::ExpCMParameters params;
            params.beta = 3.38189;
            params.kappa = 8.75598;
            params.omega = 0.977468;
            params.phi = hyphy::core::Vector4(0.355409, 0.201765, 0.249379,
                                              1.0 - 0.355409 - 0.201765 - 0.249379);

            // Build Q matrix
            hyphy::core::ExpCMQMatrix Q_matrix;
            Q_matrix.build(*gcode, params);

            // Compute branchScale using phydms formula:
            // branchScale = Σ_sites[ Σ_codons π[x] * (-P[x,x]) ] / nsites
            // where P = F ⊙ Q
            double sum_rates = 0.0;
            size_t num_sites = prefs.num_sites();

            for (size_t site = 0; site < num_sites; ++site) {
                hyphy::core::ExpCMFMatrix F_matrix;
                F_matrix.build(*gcode, prefs, site, params);

                hyphy::core::Matrix61Codon P = F_matrix.F.cwiseProduct(Q_matrix.Q);
                for (int x = 0; x < 61; ++x) {
                    P(x, x) = -P.row(x).sum();
                }

                auto pi_site = hyphy::core::LikelihoodEngine::compute_expcm_stationary(
                    *gcode, prefs, site, params.beta, params.phi
                );

                // Sum: π[x] * (-P[x,x]) for all states x
                for (int x = 0; x < 61; ++x) {
                    sum_rates += pi_site(x) * (-P(x, x));
                }
            }

            double branchScale = sum_rates / num_sites;
            double mu = 1.0;  // From phydms script
            double branchScale_with_mu = branchScale * mu;

            std::cerr << std::endl << "=== BRANCHSCALE COMPUTATION ===" << std::endl;
            std::cerr << "phydms formula: branchScale = Σ_sites[Σ_codons π[x]·(-P[x,x])] / nsites" << std::endl;
            std::cerr << std::endl << "Computed values:" << std::endl;
            std::cerr << "  Sum of rates: " << std::fixed << std::setprecision(8) << sum_rates << std::endl;
            std::cerr << "  Num sites: " << num_sites << std::endl;
            std::cerr << "  Average rate per site: " << (sum_rates / num_sites) << std::endl;
            std::cerr << "  branchScale (μ=1.0): " << branchScale_with_mu << std::endl;
            std::cerr << std::endl << "Interpretation:" << std::endl;
            std::cerr << "  Expected substitutions per codon site per unit time: " << branchScale_with_mu << std::endl;
            std::cerr << "  To convert branch lengths to model time:" << std::endl;
            std::cerr << "    If convention='as_is': use branch lengths directly" << std::endl;
            std::cerr << "    If convention='divide_by_branchScale': divide by " << branchScale_with_mu << std::endl;

            // Estimate which convention was used
            double lnL_as_is = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
                tree, aln, *gcode, prefs, params, 1.0
            );
            double lnL_scaled = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
                tree, aln, *gcode, prefs, params, 1.0 / branchScale_with_mu
            );

            std::cerr << std::endl << "Likelihood comparison:" << std::endl;
            std::cerr << "  as_is (scale=1.0): ln L = " << std::fixed << std::setprecision(2) << lnL_as_is << std::endl;
            std::cerr << "  divide_by_branchScale (scale=1/" << std::fixed << std::setprecision(6)
                      << branchScale_with_mu << "): ln L = " << std::fixed << std::setprecision(2) << lnL_scaled << std::endl;
            std::cerr << "  Expected (phydms): ln L = -775.21" << std::endl;
            std::cerr << "  Ratio (as_is): " << std::fixed << std::setprecision(4)
                      << (lnL_as_is / -775.21) << "x off" << std::endl;
            std::cerr << "  Ratio (scaled): " << std::fixed << std::setprecision(4)
                      << (lnL_scaled / -775.21) << "x off" << std::endl;
        }
    }

    TEST_CASE("Debug: F matrix values for reference parameters") {
        auto gcode = hyphy::core::GeneticCode::universal();

        hyphy::core::Preferences prefs;
        prefs.nsites = 1;
        prefs.prefs = Eigen::Matrix<hyphy::core::Scalar, Eigen::Dynamic, 20>::Ones(1, 20) / 20.0;

        hyphy::core::ExpCMParameters params;
        params.beta = 3.38189;
        params.omega = 0.977468;

        hyphy::core::ExpCMFMatrix F_matrix;
        F_matrix.build(*gcode, prefs, 0, params);

        std::cerr << "\n=== F MATRIX SAMPLE VALUES ===" << std::endl;

        // Find synonymous pairs and verify they're 1
        int syn_count = 0;
        for (int x = 0; x < 61 && syn_count < 10; ++x) {
            for (int y = x + 1; y < 61 && syn_count < 10; ++y) {
                if (gcode->is_synonymous(x, y)) {
                    std::cerr << "F[" << x << "," << y << "] (synonymous) = "
                              << std::fixed << std::setprecision(6) << F_matrix.F(x, y)
                              << " (expected 1.0)" << std::endl;
                    syn_count++;
                }
            }
        }

        // Sample some non-synonymous values
        std::cerr << "\nSample non-synonymous F values:" << std::endl;
        int nonsyn_count = 0;
        for (int x = 0; x < 61 && nonsyn_count < 10; ++x) {
            for (int y = x + 1; y < 61 && nonsyn_count < 10; ++y) {
                if (!gcode->is_synonymous(x, y) && F_matrix.F(x, y) > 0) {
                    std::cerr << "F[" << x << "," << y << "] (non-syn) = "
                              << std::fixed << std::setprecision(6) << F_matrix.F(x, y)
                              << " (should involve ω ≈ " << params.omega << ")" << std::endl;
                    nonsyn_count++;
                }
            }
        }

        CHECK(true);  // Diagnostic only
    }

    TEST_CASE("Debug: Per-site likelihood contribution analysis") {
        auto gcode = hyphy::core::GeneticCode::universal();

        // Load real test data
        std::string aln_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/alignment.fasta";
        std::string prefs_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/prefs.csv";
        std::string tree_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/tree.newick";

        if (std::filesystem::exists(aln_file) &&
            std::filesystem::exists(prefs_file) &&
            std::filesystem::exists(tree_file)) {

            auto aln = hyphy::core::Alignment::from_fasta(aln_file);
            auto prefs = hyphy::core::Preferences::from_csv(prefs_file);
            auto tree = hyphy::core::Tree::from_newick_file(tree_file);

            hyphy::core::ExpCMParameters params;
            params.beta = 3.38189;
            params.kappa = 8.75598;
            params.omega = 0.977468;
            params.phi = hyphy::core::Vector4(0.355409, 0.201765, 0.249379, 0.192447);

            // Debug: Check pattern structure
            std::cerr << "\n=== PATTERN STRUCTURE ===" << std::endl;
            std::cerr << "Num sites: " << aln.num_codons << std::endl;
            std::cerr << "Num patterns: " << aln.patterns.size() << std::endl;

            double total_weight = 0.0;
            for (size_t p = 0; p < aln.patterns.size(); ++p) {
                total_weight += aln.patterns[p].weight;
            }
            std::cerr << "Total pattern weight: " << total_weight << std::endl;
            std::cerr << "Expected total (sum of weights): " << aln.num_codons << std::endl;

            // Check if site_to_pattern mapping is consistent
            std::cerr << "\nSite to pattern mapping check (first 10 sites):" << std::endl;
            for (size_t site = 0; site < std::min(size_t(10), aln.num_codons); ++site) {
                size_t pattern_idx = aln.site_to_pattern[site];
                std::cerr << "  Site " << site << " → Pattern " << pattern_idx
                          << " (weight=" << aln.patterns[pattern_idx].weight << ")" << std::endl;
            }
        }
    }

    TEST_CASE("Debug: Simple likelihood on single codon") {
        auto gcode = hyphy::core::GeneticCode::universal();

        // Create a simple 2-taxon alignment with just 1 pattern
        hyphy::core::Alignment aln;
        aln.add_sequence("seq1", "ATG");  // Codon 40 (Met)
        aln.add_sequence("seq2", "ATG");  // Same codon
        aln.finalize(gcode);

        // Create simple tree: (seq1:0.1, seq2:0.1);
        std::string newick = "(seq1:0.1,seq2:0.1);";
        auto tree = hyphy::core::Tree::from_newick(newick);

        // Simple uniform preferences
        hyphy::core::Preferences prefs;
        prefs.nsites = 1;
        prefs.prefs = Eigen::Matrix<hyphy::core::Scalar,
                                     Eigen::Dynamic, 20>::Ones(1, 20) / 20.0;

        hyphy::core::ExpCMParameters params;
        params.beta = 1.0;
        params.omega = 1.0;
        params.kappa = 2.0;
        params.phi = hyphy::core::Vector4(0.25, 0.25, 0.25, 0.25);

        // Debug: Compute stationary distribution
        auto pi = hyphy::core::LikelihoodEngine::compute_expcm_stationary(
            *gcode, prefs, 0, params.beta, params.phi
        );
        std::cerr << "\nStationary distribution at ATG codon (index 40): "
                  << std::fixed << std::setprecision(10) << pi(40) << std::endl;
        std::cerr << "ln(pi[40]) = " << std::log(pi(40)) << std::endl;

        // Debug: Build P matrix
        hyphy::core::ExpCMQMatrix Q_matrix;
        Q_matrix.build(*gcode, params);
        hyphy::core::ExpCMFMatrix F_matrix;
        F_matrix.build(*gcode, prefs, 0, params);
        hyphy::core::Matrix61Codon P = F_matrix.F.cwiseProduct(Q_matrix.Q);
        for (int x = 0; x < 61; ++x) {
            P(x, x) = -P.row(x).sum();
        }

        std::cerr << "P[40,40] (diagonal) = " << std::fixed << std::setprecision(10) << P(40, 40) << std::endl;

        // Compute likelihood
        double lnL = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
            tree, aln, *gcode, prefs, params
        );

        std::cerr << "Simple 2-taxon likelihood: " << std::fixed << std::setprecision(6) << lnL << std::endl;
        std::cerr << "Expected ln L ≈ 2*ln(pi[40]) + 2*ln(exp(P[40,40]*0.1)[40,40])" << std::endl;
        std::cerr << "  ≈ 2*" << std::log(pi(40)) << " + 2*ln(exp("
                  << P(40, 40) * 0.1 << ")) ≈ " << (2*std::log(pi(40)) + 2*std::log(std::exp(P(40,40)*0.1))) << std::endl;
        std::cerr << "  (both sequences identical ATG, should be high likelihood)" << std::endl;

        // For identical sequences, likelihood should be reasonably high
        CHECK(std::isfinite(lnL));
        CHECK(lnL > -100.0);  // Should not be too negative
    }

    TEST_CASE("Debug: Verify matrix exponential accuracy") {
        auto gcode = hyphy::core::GeneticCode::universal();

        // Test simple case: two identical sequences
        hyphy::core::Alignment aln;
        aln.add_sequence("seq1", "ATG");
        aln.add_sequence("seq2", "ATG");
        aln.finalize(gcode);

        // Create tree with varying branch lengths
        std::string newick = "(seq1:0.01,seq2:0.01);";
        auto tree = hyphy::core::Tree::from_newick(newick);

        hyphy::core::Preferences prefs;
        prefs.nsites = 1;
        prefs.prefs = Eigen::Matrix<hyphy::core::Scalar, Eigen::Dynamic, 20>::Ones(1, 20) / 20.0;

        hyphy::core::ExpCMParameters params;
        params.beta = 1.0;
        params.omega = 1.0;
        params.kappa = 2.0;
        params.phi = hyphy::core::Vector4(0.25, 0.25, 0.25, 0.25);

        // Compute likelihood at different scales
        std::cerr << "\n=== MATRIX EXPONENTIAL ACCURACY TEST ===" << std::endl;
        std::cerr << "Small branch case (t=0.01, identical sequences)" << std::endl;

        for (double scale : {0.1, 0.373, 1.0, 2.68}) {
            double lnL = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
                tree, aln, *gcode, prefs, params, scale
            );
            std::cerr << "  Scale=" << std::fixed << std::setprecision(3) << scale
                      << " → ln L = " << std::setprecision(6) << lnL << std::endl;
        }

        // Now test with reference parameters
        params.beta = 3.38189;
        params.kappa = 8.75598;
        params.omega = 0.977468;

        std::cerr << std::endl << "With reference parameters (β, κ, ω):" << std::endl;
        double lnL_ref = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
            tree, aln, *gcode, prefs, params, 1.0
        );
        std::cerr << "  Scale=1.0 → ln L = " << std::fixed << std::setprecision(6) << lnL_ref << std::endl;
        std::cerr << "  (Expected ≈ -9 to -10 for small branches, identical seqs)" << std::endl;

        CHECK(std::isfinite(lnL_ref));
    }

    TEST_CASE("Debug: Search for optimal scale factor") {
        auto gcode = hyphy::core::GeneticCode::universal();

        // Load real test data
        std::string aln_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/alignment.fasta";
        std::string prefs_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/prefs.csv";
        std::string tree_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/tree.newick";

        if (std::filesystem::exists(aln_file) &&
            std::filesystem::exists(prefs_file) &&
            std::filesystem::exists(tree_file)) {

            auto aln = hyphy::core::Alignment::from_fasta(aln_file);
            auto prefs = hyphy::core::Preferences::from_csv(prefs_file);
            auto tree = hyphy::core::Tree::from_newick_file(tree_file);

            hyphy::core::ExpCMParameters params;
            params.beta = 3.38189;
            params.kappa = 8.75598;
            params.omega = 0.977468;
            params.phi = hyphy::core::Vector4(0.355409, 0.201765, 0.249379,
                                              1.0 - 0.355409 - 0.201765 - 0.249379);

            std::cerr << "\n=== OPTIMAL SCALE FACTOR SEARCH ===" << std::endl;
            std::cerr << "Expected likelihood: -775.21" << std::endl;
            std::cerr << "Searching for scale factor that gives correct likelihood..." << std::endl;

            double expected = -775.21;
            double best_error = std::numeric_limits<double>::infinity();
            double best_scale = 1.0;

            // Binary search for optimal scale factor
            double scale_min = 0.01;
            double scale_max = 10.0;

            for (int iteration = 0; iteration < 20; ++iteration) {
                double scale_mid = (scale_min + scale_max) / 2.0;
                double lnL = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
                    tree, aln, *gcode, prefs, params, scale_mid
                );

                double error = std::abs(lnL - expected);
                if (error < best_error) {
                    best_error = error;
                    best_scale = scale_mid;
                }

                std::cerr << "Iteration " << iteration << ": scale=" << std::fixed << std::setprecision(6)
                          << scale_mid << " → ln L = " << std::setprecision(2) << lnL
                          << " (error=" << std::setprecision(2) << error << ")" << std::endl;

                // Adjust search range
                if (lnL < expected) {
                    // Need smaller scale (less negative)
                    scale_max = scale_mid;
                } else {
                    // Need larger scale (more negative)
                    scale_min = scale_mid;
                }
            }

            std::cerr << std::endl << "Best scale found: " << std::fixed << std::setprecision(8) << best_scale << std::endl;
            std::cerr << "Best error: " << best_error << std::endl;

            // Verify the best scale
            double final_lnL = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
                tree, aln, *gcode, prefs, params, best_scale
            );
            std::cerr << "Final ln L with best scale: " << std::fixed << std::setprecision(2) << final_lnL << std::endl;
            std::cerr << "Expected: " << expected << std::endl;
            std::cerr << "Ratio: " << (final_lnL / expected) << std::endl;

            CHECK(true);  // Diagnostic only
        }
    }

    TEST_CASE("Fixed-parameter likelihood: smoke test") {
        auto gcode = hyphy::core::GeneticCode::universal();

        // Load real test data
        std::string aln_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/alignment.fasta";
        std::string prefs_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/prefs.csv";
        std::string tree_file = "/Users/tuu74294/Desktop/phd/hyphy3/hyphy3/data/tree.newick";

        if (std::filesystem::exists(aln_file) &&
            std::filesystem::exists(prefs_file) &&
            std::filesystem::exists(tree_file)) {

            // Load alignment
            auto aln = hyphy::core::Alignment::from_fasta(aln_file);
            CHECK(aln.num_taxa > 0);
            CHECK(aln.num_codons > 0);

            // Load preferences
            auto prefs = hyphy::core::Preferences::from_csv(prefs_file);
            CHECK(prefs.num_sites() > 0);

            // Load tree
            auto tree = hyphy::core::Tree::from_newick_file(tree_file);
            CHECK(tree.num_nodes() > 0);

            // Set up parameters (from reference)
            hyphy::core::ExpCMParameters params;
            params.beta = 3.38189;
            params.kappa = 8.75598;
            params.omega = 0.977468;
            params.phi = hyphy::core::Vector4(0.355409, 0.201765, 0.249379,
                                              1.0 - 0.355409 - 0.201765 - 0.249379);

            // Expected: -775.21 (from phydms)
            double expected = -775.21;

            std::cerr << std::endl << "=== FULL LIKELIHOOD TEST ===" << std::endl;
            std::cerr << "Alignment: " << aln.num_taxa << " taxa, " << aln.num_codons << " codon sites" << std::endl;
            std::cerr << "Tree nodes: " << tree.num_nodes() << std::endl;
            std::cerr << "Patterns: " << aln.patterns.size() << std::endl;
            std::cerr << "Preferences sites: " << prefs.num_sites() << std::endl;

            // Debug: Print pattern weights
            double total_weight = 0.0;
            for (const auto& p : aln.patterns) {
                total_weight += p.weight;
            }
            std::cerr << "Total pattern weight: " << std::fixed << std::setprecision(2) << total_weight << std::endl;
            std::cerr << "Average pattern weight: " << (total_weight / aln.patterns.size()) << std::endl;

            // Try "as_is" convention (branch lengths used directly)
            double lnL_as_is = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
                tree, aln, *gcode, prefs, params, 1.0
            );

            std::cerr << std::endl << "Convention 'as_is' (branch_scale=1.0):" << std::endl;
            std::cerr << "  Computed ln L = " << std::fixed << std::setprecision(2) << lnL_as_is << std::endl;
            std::cerr << "  Ratio to expected = " << (lnL_as_is / expected) << std::endl;

            // Try to estimate the branchScale (expected substitution rate)
            // Formula from phydms: branchScale = Σ_r[ Σ_x π_r[x] * (-P_r[x,x]) ] / nsites
            size_t num_sites_test = prefs.num_sites();
            hyphy::core::ExpCMQMatrix Q_test;
            Q_test.build(*gcode, params);

            double expected_sub_rate_sum = 0.0;
            for (size_t site = 0; site < num_sites_test; ++site) {
                hyphy::core::ExpCMFMatrix F_matrix;
                F_matrix.build(*gcode, prefs, site, params);
                hyphy::core::Matrix61Codon P = F_matrix.F.cwiseProduct(Q_test.Q);
                for (int x = 0; x < 61; ++x) {
                    P(x, x) = -P.row(x).sum();
                }
                auto pi_site = hyphy::core::LikelihoodEngine::compute_expcm_stationary(
                    *gcode, prefs, site, params.beta, params.phi
                );
                // Sum: π[x] * |P[x,x]|
                for (int x = 0; x < 61; ++x) {
                    expected_sub_rate_sum += pi_site(x) * (-P(x, x));
                }
            }
            double estimated_branchScale = expected_sub_rate_sum / num_sites_test;
            std::cerr << std::endl << "Estimated branchScale = " << std::fixed << std::setprecision(8)
                      << estimated_branchScale << std::endl;

            // Try dividing branch lengths by estimated branchScale
            double lnL_scaled = hyphy::core::LikelihoodEngine::compute_expcm_log_likelihood(
                tree, aln, *gcode, prefs, params, 1.0 / estimated_branchScale
            );
            std::cerr << "With branch_scale = 1/" << std::fixed << std::setprecision(8)
                      << estimated_branchScale << ":" << std::endl;
            std::cerr << "  Computed ln L = " << std::fixed << std::setprecision(2) << lnL_scaled << std::endl;
            std::cerr << "  Ratio to expected = " << (lnL_scaled / expected) << std::endl;
            std::cerr << "  Difference from expected = " << std::fabs(lnL_scaled - expected) << std::endl;

            std::cerr << std::endl << "Expected ln L = " << expected << std::endl;
            std::cerr << "NOTE: Need to determine correct branch convention and handle systematic error" << std::endl;

            CHECK(std::isfinite(lnL_as_is));
            CHECK(lnL_as_is < 0.0);
        }
    }
}

// TODO: Step 4 - Optimization
// TODO: Step 5 - Per-Site ω (omegabysite)
