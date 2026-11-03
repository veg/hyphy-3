#pragma once

#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/types.hpp"
#include <Eigen/Core>
#include <vector>
#include <cmath>

namespace hyphy::analyses {

/**
 * Per-Site Parameter Optimizer
 *
 * Optimizes individual site parameters (ω_r and optionally μ_r) while keeping
 * global parameters (β, κ, φ) fixed from global optimization.
 */
class PersiteOptimizer {
public:
    struct SiteResult {
        size_t site_id = 0;
        double omega_mle = 1.0;       // Maximum likelihood estimate
        double lnl_null = 0.0;        // Likelihood under H0 (ω = 1.0)
        double lnl_alt = 0.0;         // Likelihood under H1 (ω = ω_mle)
        double test_statistic = 0.0;  // 2(lnL_alt - lnL_null)
        double p_value = 1.0;         // Chi-squared p-value
        double q_value = 1.0;         // FDR-corrected p-value (computed later)
    };

    static SiteResult optimize_site_omega_grid(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        size_t site,
        double beta_fixed,
        double kappa_fixed,
        const hyphy::core::Vector4& phi = {0.35, 0.20, 0.25, 0.20}
    );

    static SiteResult optimize_site_omega_1d(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        size_t site,
        double beta_fixed,
        double kappa_fixed,
        const hyphy::core::Vector4& phi = {0.35, 0.20, 0.25, 0.20}
    );

    static std::vector<SiteResult> optimize_all_sites(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        double beta_fixed,
        double kappa_fixed,
        const hyphy::core::Vector4& phi = {0.35, 0.20, 0.25, 0.20}
    );

    static void apply_fdr_correction(
        std::vector<SiteResult>& results,
        double alpha = 0.05
    );

    static void write_omegabysite_file(
        const std::vector<SiteResult>& results,
        const std::string& filename,
        bool use_1based = true
    );

    /**
     * Compute likelihood for a site under specific μ_r (branch_scale) and ω_r.
     */
    static double compute_site_likelihood(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        size_t site,
        double beta,
        double mu,
        double omega,
        double kappa,
        const hyphy::core::Vector4& phi = {0.35, 0.20, 0.25, 0.20}
    );

    /**
     * Optimize μ_r (site-specific synonymous rate / branch-length scaler) for a given ω
     * using golden-section search.
     */
    static std::pair<double, double> optimize_mu_for_omega(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        size_t site,
        double beta,
        double kappa,
        double omega,
        const hyphy::core::Vector4& phi = {0.35, 0.20, 0.25, 0.20}
    );

private:
    static constexpr double OMEGA_GRID[] = {
        0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0,
        3.0, 5.0, 10.0, 100.0
    };
    static constexpr size_t OMEGA_GRID_SIZE = 20;

    static double chi2_pvalue(double x, int df = 1);
};

} // namespace hyphy::analyses
