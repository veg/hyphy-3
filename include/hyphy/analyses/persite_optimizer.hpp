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
 * Optimizes individual site parameters (currently ω) while keeping
 * global parameters (β, κ) fixed from Phase 3 optimization.
 */
class PersiteOptimizer {
public:
    // Result for single-site optimization
    struct SiteResult {
        size_t site_id;
        double omega_mle;           // Maximum likelihood estimate
        double lnl_null;            // Likelihood under H0 (ω = ω_global)
        double lnl_alt;             // Likelihood under H1 (ω = ω_mle)
        double test_statistic;      // 2(lnL_alt - lnL_null)
        double p_value;             // Chi-squared p-value
        double q_value;             // FDR-corrected p-value (computed later)
    };

    /**
     * Optimize ω for a single site using grid search
     *
     * @param tree           Phylogenetic tree
     * @param aln            Multiple sequence alignment
     * @param gcode          Genetic code
     * @param prefs          Site-specific preferences
     * @param site           Site index
     * @param beta_fixed     Fixed β parameter (from global optimization)
     * @param kappa_fixed    Fixed κ parameter (from global optimization)
     * @return               SiteResult with optimization results
     *                       Null hypothesis: ω = 1.0 (per spec)
     */
    static SiteResult optimize_site_omega_grid(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        size_t site,
        double beta_fixed,
        double kappa_fixed
    );

    /**
     * Optimize ω for a single site using 1D optimization
     * More accurate than grid search but slower
     *
     * @param tree           Phylogenetic tree
     * @param aln            Multiple sequence alignment
     * @param gcode          Genetic code
     * @param prefs          Site-specific preferences
     * @param site           Site index
     * @param beta_fixed     Fixed β parameter
     * @param kappa_fixed    Fixed κ parameter
     * @return               SiteResult with optimization results
     *                       Null hypothesis: ω = 1.0 (per spec)
     */
    static SiteResult optimize_site_omega_1d(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        size_t site,
        double beta_fixed,
        double kappa_fixed
    );

    /**
     * Optimize ω for all sites and return results
     * Uses grid search for speed (recommended)
     *
     * @param tree           Phylogenetic tree
     * @param aln            Multiple sequence alignment
     * @param gcode          Genetic code
     * @param prefs          Site-specific preferences
     * @param beta_fixed     Fixed β parameter
     * @param kappa_fixed    Fixed κ parameter
     * @return               Vector of SiteResult for all sites
     *                       Null hypothesis: ω = 1.0 (per spec)
     */
    static std::vector<SiteResult> optimize_all_sites(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        double beta_fixed,
        double kappa_fixed
    );

    /**
     * Apply Benjamini-Hochberg FDR correction to p-values
     * Modifies q_value field of each result in-place
     * Computes q-values separately for sites with ω > 1 and ω < 1 (per phydms)
     *
     * @param results        Vector of SiteResult to correct (modified in-place)
     * @param alpha          FDR level (default 0.05)
     */
    static void apply_fdr_correction(
        std::vector<SiteResult>& results,
        double alpha = 0.05
    );

    /**
     * Write results to _omegabysite.txt file in phydms format
     *
     * Format (tab-separated):
     *   site  omega  P  dLnL  Q
     * Where:
     *   - site: index (0-based or 1-based per convention)
     *   - omega: ω_r MLE value
     *   - P: p-value from LRT test
     *   - dLnL: ln L_alt - ln L_null (log-likelihood difference)
     *   - Q: FDR-corrected q-value
     *
     * @param results        Vector of SiteResult (must have q_values computed)
     * @param filename       Output file path
     * @param use_1based     If true, use 1-based indexing; else 0-based
     */
    static void write_omegabysite_file(
        const std::vector<SiteResult>& results,
        const std::string& filename,
        bool use_1based = true
    );

private:
    // Grid of ω values for grid search
    static constexpr double OMEGA_GRID[] = {
        0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0,
        3.0, 5.0, 10.0, 100.0
    };
    static constexpr size_t OMEGA_GRID_SIZE = 20;

    /**
     * Compute likelihood for a site under specific μ_r and ω
     *
     * @param mu             Synonymous rate for the site
     * @param omega          Non-synonymous/synonymous rate ratio
     * @return               Log-likelihood for the site (negated for minimization)
     */
    struct MuOptimizationContext {
        const hyphy::core::Tree* tree;
        const hyphy::core::Alignment* aln;
        const hyphy::core::GeneticCode* gcode;
        const hyphy::core::Preferences* prefs;
        size_t site;
        double beta;
        double kappa;
        double omega;
    };

    static double evaluate_likelihood_for_mu(double mu, const MuOptimizationContext& ctx);

    /**
     * Optimize μ_r for a given ω using golden section search
     *
     * @return               Optimal μ_r value and corresponding ln L
     */
    static std::pair<double, double> optimize_mu_for_omega(
        const hyphy::core::Tree& tree,
        const hyphy::core::Alignment& aln,
        const hyphy::core::GeneticCode& gcode,
        const hyphy::core::Preferences& prefs,
        size_t site,
        double beta,
        double kappa,
        double omega
    );

    /**
     * Compute likelihood for a site under specific μ, ω (internal use)
     *
     * @return               Log-likelihood for the site
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
        double kappa
    );

    /**
     * Chi-squared CDF for computing p-values
     *
     * @param x              Test statistic value
     * @param df             Degrees of freedom (typically 1)
     * @return               P-value (upper tail probability)
     */
    static double chi2_pvalue(double x, int df = 1);
};

} // namespace hyphy::analyses
