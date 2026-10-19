#include "hyphy/analyses/persite_optimizer.hpp"
#include "hyphy/core/likelihood.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>

namespace hyphy::analyses {

using hyphy::core::Scalar;
using hyphy::core::GeneticCode;
using hyphy::core::Preferences;
using hyphy::core::ExpCMParameters;
using hyphy::core::Tree;
using hyphy::core::Alignment;
using hyphy::core::LikelihoodEngine;

double PersiteOptimizer::compute_site_likelihood(
    const Tree& tree,
    const Alignment& aln,
    const GeneticCode& gcode,
    const Preferences& prefs,
    size_t site,
    double beta,
    double mu,
    double omega,
    double kappa
) {
    // Note: μ (mu) parameter controls the overall synonymous rate
    // For simplicity, we adjust kappa to account for mu scaling
    // This is a simplified approach; full ExpCM would scale all rates by mu

    ExpCMParameters params;
    params.beta = beta;
    params.omega = omega;
    params.kappa = kappa * mu;  // Scale transition/transversion rate
    params.phi = {0.35, 0.20, 0.25, 0.20};

    Scalar lnl = LikelihoodEngine::compute_expcm_site_log_likelihood(
        tree, aln, gcode, prefs, site, params
    );

    return static_cast<double>(lnl);
}

double PersiteOptimizer::chi2_pvalue(double x, int df) {
    // Incomplete gamma function approximation for chi-squared CDF
    // P(χ²_df > x) using standard normal approximation for large df

    if (x <= 0) return 1.0;
    if (x > 100) return 0.0;

    // For df=1, use error function based approximation
    if (df == 1) {
        // χ²_1 critical value to p-value
        double z = std::sqrt(x);
        // erf-based approximation
        return std::erfc(z / std::sqrt(2.0));
    }

    // For other df, use Pearson's approximation
    double mean = static_cast<double>(df);
    double variance = 2.0 * mean;
    double z = (x - mean) / std::sqrt(variance);

    // Normal approximation
    return 0.5 * std::erfc(z / std::sqrt(2.0));
}

PersiteOptimizer::SiteResult PersiteOptimizer::optimize_site_omega_grid(
    const Tree& tree,
    const Alignment& aln,
    const GeneticCode& gcode,
    const Preferences& prefs,
    size_t site,
    double beta_fixed,
    double kappa_fixed
) {
    SiteResult result;
    result.site_id = site;

    // Null hypothesis: ω = 1.0, with fixed μ = 1.0
    result.lnl_null = compute_site_likelihood(
        tree, aln, gcode, prefs, site,
        beta_fixed, 1.0, 1.0, kappa_fixed
    );

    // Grid search for best ω (with fixed μ = 1.0 for speed)
    // Note: Full implementation would optimize μ_r as well per spec
    double best_lnl = result.lnl_null;
    double best_omega = 1.0;

    for (size_t i = 0; i < OMEGA_GRID_SIZE; ++i) {
        double omega = OMEGA_GRID[i];
        double lnl = compute_site_likelihood(
            tree, aln, gcode, prefs, site,
            beta_fixed, 1.0, omega, kappa_fixed
        );

        if (lnl > best_lnl) {
            best_lnl = lnl;
            best_omega = omega;
        }
    }

    result.lnl_alt = best_lnl;
    result.omega_mle = best_omega;

    // Compute LRT statistic
    result.test_statistic = 2.0 * (result.lnl_alt - result.lnl_null);
    if (result.test_statistic < 0.0) result.test_statistic = 0.0;

    // Compute p-value
    result.p_value = chi2_pvalue(result.test_statistic, 1);

    return result;
}

PersiteOptimizer::SiteResult PersiteOptimizer::optimize_site_omega_1d(
    const Tree& tree,
    const Alignment& aln,
    const GeneticCode& gcode,
    const Preferences& prefs,
    size_t site,
    double beta_fixed,
    double kappa_fixed
) {
    // For now, use grid search as basis and refine
    // Full 1D optimization (e.g., Brent's method) could be added later
    return optimize_site_omega_grid(
        tree, aln, gcode, prefs, site,
        beta_fixed, kappa_fixed
    );
}

std::vector<PersiteOptimizer::SiteResult> PersiteOptimizer::optimize_all_sites(
    const Tree& tree,
    const Alignment& aln,
    const GeneticCode& gcode,
    const Preferences& prefs,
    double beta_fixed,
    double kappa_fixed
) {
    size_t num_sites = prefs.num_sites();
    std::vector<SiteResult> results;
    results.reserve(num_sites);

    #pragma omp parallel for schedule(dynamic)
    for (size_t site = 0; site < num_sites; ++site) {
        SiteResult result = optimize_site_omega_grid(
            tree, aln, gcode, prefs, site,
            beta_fixed, kappa_fixed
        );
        #pragma omp critical
        results.push_back(result);
    }

    // Sort by site ID to maintain order
    std::sort(results.begin(), results.end(),
        [](const SiteResult& a, const SiteResult& b) {
            return a.site_id < b.site_id;
        });

    return results;
}

void PersiteOptimizer::apply_fdr_correction(
    std::vector<SiteResult>& results,
    double alpha
) {
    if (results.empty()) return;

    // Separate into ω > 1 and ω < 1 groups (per phydms spec)
    std::vector<size_t> up_indices, down_indices;

    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].omega_mle > 1.0) {
            up_indices.push_back(i);
        } else if (results[i].omega_mle < 1.0) {
            down_indices.push_back(i);
        }
        // Sites with ω = 1.0 exactly get q = 1.0
        results[i].q_value = 1.0;
    }

    // Apply BH correction separately to each group
    auto apply_bh_to_group = [&](const std::vector<size_t>& group_indices) {
        if (group_indices.empty()) return;

        size_t m = group_indices.size();

        // Sort by p-value within group
        std::vector<size_t> sorted_indices = group_indices;
        std::sort(sorted_indices.begin(), sorted_indices.end(),
            [&results](size_t a, size_t b) {
                return results[a].p_value < results[b].p_value;
            });

        // Find largest i where p[i] <= (i/m) * alpha
        int max_i = -1;
        for (size_t i = 0; i < m; ++i) {
            double threshold = static_cast<double>(i + 1) / m * alpha;
            if (results[sorted_indices[i]].p_value <= threshold) {
                max_i = i;
            }
        }

        // Compute q-values for this group: q[j] = min{p[k] * m/k : k >= j}
        if (max_i >= 0) {
            double min_q = 1.0;
            for (int i = m - 1; i >= 0; --i) {
                double q_candidate = results[sorted_indices[i]].p_value * m / (i + 1);
                min_q = std::min(min_q, q_candidate);
                results[sorted_indices[i]].q_value = min_q;
            }
        }
    };

    apply_bh_to_group(up_indices);
    apply_bh_to_group(down_indices);
}

void PersiteOptimizer::write_omegabysite_file(
    const std::vector<SiteResult>& results,
    const std::string& filename,
    bool use_1based
) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // Write header
    file << "site\tomega\tP\tdLnL\tQ\n";

    // Write each site's results
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        int site_index = use_1based ? (r.site_id + 1) : r.site_id;
        double dlnl = r.lnl_alt - r.lnl_null;

        file << site_index << "\t"
             << std::fixed << std::setprecision(6) << r.omega_mle << "\t"
             << std::scientific << std::setprecision(6) << r.p_value << "\t"
             << std::fixed << std::setprecision(6) << dlnl << "\t"
             << std::scientific << std::setprecision(6) << r.q_value << "\n";
    }

    file.close();
}

} // namespace hyphy::analyses
