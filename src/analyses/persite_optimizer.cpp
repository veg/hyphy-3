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
using hyphy::core::Vector4;
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
    double kappa,
    const Vector4& phi
) {
    // μ_r (mu) scales the site-specific rate (i.e. branch_scale = mu), NOT kappa * mu!
    ExpCMParameters params;
    params.beta = beta;
    params.omega = omega;
    params.kappa = kappa;
    params.phi = phi;

    Scalar lnl = LikelihoodEngine::compute_expcm_site_log_likelihood(
        tree, aln, gcode, prefs, site, params, mu
    );

    return static_cast<double>(lnl);
}

std::pair<double, double> PersiteOptimizer::optimize_mu_for_omega(
    const Tree& tree,
    const Alignment& aln,
    const GeneticCode& gcode,
    const Preferences& prefs,
    size_t site,
    double beta,
    double kappa,
    double omega,
    const Vector4& phi
) {
    // Golden-section search over log(mu) in [log(0.01), log(20.0)]
    const double inv_phi = (std::sqrt(5.0) - 1.0) * 0.5;
    double a = std::log(0.01);
    double b = std::log(20.0);
    double c = b - inv_phi * (b - a);
    double d = a + inv_phi * (b - a);

    double fc = compute_site_likelihood(tree, aln, gcode, prefs, site, beta, std::exp(c), omega, kappa, phi);
    double fd = compute_site_likelihood(tree, aln, gcode, prefs, site, beta, std::exp(d), omega, kappa, phi);

    for (int iter = 0; iter < 15; ++iter) {
        if (fc > fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - inv_phi * (b - a);
            fc = compute_site_likelihood(tree, aln, gcode, prefs, site, beta, std::exp(c), omega, kappa, phi);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + inv_phi * (b - a);
            fd = compute_site_likelihood(tree, aln, gcode, prefs, site, beta, std::exp(d), omega, kappa, phi);
        }
    }

    double best_mu = std::exp(0.5 * (a + b));
    double best_lnl = compute_site_likelihood(tree, aln, gcode, prefs, site, beta, best_mu, omega, kappa, phi);
    double lnl_unit = compute_site_likelihood(tree, aln, gcode, prefs, site, beta, 1.0, omega, kappa, phi);
    if (lnl_unit >= best_lnl) {
        return {1.0, lnl_unit};
    }
    return {best_mu, best_lnl};
}

double PersiteOptimizer::chi2_pvalue(double x, int df) {
    if (x <= 0) return 1.0;
    if (x > 100) return 0.0;

    if (df == 1) {
        double z = std::sqrt(x);
        return std::erfc(z / std::sqrt(2.0));
    }

    double mean = static_cast<double>(df);
    double variance = 2.0 * mean;
    double z = (x - mean) / std::sqrt(variance);
    return 0.5 * std::erfc(z / std::sqrt(2.0));
}

PersiteOptimizer::SiteResult PersiteOptimizer::optimize_site_omega_grid(
    const Tree& tree,
    const Alignment& aln,
    const GeneticCode& gcode,
    const Preferences& prefs,
    size_t site,
    double beta_fixed,
    double kappa_fixed,
    const Vector4& phi
) {
    SiteResult result;
    result.site_id = site;

    // Null hypothesis: ω = 1.0, with fixed μ = 1.0
    result.lnl_null = compute_site_likelihood(
        tree, aln, gcode, prefs, site,
        beta_fixed, 1.0, 1.0, kappa_fixed, phi
    );

    // 1. Coarse 20-point grid search for best ω
    double best_lnl = result.lnl_null;
    double best_omega = 1.0;
    size_t best_idx = 12; // index of 1.0 in OMEGA_GRID

    for (size_t i = 0; i < OMEGA_GRID_SIZE; ++i) {
        double omega = OMEGA_GRID[i];
        double lnl = compute_site_likelihood(
            tree, aln, gcode, prefs, site,
            beta_fixed, 1.0, omega, kappa_fixed, phi
        );

        if (lnl > best_lnl + 1e-12) {
            best_lnl = lnl;
            best_omega = omega;
            best_idx = i;
        }
    }

    // 2. Continuous golden-section polish in log(ω) around the best grid point
    if (best_lnl > result.lnl_null + 1e-9) {
        double omega_low = (best_idx > 0) ? OMEGA_GRID[best_idx - 1] : OMEGA_GRID[0];
        double omega_high = (best_idx + 1 < OMEGA_GRID_SIZE) ? OMEGA_GRID[best_idx + 1] : OMEGA_GRID[OMEGA_GRID_SIZE - 1];

        const double inv_phi = (std::sqrt(5.0) - 1.0) * 0.5;
        double a = std::log(omega_low);
        double b = std::log(omega_high);
        double c = b - inv_phi * (b - a);
        double d = a + inv_phi * (b - a);

        double fc = compute_site_likelihood(tree, aln, gcode, prefs, site, beta_fixed, 1.0, std::exp(c), kappa_fixed, phi);
        double fd = compute_site_likelihood(tree, aln, gcode, prefs, site, beta_fixed, 1.0, std::exp(d), kappa_fixed, phi);

        for (int iter = 0; iter < 16; ++iter) {
            if (fc > fd) {
                b = d;
                d = c;
                fd = fc;
                c = b - inv_phi * (b - a);
                fc = compute_site_likelihood(tree, aln, gcode, prefs, site, beta_fixed, 1.0, std::exp(c), kappa_fixed, phi);
            } else {
                a = c;
                c = d;
                fc = fd;
                d = a + inv_phi * (b - a);
                fd = compute_site_likelihood(tree, aln, gcode, prefs, site, beta_fixed, 1.0, std::exp(d), kappa_fixed, phi);
            }
        }

        double polished_omega = std::clamp(std::exp(0.5 * (a + b)), OMEGA_GRID[0], OMEGA_GRID[OMEGA_GRID_SIZE - 1]);
        double polished_lnl = compute_site_likelihood(
            tree, aln, gcode, prefs, site, beta_fixed, 1.0, polished_omega, kappa_fixed, phi
        );
        if (polished_lnl > best_lnl) {
            best_lnl = polished_lnl;
            best_omega = polished_omega;
        }
    }

    result.lnl_alt = best_lnl;
    result.omega_mle = best_omega;

    result.test_statistic = 2.0 * (result.lnl_alt - result.lnl_null);
    if (result.test_statistic < 0.0) result.test_statistic = 0.0;

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
    double kappa_fixed,
    const Vector4& phi
) {
    return optimize_site_omega_grid(
        tree, aln, gcode, prefs, site,
        beta_fixed, kappa_fixed, phi
    );
}

std::vector<PersiteOptimizer::SiteResult> PersiteOptimizer::optimize_all_sites(
    const Tree& tree,
    const Alignment& aln,
    const GeneticCode& gcode,
    const Preferences& prefs,
    double beta_fixed,
    double kappa_fixed,
    const Vector4& phi
) {
    size_t num_sites = prefs.num_sites();
    std::vector<SiteResult> results(num_sites);

    #pragma omp parallel for schedule(dynamic)
    for (size_t site = 0; site < num_sites; ++site) {
        results[site] = optimize_site_omega_grid(
            tree, aln, gcode, prefs, site,
            beta_fixed, kappa_fixed, phi
        );
    }

    return results;
}

void PersiteOptimizer::apply_fdr_correction(
    std::vector<SiteResult>& results,
    double alpha
) {
    if (results.empty()) return;

    std::vector<size_t> up_indices, down_indices;

    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].omega_mle > 1.0) {
            up_indices.push_back(i);
        } else if (results[i].omega_mle < 1.0) {
            down_indices.push_back(i);
        }
        results[i].q_value = 1.0;
    }

    auto apply_bh_to_group = [&](const std::vector<size_t>& group_indices) {
        if (group_indices.empty()) return;

        size_t m = group_indices.size();

        std::vector<size_t> sorted_indices = group_indices;
        std::sort(sorted_indices.begin(), sorted_indices.end(),
            [&results](size_t a, size_t b) {
                return results[a].p_value < results[b].p_value;
            });

        int max_i = -1;
        for (size_t i = 0; i < m; ++i) {
            double threshold = static_cast<double>(i + 1) / m * alpha;
            if (results[sorted_indices[i]].p_value <= threshold) {
                max_i = static_cast<int>(i);
            }
        }

        if (max_i >= 0) {
            double min_q = 1.0;
            for (int i = static_cast<int>(m) - 1; i >= 0; --i) {
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

    file << "site\tomega\tP\tdLnL\tQ\n";

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        int site_index = use_1based ? static_cast<int>(r.site_id + 1) : static_cast<int>(r.site_id);
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
