#pragma once

#include <hyphy/core/types.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace hyphy::core::stats {

/**
 * Lower regularized incomplete gamma function P(a, z) = gamma(a, z) / Gamma(a).
 * Uses series expansion for z < a + 1 and Lentz's continued fraction for z >= a + 1.
 */
inline Scalar regularized_gamma_p(Scalar a, Scalar z) {
    if (z <= 0.0) return 0.0;
    if (a <= 0.0) return 1.0;
    if (z < a + 1.0) {
        Scalar sum = 1.0 / a;
        Scalar term = sum;
        for (int n = 1; n < 200; ++n) {
            term *= z / (a + n);
            sum += term;
            if (std::abs(term) < std::abs(sum) * 1e-15) break;
        }
        return std::clamp(sum * std::exp(-z + a * std::log(z) - std::lgamma(a)), 0.0, 1.0);
    }
    Scalar tiny = 1e-30;
    Scalar b = z + 1.0 - a;
    Scalar c = 1.0 / tiny;
    Scalar d = 1.0 / b;
    Scalar h = d;
    for (int i = 1; i < 200; ++i) {
        Scalar an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < tiny) d = tiny;
        c = b + an / c;
        if (std::abs(c) < tiny) c = tiny;
        d = 1.0 / d;
        Scalar del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < 1e-15) break;
    }
    Scalar q = std::exp(-z + a * std::log(z) - std::lgamma(a)) * h;
    return std::clamp(1.0 - q, 0.0, 1.0);
}

/**
 * Upper regularized incomplete gamma function Q(a, z) = 1 - P(a, z).
 * Uses Lentz's continued fraction directly for z >= a + 1 to avoid catastrophic cancellation in tails.
 */
inline Scalar regularized_gamma_q(Scalar a, Scalar z) {
    if (z <= 0.0) return 1.0;
    if (a <= 0.0) return 0.0;
    if (z < a + 1.0) {
        return std::clamp(1.0 - regularized_gamma_p(a, z), 0.0, 1.0);
    }
    Scalar tiny = 1e-30;
    Scalar b = z + 1.0 - a;
    Scalar c = 1.0 / tiny;
    Scalar d = 1.0 / b;
    Scalar h = d;
    for (int i = 1; i < 200; ++i) {
        Scalar an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < tiny) d = tiny;
        c = b + an / c;
        if (std::abs(c) < tiny) c = tiny;
        d = 1.0 / d;
        Scalar del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < 1e-15) break;
    }
    Scalar q = std::exp(-z + a * std::log(z) - std::lgamma(a)) * h;
    return std::clamp(q, 0.0, 1.0);
}

/**
 * Right-tail p-value for a chi-squared test statistic with `df` degrees of freedom.
 * Fast exact paths for df = 1 (erfc) and df = 2 (exp), general incomplete gamma otherwise.
 */
inline Scalar chi2_p_value(Scalar chi2, Scalar df = 1.0) {
    if (chi2 <= 0.0 || df <= 0.0) return 1.0;
    if (df == 1.0) {
        return std::clamp(std::erfc(std::sqrt(chi2 / 2.0)), 0.0, 1.0);
    }
    if (df == 2.0) {
        return std::clamp(std::exp(-chi2 / 2.0), 0.0, 1.0);
    }
    return regularized_gamma_q(df / 2.0, chi2 / 2.0);
}

/**
 * Asymptotic right-tail p-value for a boundary mixture of chi-squared distributions:
 *   w0 * chi2_0 + w1 * chi2_1 + w2 * chi2_2
 * Default weights correspond to MEME's mixture: (1/3)*chi2_0 + 0.30*chi2_1 + 0.3667*chi2_2.
 */
inline Scalar mixture_chi2_p_value(Scalar lrt, Scalar w1, Scalar w2) {
    if (lrt <= 0.0) return 1.0;
    Scalar p1 = std::erfc(std::sqrt(lrt / 2.0));
    Scalar p2 = std::exp(-lrt / 2.0);
    return std::clamp(w1 * p1 + w2 * p2, 0.0, 1.0);
}

/**
 * Small-sample corrected Akaike Information Criterion (AICc):
 *   AICc = -2 * log_l + 2 * k * n / (n - k - 1)
 * Falls back to standard AIC (-2 * log_l + 2 * k) when n <= k + 1.
 */
inline Scalar compute_aicc(Scalar log_l, size_t k_params, size_t n_obs) {
    Scalar k = static_cast<Scalar>(k_params);
    Scalar n = static_cast<Scalar>(n_obs);
    if (n > k + 1.0) {
        return -2.0 * log_l + (2.0 * k * n) / (n - k - 1.0);
    }
    return -2.0 * log_l + 2.0 * k;
}

/**
 * Benjamini-Hochberg False Discovery Rate (FDR) step-up q-values.
 * Preserves input ordering and guarantees monotonicity in sorted p-value rank.
 */
inline std::vector<Scalar> benjamini_hochberg_fdr(const std::vector<Scalar>& p_values) {
    size_t m = p_values.size();
    std::vector<Scalar> q_values(m, 1.0);
    if (m == 0) return q_values;

    std::vector<size_t> order(m);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return p_values[a] < p_values[b];
    });

    Scalar running_min_q = 1.0;
    for (int rank = static_cast<int>(m) - 1; rank >= 0; --rank) {
        size_t idx = order[rank];
        Scalar q = (p_values[idx] * static_cast<Scalar>(m)) / static_cast<Scalar>(rank + 1);
        running_min_q = std::clamp(std::min(running_min_q, q), 0.0, 1.0);
        q_values[idx] = running_min_q;
    }
    return q_values;
}

/**
 * Holm-Bonferroni step-down family-wise error rate (FWER) adjusted p-values.
 * Preserves input ordering and guarantees monotonicity in sorted p-value rank.
 */
inline std::vector<Scalar> holm_bonferroni(const std::vector<Scalar>& p_values) {
    size_t m = p_values.size();
    std::vector<Scalar> adj_p(m, 1.0);
    if (m == 0) return adj_p;

    std::vector<size_t> order(m);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return p_values[a] < p_values[b];
    });

    Scalar running_max = 0.0;
    for (size_t rank = 0; rank < m; ++rank) {
        size_t idx = order[rank];
        Scalar p = p_values[idx] * static_cast<Scalar>(m - rank);
        running_max = std::clamp(std::max(running_max, p), 0.0, 1.0);
        adj_p[idx] = running_max;
    }
    return adj_p;
}

} // namespace hyphy::core::stats
