#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/preferences.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <stdexcept>

namespace hyphy::core {

struct GTRParameters {
    Scalar theta_AC = 1.0;
    Scalar theta_AT = 1.0;
    Scalar theta_CG = 1.0;
    Scalar theta_CT = 1.0;
    Scalar theta_GT = 1.0;
    // theta_AG is fixed to 1.0 as reference
};

class GTRMatrix {
public:
    Matrix4 Q = Matrix4::Zero();
    Vector4 pi = Vector4::Zero();
    Vector4 sqrt_pi = Vector4::Zero();
    Vector4 inv_sqrt_pi = Vector4::Zero();

    // Eigendecomposition of S = diag(sqrt(pi)) * Q * diag(1/sqrt(pi))
    Vector4 eigenvalues = Vector4::Zero();
    Matrix4 eigenvectors = Matrix4::Zero(); // V
    Matrix4 inv_eigenvectors = Matrix4::Zero(); // V^T
    Matrix4 U_pi = Matrix4::Zero();
    Matrix4 V_pi = Matrix4::Zero();

    Scalar scale_factor = 1.0;

    void update(const GTRParameters& params, const Vector4& frequencies) {
        pi = frequencies;
        for (int i = 0; i < 4; ++i) {
            Scalar p = std::max(pi(i), Scalar(1e-12));
            sqrt_pi(i) = std::sqrt(p);
            inv_sqrt_pi(i) = 1.0 / sqrt_pi(i);
        }

        // Fill Q off-diagonals
        Scalar r[4][4] = {
            {0.0, params.theta_AC, 1.0, params.theta_AT},
            {params.theta_AC, 0.0, params.theta_CG, params.theta_CT},
            {1.0, params.theta_CG, 0.0, params.theta_GT},
            {params.theta_AT, params.theta_CT, params.theta_GT, 0.0}
        };

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i != j) {
                    Q(i, j) = r[i][j] * pi(j);
                } else {
                    Q(i, j) = 0.0;
                }
            }
        }

        // Set row diagonals
        for (int i = 0; i < 4; ++i) {
            Q(i, i) = -Q.row(i).sum();
        }

        // Compute scaling factor so that average rate is 1.0
        Scalar avg_rate = 0.0;
        for (int i = 0; i < 4; ++i) {
            avg_rate -= pi(i) * Q(i, i);
        }
        if (avg_rate > 0.0) {
            scale_factor = avg_rate;
            Q /= scale_factor;
        }

        // Symmetrize for eigendecomposition: S = diag(sqrt(pi)) * Q * diag(1/sqrt(pi))
        Matrix4 S = sqrt_pi.asDiagonal() * Q * inv_sqrt_pi.asDiagonal();
        // Force exact symmetry
        S = 0.5 * (S + S.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix4> solver(S);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("GTR eigendecomposition failed");
        }

        eigenvalues = solver.eigenvalues().cwiseMin(0.0);
        eigenvectors = solver.eigenvectors(); // V
        inv_eigenvectors = eigenvectors.transpose(); // V^T
        U_pi = inv_sqrt_pi.asDiagonal() * eigenvectors;
        V_pi = inv_eigenvectors * sqrt_pi.asDiagonal();
    }

    Matrix4 transition_matrix(Scalar branch_length) const {
        Vector4 exp_lambda = (eigenvalues * branch_length).array().exp();
        Matrix4 P = (U_pi * exp_lambda.asDiagonal()) * V_pi;
        // Eigendecomposition noise can yield small negative entries for
        // near-zero-probability transitions; probabilities cannot be negative.
        return P.cwiseMax(0.0);
    }
};

struct MG94Parameters {
    Scalar alpha = 1.0; // synonymous rate
    Scalar beta = 1.0;  // nonsynonymous rate (or omega = beta/alpha)
    Scalar theta_AC = 1.0;
    Scalar theta_AT = 1.0;
    Scalar theta_CG = 1.0;
    Scalar theta_CT = 1.0;
    Scalar theta_GT = 1.0;
    // theta_AG = 1.0 reference
    Scalar delta = 0.0; // Double-hit rate parameter (default 0.0 = disabled)
    Scalar psi = 0.0;   // Triple-hit rate parameter (default 0.0 = disabled)
    const std::vector<Scalar>* syn_pair_rates = nullptr; // Optional S*S flat array of relative synonymous exchangeabilities
};

class MG94Matrix {
public:
    Matrix Q;
    Vector pi;
    Vector sqrt_pi;
    Vector inv_sqrt_pi;

    Vector eigenvalues;
    Matrix eigenvectors;
    Matrix inv_eigenvectors;
    Matrix U_pi;
    Matrix V_pi;

    Scalar scale_factor = 1.0;

    struct SubstitutionRates {
        Scalar rate_single = 0.0;
        Scalar rate_double = 0.0;
        Scalar rate_triple = 0.0;
        Scalar total_rate = 0.0;
        Scalar frac_double = 0.0;
        Scalar frac_triple = 0.0;
    };

    SubstitutionRates compute_substitution_rates(const GeneticCode& code) const {
        SubstitutionRates sr;
        int S = code.num_sense_codons;
        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                Scalar flux = pi(i) * Q(i, j);
                int nd = code.diff_matrix[i][j].num_diffs;
                if (nd == 1) sr.rate_single += flux;
                else if (nd == 2) sr.rate_double += flux;
                else if (nd == 3) sr.rate_triple += flux;
            }
        }
        sr.total_rate = sr.rate_single + sr.rate_double + sr.rate_triple;
        if (sr.total_rate > 1e-12) {
            sr.frac_double = sr.rate_double / sr.total_rate;
            sr.frac_triple = sr.rate_triple / sr.total_rate;
        }
        return sr;
    }

    void update(
        const MG94Parameters& params,
        const Eigen::Matrix<Scalar, 3, 4>& pos_nuc_freqs,
        const Vector& codon_freqs,
        const GeneticCode& code = *GeneticCode::universal()
    ) {
        int S = code.num_sense_codons;
        pi = codon_freqs;
        sqrt_pi.resize(S);
        inv_sqrt_pi.resize(S);
        for (int i = 0; i < S; ++i) {
            sqrt_pi(i) = std::sqrt(std::max(pi(i), 1e-12));
            inv_sqrt_pi(i) = 1.0 / sqrt_pi(i);
        }

        Scalar r_nuc[4][4] = {
            {0.0, params.theta_AC, 1.0, params.theta_AT},
            {params.theta_AC, 0.0, params.theta_CG, params.theta_CT},
            {1.0, params.theta_CG, 0.0, params.theta_GT},
            {params.theta_AT, params.theta_CT, params.theta_GT, 0.0}
        };

        Q.resize(S, S);
        Q.setZero();

        const bool has_mss = (params.syn_pair_rates != nullptr &&
                              params.syn_pair_rates->size() >= static_cast<size_t>(S * S));

        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                const auto& diff = code.diff_matrix[i][j];
                Scalar rate_modifier = diff.synonymous ? params.alpha : params.beta;
                if (diff.num_diffs == 1) {
                    if (diff.synonymous && has_mss) {
                        rate_modifier *= (*params.syn_pair_rates)[i * S + j];
                    }
                    Scalar nuc_rate = r_nuc[diff.diff_from[0]][diff.diff_to[0]];
                    Scalar target_nuc_freq = pos_nuc_freqs(diff.diff_pos[0], diff.diff_to[0]);
                    Q(i, j) = rate_modifier * nuc_rate * target_nuc_freq;
                } else if (diff.num_diffs == 2 && params.delta > 0.0) {
                    Scalar nuc_rate = r_nuc[diff.diff_from[0]][diff.diff_to[0]] *
                                      r_nuc[diff.diff_from[1]][diff.diff_to[1]];
                    Scalar target_nuc_freq = pos_nuc_freqs(diff.diff_pos[0], diff.diff_to[0]) *
                                             pos_nuc_freqs(diff.diff_pos[1], diff.diff_to[1]);
                    Q(i, j) = rate_modifier * params.delta * nuc_rate * target_nuc_freq;
                } else if (diff.num_diffs == 3 && params.psi > 0.0) {
                    Scalar nuc_rate = r_nuc[diff.diff_from[0]][diff.diff_to[0]] *
                                      r_nuc[diff.diff_from[1]][diff.diff_to[1]] *
                                      r_nuc[diff.diff_from[2]][diff.diff_to[2]];
                    Scalar target_nuc_freq = pos_nuc_freqs(diff.diff_pos[0], diff.diff_to[0]) *
                                             pos_nuc_freqs(diff.diff_pos[1], diff.diff_to[1]) *
                                             pos_nuc_freqs(diff.diff_pos[2], diff.diff_to[2]);
                    Q(i, j) = rate_modifier * params.psi * nuc_rate * target_nuc_freq;
                }
            }
        }

        // Set row diagonals
        for (int i = 0; i < S; ++i) {
            Q(i, i) = -Q.row(i).sum();
        }

        // Average substitution rate
        Scalar avg_rate = 0.0;
        for (int i = 0; i < S; ++i) {
            avg_rate -= pi(i) * Q(i, i);
        }
        scale_factor = avg_rate;

        // Symmetrize for eigendecomposition: S = diag(sqrt(pi)) * Q * diag(1/sqrt(pi))
        Matrix S_mat = sqrt_pi.asDiagonal() * Q * inv_sqrt_pi.asDiagonal();
        S_mat = 0.5 * (S_mat + S_mat.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix> solver(S_mat);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("MG94 eigendecomposition failed");
        }

        eigenvalues = solver.eigenvalues().cwiseMin(0.0);
        eigenvectors = solver.eigenvectors();
        inv_eigenvectors = eigenvectors.transpose();
        U_pi = inv_sqrt_pi.asDiagonal() * eigenvectors;
        V_pi = inv_eigenvectors * sqrt_pi.asDiagonal();
    }

    Matrix transition_matrix(Scalar branch_length) const {
        Vector exp_lambda = (eigenvalues * branch_length).array().exp();
        Matrix P = (U_pi * exp_lambda.asDiagonal()) * V_pi;
        // Eigendecomposition noise can yield small negative entries for
        // near-zero-probability transitions; probabilities cannot be negative.
        return P.cwiseMax(0.0);
    }
};

/**
 * ExpCM (Experimentally informed Codon Model) Parameters
 *
 * Parameters for ExpCM: HKY-like mutation rates combined with
 * site-specific selection based on amino-acid preferences.
 */
struct ExpCMParameters {
    Scalar beta = 1.0;        // Stringency (preference exponent)
    Scalar omega = 0.5;       // Non-synonymous / synonymous rate ratio
    Scalar kappa = 2.0;       // Transition / transversion ratio
    Vector4 phi = Vector4::Ones() / 4.0;  // Nucleotide frequencies (must sum to 1)

    ExpCMParameters() = default;

    ExpCMParameters(Scalar b, Scalar w, Scalar k, const Vector4& p)
        : beta(b), omega(w), kappa(k), phi(p) {}

    void validate() const {
        if (beta <= 0.0 || beta > 10.0) {
            throw std::invalid_argument("beta must be in (0, 10]");
        }
        if (omega <= 0.0 || omega > 100.0) {
            throw std::invalid_argument("omega must be in (0, 100]");
        }
        if (kappa <= 0.0 || kappa > 100.0) {
            throw std::invalid_argument("kappa must be in (0, 100]");
        }
        Scalar phi_sum = phi.sum();
        if (std::fabs(phi_sum - 1.0) > 1e-6) {
            throw std::invalid_argument("phi must sum to 1.0");
        }
        if ((phi.array() <= 0.0).any()) {
            throw std::invalid_argument("phi values must be positive");
        }
    }
};

// Fixed-size 61x61 matrix and 61-vector for codon substitution rates
using Matrix61Codon = Eigen::Matrix<Scalar, 61, 61>;
using Vector61Codon = Eigen::Matrix<Scalar, 61, 1>;
static constexpr Scalar EXPCM_ALMOST_ZERO = 1e-8;

/**
 * ExpCM Q Matrix (Mutation Term)
 *
 * HKY-like substitution rates for codons. Independent of site-specific
 * preferences. Size: 61 × 61 (sense codons).
 */
class ExpCMQMatrix {
public:
    Matrix61Codon Q = Matrix61Codon::Zero();

    void build(const GeneticCode& gcode, const ExpCMParameters& params) {
        Q.setZero();

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            for (int y = 0; y < gcode.num_sense_codons; ++y) {
                if (x == y) continue;

                const auto& diff = gcode.diff_matrix[x][y];
                if (diff.num_diffs != 1) continue;

                int8_t nuc_to = diff.nuc_to;
                Q(x, y) = params.phi(nuc_to);

                int8_t nuc_from = diff.nuc_from;
                if (is_transition(nuc_from, nuc_to)) {
                    Q(x, y) *= params.kappa;
                }
            }
        }

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            Q(x, x) = -Q.row(x).sum();
        }
    }

private:
    static bool is_transition(int8_t nuc_from, int8_t nuc_to) {
        return (nuc_from == 0 && nuc_to == 2) ||
               (nuc_from == 2 && nuc_to == 0) ||
               (nuc_from == 1 && nuc_to == 3) ||
               (nuc_from == 3 && nuc_to == 1);
    }
};

/**
 * ExpCM F Matrix (Fixation Term)
 *
 * Site-specific selection based on amino-acid preferences.
 * For each site, creates a 61 × 61 matrix with selection values.
 */
class ExpCMFMatrix {
public:
    Matrix61Codon F = Matrix61Codon::Zero();

    void build(const GeneticCode& gcode, const Preferences& prefs,
               size_t site, const ExpCMParameters& params) {
        if (site >= prefs.num_sites()) {
            throw std::out_of_range("Site index out of range");
        }

        F.setOnes();

        for (int x = 0; x < gcode.num_sense_codons; ++x) {
            for (int y = 0; y < gcode.num_sense_codons; ++y) {
                if (x == y) continue;

                char aa_x = gcode.sense_codon_aa[x];
                char aa_y = gcode.sense_codon_aa[y];

                if (aa_x == aa_y) {
                    F(x, y) = 1.0;
                } else {
                    Scalar pi_x = prefs.at(site, aa_x);
                    Scalar pi_y = prefs.at(site, aa_y);

                    Scalar ratio = pi_x / pi_y;
                    Scalar ratio_beta = std::pow(ratio, params.beta);

                    if (std::fabs(1.0 - ratio_beta) < EXPCM_ALMOST_ZERO) {
                        F(x, y) = params.omega;
                    } else {
                        Scalar ln_ratio_beta = std::log(ratio_beta);
                        Scalar numerator = -ln_ratio_beta;
                        Scalar denominator = 1.0 - ratio_beta;
                        F(x, y) = params.omega * (numerator / denominator);
                    }
                }
            }
        }

        F.diagonal().setZero();
    }
};

/**
 * ExpCM P Matrix (Combined Rate Matrix)
 *
 * Final substitution rate matrix: P = F ⊙ Q (element-wise multiplication).
 * Diagonals set so that row sums = 0.
 */
class ExpCMPMatrix {
public:
    Matrix61Codon P = Matrix61Codon::Zero();

    void build(const Matrix61Codon& F, const Matrix61Codon& Q) {
        P = F.cwiseProduct(Q);
        for (int x = 0; x < 61; ++x) {
            P(x, x) = -P.row(x).sum();
        }
    }
};

/**
 * Reversible spectral decomposition of the site-specific ExpCM rate matrix P_r.
 *
 * Because P_r = F_r ⊙ Q satisfies detailed balance with respect to
 * p_{r,x} ∝ (π_{r,A(x)})^β φ_{x1} φ_{x2} φ_{x3}:
 *   S_r = diag(p_r^{1/2}) P_r diag(p_r^{-1/2})
 * is real symmetric, enabling a single O(61^3) SelfAdjointEigenSolver per site
 * and fast O(61^2) evaluation of P_{r,b}(t_b) = U_pi exp(t_b Λ_r) V_pi across all branches.
 */
struct ExpCMSiteEigen {
    Vector61Codon pi = Vector61Codon::Zero();
    Vector61Codon sqrt_pi = Vector61Codon::Zero();
    Vector61Codon inv_sqrt_pi = Vector61Codon::Zero();
    Vector61Codon eigenvalues = Vector61Codon::Zero();
    Matrix61Codon U_pi = Matrix61Codon::Zero(); // diag(p_r^{-1/2}) V
    Matrix61Codon V_pi = Matrix61Codon::Zero(); // V^T diag(p_r^{1/2})

    void decompose(const Matrix61Codon& P, const Vector61Codon& pi_site) {
        pi = pi_site;
        for (int i = 0; i < 61; ++i) {
            Scalar p = std::max(pi(i), Scalar(1e-30));
            sqrt_pi(i) = std::sqrt(p);
            inv_sqrt_pi(i) = 1.0 / sqrt_pi(i);
        }

        Matrix61Codon S = sqrt_pi.asDiagonal() * P * inv_sqrt_pi.asDiagonal();
        S = 0.5 * (S + S.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix61Codon> solver(S);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("ExpCM site eigendecomposition failed");
        }

        eigenvalues = solver.eigenvalues().cwiseMin(0.0);
        const Matrix61Codon& V = solver.eigenvectors();
        U_pi = inv_sqrt_pi.asDiagonal() * V;
        V_pi = V.transpose() * sqrt_pi.asDiagonal();
    }

    Matrix61Codon transition_matrix(Scalar branch_length) const {
        Vector61Codon exp_lambda = (eigenvalues * branch_length).array().exp();
        Matrix61Codon P_t = (U_pi * exp_lambda.asDiagonal()) * V_pi;
        return P_t.cwiseMax(0.0);
    }
};

} // namespace hyphy::core
