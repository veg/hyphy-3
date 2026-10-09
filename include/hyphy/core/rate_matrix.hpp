#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
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

    Scalar scale_factor = 1.0;

    void update(const GTRParameters& params, const Vector4& frequencies) {
        pi = frequencies;
        for (int i = 0; i < 4; ++i) {
            sqrt_pi(i) = std::sqrt(pi(i));
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

        eigenvalues = solver.eigenvalues();
        eigenvectors = solver.eigenvectors(); // V
        inv_eigenvectors = eigenvectors.transpose(); // V^T
    }

    Matrix4 transition_matrix(Scalar branch_length) const {
        Vector4 exp_lambda = (eigenvalues * branch_length).array().exp();
        Matrix4 P = inv_sqrt_pi.asDiagonal() * eigenvectors * exp_lambda.asDiagonal() * inv_eigenvectors * sqrt_pi.asDiagonal();
        return P;
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

    Scalar scale_factor = 1.0;

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

        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                const auto& diff = code.diff_matrix[i][j];
                if (diff.pos >= 0) {
                    Scalar nuc_rate = r_nuc[diff.nuc_from][diff.nuc_to];
                    Scalar rate_modifier = diff.synonymous ? params.alpha : params.beta;
                    Scalar target_nuc_freq = pos_nuc_freqs(diff.pos, diff.nuc_to);
                    Q(i, j) = rate_modifier * nuc_rate * target_nuc_freq;
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
    }

    Matrix transition_matrix(Scalar branch_length) const {
        Vector exp_lambda = (eigenvalues * branch_length).array().exp();
        Matrix P = inv_sqrt_pi.asDiagonal() * eigenvectors * exp_lambda.asDiagonal() * inv_eigenvectors * sqrt_pi.asDiagonal();
        return P;
    }
};

} // namespace hyphy::core
