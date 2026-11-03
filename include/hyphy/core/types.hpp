#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <span>
#include <array>
#include <cmath>
#include <memory>
#include <limits>
#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace hyphy::core {

using Scalar = double;
using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

using Vector4 = Eigen::Matrix<Scalar, 4, 1>;
using Matrix4 = Eigen::Matrix<Scalar, 4, 4>;

using Vector61 = Vector;
using Matrix61 = Matrix;

constexpr int32_t INVALID_INDEX = -1;
constexpr size_t NUM_NUCLEOTIDES = 4;
constexpr size_t NUM_SENSE_CODONS = 61;

// Standard Nucleotide Indices: A=0, C=1, G=2, T=3
inline constexpr int8_t char_to_nuc(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1; // gap or ambiguous
    }
}

inline constexpr char nuc_to_char(int8_t idx) {
    constexpr char nucs[] = {'A', 'C', 'G', 'T'};
    return (idx >= 0 && idx < 4) ? nucs[idx] : '-';
}

// Maps standard nucleotide and IUPAC ambiguity characters to a 4-bit mask:
// bit 0 (1): A, bit 1 (2): C, bit 2 (4): G, bit 3 (8): T/U, 0xF: N/-/?
inline constexpr uint8_t char_to_iupac_mask(char c) {
    switch (c) {
        case 'A': case 'a': return 0b0001; // A
        case 'C': case 'c': return 0b0010; // C
        case 'G': case 'g': return 0b0100; // G
        case 'T': case 't': case 'U': case 'u': return 0b1000; // T/U
        case 'R': case 'r': return 0b0101; // A or G
        case 'Y': case 'y': return 0b1010; // C or T
        case 'S': case 's': return 0b0110; // C or G
        case 'W': case 'w': return 0b1001; // A or T
        case 'K': case 'k': return 0b1100; // G or T
        case 'M': case 'm': return 0b0011; // A or C
        case 'B': case 'b': return 0b1110; // C, G, or T (not A)
        case 'D': case 'd': return 0b1101; // A, G, or T (not C)
        case 'H': case 'h': return 0b1011; // A, C, or T (not G)
        case 'V': case 'v': return 0b0111; // A, C, or G (not T)
        default:            return 0b1111; // N, -, ?, or unknown
    }
}

} // namespace hyphy::core
