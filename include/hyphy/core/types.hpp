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

// Universal Genetic Code Sense Codons in standard HyPhy alphabetical order (0..60)
inline const std::array<std::string_view, 61> SENSE_CODONS = {
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
    "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
    "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
    "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG",
    "TGT", "TTA", "TTC", "TTG", "TTT"
};

// Amino Acid translation for each of the 61 sense codons in Universal code
inline const std::array<char, 61> SENSE_CODON_AA = {
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
    'Y', 'Y', 'S', 'S', 'S', 'S', 'C', 'W',
    'C', 'L', 'F', 'L', 'F'
};

// Fast lookup from 3-nucleotide indices (0..3 each) to 61 sense codon index (-1 for stops)
struct GeneticCodeTable {
    int8_t codon_lookup[4][4][4];
    char   codon_aa[4][4][4];
    bool   is_stop[4][4][4];

    GeneticCodeTable() {
        // Initialize all to stop / invalid
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    codon_lookup[i][j][k] = -1;
                    codon_aa[i][j][k] = '*';
                    is_stop[i][j][k] = false;
                }
            }
        }
        // Universal stops: TAA (3,0,0), TAG (3,0,2), TGA (3,2,0)
        is_stop[3][0][0] = true;
        is_stop[3][0][2] = true;
        is_stop[3][2][0] = true;

        for (size_t idx = 0; idx < SENSE_CODONS.size(); ++idx) {
            auto s = SENSE_CODONS[idx];
            int8_t n1 = char_to_nuc(s[0]);
            int8_t n2 = char_to_nuc(s[1]);
            int8_t n3 = char_to_nuc(s[2]);
            codon_lookup[n1][n2][n3] = static_cast<int8_t>(idx);
            codon_aa[n1][n2][n3] = SENSE_CODON_AA[idx];
        }
    }

    static const GeneticCodeTable& universal() {
        static GeneticCodeTable table;
        return table;
    }
};

inline int8_t triplet_to_sense_codon(char c1, char c2, char c3) {
    int8_t n1 = char_to_nuc(c1);
    int8_t n2 = char_to_nuc(c2);
    int8_t n3 = char_to_nuc(c3);
    if (n1 < 0 || n2 < 0 || n3 < 0) return -1;
    return GeneticCodeTable::universal().codon_lookup[n1][n2][n3];
}

inline bool is_synonymous(int8_t codon_from, int8_t codon_to) {
    if (codon_from < 0 || codon_from >= 61 || codon_to < 0 || codon_to >= 61) return false;
    return SENSE_CODON_AA[codon_from] == SENSE_CODON_AA[codon_to];
}

inline int count_nucleotide_differences(int8_t codon_from, int8_t codon_to) {
    if (codon_from == codon_to) return 0;
    auto s1 = SENSE_CODONS[codon_from];
    auto s2 = SENSE_CODONS[codon_to];
    int diffs = 0;
    if (s1[0] != s2[0]) ++diffs;
    if (s1[1] != s2[1]) ++diffs;
    if (s1[2] != s2[2]) ++diffs;
    return diffs;
}

} // namespace hyphy::core
