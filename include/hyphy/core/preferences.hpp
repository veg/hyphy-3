#pragma once

#include "hyphy/core/types.hpp"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace hyphy::core {

/**
 * Preferences: Site-specific amino-acid preferences from deep mutational scanning.
 *
 * For ExpCM (Experimentally informed Codon Model), each site has a preference
 * vector over the 20 standard amino acids, derived from experimental data
 * (e.g., deep mutational scanning assays).
 *
 * Format: CSV with columns [site, A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y]
 *   - First column: site number (1-indexed)
 *   - Next 20 columns: preference for each amino acid
 *   - Each row (except header) sums to 1.0
 *
 * During loading:
 *   - Preferences < minpref (default 0.002) are raised to minpref
 *   - Array is re-normalized to sum to 1
 *
 * Storage: preferences[r][a] = preference of site r for amino acid a
 *   - r: 0-indexed site (0 to nsites-1)
 *   - a: 0-indexed amino acid (0 to 19, in order: A C D E F G H I K L M N P Q R S T V W Y)
 */
class Preferences {
public:
    static constexpr size_t N_AA = 20;  // Number of standard amino acids
    static constexpr Scalar MIN_PREF_DEFAULT = 0.002;
    static constexpr Scalar ALMOST_ZERO = 1e-8;

    // Amino acid letter to index (0-19)
    static const std::unordered_map<char, size_t> AA_TO_INDEX;
    static const std::vector<char> INDEX_TO_AA;

    size_t nsites = 0;
    // preferences[r][a] = preference of site r for amino acid a
    Eigen::Matrix<Scalar, Eigen::Dynamic, 20> prefs;

    /**
     * Load preferences from CSV file.
     *
     * Expected format:
     *   site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
     *   1,0.05,0.05,...,0.05
     *   2,0.10,0.05,...,0.05
     *   ...
     *
     * @param filepath  Path to CSV file
     * @param minpref   Minimum preference value (default 0.002)
     * @return          Preferences object
     */
    static Preferences from_csv(const std::string& filepath, Scalar minpref = MIN_PREF_DEFAULT) {
        Preferences prefs;

        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open preferences file: " + filepath);
        }

        std::string line;
        size_t line_num = 0;
        std::vector<std::vector<Scalar>> data;  // Temporary storage

        while (std::getline(file, line)) {
            line_num++;

            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);

            if (line.empty()) continue;
            if (line[0] == '#') continue;  // Skip comments

            // Parse CSV line
            std::vector<Scalar> row;
            std::stringstream ss(line);
            std::string field;

            while (std::getline(ss, field, ',')) {
                // Trim field
                field.erase(0, field.find_first_not_of(" \t"));
                field.erase(field.find_last_not_of(" \t") + 1);

                // First field is site number (skip validation, just parse)
                if (row.empty()) {
                    try {
                        std::stoi(field);  // Just check it parses as integer
                    } catch (...) {
                        // Header line or malformed, skip
                        if (line_num == 1) {
                            // Likely header
                            continue;
                        } else {
                            throw std::runtime_error("Line " + std::to_string(line_num) +
                                                   ": first field is not a site number: " + field);
                        }
                    }
                    row.push_back(0.0);  // Placeholder, not used
                } else {
                    // Parse preference value
                    try {
                        Scalar val = std::stod(field);
                        if (val < 0.0 || val > 1.0) {
                            throw std::runtime_error("Line " + std::to_string(line_num) +
                                                   ": preference out of range [0, 1]: " + field);
                        }
                        row.push_back(val);
                    } catch (const std::exception& e) {
                        throw std::runtime_error("Line " + std::to_string(line_num) +
                                               ": could not parse preference: " + field +
                                               " (" + e.what() + ")");
                    }
                }
            }

            // Validate: should have 21 fields (site + 20 AAs)
            if (row.size() != N_AA + 1) {
                if (line_num == 1) {
                    // Likely header, skip
                    continue;
                }
                throw std::runtime_error("Line " + std::to_string(line_num) + ": expected " +
                                       std::to_string(N_AA + 1) + " fields, got " +
                                       std::to_string(row.size()));
            }

            // Remove site number column (index 0), keep only preferences
            row.erase(row.begin());
            data.push_back(row);
        }

        if (data.empty()) {
            throw std::runtime_error("No preference data found in file: " + filepath);
        }

        // Now process: apply minpref floor and normalize
        prefs.nsites = data.size();
        prefs.prefs = Eigen::Matrix<Scalar, Eigen::Dynamic, 20>(prefs.nsites, N_AA);

        for (size_t r = 0; r < prefs.nsites; ++r) {
            // Apply minpref floor
            for (size_t a = 0; a < N_AA; ++a) {
                if (data[r][a] < minpref) {
                    data[r][a] = minpref;
                }
            }

            // Normalize to sum to 1
            Scalar sum = 0.0;
            for (size_t a = 0; a < N_AA; ++a) {
                sum += data[r][a];
            }

            if (sum < ALMOST_ZERO) {
                throw std::runtime_error("Site " + std::to_string(r + 1) +
                                       ": sum of preferences is near zero after minpref floor");
            }

            for (size_t a = 0; a < N_AA; ++a) {
                prefs.prefs(r, a) = data[r][a] / sum;
            }

            // Validate: sum should be ~1
            Scalar normalized_sum = 0.0;
            for (size_t a = 0; a < N_AA; ++a) {
                normalized_sum += prefs.prefs(r, a);
            }
            if (std::fabs(normalized_sum - 1.0) > 1e-6) {
                throw std::runtime_error("Site " + std::to_string(r + 1) +
                                       ": normalized sum != 1.0 (got " +
                                       std::to_string(normalized_sum) + ")");
            }
        }

        return prefs;
    }

    /**
     * Get preference value for site r, amino acid a (by index).
     *
     * @param r  Site index (0-based)
     * @param a  Amino acid index (0-19)
     * @return   Preference value (in [0, 1])
     */
    Scalar at(size_t r, size_t a) const {
        if (r >= nsites || a >= N_AA) {
            throw std::out_of_range("Preferences::at: index out of range");
        }
        return prefs(r, a);
    }

    /**
     * Get preference value for site r, amino acid letter.
     *
     * @param r   Site index (0-based)
     * @param aa  Amino acid letter (A-Y)
     * @return    Preference value (in [0, 1])
     */
    Scalar at(size_t r, char aa) const {
        auto it = AA_TO_INDEX.find(aa);
        if (it == AA_TO_INDEX.end()) {
            throw std::invalid_argument(std::string("Unknown amino acid: ") + aa);
        }
        return at(r, it->second);
    }

    /**
     * Get row of preferences for site r.
     *
     * @param r  Site index (0-based)
     * @return   Vector of 20 preference values
     */
    Eigen::Matrix<Scalar, 1, 20> row(size_t r) const {
        if (r >= nsites) {
            throw std::out_of_range("Preferences::row: site index out of range");
        }
        return prefs.row(r);
    }

    /**
     * Get number of sites.
     */
    size_t num_sites() const { return nsites; }

    /**
     * Validate consistency with an alignment (number of sites must match).
     *
     * @param num_alignment_sites  Number of codon sites in alignment
     * @return true if valid, false otherwise
     */
    bool validate_against_alignment(size_t num_alignment_sites) const {
        if (nsites != num_alignment_sites) {
            std::cerr << "Preferences::validate_against_alignment: "
                      << "nsites (" << nsites << ") != "
                      << "num_alignment_sites (" << num_alignment_sites << ")"
                      << std::endl;
            return false;
        }
        return true;
    }
};

// Static member definitions
const std::unordered_map<char, size_t> Preferences::AA_TO_INDEX = {
    {'A', 0},  {'C', 1},  {'D', 2},  {'E', 3},  {'F', 4},
    {'G', 5},  {'H', 6},  {'I', 7},  {'K', 8},  {'L', 9},
    {'M', 10}, {'N', 11}, {'P', 12}, {'Q', 13}, {'R', 14},
    {'S', 15}, {'T', 16}, {'V', 17}, {'W', 18}, {'Y', 19}
};

const std::vector<char> Preferences::INDEX_TO_AA = {
    'A', 'C', 'D', 'E', 'F',
    'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R',
    'S', 'T', 'V', 'W', 'Y'
};

}  // namespace hyphy::core
