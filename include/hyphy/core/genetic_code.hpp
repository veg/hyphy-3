#pragma once

#include "hyphy/core/types.hpp"
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include <stdexcept>

namespace hyphy::core {

class GeneticCode {
public:
    std::string name;
    int id = 0;
    int num_sense_codons = 61;

    std::vector<std::string> sense_codons;
    std::vector<char> sense_codon_aa;
    std::vector<std::string> stop_codons;

    int8_t codon_lookup[4][4][4];
    char codon_aa[4][4][4];
    bool is_stop[4][4][4];

    int8_t codon64_to_sense[64];
    int8_t sense_to_codon64[64];

    // Precomputed substitution graph properties between sense codons
    struct StepDiff {
        int8_t pos = -1;       // 0, 1, or 2
        int8_t nuc_from = -1;  // 0..3
        int8_t nuc_to = -1;    // 0..3
        bool synonymous = false;
    };
    std::vector<std::vector<StepDiff>> diff_matrix; // [S][S]

    GeneticCode(std::string code_name, int code_id, const std::array<int8_t, 64>& translation)
        : name(std::move(code_name)), id(code_id) {
        init(translation);
    }

    int8_t triplet_to_sense(char c1, char c2, char c3) const {
        int8_t n1 = char_to_nuc(c1);
        int8_t n2 = char_to_nuc(c2);
        int8_t n3 = char_to_nuc(c3);
        if (n1 < 0 || n2 < 0 || n3 < 0) return -1;
        return codon_lookup[n1][n2][n3];
    }

    bool is_stop_triplet(char c1, char c2, char c3) const {
        int8_t n1 = char_to_nuc(c1);
        int8_t n2 = char_to_nuc(c2);
        int8_t n3 = char_to_nuc(c3);
        if (n1 < 0 || n2 < 0 || n3 < 0) return false;
        return is_stop[n1][n2][n3];
    }

    bool is_synonymous(int8_t sense_from, int8_t sense_to) const {
        if (sense_from < 0 || sense_from >= num_sense_codons ||
            sense_to < 0 || sense_to >= num_sense_codons) return false;
        return sense_codon_aa[sense_from] == sense_codon_aa[sense_to];
    }

    int count_nucleotide_differences(int8_t sense_from, int8_t sense_to) const {
        if (sense_from == sense_to) return 0;
        const auto& s1 = sense_codons[sense_from];
        const auto& s2 = sense_codons[sense_to];
        int d = 0;
        if (s1[0] != s2[0]) ++d;
        if (s1[1] != s2[1]) ++d;
        if (s1[2] != s2[2]) ++d;
        return d;
    }

    static std::shared_ptr<const GeneticCode> from_name(std::string name);
    static std::shared_ptr<const GeneticCode> universal();
    static std::shared_ptr<const GeneticCode> vertebrate_mtdna();

private:
    void init(const std::array<int8_t, 64>& translation) {
        const std::string hyphy_aa = "FLIMVSPTAYXHQNKDECWRG";
        const char nucs[] = {'A', 'C', 'G', 'T'};

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    codon_lookup[i][j][k] = -1;
                    codon_aa[i][j][k] = '*';
                    is_stop[i][j][k] = false;
                }
            }
        }
        std::fill(std::begin(codon64_to_sense), std::end(codon64_to_sense), static_cast<int8_t>(-1));
        std::fill(std::begin(sense_to_codon64), std::end(sense_to_codon64), static_cast<int8_t>(-1));

        sense_codons.clear();
        sense_codon_aa.clear();
        stop_codons.clear();

        int8_t s_idx = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    int c64 = 16 * i + 4 * j + k;
                    int8_t aa_idx = translation[c64];
                    std::string triplet = {nucs[i], nucs[j], nucs[k]};

                    if (aa_idx == 10) { // Stop codon (X)
                        is_stop[i][j][k] = true;
                        stop_codons.push_back(triplet);
                    } else {
                        char aa = hyphy_aa[aa_idx];
                        codon_lookup[i][j][k] = s_idx;
                        codon_aa[i][j][k] = aa;
                        codon64_to_sense[c64] = s_idx;
                        sense_to_codon64[s_idx] = static_cast<int8_t>(c64);

                        sense_codons.push_back(triplet);
                        sense_codon_aa.push_back(aa);
                        s_idx++;
                    }
                }
            }
        }

        num_sense_codons = s_idx;

        // Precompute diff matrix
        diff_matrix.assign(num_sense_codons, std::vector<StepDiff>(num_sense_codons));
        for (int u = 0; u < num_sense_codons; ++u) {
            for (int v = 0; v < num_sense_codons; ++v) {
                if (u == v) continue;
                const auto& s1 = sense_codons[u];
                const auto& s2 = sense_codons[v];
                int diffs = 0;
                int diff_p = -1;
                for (int p = 0; p < 3; ++p) {
                    if (s1[p] != s2[p]) {
                        diffs++;
                        diff_p = p;
                    }
                }
                if (diffs == 1) {
                    diff_matrix[u][v].pos = diff_p;
                    diff_matrix[u][v].nuc_from = char_to_nuc(s1[diff_p]);
                    diff_matrix[u][v].nuc_to = char_to_nuc(s2[diff_p]);
                    diff_matrix[u][v].synonymous = (sense_codon_aa[u] == sense_codon_aa[v]);
                }
            }
        }
    }
};

class GeneticCodeRegistry {
public:
    static GeneticCodeRegistry& instance() {
        static GeneticCodeRegistry reg;
        return reg;
    }

    std::shared_ptr<const GeneticCode> get(const std::string& name) const {
        std::string lower_name = name;
        std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
        auto it = code_map.find(lower_name);
        if (it != code_map.end()) {
            return it->second;
        }
        throw std::runtime_error("Unknown genetic code: '" + name + "'");
    }

private:
    std::unordered_map<std::string, std::shared_ptr<const GeneticCode>> code_map;

    GeneticCodeRegistry() {
        const std::array<int8_t, 64> universal_table = {
            14, 13, 14, 13,  7,  7,  7,  7, 19,  5, 19,  5,  2,  2,  3,  2,
            12, 11, 12, 11,  6,  6,  6,  6, 19, 19, 19, 19,  1,  1,  1,  1,
            16, 15, 16, 15,  8,  8,  8,  8, 20, 20, 20, 20,  4,  4,  4,  4,
            10,  9, 10,  9,  5,  5,  5,  5, 10, 17, 18, 17,  1,  0,  1,  0
        };

        auto register_code = [&](const std::string& code_name, int code_id, const std::vector<std::pair<int, int8_t>>& overrides, const std::vector<std::string>& aliases) {
            auto tbl = universal_table;
            for (const auto& [idx, aa] : overrides) {
                tbl[idx] = aa;
            }
            auto gc = std::make_shared<GeneticCode>(code_name, code_id, tbl);
            std::string lname = code_name;
            std::transform(lname.begin(), lname.end(), lname.begin(), ::tolower);
            code_map[lname] = gc;
            for (const auto& alias : aliases) {
                std::string lalias = alias;
                std::transform(lalias.begin(), lalias.end(), lalias.begin(), ::tolower);
                code_map[lalias] = gc;
            }
        };

        // 1. Universal
        register_code("Universal", 0, {}, {"standard", "default", "transl_table=1", "1"});

        // 2. Vertebrate-mtDNA
        register_code("Vertebrate-mtDNA", 1, {
            {8, 10},   // AGA => stop
            {10, 10},  // AGG => stop
            {12, 3},   // ATA => Met
            {56, 18}   // TGA => Trp
        }, {"vertebrate_mtdna", "vertebrate-mtdna", "transl_table=2", "2"});

        // 3. Yeast-mtDNA
        register_code("Yeast-mtDNA", 2, {
            {12, 3}, {28, 7}, {29, 7}, {30, 7}, {31, 7}, {56, 18}
        }, {"yeast_mtdna", "transl_table=3", "3"});

        // 4. Mold-Protozoan-mtDNA
        register_code("Mold-Protozoan-mtDNA", 3, {
            {56, 18}
        }, {"mold_protozoan_mtdna", "mold-mtdna", "transl_table=4", "4"});

        // 5. Invertebrate-mtDNA
        register_code("Invertebrate-mtDNA", 4, {
            {8, 5}, {10, 5}, {12, 3}, {56, 18}
        }, {"invertebrate_mtdna", "transl_table=5", "5"});

        // 6. Ciliate-Nuclear
        register_code("Ciliate-Nuclear", 5, {
            {48, 12}, {50, 12}
        }, {"ciliate_nuclear", "ciliate", "transl_table=6", "6"});

        // 7. Echinoderm-mtDNA
        register_code("Echinoderm-mtDNA", 6, {
            {0, 13}, {8, 5}, {10, 5}, {56, 18}
        }, {"echinoderm_mtdna", "transl_table=9", "9"});

        // 8. Euplotid-Nuclear
        register_code("Euplotid-Nuclear", 7, {
            {56, 17}
        }, {"euplotid_nuclear", "euplotid", "transl_table=10", "10"});

        // 9. Alt-Yeast-Nuclear
        register_code("Alt-Yeast-Nuclear", 8, {
            {30, 5}
        }, {"alt_yeast_nuclear", "transl_table=12", "12"});

        // 10. Ascidian-mtDNA
        register_code("Ascidian-mtDNA", 9, {
            {8, 20}, {10, 20}, {12, 3}, {56, 18}
        }, {"ascidian_mtdna", "transl_table=13", "13"});

        // 11. Flatworm-mtDNA
        register_code("Flatworm-mtDNA", 10, {
            {0, 13}, {8, 5}, {10, 5}, {48, 9}, {56, 18}
        }, {"flatworm_mtdna", "transl_table=14", "14"});

        // 12. Blepharisma-Nuclear
        register_code("Blepharisma-Nuclear", 11, {
            {50, 12}
        }, {"blepharisma_nuclear", "transl_table=15", "15"});
    }
};

inline std::shared_ptr<const GeneticCode> GeneticCode::from_name(std::string name) {
    return GeneticCodeRegistry::instance().get(name);
}

inline std::shared_ptr<const GeneticCode> GeneticCode::universal() {
    return GeneticCodeRegistry::instance().get("Universal");
}

inline std::shared_ptr<const GeneticCode> GeneticCode::vertebrate_mtdna() {
    return GeneticCodeRegistry::instance().get("Vertebrate-mtDNA");
}

} // namespace hyphy::core
