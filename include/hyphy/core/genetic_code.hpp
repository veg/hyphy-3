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
        int8_t pos = -1;       // 0, 1, or 2 (backwards-compatible for 1-step)
        int8_t nuc_from = -1;  // 0..3
        int8_t nuc_to = -1;    // 0..3
        bool synonymous = false;
        int8_t num_diffs = 0;  // 0, 1, 2, or 3
        int8_t diff_pos[3] = {-1, -1, -1};
        int8_t diff_from[3] = {-1, -1, -1};
        int8_t diff_to[3] = {-1, -1, -1};
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

    uint64_t resolve_triplet_mask(char c1, char c2, char c3) const {
        uint8_t m1 = char_to_iupac_mask(c1);
        uint8_t m2 = char_to_iupac_mask(c2);
        uint8_t m3 = char_to_iupac_mask(c3);
        uint64_t all_sense = (num_sense_codons >= 64) ? ~uint64_t(0) : ((uint64_t(1) << num_sense_codons) - 1);
        if (m1 == 0xF && m2 == 0xF && m3 == 0xF) {
            return all_sense;
        }
        uint64_t mask = 0;
        for (int i = 0; i < 4; ++i) {
            if (!(m1 & (1 << i))) continue;
            for (int j = 0; j < 4; ++j) {
                if (!(m2 & (1 << j))) continue;
                for (int k = 0; k < 4; ++k) {
                    if (!(m3 & (1 << k))) continue;
                    int8_t s = codon_lookup[i][j][k];
                    if (s >= 0) {
                        mask |= (uint64_t(1) << s);
                    }
                }
            }
        }
        return (mask != 0) ? mask : all_sense;
    }

    std::vector<int> resolve_ambiguous_triplet(char c1, char c2, char c3) const {
        uint64_t mask = resolve_triplet_mask(c1, c2, c3);
        std::vector<int> res;
        for (int s = 0; s < num_sense_codons; ++s) {
            if (mask & (uint64_t(1) << s)) {
                res.push_back(s);
            }
        }
        return res;
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
                int8_t diffs = 0;
                for (int p = 0; p < 3; ++p) {
                    if (s1[p] != s2[p]) {
                        diff_matrix[u][v].diff_pos[diffs] = static_cast<int8_t>(p);
                        diff_matrix[u][v].diff_from[diffs] = char_to_nuc(s1[p]);
                        diff_matrix[u][v].diff_to[diffs] = char_to_nuc(s2[p]);
                        diffs++;
                    }
                }
                diff_matrix[u][v].num_diffs = diffs;
                diff_matrix[u][v].synonymous = (sense_codon_aa[u] == sense_codon_aa[v]);
                if (diffs == 1) {
                    diff_matrix[u][v].pos = diff_matrix[u][v].diff_pos[0];
                    diff_matrix[u][v].nuc_from = diff_matrix[u][v].diff_from[0];
                    diff_matrix[u][v].nuc_to = diff_matrix[u][v].diff_to[0];
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

        // 13. Chlorophycean-mtDNA
        register_code("Chlorophycean-mtDNA", 12, {
            {50, 1}    // TAG => Leu
        }, {"chlorophycean_mtdna", "chlorophycean", "transl_table=16", "16"});

        // 14. Trematode-mtDNA
        register_code("Trematode-mtDNA", 13, {
            {0, 13},   // AAA => Asn
            {8, 5},    // AGA => Ser
            {10, 5},   // AGG => Ser
            {12, 3},   // ATA => Met
            {56, 18}   // TGA => Trp
        }, {"trematode_mtdna", "trematode", "transl_table=21", "21"});

        // 15. Scenedesmus-obliquus-mtDNA
        register_code("Scenedesmus-obliquus-mtDNA", 14, {
            {50, 1},   // TAG => Leu
            {52, 10}   // TCA => Stop
        }, {"scenedesmus_obliquus_mtdna", "scenedesmus-mtdna", "scenedesmus_mtdna", "transl_table=22", "22"});

        // 16. Thraustochytrium-mtDNA
        register_code("Thraustochytrium-mtDNA", 15, {
            {60, 10}   // TTA => Stop
        }, {"thraustochytrium_mtdna", "thraustochytrium", "transl_table=23", "23"});

        // 17. Pterobranchia-mtDNA
        register_code("Pterobranchia-mtDNA", 16, {
            {8, 5},    // AGA => Ser
            {10, 14},  // AGG => Lys
            {56, 18}   // TGA => Trp
        }, {"pterobranchia_mtdna", "pterobranchia", "transl_table=24", "24"});

        // 18. SR1-and-Gracilibacteria
        register_code("SR1-and-Gracilibacteria", 17, {
            {56, 20}   // TGA => Gly
        }, {"sr1_and_gracilibacteria", "sr1", "gracilibacteria", "transl_table=25", "25"});

        // 19. Pachysolen-Nuclear
        register_code("Pachysolen-Nuclear", 18, {
            {30, 8}    // CTG => Ala
        }, {"pachysolen_nuclear", "pachysolen", "transl_table=26", "26"});

        // 20. Mesodinium-Nuclear
        register_code("Mesodinium-Nuclear", 19, {
            {48, 9},   // TAA => Tyr
            {50, 9}    // TAG => Tyr
        }, {"mesodinium_nuclear", "mesodinium", "transl_table=29", "29"});

        // 21. Peritrich-Nuclear
        register_code("Peritrich-Nuclear", 20, {
            {48, 16},  // TAA => Glu
            {50, 16}   // TAG => Glu
        }, {"peritrich_nuclear", "peritrich", "transl_table=30", "30"});

        // 22. Cephalodiscidae-mtDNA
        register_code("Cephalodiscidae-mtDNA", 21, {
            {8, 5},    // AGA => Ser
            {10, 14},  // AGG => Lys
            {48, 9},   // TAA => Tyr
            {56, 18}   // TGA => Trp
        }, {"cephalodiscidae_mtdna", "cephalodiscidae", "transl_table=33", "33"});
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
