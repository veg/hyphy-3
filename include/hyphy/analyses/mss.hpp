#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/stats.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "nlohmann/json.hpp"

#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

enum class MSSNormalization {
    NeutralRef,  // Explicit neutral reference class fixed to 1.0
    Flux,        // Equilibrium synonymous flux-weighted mean = 1.0
    Unweighted   // Unweighted arithmetic mean across classes = 1.0 (matches HyPhy 2.5 SynREV/SynREVCodon)
};

struct MSSSynPair {
    int u = -1;                 // Sense codon index 1 (u < v)
    int v = -1;                 // Sense codon index 2
    std::string codon1;         // 3-letter codon 1 (lexicographically smaller)
    std::string codon2;         // 3-letter codon 2
    char aa = '*';              // Single-letter amino acid
    int pos = -1;               // Mutated codon position (0, 1, or 2)
    int nuc_from = -1;          // 0..3
    int nuc_to = -1;            // 0..3
    bool is_transition = false; // A<->G or C<->T
    int class_id = 0;           // Assigned MSS rate class index (0 .. K-1)
};

struct MSSPartition {
    std::string name = "SynREV";
    std::vector<MSSSynPair> pairs;             // All 68 single-step synonymous pairs (for Universal code)
    std::vector<std::string> class_names;      // Size K
    int neutral_class_id = -1;                 // >= 0 if a neutral reference class is designated
    bool is_empirical_fixed = false;           // True for Empirical TSV with fixed pair rates
    std::vector<Scalar> empirical_pair_rates;  // Size = pairs.size() if is_empirical_fixed
    std::vector<int> pair_lookup;              // S * S flat lookup -> index in `pairs` (-1 if not 1-step syn)

    size_t num_classes() const { return class_names.size(); }
    size_t num_pairs() const { return pairs.size(); }

    // Enumerate all unordered single-step synonymous sense codon pairs (u < v)
    static std::vector<MSSSynPair> enumerate_synonymous_pairs(const GeneticCode& code, std::vector<int>* out_lookup = nullptr) {
        int S = code.num_sense_codons;
        if (out_lookup) {
            out_lookup->assign(S * S, -1);
        }
        // Sort sense codons lexicographically so pair ordering is canonical
        std::vector<int> order(S);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return code.sense_codons[a] < code.sense_codons[b];
        });

        std::vector<MSSSynPair> list;
        for (int ia = 0; ia < S; ++ia) {
            int u = order[ia];
            for (int ib = ia + 1; ib < S; ++ib) {
                int v = order[ib];
                const auto& diff = code.diff_matrix[u][v];
                if (diff.num_diffs == 1 && diff.synonymous) {
                    MSSSynPair sp;
                    sp.u = u;
                    sp.v = v;
                    sp.codon1 = code.sense_codons[u];
                    sp.codon2 = code.sense_codons[v];
                    sp.aa = code.sense_codon_aa[u];
                    sp.pos = diff.pos;
                    sp.nuc_from = diff.nuc_from;
                    sp.nuc_to = diff.nuc_to;
                    // Transition: A(0)<->G(2) or C(1)<->T(3)
                    sp.is_transition = ((sp.nuc_from == 0 && sp.nuc_to == 2) ||
                                        (sp.nuc_from == 2 && sp.nuc_to == 0) ||
                                        (sp.nuc_from == 1 && sp.nuc_to == 3) ||
                                        (sp.nuc_from == 3 && sp.nuc_to == 1));
                    int idx = static_cast<int>(list.size());
                    list.push_back(sp);
                    if (out_lookup) {
                        (*out_lookup)[u * S + v] = idx;
                        (*out_lookup)[v * S + u] = idx;
                    }
                }
            }
        }
        return list;
    }

    // Build flat S*S synonymous rate multiplier table from class rates (length K)
    std::vector<Scalar> expand_to_matrix(int S, const std::vector<Scalar>& class_rates) const {
        std::vector<Scalar> mat(S * S, 1.0);
        if (is_empirical_fixed) {
            for (size_t p = 0; p < pairs.size(); ++p) {
                Scalar r = (p < empirical_pair_rates.size()) ? empirical_pair_rates[p] : 1.0;
                mat[pairs[p].u * S + pairs[p].v] = r;
                mat[pairs[p].v * S + pairs[p].u] = r;
            }
            return mat;
        }
        for (const auto& sp : pairs) {
            Scalar r = (sp.class_id >= 0 && static_cast<size_t>(sp.class_id) < class_rates.size())
                           ? class_rates[sp.class_id]
                           : 1.0;
            mat[sp.u * S + sp.v] = r;
            mat[sp.v * S + sp.u] = r;
        }
        return mat;
    }

    static std::string trim(const std::string& s) {
        size_t a = 0, b = s.size();
        while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
        while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
        return s.substr(a, b - a);
    }

    static std::string to_upper(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        return s;
    }

    // Construct a preset or TSV-defined MSS partition
    static MSSPartition build(
        const std::string& spec,
        const GeneticCode& code = *GeneticCode::universal(),
        const std::string& neutral_ref_name = "",
        bool shared_between_class_rate = false
    ) {
        MSSPartition part;
        part.pairs = enumerate_synonymous_pairs(code, &part.pair_lookup);
        std::string up = to_upper(trim(spec));

        if (up == "SYNREV" || up == "AA" || up == "AMINO-ACID") {
            part.name = "SynREV";
            std::map<char, int> aa_to_class;
            for (auto& sp : part.pairs) {
                auto it = aa_to_class.find(sp.aa);
                if (it == aa_to_class.end()) {
                    int cid = static_cast<int>(part.class_names.size());
                    aa_to_class[sp.aa] = cid;
                    part.class_names.push_back(std::string("AA_") + sp.aa);
                    sp.class_id = cid;
                } else {
                    sp.class_id = it->second;
                }
            }
        } else if (up == "SYNREVCODON" || up == "CODON" || up == "PAIRS" || up == "PAIR") {
            part.name = "SynREVCodon";
            for (size_t i = 0; i < part.pairs.size(); ++i) {
                part.pairs[i].class_id = static_cast<int>(i);
                part.class_names.push_back(std::string(1, part.pairs[i].aa) + "_" + part.pairs[i].codon1 + "_" + part.pairs[i].codon2);
            }
        } else if (up == "SYNREV2" || up == "SYNREV2G") {
            bool group_between = (up == "SYNREV2G") || shared_between_class_rate;
            part.name = group_between ? "SynREV2g" : "SynREV2";
            // Partition codons within each AA into transition-connected pairs (pos 0 & 1 equal, pos 2 A<->G or C<->T)
            int S = code.num_sense_codons;
            std::unordered_map<std::string, std::string> codon_subtag;
            std::map<char, std::vector<std::string>> aa_codons;
            for (int i = 0; i < S; ++i) {
                aa_codons[code.sense_codon_aa[i]].push_back(code.sense_codons[i]);
            }
            for (auto& [aa, c_list] : aa_codons) {
                std::sort(c_list.begin(), c_list.end());
                std::vector<std::vector<std::string>> subsets;
                for (const auto& c : c_list) {
                    bool matched = false;
                    for (auto& sub : subsets) {
                        if (sub.size() == 1 && sub[0][0] == c[0] && sub[0][1] == c[1]) {
                            char l3 = sub[0][2], r3 = c[2];
                            if ((l3 == 'A' && r3 == 'G') || (l3 == 'G' && r3 == 'A') ||
                                (l3 == 'C' && r3 == 'T') || (l3 == 'T' && r3 == 'C')) {
                                sub.push_back(c);
                                matched = true;
                                break;
                            }
                        }
                    }
                    if (!matched) {
                        subsets.push_back({c});
                    }
                }
                for (const auto& sub : subsets) {
                    std::string tag(1, aa);
                    for (const auto& c : sub) tag += "_" + c;
                    for (const auto& c : sub) codon_subtag[c] = tag;
                }
            }
            std::map<std::string, int> cname_to_id;
            for (auto& sp : part.pairs) {
                std::string t1 = codon_subtag[sp.codon1];
                std::string t2 = codon_subtag[sp.codon2];
                std::string cname;
                if (t1 == t2) {
                    cname = t1;
                } else if (group_between) {
                    cname = "BETWEEN_CLASS";
                } else {
                    cname = (t1 < t2) ? (t1 + "__" + t2) : (t2 + "__" + t1);
                }
                auto it = cname_to_id.find(cname);
                if (it == cname_to_id.end()) {
                    int cid = static_cast<int>(part.class_names.size());
                    cname_to_id[cname] = cid;
                    part.class_names.push_back(cname);
                    sp.class_id = cid;
                } else {
                    sp.class_id = it->second;
                }
            }
        } else {
            // Load from TSV file (auto-detects Codon TSV, 3-col Pair TSV, 5-col processor-codon TSV, or Empirical rate TSV)
            part = load_from_tsv(spec, code, neutral_ref_name, shared_between_class_rate);
        }

        // Resolve neutral reference class if requested
        if (!neutral_ref_name.empty()) {
            std::string target = to_upper(trim(neutral_ref_name));
            for (size_t k = 0; k < part.class_names.size(); ++k) {
                std::string ck = to_upper(part.class_names[k]);
                if (ck == target || ck == "AA_" + target || ck == "WITHIN_" + target) {
                    part.neutral_class_id = static_cast<int>(k);
                    break;
                }
            }
            if (part.neutral_class_id < 0) {
                throw std::runtime_error("Requested MSS neutral reference class '" + neutral_ref_name + "' not found in partition '" + part.name + "'.");
            }
        } else if (part.neutral_class_id < 0) {
            // Auto-detect if a class is named NEUTRAL or WITHIN_NEUTRAL
            for (size_t k = 0; k < part.class_names.size(); ++k) {
                std::string ck = to_upper(part.class_names[k]);
                if (ck == "NEUTRAL" || ck == "WITHIN_NEUTRAL") {
                    part.neutral_class_id = static_cast<int>(k);
                    break;
                }
            }
        }

        return part;
    }

    static MSSPartition load_from_tsv(
        const std::string& filepath,
        const GeneticCode& code,
        [[maybe_unused]] const std::string& neutral_ref_name = "",
        bool shared_between_class_rate = false
    ) {
        std::ifstream in(filepath);
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open MSS partition file: '" + filepath + "'. Also note that deprecated HyPhy 2.5 modes ('Full', 'Random') have been removed; use 'SynREV', 'SynREVCodon', 'SynREV2', 'SynREV2g', or a valid TSV file.");
        }

        MSSPartition part;
        part.name = filepath;
        part.pairs = enumerate_synonymous_pairs(code, &part.pair_lookup);

        std::vector<std::vector<std::string>> rows;
        std::string line;
        while (std::getline(in, line)) {
            std::string tline = trim(line);
            if (tline.empty() || tline[0] == '#') continue;
            std::vector<std::string> cols;
            std::stringstream ss(tline);
            std::string cell;
            if (tline.find('\t') != std::string::npos) {
                while (std::getline(ss, cell, '\t')) cols.push_back(trim(cell));
            } else {
                while (ss >> cell) cols.push_back(trim(cell));
            }
            if (cols.size() >= 3) {
                rows.push_back(cols);
            }
        }
        if (rows.empty()) {
            throw std::runtime_error("MSS TSV file '" + filepath + "' contains no valid rows.");
        }

        auto is_codon3 = [&](const std::string& s) {
            if (s.size() != 3) return false;
            std::string u = to_upper(s);
            return code.triplet_to_sense(u[0], u[1], u[2]) >= 0;
        };
        auto is_pair6 = [&](const std::string& s) {
            if (s.size() != 6) return false;
            return is_codon3(s.substr(0, 3)) && is_codon3(s.substr(3, 3));
        };

        // Check if row 0 is a header
        size_t start_row = 0;
        const auto& r0 = rows[0];
        if (!is_codon3(r0[0]) && !is_codon3(r0[1]) && !is_pair6(r0[1])) {
            start_row = 1;
        }
        if (start_row >= rows.size()) {
            throw std::runtime_error("MSS TSV file '" + filepath + "' has header only and no data rows.");
        }

        const auto& first = rows[start_row];
        // Format detection:
        // Case A: 4+ cols where col[1] and col[2] are 3-letter codons -> (AA, CODON1, CODON2, GROUP [, SUPPORT])
        // Case B: 3 cols where col[0] and col[1] are 3-letter codons -> (CODON1, CODON2, GROUP_OR_RATE)
        // Case C: 3 cols where col[1] is 6-letter codon pair -> (AA, CODON1CODON2, GROUP)
        // Case D: 3 cols where col[1] is 3-letter codon -> (AA, CODON, GROUP)

        if (first.size() >= 4 && is_codon3(first[1]) && is_codon3(first[2])) {
            // Case A: processor-codon.bf 5-col TSV (AA, CODON1, CODON2, GROUP, SUPPORT)
            std::unordered_map<std::string, std::string> pair_to_group;
            for (size_t r = start_row; r < rows.size(); ++r) {
                if (rows[r].size() < 4) continue;
                std::string c1 = to_upper(rows[r][1]);
                std::string c2 = to_upper(rows[r][2]);
                if (c1 > c2) std::swap(c1, c2);
                pair_to_group[c1 + "|" + c2] = rows[r][3];
            }
            assign_pair_classes(part, pair_to_group);
        } else if (is_codon3(first[0]) && is_codon3(first[1])) {
            // Case B: (CODON1, CODON2, CLASS_OR_EMPIRICAL_RATE)
            bool all_numeric = true;
            for (size_t r = start_row; r < rows.size(); ++r) {
                char* endptr = nullptr;
                std::strtod(rows[r][2].c_str(), &endptr);
                if (!endptr || *endptr != '\0') {
                    all_numeric = false;
                    break;
                }
            }
            if (all_numeric) {
                part.is_empirical_fixed = true;
                part.class_names = {"Empirical"};
                part.empirical_pair_rates.assign(part.pairs.size(), 1.0);
                std::unordered_map<std::string, Scalar> pair_rate;
                for (size_t r = start_row; r < rows.size(); ++r) {
                    std::string c1 = to_upper(rows[r][0]);
                    std::string c2 = to_upper(rows[r][1]);
                    if (c1 > c2) std::swap(c1, c2);
                    pair_rate[c1 + "|" + c2] = std::stod(rows[r][2]);
                }
                for (size_t p = 0; p < part.pairs.size(); ++p) {
                    std::string key = part.pairs[p].codon1 + "|" + part.pairs[p].codon2;
                    auto it = pair_rate.find(key);
                    if (it == pair_rate.end()) {
                        throw std::runtime_error("Missing empirical synonymous rate for codon pair " + key + " in '" + filepath + "'.");
                    }
                    part.empirical_pair_rates[p] = it->second;
                    part.pairs[p].class_id = 0;
                }
            } else {
                std::unordered_map<std::string, std::string> pair_to_group;
                for (size_t r = start_row; r < rows.size(); ++r) {
                    std::string c1 = to_upper(rows[r][0]);
                    std::string c2 = to_upper(rows[r][1]);
                    if (c1 > c2) std::swap(c1, c2);
                    pair_to_group[c1 + "|" + c2] = rows[r][2];
                }
                assign_pair_classes(part, pair_to_group);
            }
        } else if (is_pair6(first[1])) {
            // Case C: (AA, CODON1CODON2, CLASS) like lib/pairs.tsv
            std::unordered_map<std::string, std::string> pair_to_group;
            for (size_t r = start_row; r < rows.size(); ++r) {
                std::string p6 = to_upper(rows[r][1]);
                if (p6.size() != 6) continue;
                std::string c1 = p6.substr(0, 3);
                std::string c2 = p6.substr(3, 3);
                if (c1 > c2) std::swap(c1, c2);
                pair_to_group[c1 + "|" + c2] = rows[r][2];
            }
            assign_pair_classes(part, pair_to_group);
        } else if (is_codon3(first[1])) {
            // Case D: (AA, CODON, CLASS) like codons.tsv -> induces within-class and between-class pair rates
            std::unordered_map<std::string, std::string> codon_to_group;
            for (size_t r = start_row; r < rows.size(); ++r) {
                std::string c = to_upper(rows[r][1]);
                codon_to_group[c] = rows[r][2];
            }
            std::map<std::string, int> cname_to_id;
            for (auto& sp : part.pairs) {
                auto it1 = codon_to_group.find(sp.codon1);
                auto it2 = codon_to_group.find(sp.codon2);
                if (it1 == codon_to_group.end() || it2 == codon_to_group.end()) {
                    throw std::runtime_error("Codon '" + sp.codon1 + "' or '" + sp.codon2 + "' missing from MSS codon class TSV '" + filepath + "'.");
                }
                std::string g1 = it1->second;
                std::string g2 = it2->second;
                std::string cname;
                if (g1 == g2) {
                    cname = g1;
                } else if (shared_between_class_rate) {
                    cname = "BETWEEN_CLASSES";
                } else {
                    cname = (g1 < g2) ? (g1 + "__" + g2) : (g2 + "__" + g1);
                }
                auto itc = cname_to_id.find(cname);
                if (itc == cname_to_id.end()) {
                    int cid = static_cast<int>(part.class_names.size());
                    cname_to_id[cname] = cid;
                    part.class_names.push_back(cname);
                    sp.class_id = cid;
                } else {
                    sp.class_id = itc->second;
                }
            }
        } else {
            throw std::runtime_error("Unrecognized MSS TSV format in '" + filepath + "'. Expected (AA, Codon, Class), (Codon1, Codon2, Class), (AA, Codon1Codon2, Class), or (AA, Codon1, Codon2, Class, Support).");
        }

        return part;
    }

private:
    static void assign_pair_classes(MSSPartition& part, const std::unordered_map<std::string, std::string>& pair_to_group) {
        std::map<std::string, int> cname_to_id;
        for (auto& sp : part.pairs) {
            std::string key = sp.codon1 + "|" + sp.codon2;
            auto it = pair_to_group.find(key);
            std::string cname = (it != pair_to_group.end()) ? it->second : "NEUTRAL";
            auto itc = cname_to_id.find(cname);
            if (itc == cname_to_id.end()) {
                int cid = static_cast<int>(part.class_names.size());
                cname_to_id[cname] = cid;
                part.class_names.push_back(cname);
                sp.class_id = cid;
            } else {
                sp.class_id = itc->second;
            }
        }
    }
};

struct MSSRateEstimate {
    std::string class_name;
    Scalar rate = 1.0;
    Scalar std_error = 0.0;
    Scalar ci_low = 1.0;
    Scalar ci_high = 1.0;
    Scalar lrt_vs_one = 0.0;
    Scalar p_value_uncorrected = 1.0;
    Scalar p_value_holm = 1.0;
    size_t num_pairs = 0;
    Scalar equilibrium_flux_share = 0.0;
    bool is_reference = false;
};

struct MSSGeneFitSummary {
    std::string filepath;
    size_t num_taxa = 0;
    size_t num_codons = 0;
    Scalar baseline_log_l = 0.0;
    Scalar mss_log_l = 0.0;
    Scalar baseline_omega = 1.0;
    Scalar mss_omega = 1.0;
    Scalar tree_scaler = 1.0;
    Scalar baseline_tree_length = 0.0;
    Scalar synonymous_tree_length = 0.0;
    Scalar nonsynonymous_tree_length = 0.0;
};

struct MSSPairSelectionEntry {
    char aa = '*';
    std::string codon1;
    std::string codon2;
    Scalar mle_rate = 1.0;
    Scalar curvature = 1.0;
    int assigned_class = 0;
    std::string assigned_label;
    Scalar model_averaged_rate = 1.0;
    std::vector<Scalar> class_support; // Akaike posterior support across K classes
    Scalar max_support = 1.0;
};

struct MSSSelectionResult {
    bool enabled = false;
    size_t num_classes = 2;
    std::string ic_type = "BIC";
    Scalar synrev_codon_log_l = 0.0;
    Scalar synrev_codon_ic = 0.0;
    Scalar best_partition_log_l = 0.0;
    Scalar best_partition_ic = 0.0;
    Scalar null_mg94_ic = 0.0;
    size_t models_evaluated = 0;
    std::vector<Scalar> class_rates;
    std::vector<std::string> class_labels;
    std::vector<MSSPairSelectionEntry> pair_entries;
};

struct MSSSettings {
    std::string model_type = "SynREV"; // "SynREV", "SynREVCodon", "SynREV2", "SynREV2g", or TSV path
    std::string neutral_ref = "";      // Optional neutral reference class name
    MSSNormalization normalization = MSSNormalization::Flux;
    bool normalization_explicit = false;
    bool shared_between_rate = false;  // Share a single between-class rate for Codon TSV partitions
    bool fix_omega = false;            // Fix per-gene omega at MG94 baseline (default false for 1 gene, configurable for multi-gene)
    bool compute_ci = true;            // Compute 95% profile CIs for synonymous rates
    bool compute_lrt = true;           // Perform LRTs (rate == 1) with Holm-Bonferroni correction
    bool use_f81 = false;              // Constrain nucleotide theta_ij = 1 (F81 nucleotide component)
    bool select_partition = false;     // Run fast K-class partition selector on the 68 codon pairs
    size_t select_classes = 2;         // Number of synonymous rate classes K for partition selection (or max K if auto_k)
    bool select_auto_k = false;        // Step-up K = 2..max_k by BIC/AICc
    std::string ic_type = "BIC";       // "BIC" or "AIC-c"
    std::string output_tsv = "";       // Optional path to write selected partition TSV (compatible with --mss <file>)
};

struct MSSResult {
    MSSSettings settings;
    MSSPartition partition;
    Scalar baseline_log_l = 0.0;
    Scalar baseline_aicc = 0.0;
    Scalar baseline_bic = 0.0;
    Scalar mss_log_l = 0.0;
    Scalar mss_aicc = 0.0;
    Scalar mss_bic = 0.0;
    Scalar omnibus_lrt = 0.0;
    Scalar omnibus_p_value = 1.0;
    size_t omnibus_df = 0;
    size_t total_parameters = 0;
    size_t total_sample_size = 0;
    std::vector<MSSRateEstimate> rates;
    std::vector<MSSGeneFitSummary> genes;
    MSSSelectionResult selection;
    double runtime_seconds = 0.0;
};

class MSSAnalyzer {
public:
    struct GeneContext {
        std::string filepath;
        Alignment aln;
        Tree tree;
        std::vector<size_t> leaf_to_taxon;
        GTRParameters gtr_params;
        Scalar gtr_log_l = 0.0;
        MG94Parameters mg94_params;
        Scalar mg94_log_l = 0.0;
        Scalar mg94_omega = 0.2;
        Scalar mg94_scale_factor = 1.0;
        std::vector<Scalar> branch_tau; // Unnormalized time multipliers tau_b = bl_b * (3.0 / scale_factor)
        Scalar tree_scaler = 1.0;
        Scalar omega = 0.2;
        std::vector<Scalar> pair_flux;  // Equilibrium baseline flux per synonymous pair
    };

    std::shared_ptr<const GeneticCode> code;
    std::vector<GeneContext> genes;

    explicit MSSAnalyzer(std::shared_ptr<const GeneticCode> gcode = GeneticCode::universal())
        : code(std::move(gcode)) {}

    void add_dataset(const std::string& path, Tree t, Alignment a, bool use_f81 = false) {
        GeneContext ctx;
        ctx.filepath = path;
        ctx.aln = std::move(a);
        ctx.tree = std::move(t);
        ctx.leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(ctx.tree, ctx.aln);

        // 1. Fit GTR (or F81 if requested)
        GTRFitter gtr_fitter(ctx.tree, ctx.aln);
        auto gtr_res = gtr_fitter.fit();
        ctx.gtr_log_l = gtr_res.log_likelihood;
        ctx.gtr_params = gtr_res.params;
        ctx.tree = gtr_res.tree;

        if (use_f81) {
            ctx.gtr_params.theta_AC = 1.0;
            ctx.gtr_params.theta_AT = 1.0;
            ctx.gtr_params.theta_CG = 1.0;
            ctx.gtr_params.theta_CT = 1.0;
            ctx.gtr_params.theta_GT = 1.0;
        }

        // 2. Fit baseline MG94 (global omega + nucleotide exchangeabilities + branch lengths)
        MG94Parameters mg_p;
        mg_p.alpha = 1.0;
        mg_p.beta = 0.5;
        mg_p.theta_AC = ctx.gtr_params.theta_AC;
        mg_p.theta_AT = ctx.gtr_params.theta_AT;
        mg_p.theta_CG = ctx.gtr_params.theta_CG;
        mg_p.theta_CT = ctx.gtr_params.theta_CT;
        mg_p.theta_GT = ctx.gtr_params.theta_GT;

        if (!use_f81) {
            MG94Fitter mg_fitter(ctx.tree, ctx.aln);
            auto mg_res = mg_fitter.fit_full_model(1.0, mg_p, nullptr, 35, 1e-4);
            mg_p = mg_res.params;
            ctx.tree = mg_res.tree;
            ctx.gtr_params.theta_AC = mg_p.theta_AC;
            ctx.gtr_params.theta_AT = mg_p.theta_AT;
            ctx.gtr_params.theta_CG = mg_p.theta_CG;
            ctx.gtr_params.theta_CT = mg_p.theta_CT;
            ctx.gtr_params.theta_GT = mg_p.theta_GT;
        } else {
            auto opt_mg = NelderMead2D::minimize([&](Scalar w, Scalar s_scale) {
                if (w < 1e-5 || w > 50.0 || s_scale < 0.05 || s_scale > 50.0) return 1e20;
                MG94Parameters p = mg_p;
                p.beta = w;
                MG94Matrix m;
                m.update(p, ctx.aln.pos_nuc_frequencies, ctx.aln.codon_frequencies_f3x4, *code);
                Scalar conv = (m.scale_factor > 0.0) ? (3.0 / m.scale_factor) : 3.0;
                return -LikelihoodEngine::compute_mg94_log_likelihood(ctx.tree, ctx.aln, m, &ctx.leaf_to_taxon, s_scale * conv);
            }, 0.25, 1.0, 1e-4, 20.0, 1e-4, 50);

            mg_p.beta = std::max(1e-4, opt_mg.a);
            Scalar init_scale = std::max(0.05, opt_mg.b);
            for (auto& node : ctx.tree.nodes) {
                if (node.id != ctx.tree.root_id) {
                    node.branch_length *= init_scale;
                }
            }
            MG94Matrix mg_mat_f81;
            mg_mat_f81.update(mg_p, ctx.aln.pos_nuc_frequencies, ctx.aln.codon_frequencies_f3x4, *code);
            Scalar conv_f81 = (mg_mat_f81.scale_factor > 0.0) ? (3.0 / mg_mat_f81.scale_factor) : 3.0;
            Scalar prev_sweep_ll = LikelihoodEngine::compute_mg94_log_likelihood(ctx.tree, ctx.aln, mg_mat_f81, &ctx.leaf_to_taxon, conv_f81);
            for (int sweep = 0; sweep < 5; ++sweep) {
                sweep_codon_branches_cavity(ctx.tree, ctx.aln, mg_mat_f81, ctx.leaf_to_taxon, conv_f81, 1e-4, 25);
                Scalar cur_ll = LikelihoodEngine::compute_mg94_log_likelihood(ctx.tree, ctx.aln, mg_mat_f81, &ctx.leaf_to_taxon, conv_f81);
                if (std::abs(cur_ll - prev_sweep_ll) < 1e-3) {
                    break;
                }
                prev_sweep_ll = cur_ll;
            }
        }

        MG94Matrix mg_mat;
        mg_mat.update(mg_p, ctx.aln.pos_nuc_frequencies, ctx.aln.codon_frequencies_f3x4, *code);
        Scalar conv = (mg_mat.scale_factor > 0.0) ? (3.0 / mg_mat.scale_factor) : 3.0;
        ctx.mg94_log_l = LikelihoodEngine::compute_mg94_log_likelihood(ctx.tree, ctx.aln, mg_mat, &ctx.leaf_to_taxon, conv);
        ctx.mg94_params = mg_p;
        ctx.mg94_omega = mg_p.beta;
        ctx.mg94_scale_factor = mg_mat.scale_factor;
        ctx.omega = mg_p.beta;
        ctx.tree_scaler = 1.0;

        size_t num_nodes = ctx.tree.num_nodes();
        ctx.branch_tau.assign(num_nodes, 0.0);
        for (const auto& node : ctx.tree.nodes) {
            if (node.id != ctx.tree.root_id) {
                ctx.branch_tau[node.id] = node.branch_length * conv;
            }
        }

        genes.push_back(std::move(ctx));
    }

    // Fast spectral-cavity branch length sweep on the codon tree
    static void sweep_codon_branches_cavity(
        Tree& curr_tree,
        const Alignment& aln,
        const MG94Matrix& mg,
        const std::vector<size_t>& leaf_to_taxon,
        Scalar conv,
        Scalar tol = 1e-4,
        int max_iter = 25
    ) {
        size_t num_nodes = curr_tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = mg.pi.size();
        if (num_nodes == 0 || num_patterns == 0 || !curr_tree.bridges.empty()) return;

        std::vector<std::vector<Vector>> D(num_nodes, std::vector<Vector>(num_patterns, Vector::Ones(S)));
        std::vector<std::vector<Vector>> P_D(num_nodes, std::vector<Vector>(num_patterns, Vector::Ones(S)));
        std::vector<std::vector<Vector>> U(num_nodes, std::vector<Vector>(num_patterns, Vector::Ones(S)));
        std::vector<Matrix> P_branch(num_nodes, Matrix::Identity(S, S));

        for (const auto& node : curr_tree.nodes) {
            if (node.id != curr_tree.root_id && (node.is_leaf || !node.children.empty())) {
                P_branch[node.id] = mg.transition_matrix(node.branch_length * conv);
            }
        }

        for (int32_t nid : curr_tree.post_order) {
            const auto& node = curr_tree.nodes[nid];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                for (size_t p = 0; p < num_patterns; ++p) {
                    int8_t state = (t_idx != static_cast<size_t>(-1)) ? aln.patterns[p].states[t_idx] : -1;
                    if (state >= 0 && state < S) {
                        D[nid][p].setZero();
                        D[nid][p](state) = 1.0;
                    } else {
                        D[nid][p].setOnes();
                    }
                    if (nid != curr_tree.root_id) {
                        P_D[nid][p] = P_branch[nid] * D[nid][p];
                    }
                }
            } else {
                for (size_t p = 0; p < num_patterns; ++p) {
                    Vector d_vec = Vector::Ones(S);
                    for (int32_t cid : node.children) {
                        d_vec.array() *= P_D[cid][p].array();
                        Scalar mx = d_vec.maxCoeff();
                        if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                    }
                    D[nid][p] = d_vec;
                    if (nid != curr_tree.root_id) {
                        P_D[nid][p] = P_branch[nid] * d_vec;
                    }
                }
            }
        }

        for (size_t p = 0; p < num_patterns; ++p) {
            U[curr_tree.root_id][p] = mg.pi;
        }

        std::vector<Vector> c_spec(num_patterns, Vector::Zero(S));
        std::function<void(int32_t)> dfs = [&](int32_t uid) {
            const auto& node = curr_tree.nodes[uid];
            if (uid != curr_tree.root_id && (node.is_leaf || !node.children.empty())) {
                for (size_t p = 0; p < num_patterns; ++p) {
                    Vector u_hat = mg.U_pi.transpose() * U[uid][p];
                    Vector d_hat = mg.V_pi * D[uid][p];
                    c_spec[p] = u_hat.cwiseProduct(d_hat);
                }
                Scalar init_bl = curr_tree.nodes[uid].branch_length;
                auto bl_obj = [&](Scalar bl) -> Scalar {
                    Vector exp_lam = (mg.eigenvalues * (bl * conv)).array().exp();
                    Scalar lnl = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Scalar lk = c_spec[p].dot(exp_lam);
                        lnl += aln.patterns[p].weight * ((lk > 1e-300) ? std::log(lk) : -690.7755);
                    }
                    return -lnl;
                };
                auto [best_bl, _] = Brent1D::minimize(bl_obj, 1e-6, std::max<Scalar>(init_bl, 1e-4), 10.0, tol, max_iter);
                curr_tree.nodes[uid].branch_length = best_bl;
                P_branch[uid] = mg.transition_matrix(best_bl * conv);
                for (size_t p = 0; p < num_patterns; ++p) {
                    P_D[uid][p] = P_branch[uid] * D[uid][p];
                }
            }

            if (!node.is_leaf) {
                int deg = static_cast<int>(node.children.size());
                for (int idx = 0; idx < deg; ++idx) {
                    int32_t cid = node.children[idx];
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Vector cav = (uid == curr_tree.root_id) ? U[uid][p] : (P_branch[uid].transpose() * U[uid][p]);
                        for (int j = 0; j < deg; ++j) {
                            if (j == idx) continue;
                            cav.array() *= P_D[node.children[j]][p].array();
                            Scalar mx = cav.maxCoeff();
                            if (mx > 0.0 && mx < 1e-100) cav /= mx;
                        }
                        Scalar mx = cav.maxCoeff();
                        if (mx > 0.0) cav /= mx;
                        U[cid][p] = cav;
                    }
                    dfs(cid);
                }
                for (size_t p = 0; p < num_patterns; ++p) {
                    Vector d_vec = Vector::Ones(S);
                    for (int32_t cid : node.children) {
                        d_vec.array() *= P_D[cid][p].array();
                        Scalar mx = d_vec.maxCoeff();
                        if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                    }
                    D[uid][p] = d_vec;
                    if (uid != curr_tree.root_id) {
                        P_D[uid][p] = P_branch[uid] * d_vec;
                    }
                }
            }
        };
        dfs(curr_tree.root_id);
    }

    // Compute equilibrium synonymous flux weights for each of the 68 pairs across all loaded genes
    std::vector<Scalar> compute_pair_flux_weights(const MSSPartition& part) const {
        size_t P = part.num_pairs();
        std::vector<Scalar> pair_flux(P, 0.0);
        for (const auto& g : genes) {
            Scalar r_nuc[4][4] = {
                {0.0, g.gtr_params.theta_AC, 1.0, g.gtr_params.theta_AT},
                {g.gtr_params.theta_AC, 0.0, g.gtr_params.theta_CG, g.gtr_params.theta_CT},
                {1.0, g.gtr_params.theta_CG, 0.0, g.gtr_params.theta_GT},
                {g.gtr_params.theta_AT, g.gtr_params.theta_CT, g.gtr_params.theta_GT, 0.0}
            };
            Scalar w_gene = static_cast<Scalar>(g.aln.num_codons);
            for (size_t p = 0; p < P; ++p) {
                const auto& sp = part.pairs[p];
                Scalar theta = r_nuc[sp.nuc_from][sp.nuc_to];
                Scalar f_uv = g.aln.codon_frequencies_f3x4(sp.u) * theta * g.aln.pos_nuc_frequencies(sp.pos, sp.nuc_to);
                Scalar f_vu = g.aln.codon_frequencies_f3x4(sp.v) * theta * g.aln.pos_nuc_frequencies(sp.pos, sp.nuc_from);
                pair_flux[p] += w_gene * (f_uv + f_vu);
            }
        }
        Scalar sum_f = std::accumulate(pair_flux.begin(), pair_flux.end(), 0.0);
        if (sum_f > 0.0) {
            for (auto& v : pair_flux) v /= sum_f;
        } else if (P > 0) {
            std::fill(pair_flux.begin(), pair_flux.end(), 1.0 / static_cast<Scalar>(P));
        }
        return pair_flux;
    }

    // Normalize class rates according to normalization scheme
    static void normalize_class_rates(
        std::vector<Scalar>& rates,
        MSSNormalization norm,
        int neutral_class_id,
        const std::vector<Scalar>& class_flux_weights
    ) {
        size_t K = rates.size();
        if (K == 0) return;
        for (auto& r : rates) r = std::clamp(r, 1e-5, 100.0);

        if (norm == MSSNormalization::NeutralRef && neutral_class_id >= 0 && static_cast<size_t>(neutral_class_id) < K) {
            Scalar ref_val = rates[neutral_class_id];
            if (ref_val > 1e-9) {
                for (auto& r : rates) r /= ref_val;
            }
            rates[neutral_class_id] = 1.0;
        } else if (norm == MSSNormalization::Unweighted) {
            Scalar mean_r = std::accumulate(rates.begin(), rates.end(), 0.0) / static_cast<Scalar>(K);
            if (mean_r > 1e-9) {
                for (auto& r : rates) r /= mean_r;
            }
        } else {
            // Flux-weighted mean = 1.0
            Scalar mean_r = 0.0;
            for (size_t k = 0; k < K; ++k) {
                Scalar w = (k < class_flux_weights.size()) ? class_flux_weights[k] : (1.0 / static_cast<Scalar>(K));
                mean_r += w * rates[k];
            }
            if (mean_r > 1e-9) {
                for (auto& r : rates) r /= mean_r;
            }
        }
    }

    Scalar evaluate_gene_log_l(
        const GeneContext& g,
        const MG94Matrix& m,
        Scalar sc,
        int S
    ) const {
        size_t nn = g.tree.num_nodes();
        std::vector<Matrix> P_br(nn);
        for (const auto& node : g.tree.nodes) {
            if (node.id != g.tree.root_id) {
                P_br[node.id] = m.transition_matrix(g.branch_tau[node.id] * sc);
            }
        }
        std::vector<std::array<Matrix, 3>> P_bridges;
        const std::vector<std::array<Matrix, 3>>* p_bridges_ptr = nullptr;
        if (!g.tree.bridges.empty()) {
            Scalar conv = (g.mg94_scale_factor > 0.0) ? (3.0 / g.mg94_scale_factor) : 3.0;
            P_bridges.resize(g.tree.bridges.size());
            for (size_t k = 0; k < g.tree.bridges.size(); ++k) {
                for (int r = 0; r < 3; ++r) {
                    if (g.tree.bridges[k].weights[r] > 0.0 || r == 0) {
                        P_bridges[k][r] = m.transition_matrix(g.tree.bridges[k].branch_lengths[r] * conv * sc);
                    } else {
                        P_bridges[k][r] = Matrix::Identity(S, S);
                    }
                }
            }
            p_bridges_ptr = &P_bridges;
        }
        return LikelihoodEngine::compute_tree_log_likelihood_from_P(
            g.tree, g.aln, P_br, m.pi, S, &g.leaf_to_taxon, p_bridges_ptr
        );
    }

    // Evaluate composite log-likelihood across all genes for a given MSSPartition + class_rates
    Scalar evaluate_composite_log_l(
        const MSSPartition& part,
        const std::vector<Scalar>& class_rates,
        const std::vector<Scalar>* custom_gene_scalers = nullptr,
        const std::vector<Scalar>* custom_gene_omegas = nullptr
    ) const {
        int S = code->num_sense_codons;
        std::vector<Scalar> syn_mat = part.expand_to_matrix(S, class_rates);
        Scalar total_ll = 0.0;

        for (size_t g_idx = 0; g_idx < genes.size(); ++g_idx) {
            const auto& g = genes[g_idx];
            MG94Parameters p = g.mg94_params;
            p.alpha = 1.0;
            p.beta = custom_gene_omegas ? (*custom_gene_omegas)[g_idx] : g.omega;
            p.syn_pair_rates = &syn_mat;

            MG94Matrix m;
            m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);

            Scalar g_scale = custom_gene_scalers ? (*custom_gene_scalers)[g_idx] : g.tree_scaler;
            total_ll += evaluate_gene_log_l(g, m, g_scale, S);
        }
        return total_ll;
    }

    // Optimize class_rates (and per-gene tree_scaler + optional per-gene omega) for a given MSSPartition
    Scalar fit_partition(
        const MSSPartition& part,
        MSSNormalization norm,
        std::vector<Scalar>& class_rates,
        std::vector<Scalar>& gene_scalers,
        std::vector<Scalar>& gene_omegas,
        bool fix_omega,
        int max_cycles = 8
    ) const {
        size_t K = part.num_classes();
        size_t G = genes.size();
        if (class_rates.size() != K) class_rates.assign(K, 1.0);
        if (gene_scalers.size() != G) {
            gene_scalers.resize(G);
            for (size_t g = 0; g < G; ++g) gene_scalers[g] = genes[g].tree_scaler;
        }
        if (gene_omegas.size() != G) {
            gene_omegas.resize(G);
            for (size_t g = 0; g < G; ++g) gene_omegas[g] = genes[g].mg94_omega;
        }

        // Compute class flux weights from pair flux weights
        std::vector<Scalar> pair_flux = compute_pair_flux_weights(part);
        std::vector<Scalar> class_flux(K, 0.0);
        for (size_t p = 0; p < part.pairs.size(); ++p) {
            int cid = part.pairs[p].class_id;
            if (cid >= 0 && static_cast<size_t>(cid) < K) {
                class_flux[cid] += pair_flux[p];
            }
        }

        normalize_class_rates(class_rates, norm, part.neutral_class_id, class_flux);

        if (part.is_empirical_fixed || K <= 1) {
            // Only optimize gene_scalers and (if !fix_omega) gene_omegas
            int S = code->num_sense_codons;
            std::vector<Scalar> syn_mat = part.expand_to_matrix(S, class_rates);
            for (size_t g_idx = 0; g_idx < G; ++g_idx) {
                const auto& g = genes[g_idx];
                if (fix_omega) {
                    MG94Parameters p = g.mg94_params;
                    p.alpha = 1.0;
                    p.beta = gene_omegas[g_idx];
                    p.syn_pair_rates = &syn_mat;
                    MG94Matrix m;
                    m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);
                    auto [best_s, _] = Brent1D::minimize([&](Scalar sc) {
                        return -evaluate_gene_log_l(g, m, sc, S);
                    }, 1e-3, gene_scalers[g_idx], 1000.0, 1e-4, 40);
                    gene_scalers[g_idx] = best_s;
                } else {
                    auto opt2 = NelderMead2D::minimize([&](Scalar w, Scalar sc) {
                        if (w < 1e-5 || w > 50.0 || sc < 1e-3 || sc > 1000.0) return 1e20;
                        MG94Parameters p = g.mg94_params;
                        p.alpha = 1.0;
                        p.beta = w;
                        p.syn_pair_rates = &syn_mat;
                        MG94Matrix m;
                        m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);
                        return -evaluate_gene_log_l(g, m, sc, S);
                    }, gene_omegas[g_idx], gene_scalers[g_idx], 1e-5, 1000.0, 1e-4, 40);
                    gene_omegas[g_idx] = std::clamp(opt2.a, 1e-5, 50.0);
                    gene_scalers[g_idx] = std::clamp(opt2.b, 1e-3, 1000.0);
                }
            }
            return evaluate_composite_log_l(part, class_rates, &gene_scalers, &gene_omegas);
        }

        // Determine reference class index to hold fixed during coordinate/L-BFGS updates
        int ref_k = part.neutral_class_id;
        if (norm != MSSNormalization::NeutralRef || ref_k < 0 || static_cast<size_t>(ref_k) >= K) {
            // Pick class with largest flux weight as pivot during optimization, then re-normalize
            ref_k = static_cast<int>(std::max_element(class_flux.begin(), class_flux.end()) - class_flux.begin());
        }

        Scalar prev_ll = evaluate_composite_log_l(part, class_rates, &gene_scalers, &gene_omegas);

        for (int cycle = 0; cycle < max_cycles; ++cycle) {
            // Step A: Optimize each non-pivot class rate via 1D Brent
            for (size_t k = 0; k < K; ++k) {
                if (static_cast<int>(k) == ref_k) continue;
                if (class_flux[k] <= 1e-12) continue;

                Scalar init_r = class_rates[k];
                auto obj_k = [&](Scalar r_val) -> Scalar {
                    std::vector<Scalar> cand = class_rates;
                    cand[k] = r_val;
                    if (norm != MSSNormalization::NeutralRef) {
                        normalize_class_rates(cand, norm, part.neutral_class_id, class_flux);
                    }
                    return -evaluate_composite_log_l(part, cand, &gene_scalers, &gene_omegas);
                };
                Scalar lo = std::max(1e-3, init_r * 0.15);
                Scalar hi = std::min(100.0, std::max(10.0, init_r * 6.0));
                auto [best_r, _] = Brent1D::minimize(obj_k, lo, init_r, hi, 1e-4, 40);
                class_rates[k] = best_r;
                if (norm != MSSNormalization::NeutralRef) {
                    normalize_class_rates(class_rates, norm, part.neutral_class_id, class_flux);
                }
            }

            // Step B: Update per-gene tree scaler and (optionally) omega
            int S = code->num_sense_codons;
            std::vector<Scalar> syn_mat = part.expand_to_matrix(S, class_rates);
            #pragma omp parallel for schedule(dynamic) if(G > 1)
            for (size_t g_idx = 0; g_idx < G; ++g_idx) {
                const auto& g = genes[g_idx];
                if (fix_omega) {
                    MG94Parameters p = g.mg94_params;
                    p.alpha = 1.0;
                    p.beta = gene_omegas[g_idx];
                    p.syn_pair_rates = &syn_mat;
                    MG94Matrix m;
                    m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);
                    auto [best_s, _] = Brent1D::minimize([&](Scalar sc) {
                        return -evaluate_gene_log_l(g, m, sc, S);
                    }, 1e-3, gene_scalers[g_idx], 1000.0, 1e-4, 40);
                    gene_scalers[g_idx] = best_s;
                } else {
                    auto opt2 = NelderMead2D::minimize([&](Scalar w, Scalar sc) {
                        if (w < 1e-5 || w > 50.0 || sc < 1e-3 || sc > 1000.0) return 1e20;
                        MG94Parameters p = g.mg94_params;
                        p.alpha = 1.0;
                        p.beta = w;
                        p.syn_pair_rates = &syn_mat;
                        MG94Matrix m;
                        m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);
                        return -evaluate_gene_log_l(g, m, sc, S);
                    }, gene_omegas[g_idx], gene_scalers[g_idx], 1e-5, 1000.0, 1e-4, 40);
                    gene_omegas[g_idx] = std::clamp(opt2.a, 1e-5, 50.0);
                    gene_scalers[g_idx] = std::clamp(opt2.b, 1e-3, 1000.0);
                }
            }

            normalize_class_rates(class_rates, norm, part.neutral_class_id, class_flux);
            Scalar cur_ll = evaluate_composite_log_l(part, class_rates, &gene_scalers, &gene_omegas);
            Scalar gain = cur_ll - prev_ll;
            if (cycle >= 2 && gain < 1e-4) {
                prev_ll = cur_ll;
                break;
            }
            prev_ll = cur_ll;
        }

        // Single-gene 1D Brent scaler polish if G == 1 (F-061)
        if (G == 1 && !fix_omega) {
            int S = code->num_sense_codons;
            std::vector<Scalar> syn_mat = part.expand_to_matrix(S, class_rates);
            const auto& g = genes[0];
            MG94Parameters p = g.mg94_params;
            p.alpha = 1.0;
            p.beta = gene_omegas[0];
            p.syn_pair_rates = &syn_mat;
            MG94Matrix m;
            m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);
            auto [best_s, best_nll] = Brent1D::minimize([&](Scalar sc) {
                return -evaluate_gene_log_l(g, m, sc, S);
            }, 1e-3, gene_scalers[0], 1000.0, 1e-4, 40);
            if (-best_nll > prev_ll) {
                gene_scalers[0] = best_s;
                prev_ll = -best_nll;
            }
        }

        return prev_ll;
    }

    // Fast K-class Partition Selector via 1D Profile-Curvature Dynamic Programming + Local Swap Polish + Akaike Model Averaging
    MSSSelectionResult run_partition_selector(
        const MSSSettings& settings,
        size_t baseline_params,
        size_t total_samples
    ) const {
        MSSSelectionResult sel;
        sel.enabled = true;
        sel.num_classes = std::clamp<size_t>(settings.select_classes, 2, 10);
        sel.ic_type = settings.ic_type;
        bool use_bic = (MSSPartition::to_upper(settings.ic_type) == "BIC");

        auto compute_ic = [&](Scalar ll, size_t k_extra) -> Scalar {
            size_t k_tot = baseline_params + k_extra;
            if (use_bic) {
                return -2.0 * ll + std::log(static_cast<Scalar>(std::max<size_t>(total_samples, 2))) * static_cast<Scalar>(k_tot);
            }
            return stats::compute_aicc(ll, k_tot, total_samples);
        };

        Scalar baseline_ll = 0.0;
        for (const auto& g : genes) baseline_ll += g.mg94_log_l;
        sel.null_mg94_ic = compute_ic(baseline_ll, 0);

        // Step 1: Fit SynREVCodon (68 pair rates) with 1 cycle of coordinate-wise Brent + curvature estimation
        MSSPartition codon_part = MSSPartition::build("SynREVCodon", *code);
        size_t P = codon_part.num_pairs();
        std::vector<Scalar> pair_flux = compute_pair_flux_weights(codon_part);

        std::vector<Scalar> pair_mles(P, 1.0);
        std::vector<Scalar> pair_curv(P, 1.0);
        std::vector<Scalar> g_scalers(genes.size(), 1.0);
        std::vector<Scalar> g_omegas(genes.size());
        for (size_t g = 0; g < genes.size(); ++g) g_omegas[g] = genes[g].mg94_omega;

        // Fast per-pair 1D profile scan around 1.0 to get both MLE alpha_p and quadratic curvature w_p = -d^2 lnL / d(alpha_p)^2
        for (size_t p = 0; p < P; ++p) {
            if (pair_flux[p] <= 1e-12) {
                pair_mles[p] = 1.0;
                pair_curv[p] = 1e-4;
                continue;
            }
            auto obj_p = [&](Scalar r_val) -> Scalar {
                std::vector<Scalar> cand = pair_mles;
                cand[p] = r_val;
                return -evaluate_composite_log_l(codon_part, cand, &g_scalers, &g_omegas);
            };
            auto [best_r, best_nll] = Brent1D::minimize(obj_p, 0.05, 1.0, 5.0, 2e-3, 12);
            pair_mles[p] = best_r;

            // Estimate local quadratic curvature via 3-point finite difference
            Scalar h = std::max(0.05, 0.15 * best_r);
            Scalar r_lo = std::max(0.01, best_r - h);
            Scalar r_hi = best_r + h;
            Scalar ll_mid = -best_nll;
            Scalar ll_lo = -obj_p(r_lo);
            Scalar ll_hi = -obj_p(r_hi);
            Scalar d2 = (ll_hi - 2.0 * ll_mid + ll_lo) / (h * h);
            pair_curv[p] = std::max(1e-3, -d2);
        }

        normalize_class_rates(pair_mles, MSSNormalization::Flux, -1, pair_flux);
        sel.synrev_codon_log_l = evaluate_composite_log_l(codon_part, pair_mles, &g_scalers, &g_omegas);
        sel.synrev_codon_ic = compute_ic(sel.synrev_codon_log_l, P - 1);

        // Sort the P pairs by MLE rate: alpha_{(0)} <= alpha_{(1)} <= ... <= alpha_{(P-1)}
        std::vector<size_t> sorted_idx(P);
        std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
        std::sort(sorted_idx.begin(), sorted_idx.end(), [&](size_t a, size_t b) {
            return pair_mles[a] < pair_mles[b];
        });

        // Candidate model archive for Akaike model averaging: key = class assignment vector -> (IC, logL, class_rates)
        struct ArchivedModel {
            std::vector<int> assignment; // length P, sorted by class rate
            Scalar ic = 1e30;
            Scalar log_l = -1e30;
            std::vector<Scalar> class_rates;
        };
        std::vector<ArchivedModel> archive;

        auto evaluate_and_archive = [&](const std::vector<int>& assign, size_t K_cur, int cycles = 1) -> ArchivedModel {
            MSSPartition cand_part = codon_part;
            cand_part.class_names.resize(K_cur);
            for (size_t k = 0; k < K_cur; ++k) cand_part.class_names[k] = "CLASS_" + std::to_string(k);
            for (size_t p = 0; p < P; ++p) cand_part.pairs[p].class_id = assign[p];

            // Initial class rates from curvature-weighted mean of pair_mles
            std::vector<Scalar> init_rates(K_cur, 1.0);
            for (size_t k = 0; k < K_cur; ++k) {
                Scalar num = 0.0, den = 0.0;
                for (size_t p = 0; p < P; ++p) {
                    if (assign[p] == static_cast<int>(k)) {
                        num += pair_curv[p] * pair_mles[p];
                        den += pair_curv[p];
                    }
                }
                if (den > 0.0) init_rates[k] = num / den;
            }
            std::vector<Scalar> gs = g_scalers;
            std::vector<Scalar> go = g_omegas;
            Scalar ll = fit_partition(cand_part, MSSNormalization::Flux, init_rates, gs, go, settings.fix_omega, cycles);

            // Sort classes so class 0 has smallest rate (SELECTED) and class K-1 has largest rate (NEUTRAL)
            std::vector<size_t> c_order(K_cur);
            std::iota(c_order.begin(), c_order.end(), 0);
            std::sort(c_order.begin(), c_order.end(), [&](size_t a, size_t b) {
                return init_rates[a] < init_rates[b];
            });
            std::vector<int> inv_order(K_cur);
            std::vector<Scalar> sorted_rates(K_cur);
            for (size_t rank = 0; rank < K_cur; ++rank) {
                inv_order[c_order[rank]] = static_cast<int>(rank);
                sorted_rates[rank] = init_rates[c_order[rank]];
            }
            // Normalize relative to highest rate (NEUTRAL = 1.0) for interpretability
            Scalar max_r = sorted_rates.back();
            if (max_r > 1e-9) {
                for (auto& r : sorted_rates) r /= max_r;
            }
            std::vector<int> canon_assign(P);
            for (size_t p = 0; p < P; ++p) canon_assign[p] = inv_order[assign[p]];

            ArchivedModel am;
            am.assignment = std::move(canon_assign);
            am.log_l = ll;
            am.ic = compute_ic(ll, K_cur - 1 + (settings.fix_omega ? 0 : genes.size()));
            am.class_rates = std::move(sorted_rates);
            archive.push_back(am);
            sel.models_evaluated++;
            return am;
        };

        // Step 2: Exact O(P^2 * K) 1D Weighted k-Means (Fisher-Jenks) Dynamic Programming on sorted pair MLEs
        auto solve_dp_partition = [&](size_t K_cur) -> std::vector<int> {
            // Prefix sums of w, w*x, w*x^2
            std::vector<Scalar> pref_w(P + 1, 0.0), pref_wx(P + 1, 0.0), pref_wx2(P + 1, 0.0);
            for (size_t i = 0; i < P; ++i) {
                size_t pid = sorted_idx[i];
                Scalar w = pair_curv[pid];
                Scalar x = pair_mles[pid];
                pref_w[i + 1] = pref_w[i] + w;
                pref_wx[i + 1] = pref_wx[i] + w * x;
                pref_wx2[i + 1] = pref_wx2[i] + w * x * x;
            }
            auto seg_cost = [&](size_t l, size_t r) -> Scalar {
                // inclusive [l, r] in 0..P-1
                Scalar w = pref_w[r + 1] - pref_w[l];
                if (w <= 1e-12) return 0.0;
                Scalar wx = pref_wx[r + 1] - pref_wx[l];
                Scalar wx2 = pref_wx2[r + 1] - pref_wx2[l];
                return std::max(0.0, wx2 - (wx * wx) / w);
            };

            std::vector<std::vector<Scalar>> dp(K_cur + 1, std::vector<Scalar>(P + 1, 1e30));
            std::vector<std::vector<size_t>> prev_cut(K_cur + 1, std::vector<size_t>(P + 1, 0));
            dp[0][0] = 0.0;

            for (size_t k = 1; k <= K_cur; ++k) {
                for (size_t i = k; i <= P; ++i) {
                    for (size_t j = k - 1; j < i; ++j) {
                        Scalar c = dp[k - 1][j] + seg_cost(j, i - 1);
                        if (c < dp[k][i]) {
                            dp[k][i] = c;
                            prev_cut[k][i] = j;
                        }
                    }
                }
            }

            std::vector<int> assign(P, 0);
            size_t cur = P;
            for (size_t k = K_cur; k >= 1; --k) {
                size_t prv = prev_cut[k][cur];
                for (size_t idx = prv; idx < cur; ++idx) {
                    assign[sorted_idx[idx]] = static_cast<int>(k - 1);
                }
                cur = prv;
            }
            return assign;
        };

        size_t min_k = settings.select_auto_k ? 2 : sel.num_classes;
        size_t max_k = settings.select_auto_k ? std::max<size_t>(2, sel.num_classes) : sel.num_classes;

        ArchivedModel global_best;
        global_best.ic = 1e30;
        size_t best_k = min_k;

        for (size_t K_cur = min_k; K_cur <= max_k; ++K_cur) {
            std::vector<int> base_assign = solve_dp_partition(K_cur);
            ArchivedModel best_k_model = evaluate_and_archive(base_assign, K_cur, 2);

            // Step 3: Explore boundary perturbations along the sorted 1D rate order
            std::vector<size_t> cuts;
            for (size_t i = 1; i < P; ++i) {
                if (best_k_model.assignment[sorted_idx[i]] != best_k_model.assignment[sorted_idx[i - 1]]) {
                    cuts.push_back(i);
                }
            }

            // Local cut-shift neighborhood search around DP optimum
            bool improved = true;
            int polish_rounds = 0;
            while (improved && polish_rounds < 2) {
                improved = false;
                polish_rounds++;
                for (size_t c_idx = 0; c_idx < cuts.size(); ++c_idx) {
                    for (int delta_cut : {-2, -1, 1, 2}) {
                        long new_cut = static_cast<long>(cuts[c_idx]) + delta_cut;
                        long lo_bound = (c_idx == 0) ? 1 : static_cast<long>(cuts[c_idx - 1]) + 1;
                        long hi_bound = (c_idx + 1 == cuts.size()) ? static_cast<long>(P) - 1 : static_cast<long>(cuts[c_idx + 1]) - 1;
                        if (new_cut < lo_bound || new_cut > hi_bound) continue;

                        auto trial_cuts = cuts;
                        trial_cuts[c_idx] = static_cast<size_t>(new_cut);
                        std::vector<int> trial_assign(P, 0);
                        size_t seg_start = 0;
                        for (size_t seg = 0; seg <= trial_cuts.size(); ++seg) {
                            size_t seg_end = (seg < trial_cuts.size()) ? trial_cuts[seg] : P;
                            for (size_t idx = seg_start; idx < seg_end; ++idx) {
                                trial_assign[sorted_idx[idx]] = static_cast<int>(seg);
                            }
                            seg_start = seg_end;
                        }
                        ArchivedModel cand = evaluate_and_archive(trial_assign, K_cur, 1);
                        if (cand.ic + 1e-6 < best_k_model.ic) {
                            best_k_model = evaluate_and_archive(trial_assign, K_cur, 2);
                            cuts = trial_cuts;
                            improved = true;
                        }
                    }
                }
            }

            // Step 4: Single-pair greedy boundary flip polish for pairs at class boundaries
            for (size_t i = 0; i < P; ++i) {
                size_t pid = sorted_idx[i];
                int cur_c = best_k_model.assignment[pid];
                bool near_cut = false;
                for (size_t cut_pos : cuts) {
                    if (std::abs(static_cast<long>(i) - static_cast<long>(cut_pos)) <= 1) {
                        near_cut = true;
                        break;
                    }
                }
                if (!near_cut) continue;
                for (int alt_c : {cur_c - 1, cur_c + 1}) {
                    if (alt_c < 0 || static_cast<size_t>(alt_c) >= K_cur) continue;
                    auto trial_assign = best_k_model.assignment;
                    trial_assign[pid] = alt_c;
                    std::vector<int> counts(K_cur, 0);
                    for (int c : trial_assign) counts[c]++;
                    if (*std::min_element(counts.begin(), counts.end()) == 0) continue;

                    ArchivedModel cand = evaluate_and_archive(trial_assign, K_cur, 1);
                    if (cand.ic + 1e-6 < best_k_model.ic) {
                        best_k_model = evaluate_and_archive(trial_assign, K_cur, 2);
                    }
                }
            }

            if (best_k_model.ic < global_best.ic) {
                global_best = best_k_model;
                best_k = K_cur;
            } else if (settings.select_auto_k) {
                // Stop step-up if IC stopped improving
                break;
            }
        }

        sel.num_classes = best_k;
        sel.best_partition_log_l = global_best.log_l;
        sel.best_partition_ic = global_best.ic;
        sel.class_rates = global_best.class_rates;

        // Assign human-readable labels: SELECTED (lowest rate), INTERMEDIATE_i, NEUTRAL (highest rate = 1.0)
        sel.class_labels.resize(best_k);
        sel.class_labels[0] = "SELECTED";
        sel.class_labels[best_k - 1] = "NEUTRAL";
        for (size_t k = 1; k + 1 < best_k; ++k) {
            sel.class_labels[k] = (best_k == 3) ? "INTERMEDIATE" : ("INTERMEDIATE_" + std::to_string(k));
        }

        // Compute Akaike/BIC weights across all archived models with K == best_k
        Scalar min_ic = global_best.ic;
        Scalar sum_w = 0.0;
        std::vector<std::pair<const ArchivedModel*, Scalar>> weighted_models;
        for (const auto& m : archive) {
            if (m.class_rates.size() != best_k) continue;
            Scalar w = std::exp(-0.5 * std::max(0.0, m.ic - min_ic));
            if (w > 1e-6) {
                weighted_models.push_back({&m, w});
                sum_w += w;
            }
        }
        if (sum_w <= 0.0) sum_w = 1.0;

        sel.pair_entries.resize(P);
        for (size_t p = 0; p < P; ++p) {
            MSSPairSelectionEntry pe;
            pe.aa = codon_part.pairs[p].aa;
            pe.codon1 = codon_part.pairs[p].codon1;
            pe.codon2 = codon_part.pairs[p].codon2;
            pe.mle_rate = pair_mles[p];
            pe.curvature = pair_curv[p];
            pe.assigned_class = global_best.assignment[p];
            pe.assigned_label = sel.class_labels[pe.assigned_class];
            pe.class_support.assign(best_k, 0.0);
            pe.model_averaged_rate = 0.0;

            for (const auto& [m_ptr, w_raw] : weighted_models) {
                Scalar w = w_raw / sum_w;
                int c = m_ptr->assignment[p];
                pe.class_support[c] += w;
                pe.model_averaged_rate += w * m_ptr->class_rates[c];
            }
            pe.max_support = *std::max_element(pe.class_support.begin(), pe.class_support.end());
            // Match processor-codon.bf convention: if max support < 0.5, label as AMBIGUOUS
            if (pe.max_support < 0.5) {
                pe.assigned_label = "AMBIGUOUS";
            }
            sel.pair_entries[p] = std::move(pe);
        }

        if (!settings.output_tsv.empty()) {
            std::ofstream out(settings.output_tsv);
            if (out.is_open()) {
                out << "AminoAcid\tCODON1\tCODON2\tGROUP\tSUPPORT\tMODEL_AVG_RATE\tPAIR_MLE\n";
                for (const auto& pe : sel.pair_entries) {
                    out << pe.aa << "\t" << pe.codon1 << "\t" << pe.codon2 << "\t"
                        << pe.assigned_label << "\t" << std::fixed << std::setprecision(4) << pe.max_support << "\t"
                        << pe.model_averaged_rate << "\t" << pe.mle_rate << "\n";
                }
            }
        }

        return sel;
    }

    MSSResult run(const MSSSettings& in_settings) {
        auto t_start = std::chrono::high_resolution_clock::now();
        MSSResult res;
        res.settings = in_settings;

        if (genes.empty()) {
            throw std::runtime_error("MSSAnalyzer::run called with no loaded datasets.");
        }

        // Build partition
        res.partition = MSSPartition::build(
            in_settings.model_type,
            *code,
            in_settings.neutral_ref,
            in_settings.shared_between_rate
        );

        // Determine effective normalization (F-030)
        MSSNormalization norm = in_settings.normalization;
        if (!in_settings.normalization_explicit) {
            if (res.partition.neutral_class_id >= 0) {
                norm = MSSNormalization::NeutralRef;
            } else if (in_settings.fix_omega) {
                norm = MSSNormalization::Unweighted;
            } else {
                norm = MSSNormalization::Flux;
            }
        }
        res.settings.normalization = norm;

        // Baseline summary
        size_t baseline_params = 0;
        size_t total_samples = 0;
        res.baseline_log_l = 0.0;
        for (const auto& g : genes) {
            res.baseline_log_l += g.mg94_log_l;
            size_t n_branches = (g.tree.num_nodes() > 0) ? (g.tree.num_nodes() - 1) : 0;
            baseline_params += (in_settings.use_f81 ? 0 : 5) + 9 + 1 + n_branches;
            total_samples += g.aln.num_taxa * g.aln.num_codons;
        }
        res.baseline_aicc = stats::compute_aicc(res.baseline_log_l, baseline_params, total_samples);
        res.baseline_bic = -2.0 * res.baseline_log_l + std::log(static_cast<Scalar>(std::max<size_t>(total_samples, 2))) * static_cast<Scalar>(baseline_params);

        size_t K = res.partition.num_classes();
        std::vector<Scalar> class_rates(K, 1.0);
        std::vector<Scalar> gene_scalers(genes.size(), 1.0);
        std::vector<Scalar> gene_omegas(genes.size());
        for (size_t g = 0; g < genes.size(); ++g) gene_omegas[g] = genes[g].mg94_omega;

        res.mss_log_l = fit_partition(
            res.partition,
            norm,
            class_rates,
            gene_scalers,
            gene_omegas,
            in_settings.fix_omega,
            8
        );

        size_t extra_params = (!res.partition.is_empirical_fixed && K > 1) ? (K - 1) : 0;
        res.total_parameters = baseline_params + extra_params;
        res.total_sample_size = total_samples;
        res.mss_aicc = stats::compute_aicc(res.mss_log_l, res.total_parameters, total_samples);
        res.mss_bic = -2.0 * res.mss_log_l + std::log(static_cast<Scalar>(std::max<size_t>(total_samples, 2))) * static_cast<Scalar>(res.total_parameters);

        res.omnibus_df = extra_params;
        res.omnibus_lrt = std::max(0.0, 2.0 * (res.mss_log_l - res.baseline_log_l));
        res.omnibus_p_value = (extra_params > 0)
                                  ? stats::chi2_p_value(res.omnibus_lrt, static_cast<Scalar>(extra_params))
                                  : 1.0;

        // Compute per-class pair counts and equilibrium flux shares
        std::vector<Scalar> pair_flux = compute_pair_flux_weights(res.partition);
        std::vector<Scalar> class_flux(K, 0.0);
        std::vector<size_t> class_pair_counts(K, 0);
        for (size_t p = 0; p < res.partition.pairs.size(); ++p) {
            int cid = res.partition.pairs[p].class_id;
            if (cid >= 0 && static_cast<size_t>(cid) < K) {
                class_flux[cid] += pair_flux[p];
                class_pair_counts[cid]++;
            }
        }

        // Build per-class rate estimates, profile CIs, and LRTs (alpha_k == 1)
        res.rates.resize(K);
        std::vector<Scalar> raw_p_vals;
        std::vector<size_t> tested_indices;

        for (size_t k = 0; k < K; ++k) {
            MSSRateEstimate est;
            est.class_name = res.partition.class_names[k];
            est.rate = class_rates[k];
            est.num_pairs = class_pair_counts[k];
            est.equilibrium_flux_share = class_flux[k];
            est.is_reference = (norm == MSSNormalization::NeutralRef && static_cast<int>(k) == res.partition.neutral_class_id);

            if (est.is_reference || res.partition.is_empirical_fixed || K <= 1 || class_flux[k] <= 1e-12) {
                est.ci_low = est.rate;
                est.ci_high = est.rate;
                est.std_error = 0.0;
                est.lrt_vs_one = 0.0;
                est.p_value_uncorrected = 1.0;
                est.p_value_holm = 1.0;
            } else {
                // Evaluate profile log-likelihood at target_r (with per-gene scalers & omegas re-optimized, F-027)
                auto eval_profile_at = [&](Scalar target_r) -> Scalar {
                    std::vector<Scalar> cand = class_rates;
                    if (norm == MSSNormalization::NeutralRef) {
                        cand[k] = target_r;
                    } else {
                        // Scale remaining classes proportionally to maintain mean = 1
                        cand[k] = target_r;
                        Scalar w_k = (norm == MSSNormalization::Unweighted)
                                         ? (1.0 / static_cast<Scalar>(K))
                                         : class_flux[k];
                        Scalar rem_target = std::max(1e-4, 1.0 - w_k * target_r);
                        Scalar rem_cur = 0.0;
                        for (size_t j = 0; j < K; ++j) {
                            if (j == k) continue;
                            Scalar w_j = (norm == MSSNormalization::Unweighted)
                                             ? (1.0 / static_cast<Scalar>(K))
                                             : class_flux[j];
                            rem_cur += w_j * cand[j];
                        }
                        if (rem_cur > 1e-9) {
                            Scalar factor = rem_target / rem_cur;
                            for (size_t j = 0; j < K; ++j) {
                                if (j != k) cand[j] *= factor;
                            }
                        }
                    }

                    std::vector<Scalar> local_scalers = gene_scalers;
                    std::vector<Scalar> local_omegas = gene_omegas;
                    int S_local = code->num_sense_codons;
                    std::vector<Scalar> syn_mat_local = res.partition.expand_to_matrix(S_local, cand);
                    size_t G_local = genes.size();
                    for (size_t g_idx = 0; g_idx < G_local; ++g_idx) {
                        const auto& g = genes[g_idx];
                        if (in_settings.fix_omega) {
                            MG94Parameters p = g.mg94_params;
                            p.alpha = 1.0;
                            p.beta = local_omegas[g_idx];
                            p.syn_pair_rates = &syn_mat_local;
                            MG94Matrix m;
                            m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);
                            auto [best_s, _] = Brent1D::minimize([&](Scalar sc) {
                                return -evaluate_gene_log_l(g, m, sc, S_local);
                            }, 1e-3, local_scalers[g_idx], 1000.0, 1e-4, 25);
                            local_scalers[g_idx] = best_s;
                        } else {
                            auto opt2 = NelderMead2D::minimize([&](Scalar w, Scalar sc) {
                                if (w < 1e-5 || w > 50.0 || sc < 1e-3 || sc > 1000.0) return 1e20;
                                MG94Parameters p = g.mg94_params;
                                p.alpha = 1.0;
                                p.beta = w;
                                p.syn_pair_rates = &syn_mat_local;
                                MG94Matrix m;
                                m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);
                                return -evaluate_gene_log_l(g, m, sc, S_local);
                            }, local_omegas[g_idx], local_scalers[g_idx], 1e-5, 1000.0, 1e-4, 25);
                            local_omegas[g_idx] = std::clamp(opt2.a, 1e-5, 50.0);
                            local_scalers[g_idx] = std::clamp(opt2.b, 1e-3, 1000.0);
                        }
                    }

                    return evaluate_composite_log_l(res.partition, cand, &local_scalers, &local_omegas);
                };

                // Local quadratic curvature for SE and 95% profile CI
                Scalar h = std::max(0.03, 0.10 * est.rate);
                Scalar r_lo = std::max(1e-3, est.rate - h);
                Scalar r_hi = est.rate + h;
                Scalar ll_lo = eval_profile_at(r_lo);
                Scalar ll_hi = eval_profile_at(r_hi);
                Scalar d2 = (ll_hi - 2.0 * res.mss_log_l + ll_lo) / (h * h);
                Scalar info = std::max(1e-2, -d2);
                est.std_error = 1.0 / std::sqrt(info);

                if (in_settings.compute_ci) {
                    // 95% profile likelihood CI (delta lnL = 1.9207) via secant/bisection bracketed by Wald estimate
                    Scalar target_ll = res.mss_log_l - 1.920729;
                    Scalar lo_bracket = std::max(1e-4, est.rate - 2.5 * est.std_error);
                    if (eval_profile_at(lo_bracket) > target_ll) {
                        lo_bracket = std::max(1e-4, est.rate * 0.1);
                    }
                    Scalar a = lo_bracket, b = est.rate;
                    for (int it = 0; it < 10; ++it) {
                        Scalar mid = 0.5 * (a + b);
                        if (eval_profile_at(mid) < target_ll) a = mid;
                        else b = mid;
                    }
                    est.ci_low = 0.5 * (a + b);

                    Scalar hi_bracket = est.rate + 2.5 * est.std_error;
                    if (eval_profile_at(hi_bracket) > target_ll) {
                        hi_bracket = std::max(est.rate + 1.0, est.rate * 3.0);
                    }
                    a = est.rate;
                    b = hi_bracket;
                    for (int it = 0; it < 10; ++it) {
                        Scalar mid = 0.5 * (a + b);
                        if (eval_profile_at(mid) > target_ll) a = mid;
                        else b = mid;
                    }
                    est.ci_high = 0.5 * (a + b);
                } else {
                    est.ci_low = std::max(0.0, est.rate - 1.96 * est.std_error);
                    est.ci_high = est.rate + 1.96 * est.std_error;
                }

                if (in_settings.compute_lrt) {
                    Scalar ll_one = eval_profile_at(1.0);
                    est.lrt_vs_one = std::max(0.0, 2.0 * (res.mss_log_l - ll_one));
                    est.p_value_uncorrected = stats::chi2_p_value(est.lrt_vs_one, 1.0);
                    raw_p_vals.push_back(est.p_value_uncorrected);
                    tested_indices.push_back(k);
                }
            }
            res.rates[k] = std::move(est);
        }

        if (!raw_p_vals.empty()) {
            auto holm_p = stats::holm_bonferroni(raw_p_vals);
            for (size_t i = 0; i < tested_indices.size(); ++i) {
                res.rates[tested_indices[i]].p_value_holm = holm_p[i];
            }
        }

        // Populate per-gene summaries (including synonymous and nonsynonymous component tree lengths)
        int S = code->num_sense_codons;
        std::vector<Scalar> syn_mat = res.partition.expand_to_matrix(S, class_rates);
        res.genes.resize(genes.size());
        for (size_t g_idx = 0; g_idx < genes.size(); ++g_idx) {
            const auto& g = genes[g_idx];
            MSSGeneFitSummary gs;
            gs.filepath = g.filepath;
            gs.num_taxa = g.aln.num_taxa;
            gs.num_codons = g.aln.num_codons;
            gs.baseline_log_l = g.mg94_log_l;
            gs.baseline_omega = g.mg94_omega;
            gs.mss_omega = gene_omegas[g_idx];
            gs.tree_scaler = gene_scalers[g_idx];

            MG94Parameters p = g.mg94_params;
            p.alpha = 1.0;
            p.beta = gene_omegas[g_idx];
            p.syn_pair_rates = &syn_mat;
            MG94Matrix m;
            m.update(p, g.aln.pos_nuc_frequencies, g.aln.codon_frequencies_f3x4, *code);

            // Compute synonymous vs nonsynonymous flux fractions at equilibrium
            Scalar syn_flux = 0.0, nonsyn_flux = 0.0;
            for (int u = 0; u < S; ++u) {
                for (int v = 0; v < S; ++v) {
                    if (u == v) continue;
                    Scalar f = m.pi(u) * m.Q(u, v);
                    if (code->diff_matrix[u][v].synonymous) syn_flux += f;
                    else nonsyn_flux += f;
                }
            }
            Scalar tot_flux = syn_flux + nonsyn_flux;
            Scalar syn_frac = (tot_flux > 0.0) ? (syn_flux / tot_flux) : 0.0;
            Scalar ns_frac = (tot_flux > 0.0) ? (nonsyn_flux / tot_flux) : 0.0;

            Scalar base_T = 0.0, mss_T = 0.0;
            for (const auto& node : g.tree.nodes) {
                if (node.id != g.tree.root_id) {
                    base_T += node.branch_length;
                    Scalar tau_scaled = g.branch_tau[node.id] * gs.tree_scaler;
                    mss_T += tau_scaled * (m.scale_factor / 3.0);
                }
            }
            gs.baseline_tree_length = base_T;
            gs.synonymous_tree_length = mss_T * syn_frac;
            gs.nonsynonymous_tree_length = mss_T * ns_frac;
            gs.mss_log_l = evaluate_gene_log_l(g, m, gs.tree_scaler, S);
            res.genes[g_idx] = std::move(gs);
        }

        // Run fast partition selector if requested
        if (in_settings.select_partition) {
            res.selection = run_partition_selector(in_settings, baseline_params, total_samples);
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t_end - t_start).count();
        return res;
    }

    static nlohmann::json to_json(const MSSResult& res);

    static void save_json(const std::string& filepath, const MSSResult& res) {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file for writing: " + filepath);
        }
        out << to_json(res).dump(2) << "\n";
    }
};

} // namespace hyphy::analyses
