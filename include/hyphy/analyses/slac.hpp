#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/analyses/site_visualizer.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "nlohmann/json.hpp"

#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <regex>
#include <unordered_set>
#include <unordered_map>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

/**
 * Extended Binomial Tail Probability (Left-tail)
 *
 * Computes P(X <= ebx | ebn, ebp) for extended (possibly continuous / fractional)
 * observations under the binomial distribution.
 * Matches legacy HyPhy SLAC extendedBinTail() in SelectionAnalyses/SLAC.bf.
 */
inline double extended_bin_tail(double ebn, double ebp, double ebx) {
    if (ebp <= 0.0) {
        return 0.0;
    }
    if (ebp >= 1.0) {
        return 1.0;
    }
    if (ebx < 0.0) {
        return 0.0;
    }
    if (ebx >= ebn) {
        return 1.0;
    }

    double ebr = std::floor(ebx);
    double current_bin_coeff = std::pow(1.0 - ebp, ebn);
    double bin_head = 0.0;

    for (int ebk = 0; ebk <= static_cast<int>(ebr); ++ebk) {
        bin_head += current_bin_coeff;
        current_bin_coeff = current_bin_coeff * (ebn - ebk) / (ebk + 1.0) * ebp / (1.0 - ebp);
    }

    if (ebx <= std::floor(ebn)) {
        bin_head += current_bin_coeff * (ebx - ebr);
    } else {
        double denom = ebn - std::floor(ebn);
        if (denom > 1e-12) {
            bin_head += (1.0 - bin_head) * (ebx - ebr) / denom;
        }
    }

    return std::min(1.0, std::max(0.0, bin_head));
}

/**
 * Decode standard IUPAC nucleotide character to possible base indices {0=A, 1=C, 2=G, 3=T}
 */
inline std::vector<int8_t> decode_iupac_nucleotide(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return {0};
        case 'C': return {1};
        case 'G': return {2};
        case 'T': case 'U': return {3};
        case 'R': return {0, 2};             // A or G
        case 'Y': return {1, 3};             // C or T
        case 'S': return {1, 2};             // C or G
        case 'W': return {0, 3};             // A or T
        case 'K': return {2, 3};             // G or T
        case 'M': return {0, 1};             // A or C
        case 'B': return {1, 2, 3};          // C, G, T
        case 'D': return {0, 2, 3};          // A, G, T
        case 'H': return {0, 1, 3};          // A, C, T
        case 'V': return {0, 1, 2};          // A, C, G
        case 'N': case '?': case '.': return {0, 1, 2, 3};
        default: return {}; // gap '-' or invalid
    }
}

/**
 * Resolve a codon triplet into all compatible sense codon indices (0..S-1) under genetic code.
 * If gap/missing or stop-codon only, returns an empty vector.
 */
inline std::vector<int8_t> resolve_codon_to_sense(char c1, char c2, char c3, const GeneticCode& code) {
    auto n1s = decode_iupac_nucleotide(c1);
    auto n2s = decode_iupac_nucleotide(c2);
    auto n3s = decode_iupac_nucleotide(c3);
    if (n1s.empty() || n2s.empty() || n3s.empty()) return {};

    std::vector<int8_t> candidates;
    for (int8_t n1 : n1s) {
        for (int8_t n2 : n2s) {
            for (int8_t n3 : n3s) {
                int8_t s = code.codon_lookup[n1][n2][n3];
                if (s >= 0) {
                    candidates.push_back(s);
                }
            }
        }
    }
    return candidates;
}

/**
 * Precomputed SLAC codon metrics:
 * Expected site potentials (SS, NS) and pairwise distance/potential matrices (EPS, EPN, OPS, OPN).
 */
struct SLACCodonMetric {
    int S = 61;
    std::vector<Scalar> SS_64; // Raw 64 codon index -> synonymous site potential [0..3]
    std::vector<Scalar> NS_64; // Raw 64 codon index -> non-synonymous site potential [0..3]
    std::vector<Scalar> SS;    // Sense codon index (0..S-1) -> SS
    std::vector<Scalar> NS;    // Sense codon index (0..S-1) -> NS

    Matrix EPS; // Expected synonymous sites for all sense codon pairs [S x S]
    Matrix EPN; // Expected non-synonymous sites for all sense codon pairs [S x S]
    Matrix OPS; // Inferred synonymous substitutions for all sense codon pairs [S x S]
    Matrix OPN; // Inferred non-synonymous substitutions for all sense codon pairs [S x S]

    static SLACCodonMetric compute(
        const GeneticCode& code,
        const Vector4& nuc_freqs,
        const GTRParameters& gtr_params
    ) {
        SLACCodonMetric m;
        m.S = code.num_sense_codons;
        m.SS_64.assign(64, 0.0);
        m.NS_64.assign(64, 0.0);
        m.SS.assign(m.S, 0.0);
        m.NS.assign(m.S, 0.0);
        m.EPS = Matrix::Zero(m.S, m.S);
        m.EPN = Matrix::Zero(m.S, m.S);
        m.OPS = Matrix::Zero(m.S, m.S);
        m.OPN = Matrix::Zero(m.S, m.S);

        // 1. Build counting bias matrix W (4x4)
        // W[i][j] = theta[i][j] * pi[j] for j != i
        Matrix4 W = Matrix4::Ones();
        Scalar theta[4][4];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                theta[i][j] = 1.0;
            }
        }
        theta[0][1] = theta[1][0] = gtr_params.theta_AC;
        theta[0][2] = theta[2][0] = 1.0; // AG
        theta[0][3] = theta[3][0] = gtr_params.theta_AT;
        theta[1][2] = theta[2][1] = gtr_params.theta_CG;
        theta[1][3] = theta[3][1] = gtr_params.theta_CT;
        theta[2][3] = theta[3][2] = gtr_params.theta_GT;

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i != j) {
                    W(i, j) = theta[i][j] * nuc_freqs(j);
                } else {
                    W(i, j) = 1.0;
                }
            }
        }

        // 2. Compute SS_64 and NS_64 for all 64 codons
        for (int c64 = 0; c64 < 64; ++c64) {
            int n0 = c64 / 16;
            int n1 = (c64 / 4) % 4;
            int n2 = c64 % 4;

            if (code.is_stop[n0][n1][n2]) {
                continue;
            }

            char aa = code.codon_aa[n0][n1][n2];
            int nucs[3] = {n0, n1, n2};

            for (int pos = 0; pos < 3; ++pos) {
                Scalar norm_factor = 0.0;
                Scalar s_sites = 0.0;
                Scalar ns_sites = 0.0;

                for (int new_nuc = 0; new_nuc < 4; ++new_nuc) {
                    if (new_nuc != nucs[pos]) {
                        int copy_nucs[3] = {n0, n1, n2};
                        copy_nucs[pos] = new_nuc;
                        Scalar w = W(nucs[pos], new_nuc);

                        if (!code.is_stop[copy_nucs[0]][copy_nucs[1]][copy_nucs[2]]) {
                            char new_aa = code.codon_aa[copy_nucs[0]][copy_nucs[1]][copy_nucs[2]];
                            if (new_aa == aa) {
                                s_sites += w;
                            } else {
                                ns_sites += w;
                            }
                        }
                        norm_factor += w;
                    }
                }

                if (norm_factor > 0.0) {
                    m.SS_64[c64] += s_sites / norm_factor;
                    m.NS_64[c64] += ns_sites / norm_factor;
                }
            }
        }

        // Map SS/NS to sense codon indices
        for (int s = 0; s < m.S; ++s) {
            int c64 = code.sense_to_codon64[s];
            m.SS[s] = m.SS_64[c64];
            m.NS[s] = m.NS_64[c64];
        }

        // 3. Compute pairwise metrics for all sense codon pairs
        const int all_perms[6][3] = {
            {0, 1, 2},
            {0, 2, 1},
            {1, 0, 2},
            {1, 2, 0},
            {2, 0, 1},
            {2, 1, 0}
        };

        for (int s1 = 0; s1 < m.S; ++s1) {
            int c64_1 = code.sense_to_codon64[s1];
            int n1[3] = {c64_1 / 16, (c64_1 / 4) % 4, c64_1 % 4};
            char aa1 = code.codon_aa[n1[0]][n1[1]][n1[2]];

            m.EPS(s1, s1) = m.SS[s1];
            m.EPN(s1, s1) = m.NS[s1];
            m.OPS(s1, s1) = 0.0;
            m.OPN(s1, s1) = 0.0;

            for (int s2 = s1 + 1; s2 < m.S; ++s2) {
                int c64_2 = code.sense_to_codon64[s2];
                int n2[3] = {c64_2 / 16, (c64_2 / 4) % 4, c64_2 % 4};

                int path_count = 0;
                Scalar eps_total = 0.0;
                Scalar epn_total = 0.0;
                Scalar ops_total = 0.0;
                Scalar opn_total = 0.0;

                for (int path = 0; path < 6; ++path) {
                    int curr_codon[3] = {n1[0], n1[1], n1[2]};
                    char curr_aa = aa1;
                    std::vector<int> codon_seq;
                    codon_seq.push_back(c64_1);

                    int ps = 0;
                    int pn = 0;
                    int step = 0;

                    for (step = 0; step < 3; ++step) {
                        int pos = all_perms[path][step];
                        if (curr_codon[pos] != n2[pos]) {
                            curr_codon[pos] = n2[pos];
                            int next_c64 = 16 * curr_codon[0] + 4 * curr_codon[1] + curr_codon[2];
                            if (code.is_stop[curr_codon[0]][curr_codon[1]][curr_codon[2]]) {
                                break; // invalid path through stop codon
                            }
                            char next_aa = code.codon_aa[curr_codon[0]][curr_codon[1]][curr_codon[2]];
                            codon_seq.push_back(next_c64);
                            if (curr_aa == next_aa) {
                                ps += 1;
                            } else {
                                pn += 1;
                            }
                            curr_aa = next_aa;
                        }
                    }

                    if (step == 3) {
                        path_count += 1;
                        size_t path_len = codon_seq.size();
                        Scalar pes = 0.0;
                        Scalar pns = 0.0;
                        for (int c : codon_seq) {
                            pes += m.SS_64[c];
                            pns += m.NS_64[c];
                        }
                        eps_total += pes / static_cast<Scalar>(path_len);
                        epn_total += pns / static_cast<Scalar>(path_len);
                        ops_total += ps;
                        opn_total += pn;
                    }
                }

                if (path_count > 0) {
                    Scalar eps_val = eps_total / path_count;
                    Scalar epn_val = epn_total / path_count;
                    Scalar ops_val = ops_total / path_count;
                    Scalar opn_val = opn_total / path_count;

                    m.EPS(s1, s2) = m.EPS(s2, s1) = eps_val;
                    m.EPN(s1, s2) = m.EPN(s2, s1) = epn_val;
                    m.OPS(s1, s2) = m.OPS(s2, s1) = ops_val;
                    m.OPN(s1, s2) = m.OPN(s2, s1) = opn_val;
                }
            }
        }

        return m;
    }
};

/**
 * Site-level statistical results for SLAC (11 standard columns)
 */
struct SLACSiteResult {
    size_t site_index = 0;
    Scalar ES = 0.0;
    Scalar EN = 0.0;
    Scalar S = 0.0;
    Scalar N = 0.0;
    Scalar P_S = 0.0;
    Scalar dS = 0.0;
    Scalar dN = 0.0;
    Scalar dN_minus_dS = 0.0;
    Scalar p_positive = 1.0;
    Scalar p_negative = 1.0;
    Scalar total_branch_length = 0.0;

    bool has_dS = false;
    bool has_dN = false;
    bool has_dN_minus_dS = false;

    // Resampling confidence intervals (when samples > 0)
    std::vector<Scalar> quantiles_2_5;
    std::vector<Scalar> quantiles_median;
    std::vector<Scalar> quantiles_97_5;

    std::vector<nlohmann::json> to_legacy_row() const {
        std::vector<nlohmann::json> row(11);
        row[0] = ES;
        row[1] = EN;
        row[2] = S;
        row[3] = N;
        row[4] = P_S;
        row[5] = has_dS ? nlohmann::json(dS) : nullptr;
        row[6] = has_dN ? nlohmann::json(dN) : nullptr;
        row[7] = has_dN_minus_dS ? nlohmann::json(dN_minus_dS) : nullptr;
        row[8] = p_positive;
        row[9] = p_negative;
        row[10] = total_branch_length;
        return row;
    }
};

/**
 * Branch-level statistical results for SLAC (11 standard columns)
 */
struct SLACBranchResult {
    std::string branch_name;
    Scalar ES = 0.0;
    Scalar EN = 0.0;
    Scalar S = 0.0;
    Scalar N = 0.0;
    Scalar P_S = 0.0;
    Scalar dS = 0.0;
    Scalar dN = 0.0;
    Scalar dN_minus_dS = 0.0;
    Scalar p_positive = 1.0;
    Scalar p_negative = 1.0;
    Scalar branch_length = 0.0;

    bool has_dS = false;
    bool has_dN = false;
    bool has_dN_minus_dS = false;

    std::vector<nlohmann::json> to_legacy_row() const {
        std::vector<nlohmann::json> row(11);
        row[0] = ES;
        row[1] = EN;
        row[2] = S;
        row[3] = N;
        row[4] = P_S;
        row[5] = has_dS ? nlohmann::json(dS) : nullptr;
        row[6] = has_dN ? nlohmann::json(dN) : nullptr;
        row[7] = has_dN_minus_dS ? nlohmann::json(dN_minus_dS) : nullptr;
        row[8] = p_positive;
        row[9] = p_negative;
        row[10] = branch_length;
        return row;
    }
};

/**
 * Branch-site reconstructed sequences and substitution counts
 */
struct SLACBranchSiteAttributes {
    std::vector<std::string> codons;
    std::vector<std::string> amino_acids;
    std::vector<Scalar> syn_subs;
    std::vector<Scalar> nonsyn_subs;
};

/**
 * SLAC Analysis Configuration Settings
 */
struct SLACSettings {
    Scalar p_value_threshold = 0.1;
    size_t samples = 100;
    std::string test_branches = "All";
    std::string test_branch_regex = "";
    std::vector<std::string> test_branch_names = {};
    bool full_model = true;
    bool quiet = false;
    bool all_sites = false;
    unsigned int random_seed = 42;
};

/**
 * Complete SLAC Execution Result
 */
struct SLACResult {
    SLACSettings settings;
    Scalar baseline_log_l = 0.0;
    Scalar baseline_aicc = 0.0;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    MG94Parameters mg94_params;
    Scalar global_omega = 1.0;

    Scalar total_tree_length = 0.0;
    Scalar total_tested_tree_length = 0.0;
    size_t num_tested_branches = 0;
    size_t num_background_branches = 0;
    std::vector<std::string> tested_branch_names;
    std::vector<std::string> background_branch_names;

    // Site results for RESOLVED and AVERAGED modes
    std::vector<SLACSiteResult> resolved_sites;
    std::vector<SLACSiteResult> averaged_sites;

    // Branch results for RESOLVED and AVERAGED modes
    std::vector<SLACBranchResult> resolved_branches;
    std::vector<SLACBranchResult> averaged_branches;

    // Branch site attributes: [branch_name] -> attributes
    std::unordered_map<std::string, SLACBranchSiteAttributes> branch_site_attrs;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    size_t total_sites = 0;
    size_t diversifying_sites = 0;
    size_t purifying_sites = 0;
    size_t neutral_sites = 0;
    double runtime_seconds = 0.0;

    Tree tree;
    Alignment aln;

nlohmann::json to_legacy_json(
        const std::string& input_filepath = "",
        const std::string& tree_string = ""
    ) const;

nlohmann::json to_modern_json(
        const std::string& input_filepath = "",
        const std::string& tree_string = "",
        const Provenance& prov = {}
    ) const;
};

/**
 * SLAC Analysis Engine
 */
class SLACAnalyzer {
public:
    Tree tree;
    Alignment aln;
    MG94Parameters base_params;
    SLACSettings settings;

    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    Scalar global_log_l = 0.0;
    Scalar global_aicc = 0.0;
    Scalar conversion_factor = 1.0;

    std::vector<bool> is_test_branch;
    std::vector<std::string> test_branch_names;
    std::vector<std::string> bg_branch_names;
    size_t num_tested = 0;
    size_t num_background = 0;
    Scalar total_tested_length = 0.0;

    MG94Matrix base_mg94_matrix;
    std::vector<Matrix> branch_transition_matrices;
    std::vector<size_t> leaf_to_taxon;

    SLACAnalyzer(Tree t, Alignment a, MG94Parameters params = MG94Parameters{}, SLACSettings s = SLACSettings{})
        : tree(std::move(t)), aln(std::move(a)), base_params(params), settings(std::move(s)) {
        init_leaf_mapping();
        partition_branches();
        prepare_branch_lengths();
    }

    static SLACAnalyzer create_and_fit(
        Tree input_tree,
        Alignment aln,
        SLACSettings settings = SLACSettings{},
        std::function<void(const std::string&, double)> progress_cb = nullptr
    ) {
        if (progress_cb) progress_cb("Phase 1: Fitting Nucleotide GTR Model", 0.1);
        GTRFitter gtr_fitter(input_tree, aln);
        auto gtr_res = gtr_fitter.fit();

        if (progress_cb) progress_cb("Phase 2: Refining under Global MG94 Model", 0.3);
        MG94Parameters base_p;
        base_p.theta_AC = gtr_res.params.theta_AC;
        base_p.theta_AT = gtr_res.params.theta_AT;
        base_p.theta_CG = gtr_res.params.theta_CG;
        base_p.theta_CT = gtr_res.params.theta_CT;
        base_p.theta_GT = gtr_res.params.theta_GT;

        MG94Fitter mg_fitter(gtr_res.tree, aln);
        FitResult mg_res;
        if (settings.full_model) {
            mg_res = mg_fitter.fit_full_model(1.0, base_p, progress_cb);
            base_p = mg_res.params;
        } else {
            mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
            base_p.alpha = 1.0;
            base_p.beta = mg_res.x_opt(0);
        }

        SLACAnalyzer analyzer(mg_res.tree, std::move(aln), base_p, settings);
        analyzer.global_log_l = mg_res.log_likelihood;
        size_t K_mg = 5 + 1 + (mg_res.tree.num_nodes() - 1);
        size_t N_codons = analyzer.aln.num_codons;
        analyzer.global_aicc = -2.0 * mg_res.log_likelihood + 2.0 * K_mg * N_codons / (N_codons - K_mg - 1.0);

        analyzer.has_gtr_fit = true;
        analyzer.gtr_log_l = gtr_res.log_likelihood;
        analyzer.gtr_aicc = gtr_res.aicc;
        analyzer.gtr_rates = gtr_res.params;
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                analyzer.gtr_branch_lengths[node.name] = node.branch_length;
            }
        }

        return analyzer;
    }

    void init_leaf_mapping() {
        leaf_to_taxon.assign(tree.num_nodes(), static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                auto it = aln.taxon_to_index.find(node.name);
                if (it != aln.taxon_to_index.end()) {
                    leaf_to_taxon[node.id] = it->second;
                } else {
                    std::string lower = node.name;
                    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                    auto it_lower = aln.lower_taxon_to_index.find(lower);
                    if (it_lower != aln.lower_taxon_to_index.end()) {
                        leaf_to_taxon[node.id] = it_lower->second;
                    }
                }
            }
        }
    }

    void partition_branches() {
        auto sel = tree.select_branches(
            settings.test_branches,
            settings.test_branch_regex,
            settings.test_branch_names,
            true
        );
        is_test_branch = std::move(sel.is_selected);
        test_branch_names = std::move(sel.selected_names);
        bg_branch_names = std::move(sel.unselected_names);
        num_tested = sel.num_selected;
        num_background = sel.num_unselected;
    }

    void prepare_branch_lengths() {
        base_mg94_matrix.update(base_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *(aln.code ? aln.code : GeneticCode::universal()));
        Scalar Q_scale = base_mg94_matrix.scale_factor;
        conversion_factor = (Q_scale > 1e-12) ? (3.0 / Q_scale) : 1.0;

        total_tested_length = 0.0;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                if (is_test_branch[node.id]) {
                    total_tested_length += node.branch_length;
                }
            }
        }

        // Branch transition matrices for ancestral reconstruction
        branch_transition_matrices.resize(tree.num_nodes());
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar bl_scaled = node.branch_length * conversion_factor;
                branch_transition_matrices[node.id] = base_mg94_matrix.transition_matrix(bl_scaled);
            }
        }
    }

    /**
     * Joint ML Ancestral State Reconstruction for a single pattern
     */
    std::vector<int8_t> reconstruct_joint_ml_pattern(size_t p) const {
        const auto& code = *(aln.code ? aln.code : GeneticCode::universal());
        int S = code.num_sense_codons;
        size_t num_nodes = tree.num_nodes();
        const auto& pattern = aln.patterns[p];

        std::vector<std::vector<Scalar>> M(num_nodes, std::vector<Scalar>(S, 1.0));
        std::vector<std::vector<int8_t>> ptr(num_nodes, std::vector<int8_t>(S, 0));

        // Bottom-up pass
        for (int32_t node_id : tree.post_order) {
            const auto& node = tree.nodes[node_id];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[node_id];
                int8_t st = (t_idx != static_cast<size_t>(-1)) ? pattern.states[t_idx] : -1;
                if (st >= 0 && st < S) {
                    for (int c = 0; c < S; ++c) {
                        M[node_id][c] = (c == st) ? 1.0 : 0.0;
                    }
                } else {
                    for (int c = 0; c < S; ++c) {
                        M[node_id][c] = 1.0;
                    }
                }
            }

            if (node_id != tree.root_id) {
                const Matrix& P = branch_transition_matrices[node_id];
                Scalar max_val = 0.0;
                std::vector<Scalar> val(S, 0.0);

                for (int p_state = 0; p_state < S; ++p_state) {
                    Scalar best_score = -1.0;
                    int8_t best_c = 0;
                    for (int c = 0; c < S; ++c) {
                        Scalar score = P(p_state, c) * M[node_id][c];
                        if (score > best_score) {
                            best_score = score;
                            best_c = static_cast<int8_t>(c);
                        }
                    }
                    ptr[node_id][p_state] = best_c;
                    val[p_state] = std::max(best_score, 0.0);
                    if (val[p_state] > max_val) {
                        max_val = val[p_state];
                    }
                }

                if (max_val > 0.0) {
                    Scalar scale = 1.0 / max_val;
                    for (int p_state = 0; p_state < S; ++p_state) {
                        val[p_state] *= scale;
                    }
                }

                int32_t parent_id = node.parent_id;
                Scalar p_max = 0.0;
                for (int p_state = 0; p_state < S; ++p_state) {
                    M[parent_id][p_state] *= val[p_state];
                    if (M[parent_id][p_state] > p_max) p_max = M[parent_id][p_state];
                }
                if (p_max > 0.0 && p_max < 1e-100) {
                    Scalar p_scale = 1.0 / p_max;
                    for (int p_state = 0; p_state < S; ++p_state) {
                        M[parent_id][p_state] *= p_scale;
                    }
                }
            }
        }

        // Root selection
        int8_t root_state = 0;
        Scalar best_root_score = -1.0;
        for (int c = 0; c < S; ++c) {
            Scalar score = base_mg94_matrix.pi(c) * M[tree.root_id][c];
            if (score > best_root_score) {
                best_root_score = score;
                root_state = static_cast<int8_t>(c);
            }
        }

        // Top-down traceback
        std::vector<int8_t> states(num_nodes, -1);
        states[tree.root_id] = root_state;

        std::vector<int32_t> stack;
        stack.push_back(tree.root_id);
        while (!stack.empty()) {
            int32_t curr = stack.back();
            stack.pop_back();

            for (int32_t child_id : tree.nodes[curr].children) {
                int8_t parent_st = states[curr];
                states[child_id] = ptr[child_id][parent_st];
                stack.push_back(child_id);
            }
        }

        return states;
    }

    /**
     * Compute full SLAC analysis
     */
    SLACResult run(
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        bool show_progress = false,
        bool force_progress = false
    ) {
        auto start_time = std::chrono::high_resolution_clock::now();
        SLACResult res;
        res.settings = settings;
        res.baseline_log_l = global_log_l;
        res.baseline_aicc = global_aicc;
        res.gtr_log_l = gtr_log_l;
        res.gtr_aicc = gtr_aicc;
        res.gtr_rates = gtr_rates;
        res.mg94_params = base_params;
        res.global_omega = base_params.beta;
        res.total_tested_tree_length = total_tested_length;
        res.num_tested_branches = num_tested;
        res.num_background_branches = num_background;
        res.tested_branch_names = test_branch_names;
        res.background_branch_names = bg_branch_names;
        res.gtr_branch_lengths = gtr_branch_lengths;
        res.total_sites = aln.num_codons;
        res.tree = tree;
        res.aln = aln;

        // Tree length
        res.total_tree_length = 0.0;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                res.total_tree_length += node.branch_length;
            }
        }

        if (total_tested_length <= 0.0) {
            throw std::runtime_error("SLAC cannot be applied to a branch selection with total zero branch length.");
        }

        const auto& code = *(aln.code ? aln.code : GeneticCode::universal());
        int S = code.num_sense_codons;
        size_t num_sites = aln.num_codons;
        size_t num_nodes = tree.num_nodes();

        // 1. Precompute codon differences & potential matrices
        if (progress_cb) progress_cb("Precomputing codon pairwise difference and potential matrices", 0.35);
        SLACCodonMetric metric = SLACCodonMetric::compute(code, aln.nuc_frequencies, gtr_rates);

        // 2. Joint ML Ancestral State Reconstruction for each pattern
        if (progress_cb) progress_cb("Performing joint maximum likelihood ancestral state reconstruction", 0.45);
        size_t num_patterns = aln.patterns.size();
        std::vector<std::vector<int8_t>> pattern_ancestors(num_patterns);

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            pattern_ancestors[p] = reconstruct_joint_ml_pattern(p);
        }

        // Map reconstructed states to sites
        std::vector<std::vector<int8_t>> site_ancestors(num_sites);
        for (size_t s = 0; s < num_sites; ++s) {
            site_ancestors[s] = pattern_ancestors[aln.site_to_pattern[s]];
        }

        // 3. Substitution counting on tested branches
        if (progress_cb) progress_cb("Counting substitutions and evaluating binomial selection tests", 0.6);

        // Accumulators for by-site (RESOLVED and AVERAGED)
        // 0: ES, 1: EN, 2: S, 3: N, 4: scaler
        std::vector<std::array<Scalar, 5>> site_acc_res(num_sites, {0.0, 0.0, 0.0, 0.0, total_tested_length});
        std::vector<std::array<Scalar, 5>> site_acc_avg(num_sites, {0.0, 0.0, 0.0, 0.0, total_tested_length});

        // Accumulators for by-branch (RESOLVED and AVERAGED)
        std::vector<std::array<Scalar, 4>> branch_acc_res(num_tested, {0.0, 0.0, 0.0, 0.0});
        std::vector<std::array<Scalar, 4>> branch_acc_avg(num_tested, {0.0, 0.0, 0.0, 0.0});
        std::vector<Scalar> branch_tested_lengths(num_tested, 0.0);

        // Map node_id -> tested branch index (or -1)
        std::vector<int> node_to_tested_idx(num_nodes, -1);
        size_t t_cnt = 0;
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id) && is_test_branch[i]) {
                node_to_tested_idx[i] = t_cnt;
                branch_tested_lengths[t_cnt] = tree.nodes[i].branch_length;
                t_cnt++;
            }
        }

        // Pre-resolve tip states for every site
        // [taxon_idx][site] -> candidates
        std::vector<std::vector<std::vector<int8_t>>> tip_candidates(aln.num_taxa, std::vector<std::vector<int8_t>>(num_sites));
        for (size_t t = 0; t < aln.num_taxa; ++t) {
            const auto& raw_seq = aln.sequences[t].raw_seq;
            for (size_t s = 0; s < num_sites; ++s) {
                tip_candidates[t][s] = resolve_codon_to_sense(raw_seq[3 * s], raw_seq[3 * s + 1], raw_seq[3 * s + 2], code);
            }
        }

        // Counts of observed codons across taxa at each site
        std::vector<std::vector<int>> site_codon_counts(num_sites, std::vector<int>(S, 0));
        for (size_t s = 0; s < num_sites; ++s) {
            for (size_t t = 0; t < aln.num_taxa; ++t) {
                if (tip_candidates[t][s].size() == 1) {
                    site_codon_counts[s][tip_candidates[t][s][0]]++;
                }
            }
        }

        // Substitution mapper structures for branch attributes
        std::unordered_map<std::string, SLACBranchSiteAttributes> branch_attrs;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                SLACBranchSiteAttributes bsa;
                bsa.codons.resize(num_sites, "---");
                bsa.amino_acids.resize(num_sites, "-");
                bsa.syn_subs.resize(num_sites, 0.0);
                bsa.nonsyn_subs.resize(num_sites, 0.0);
                branch_attrs[node.name] = std::move(bsa);
            }
        }

        // Main substitution counting loop
        for (size_t node_id = 0; node_id < num_nodes; ++node_id) {
            if (node_id == static_cast<size_t>(tree.root_id)) continue;
            const auto& node = tree.nodes[node_id];
            int32_t parent_id = node.parent_id;
            if (parent_id < 0) continue;
            bool is_tested = is_test_branch[node_id];
            int t_idx = node_to_tested_idx[node_id];
            Scalar bl = node.branch_length;
            Scalar rel_bl = bl / total_tested_length;
            auto& bsa = branch_attrs[node.name];

            for (size_t s = 0; s < num_sites; ++s) {
                int8_t parent_st = site_ancestors[s][parent_id];
                if (parent_st < 0 || parent_st >= S) continue;

                int8_t child_st = -1;
                std::vector<int8_t> cand;

                if (node.is_leaf) {
                    size_t taxon_idx = leaf_to_taxon[node_id];
                    if (taxon_idx != static_cast<size_t>(-1)) {
                        cand = tip_candidates[taxon_idx][s];
                    }
                    if (cand.size() == 1) {
                        child_st = cand[0];
                    } else if (cand.empty()) {
                        child_st = -1; // missing/gap
                    } else {
                        child_st = -2; // ambiguous
                    }
                } else {
                    child_st = site_ancestors[s][node_id];
                }

                // Record branch attributes (codons, amino-acids, subs)
                if (child_st >= 0 && child_st < S) {
                    bsa.codons[s] = code.sense_codons[child_st];
                    bsa.amino_acids[s] = std::string(1, code.sense_codon_aa[child_st]);
                    bsa.syn_subs[s] = metric.OPS(parent_st, child_st);
                    bsa.nonsyn_subs[s] = metric.OPN(parent_st, child_st);
                } else if (child_st == -2 && !cand.empty()) {
                    bsa.codons[s] = code.sense_codons[cand[0]];
                    std::string aa_set;
                    for (int8_t c : cand) {
                        char a = code.sense_codon_aa[c];
                        if (aa_set.find(a) == std::string::npos) aa_set += a;
                    }
                    bsa.amino_acids[s] = aa_set;
                    Scalar s_sum = 0.0, n_sum = 0.0;
                    for (int8_t c : cand) {
                        s_sum += metric.OPS(parent_st, c);
                        n_sum += metric.OPN(parent_st, c);
                    }
                    bsa.syn_subs[s] = s_sum / cand.size();
                    bsa.nonsyn_subs[s] = n_sum / cand.size();
                } else {
                    bsa.codons[s] = "---";
                    bsa.amino_acids[s] = "-";
                    bsa.syn_subs[s] = 0.0;
                    bsa.nonsyn_subs[s] = 0.0;
                }

                // Now evaluate SLAC site and branch counts
                if (!is_tested) continue;

                if (child_st >= 0 && child_st < S) {
                    // Fully resolved
                    Scalar eps = metric.EPS(parent_st, child_st);
                    Scalar epn = metric.EPN(parent_st, child_st);
                    Scalar ops = metric.OPS(parent_st, child_st);
                    Scalar opn = metric.OPN(parent_st, child_st);

                    Scalar es_site = eps * rel_bl;
                    Scalar en_site = epn * rel_bl;

                    site_acc_res[s][0] += es_site;
                    site_acc_res[s][1] += en_site;
                    site_acc_res[s][2] += ops;
                    site_acc_res[s][3] += opn;

                    site_acc_avg[s][0] += es_site;
                    site_acc_avg[s][1] += en_site;
                    site_acc_avg[s][2] += ops;
                    site_acc_avg[s][3] += opn;

                    branch_acc_res[t_idx][0] += eps;
                    branch_acc_res[t_idx][1] += epn;
                    branch_acc_res[t_idx][2] += ops;
                    branch_acc_res[t_idx][3] += opn;

                    branch_acc_avg[t_idx][0] += eps;
                    branch_acc_avg[t_idx][1] += epn;
                    branch_acc_avg[t_idx][2] += ops;
                    branch_acc_avg[t_idx][3] += opn;
                } else if (child_st == -1) {
                    // Missing / gap: child missing, reduce site scaler by branch length
                    site_acc_res[s][4] -= bl;
                    site_acc_avg[s][4] -= bl;

                    // Branch still accumulates self potentials
                    Scalar eps = metric.EPS(parent_st, parent_st);
                    Scalar epn = metric.EPN(parent_st, parent_st);
                    branch_acc_res[t_idx][0] += eps;
                    branch_acc_res[t_idx][1] += epn;
                    branch_acc_avg[t_idx][0] += eps;
                    branch_acc_avg[t_idx][1] += epn;
                } else if (child_st == -2 && !cand.empty()) {
                    // Ambiguous character: AVERAGED
                    Scalar eps_avg = 0.0, epn_avg = 0.0, ops_avg = 0.0, opn_avg = 0.0;
                    for (int8_t c : cand) {
                        eps_avg += metric.EPS(parent_st, c);
                        epn_avg += metric.EPN(parent_st, c);
                        ops_avg += metric.OPS(parent_st, c);
                        opn_avg += metric.OPN(parent_st, c);
                    }
                    Scalar k_cand = static_cast<Scalar>(cand.size());
                    eps_avg /= k_cand;
                    epn_avg /= k_cand;
                    ops_avg /= k_cand;
                    opn_avg /= k_cand;

                    site_acc_avg[s][0] += eps_avg * rel_bl;
                    site_acc_avg[s][1] += epn_avg * rel_bl;
                    site_acc_avg[s][2] += ops_avg;
                    site_acc_avg[s][3] += opn_avg;

                    branch_acc_avg[t_idx][0] += eps_avg;
                    branch_acc_avg[t_idx][1] += epn_avg;
                    branch_acc_avg[t_idx][2] += ops_avg;
                    branch_acc_avg[t_idx][3] += opn_avg;

                    // RESOLVED: filter to most frequent character at this site
                    int max_cnt = 0;
                    for (int8_t c : cand) {
                        if (site_codon_counts[s][c] > max_cnt) {
                            max_cnt = site_codon_counts[s][c];
                        }
                    }

                    std::vector<int8_t> resolved_cand;
                    if (max_cnt > 0) {
                        for (int8_t c : cand) {
                            if (site_codon_counts[s][c] == max_cnt) {
                                resolved_cand.push_back(c);
                            }
                        }
                    } else {
                        resolved_cand = cand;
                    }

                    Scalar eps_res = 0.0, epn_res = 0.0, ops_res = 0.0, opn_res = 0.0;
                    for (int8_t c : resolved_cand) {
                        eps_res += metric.EPS(parent_st, c);
                        epn_res += metric.EPN(parent_st, c);
                        ops_res += metric.OPS(parent_st, c);
                        opn_res += metric.OPN(parent_st, c);
                    }
                    Scalar k_res = static_cast<Scalar>(resolved_cand.size());
                    eps_res /= k_res;
                    epn_res /= k_res;
                    ops_res /= k_res;
                    opn_res /= k_res;

                    site_acc_res[s][0] += eps_res * rel_bl;
                    site_acc_res[s][1] += epn_res * rel_bl;
                    site_acc_res[s][2] += ops_res;
                    site_acc_res[s][3] += opn_res;

                    branch_acc_res[t_idx][0] += eps_res;
                    branch_acc_res[t_idx][1] += epn_res;
                    branch_acc_res[t_idx][2] += ops_res;
                    branch_acc_res[t_idx][3] += opn_res;
                }
            }
        }

        // 4. Normalize and finalize site-level statistics
        res.resolved_sites.resize(num_sites);
        res.averaged_sites.resize(num_sites);

        auto finalize_site_stats = [&](const std::array<Scalar, 5>& acc, size_t s) -> SLACSiteResult {
            SLACSiteResult sr;
            sr.site_index = s;
            Scalar k = acc[4]; // effective branch length
            sr.total_branch_length = k;

            sr.ES = acc[0];
            sr.EN = acc[1];
            sr.S = acc[2];
            sr.N = acc[3];

            if (k > 0.0) {
                Scalar sc = total_tested_length / k;
                sr.ES *= sc;
                sr.EN *= sc;
            }

            Scalar total_exp = sr.ES + sr.EN;
            sr.P_S = (total_exp > 0.0) ? (sr.ES / total_exp) : 0.0;

            if (sr.ES > 0.0) {
                sr.dS = sr.S / sr.ES;
                sr.has_dS = true;
            } else {
                sr.dS = 0.0;
                sr.has_dS = false;
            }

            if (sr.EN > 0.0) {
                sr.dN = sr.N / sr.EN;
                sr.has_dN = true;
            } else {
                sr.dN = 0.0;
                sr.has_dN = false;
            }

            if (sr.has_dS && sr.has_dN && k > 0.0) {
                sr.dN_minus_dS = (sr.dN - sr.dS) / k;
                sr.has_dN_minus_dS = true;
            } else {
                sr.dN_minus_dS = 0.0;
                sr.has_dN_minus_dS = false;
            }

            Scalar total_subs = sr.S + sr.N;
            if (total_subs > 0.0) {
                sr.p_positive = extended_bin_tail(total_subs, sr.P_S, sr.S);
                if (sr.S <= 0.0) {
                    sr.p_negative = 1.0;
                } else {
                    sr.p_negative = 1.0 - extended_bin_tail(total_subs, sr.P_S, std::max(0.0, sr.S - 1.0));
                }
            } else {
                sr.p_positive = 1.0;
                sr.p_negative = 1.0;
            }

            return sr;
        };

        SiteVisualizer viz(num_sites, num_sites, SiteVisualizerConfig{
            .method_name = "SLAC",
            .p_value_threshold = settings.p_value_threshold,
            .pos_label = "positive",
            .neg_label = "negative",
            .sig_label = "significant",
            .quiet = !show_progress,
            .force_tty = force_progress,
            .enable_chart = true
        });

        for (size_t s = 0; s < num_sites; ++s) {
            res.resolved_sites[s] = finalize_site_stats(site_acc_res[s], s);
            res.averaged_sites[s] = finalize_site_stats(site_acc_avg[s], s);

            Scalar p_pos = std::min(res.resolved_sites[s].p_positive, res.averaged_sites[s].p_positive);
            Scalar p_neg = std::min(res.resolved_sites[s].p_negative, res.averaged_sites[s].p_negative);
            Scalar p_min = std::min(p_pos, p_neg);

            SiteSelectionType type = SiteSelectionType::Neutral;
            std::string detail;
            Scalar metric = res.resolved_sites[s].dN_minus_dS;

            if (p_pos <= settings.p_value_threshold) {
                res.diversifying_sites++;
                type = SiteSelectionType::Positive;
                detail = "▲ dN-dS=" + SiteVisualizer::format_scalar(metric, 2) + ", p=" + SiteVisualizer::format_pvalue(p_pos);
            } else if (p_neg <= settings.p_value_threshold) {
                res.purifying_sites++;
                type = SiteSelectionType::Negative;
                detail = "▼ dN-dS=" + SiteVisualizer::format_scalar(metric, 2) + ", p=" + SiteVisualizer::format_pvalue(p_neg);
            } else {
                res.neutral_sites++;
            }

            viz.record_site(s, type, p_min, metric, detail);
        }

        viz.finish();

        // 5. Finalize branch-level statistics
        res.resolved_branches.resize(num_tested);
        res.averaged_branches.resize(num_tested);

        auto finalize_branch_stats = [&](const std::array<Scalar, 4>& acc, const std::string& bname, Scalar bl) -> SLACBranchResult {
            SLACBranchResult br;
            br.branch_name = bname;
            br.branch_length = bl;
            br.ES = acc[0];
            br.EN = acc[1];
            br.S = acc[2];
            br.N = acc[3];

            Scalar total_exp = br.ES + br.EN;
            br.P_S = (total_exp > 0.0) ? (br.ES / total_exp) : 0.0;

            if (br.ES > 0.0) {
                br.dS = br.S / br.ES;
                br.has_dS = true;
            }
            if (br.EN > 0.0) {
                br.dN = br.N / br.EN;
                br.has_dN = true;
            }
            if (br.has_dS && br.has_dN && bl > 0.0) {
                br.dN_minus_dS = (br.dN - br.dS) / bl;
                br.has_dN_minus_dS = true;
            }

            Scalar total_subs = br.S + br.N;
            if (total_subs > 0.0) {
                br.p_positive = extended_bin_tail(total_subs, br.P_S, br.S);
                if (br.S <= 0.0) {
                    br.p_negative = 1.0;
                } else {
                    br.p_negative = 1.0 - extended_bin_tail(total_subs, br.P_S, std::max(0.0, br.S - 1.0));
                }
            } else {
                br.p_positive = 1.0;
                br.p_negative = 1.0;
            }

            return br;
        };

        for (size_t b = 0; b < num_tested; ++b) {
            res.resolved_branches[b] = finalize_branch_stats(branch_acc_res[b], test_branch_names[b], branch_tested_lengths[b]);
            res.averaged_branches[b] = finalize_branch_stats(branch_acc_avg[b], test_branch_names[b], branch_tested_lengths[b]);
        }

        res.branch_site_attrs = std::move(branch_attrs);

        // 6. Ancestral State Resampling (Confidence Intervals) when samples > 0
        if (settings.samples > 0) {
            if (progress_cb) progress_cb("Generating ancestral sequence samples for confidence intervals", 0.75);

            // Precompute conditional likelihood prunings L_node for all patterns
            // L_node[p][node_id][c] = P(subtree under node_id | node_id = c)
            std::vector<std::vector<std::vector<Scalar>>> pattern_L(num_patterns, std::vector<std::vector<Scalar>>(num_nodes, std::vector<Scalar>(S, 1.0)));

            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                for (int32_t node_id : tree.post_order) {
                    const auto& node = tree.nodes[node_id];
                    if (node.is_leaf) {
                        size_t t_idx = leaf_to_taxon[node_id];
                        int8_t st = (t_idx != static_cast<size_t>(-1)) ? pattern.states[t_idx] : -1;
                        if (st >= 0 && st < S) {
                            for (int c = 0; c < S; ++c) {
                                pattern_L[p][node_id][c] = (c == st) ? 1.0 : 0.0;
                            }
                        } else {
                            for (int c = 0; c < S; ++c) {
                                pattern_L[p][node_id][c] = 1.0;
                            }
                        }
                    }

                    if (node_id != tree.root_id) {
                        const Matrix& P = branch_transition_matrices[node_id];
                        std::vector<Scalar> val(S, 0.0);
                        Scalar max_v = 0.0;

                        for (int p_state = 0; p_state < S; ++p_state) {
                            Scalar sum = 0.0;
                            for (int c = 0; c < S; ++c) {
                                sum += P(p_state, c) * pattern_L[p][node_id][c];
                            }
                            val[p_state] = sum;
                            if (sum > max_v) max_v = sum;
                        }

                        if (max_v > 0.0) {
                            Scalar scale = 1.0 / max_v;
                            for (int p_state = 0; p_state < S; ++p_state) {
                                val[p_state] *= scale;
                            }
                        }

                        int32_t parent_id = node.parent_id;
                        Scalar p_max = 0.0;
                        for (int p_state = 0; p_state < S; ++p_state) {
                            pattern_L[p][parent_id][p_state] *= val[p_state];
                            if (pattern_L[p][parent_id][p_state] > p_max) p_max = pattern_L[p][parent_id][p_state];
                        }
                        if (p_max > 0.0 && p_max < 1e-100) {
                            Scalar p_scale = 1.0 / p_max;
                            for (int p_state = 0; p_state < S; ++p_state) {
                                pattern_L[p][parent_id][p_state] *= p_scale;
                            }
                        }
                    }
                }
            }

            // Storage for sampled values across all replicates:
            // [site][column] -> vector of sampled values (size = samples)
            std::vector<std::vector<std::vector<Scalar>>> sample_vals_res(num_sites, std::vector<std::vector<Scalar>>(11, std::vector<Scalar>(settings.samples, 0.0)));
            std::vector<std::vector<std::vector<Scalar>>> sample_vals_avg(num_sites, std::vector<std::vector<Scalar>>(11, std::vector<Scalar>(settings.samples, 0.0)));

            #pragma omp parallel
            {
                std::mt19937_64 rng(settings.random_seed + omp_get_thread_num());

                #pragma omp for schedule(dynamic)
                for (size_t rep = 0; rep < settings.samples; ++rep) {
                    // Sample states for all patterns
                    std::vector<std::vector<int8_t>> sampled_states_pattern(num_patterns, std::vector<int8_t>(num_nodes, -1));

                    for (size_t p = 0; p < num_patterns; ++p) {
                        // Sample root
                        std::vector<Scalar> root_weights(S);
                        for (int c = 0; c < S; ++c) {
                            root_weights[c] = base_mg94_matrix.pi(c) * pattern_L[p][tree.root_id][c];
                        }
                        std::discrete_distribution<int> root_dist(root_weights.begin(), root_weights.end());
                        sampled_states_pattern[p][tree.root_id] = static_cast<int8_t>(root_dist(rng));

                        // Top-down sample for children
                        std::vector<int32_t> stack;
                        stack.push_back(tree.root_id);
                        while (!stack.empty()) {
                            int32_t curr = stack.back();
                            stack.pop_back();

                            for (int32_t child_id : tree.nodes[curr].children) {
                                const auto& child_node = tree.nodes[child_id];
                                int8_t parent_st = sampled_states_pattern[p][curr];

                                if (child_node.is_leaf) {
                                    size_t t_idx = leaf_to_taxon[child_id];
                                    int8_t st = (t_idx != static_cast<size_t>(-1)) ? aln.patterns[p].states[t_idx] : -1;
                                    sampled_states_pattern[p][child_id] = st;
                                } else {
                                    const Matrix& P = branch_transition_matrices[child_id];
                                    std::vector<Scalar> child_weights(S);
                                    for (int c = 0; c < S; ++c) {
                                        child_weights[c] = P(parent_st, c) * pattern_L[p][child_id][c];
                                    }
                                    std::discrete_distribution<int> child_dist(child_weights.begin(), child_weights.end());
                                    sampled_states_pattern[p][child_id] = static_cast<int8_t>(child_dist(rng));
                                    stack.push_back(child_id);
                                }
                            }
                        }
                    }

                    // Count substitution events for this sample
                    std::vector<std::array<Scalar, 5>> s_acc_res(num_sites, {0.0, 0.0, 0.0, 0.0, total_tested_length});
                    std::vector<std::array<Scalar, 5>> s_acc_avg(num_sites, {0.0, 0.0, 0.0, 0.0, total_tested_length});

                    for (size_t node_id = 0; node_id < num_nodes; ++node_id) {
                        if (node_id == static_cast<size_t>(tree.root_id) || !is_test_branch[node_id]) continue;
                        const auto& node = tree.nodes[node_id];
                        int32_t parent_id = node.parent_id;
                        if (parent_id < 0) continue;
                        Scalar bl = node.branch_length;
                        Scalar rel_bl = bl / total_tested_length;

                        for (size_t s = 0; s < num_sites; ++s) {
                            size_t p = aln.site_to_pattern[s];
                            int8_t parent_st = sampled_states_pattern[p][parent_id];
                            int8_t child_st = sampled_states_pattern[p][node_id];

                            if (parent_st < 0 || parent_st >= S) continue;

                            if (node.is_leaf) {
                                size_t t_idx = leaf_to_taxon[node_id];
                                const auto& cand = (t_idx != static_cast<size_t>(-1)) ? tip_candidates[t_idx][s] : std::vector<int8_t>{};
                                if (cand.empty()) {
                                    child_st = -1;
                                } else if (cand.size() == 1) {
                                    child_st = cand[0];
                                } else {
                                    child_st = -2;
                                }

                                if (child_st >= 0 && child_st < S) {
                                    Scalar eps = metric.EPS(parent_st, child_st);
                                    Scalar epn = metric.EPN(parent_st, child_st);
                                    Scalar ops = metric.OPS(parent_st, child_st);
                                    Scalar opn = metric.OPN(parent_st, child_st);

                                    s_acc_res[s][0] += eps * rel_bl;
                                    s_acc_res[s][1] += epn * rel_bl;
                                    s_acc_res[s][2] += ops;
                                    s_acc_res[s][3] += opn;

                                    s_acc_avg[s][0] += eps * rel_bl;
                                    s_acc_avg[s][1] += epn * rel_bl;
                                    s_acc_avg[s][2] += ops;
                                    s_acc_avg[s][3] += opn;
                                } else if (child_st == -1) {
                                    s_acc_res[s][4] -= bl;
                                    s_acc_avg[s][4] -= bl;
                                } else {
                                    // Ambiguous tip in sample
                                    Scalar eps_avg = 0.0, epn_avg = 0.0, ops_avg = 0.0, opn_avg = 0.0;
                                    for (int8_t c : cand) {
                                        eps_avg += metric.EPS(parent_st, c);
                                        epn_avg += metric.EPN(parent_st, c);
                                        ops_avg += metric.OPS(parent_st, c);
                                        opn_avg += metric.OPN(parent_st, c);
                                    }
                                    eps_avg /= cand.size();
                                    epn_avg /= cand.size();
                                    ops_avg /= cand.size();
                                    opn_avg /= cand.size();

                                    s_acc_avg[s][0] += eps_avg * rel_bl;
                                    s_acc_avg[s][1] += epn_avg * rel_bl;
                                    s_acc_avg[s][2] += ops_avg;
                                    s_acc_avg[s][3] += opn_avg;

                                    s_acc_res[s][0] += eps_avg * rel_bl;
                                    s_acc_res[s][1] += epn_avg * rel_bl;
                                    s_acc_res[s][2] += ops_avg;
                                    s_acc_res[s][3] += opn_avg;
                                }
                            } else {
                                if (child_st >= 0 && child_st < S) {
                                    Scalar eps = metric.EPS(parent_st, child_st);
                                    Scalar epn = metric.EPN(parent_st, child_st);
                                    Scalar ops = metric.OPS(parent_st, child_st);
                                    Scalar opn = metric.OPN(parent_st, child_st);

                                    s_acc_res[s][0] += eps * rel_bl;
                                    s_acc_res[s][1] += epn * rel_bl;
                                    s_acc_res[s][2] += ops;
                                    s_acc_res[s][3] += opn;

                                    s_acc_avg[s][0] += eps * rel_bl;
                                    s_acc_avg[s][1] += epn * rel_bl;
                                    s_acc_avg[s][2] += ops;
                                    s_acc_avg[s][3] += opn;
                                }
                            }
                        }
                    }

                    // Evaluate statistics for this sample
                    for (size_t s = 0; s < num_sites; ++s) {
                        auto r_res = finalize_site_stats(s_acc_res[s], s);
                        auto r_avg = finalize_site_stats(s_acc_avg[s], s);

                        sample_vals_res[s][0][rep] = r_res.ES;
                        sample_vals_res[s][1][rep] = r_res.EN;
                        sample_vals_res[s][2][rep] = r_res.S;
                        sample_vals_res[s][3][rep] = r_res.N;
                        sample_vals_res[s][4][rep] = r_res.P_S;
                        sample_vals_res[s][5][rep] = r_res.dS;
                        sample_vals_res[s][6][rep] = r_res.dN;
                        sample_vals_res[s][7][rep] = r_res.dN_minus_dS;
                        sample_vals_res[s][8][rep] = r_res.p_positive;
                        sample_vals_res[s][9][rep] = r_res.p_negative;
                        sample_vals_res[s][10][rep] = r_res.total_branch_length;

                        sample_vals_avg[s][0][rep] = r_avg.ES;
                        sample_vals_avg[s][1][rep] = r_avg.EN;
                        sample_vals_avg[s][2][rep] = r_avg.S;
                        sample_vals_avg[s][3][rep] = r_avg.N;
                        sample_vals_avg[s][4][rep] = r_avg.P_S;
                        sample_vals_avg[s][5][rep] = r_avg.dS;
                        sample_vals_avg[s][6][rep] = r_avg.dN;
                        sample_vals_avg[s][7][rep] = r_avg.dN_minus_dS;
                        sample_vals_avg[s][8][rep] = r_avg.p_positive;
                        sample_vals_avg[s][9][rep] = r_avg.p_negative;
                        sample_vals_avg[s][10][rep] = r_avg.total_branch_length;
                    }
                }
            }

            // Continuous quantile calculation (Type 7 / HyPhy stats.Quantile)
            auto calc_quantile = [](std::vector<Scalar>& v, double q) -> Scalar {
                if (v.empty()) return 0.0;
                std::sort(v.begin(), v.end());
                if (v.size() == 1) return v[0];
                double idx = q * (v.size() - 1);
                size_t lo = static_cast<size_t>(std::floor(idx));
                size_t hi = static_cast<size_t>(std::ceil(idx));
                double w = idx - lo;
                return (1.0 - w) * v[lo] + w * v[hi];
            };

            for (size_t s = 0; s < num_sites; ++s) {
                res.resolved_sites[s].quantiles_2_5.resize(11);
                res.resolved_sites[s].quantiles_median.resize(11);
                res.resolved_sites[s].quantiles_97_5.resize(11);

                res.averaged_sites[s].quantiles_2_5.resize(11);
                res.averaged_sites[s].quantiles_median.resize(11);
                res.averaged_sites[s].quantiles_97_5.resize(11);

                for (size_t col = 0; col < 11; ++col) {
                    res.resolved_sites[s].quantiles_2_5[col] = calc_quantile(sample_vals_res[s][col], 0.025);
                    res.resolved_sites[s].quantiles_median[col] = calc_quantile(sample_vals_res[s][col], 0.50);
                    res.resolved_sites[s].quantiles_97_5[col] = calc_quantile(sample_vals_res[s][col], 0.975);

                    res.averaged_sites[s].quantiles_2_5[col] = calc_quantile(sample_vals_avg[s][col], 0.025);
                    res.averaged_sites[s].quantiles_median[col] = calc_quantile(sample_vals_avg[s][col], 0.50);
                    res.averaged_sites[s].quantiles_97_5[col] = calc_quantile(sample_vals_avg[s][col], 0.975);
                }
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(end_time - start_time).count();

        if (progress_cb) progress_cb("SLAC analysis complete", 1.0);
        return res;
    }
};

} // namespace hyphy::analyses
