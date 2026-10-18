#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/stats.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/analyses/site_visualizer.hpp"
#include "nlohmann/json.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <regex>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

// ============================================================================
// Amino Acid & Property Definitions
// ============================================================================

inline int aa_char_to_index(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'D': return 2;
        case 'E': return 3;
        case 'F': return 4;
        case 'G': return 5;
        case 'H': return 6;
        case 'I': return 7;
        case 'K': return 8;
        case 'L': return 9;
        case 'M': return 10;
        case 'N': return 11;
        case 'P': return 12;
        case 'Q': return 13;
        case 'R': return 14;
        case 'S': return 15;
        case 'T': return 16;
        case 'V': return 17;
        case 'W': return 18;
        case 'Y': return 19;
        default: return -1;
    }
}

inline char aa_index_to_char(int idx) {
    static const char aa_list[20] = {
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    };
    if (idx >= 0 && idx < 20) return aa_list[idx];
    return '?';
}

struct AminoAcidProperty {
    std::string name;
    std::string description;
    std::array<Scalar, 20> raw_values;
    std::array<Scalar, 20> z_scores;
    Scalar mean = 0.0;
    Scalar std_dev = 1.0;

    void standardize() {
        Scalar sum = 0.0;
        Scalar sum_sq = 0.0;
        for (int i = 0; i < 20; ++i) {
            sum += raw_values[i];
            sum_sq += raw_values[i] * raw_values[i];
        }
        mean = sum / 20.0;
        Scalar var = (sum_sq / 20.0) - (mean * mean);
        std_dev = (var > 1e-12) ? std::sqrt(var) : 1.0;
        for (int i = 0; i < 20; ++i) {
            z_scores[i] = (raw_values[i] - mean) / std_dev;
        }
    }
};

class PropertyDatabase {
public:
    static std::vector<AminoAcidProperty> get_predefined_properties(const std::string& set_name) {
        std::string lower_name = set_name;
        std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);

        if (lower_name == "atchley") {
            return get_atchley_properties();
        } else if (lower_name == "2prop") {
            return get_2prop_properties();
        } else if (lower_name == "3prop") {
            return get_3prop_properties();
        } else if (lower_name == "4prop") {
            return get_4prop_properties();
        } else {
            // Default: 5PROP
            return get_5prop_properties();
        }
    }

    static std::vector<AminoAcidProperty> get_5prop_properties() {
        std::vector<AminoAcidProperty> props;

        // 1. AlphaHelix_Propensity_ChouFasman
        // Alphabetical ordering matches HyPhy 2.5 table indexing
        props.push_back({
            "AlphaHelix_Propensity_ChouFasman",
            "Chou-Fasman Alpha Helix Propensity",
            // A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
            {1.42, 0.70, 1.01, 1.51, 1.13, 0.57, 1.00, 1.08, 1.16, 1.21,
             1.45, 0.67, 0.57, 1.11, 0.98, 0.77, 0.83, 1.06, 1.08, 0.69},
            {}, 0.0, 1.0
        });

        // 2. BetaSheet_Propensity_ChouFasman
        props.push_back({
            "BetaSheet_Propensity_ChouFasman",
            "Chou-Fasman Beta Sheet Propensity",
            {0.83, 1.19, 0.54, 0.37, 1.38, 0.75, 0.87, 1.60, 0.74, 1.30,
             1.05, 0.89, 0.55, 1.10, 0.93, 0.75, 1.19, 1.70, 1.37, 1.47},
            {}, 0.0, 1.0
        });

        // 3. Hydrophobicity_KyteDoolittle
        props.push_back({
            "Hydrophobicity_KyteDoolittle",
            "Kyte-Doolittle Hydropathy Index",
            {1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
             1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3},
            {}, 0.0, 1.0
        });

        // 4. Isoelectric_Point_pI
        props.push_back({
            "Isoelectric_Point_pI",
            "Isoelectric Point (pI)",
            {6.00, 5.07, 2.77, 3.22, 5.48, 5.97, 7.59, 6.02, 9.74, 5.98,
             5.74, 5.41, 6.30, 5.65, 10.76, 5.68, 5.60, 5.96, 5.89, 5.66},
            {}, 0.0, 1.0
        });

        // 5. Volume_Angstrom3
        props.push_back({
            "Volume_Angstrom3",
            "Residue Molecular Volume (Angstrom^3)",
            {88.6, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, 168.6, 166.7,
             162.9, 114.1, 112.7, 143.8, 173.4, 89.0, 116.1, 140.0, 227.8, 193.6},
            {}, 0.0, 1.0
        });

        for (auto& p : props) p.standardize();
        return props;
    }

    static std::vector<AminoAcidProperty> get_3prop_properties() {
        std::vector<AminoAcidProperty> props;
        props.push_back({
            "Hydrophobicity_KyteDoolittle",
            "Kyte-Doolittle Hydropathy Index",
            {1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
             1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3},
            {}, 0.0, 1.0
        });
        props.push_back({
            "Isoelectric_Point_pI",
            "Isoelectric Point (pI)",
            {6.00, 5.07, 2.77, 3.22, 5.48, 5.97, 7.59, 6.02, 9.74, 5.98,
             5.74, 5.41, 6.30, 5.65, 10.76, 5.68, 5.60, 5.96, 5.89, 5.66},
            {}, 0.0, 1.0
        });
        props.push_back({
            "Volume_Angstrom3",
            "Residue Molecular Volume (Angstrom^3)",
            {88.6, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, 168.6, 166.7,
             162.9, 114.1, 112.7, 143.8, 173.4, 89.0, 116.1, 140.0, 227.8, 193.6},
            {}, 0.0, 1.0
        });
        for (auto& p : props) p.standardize();
        return props;
    }

    static std::vector<AminoAcidProperty> get_4prop_properties() {
        auto props = get_5prop_properties();
        // Keep AlphaHelix, Hydrophobicity, pI, Volume
        props.erase(props.begin() + 1); // remove BetaSheet
        return props;
    }

    static std::vector<AminoAcidProperty> get_2prop_properties() {
        std::vector<AminoAcidProperty> props;
        props.push_back({
            "Hydrophobicity_KyteDoolittle",
            "Kyte-Doolittle Hydropathy Index",
            {1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
             1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3},
            {}, 0.0, 1.0
        });
        props.push_back({
            "Volume_Angstrom3",
            "Residue Molecular Volume (Angstrom^3)",
            {88.6, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, 168.6, 166.7,
             162.9, 114.1, 112.7, 143.8, 173.4, 89.0, 116.1, 140.0, 227.8, 193.6},
            {}, 0.0, 1.0
        });
        for (auto& p : props) p.standardize();
        return props;
    }

    static std::vector<AminoAcidProperty> get_atchley_properties() {
        std::vector<AminoAcidProperty> props;
        props.push_back({
            "Factor I bipolar",
            "Atchley Factor I: Bipolar (Hydrophobicity / Polarity)",
            {-0.591, -1.343, 1.050, 1.357, -1.006, -0.384, 0.336, -1.239, 1.831, -1.019,
             -0.663, 0.945, 0.189, 0.931, 1.538, -0.228, -0.032, -1.337, -0.595, 0.260},
            {}, 0.0, 1.0
        });
        props.push_back({
            "Factor II secondary structure",
            "Atchley Factor II: Secondary Structure (Helix / Beta-turn)",
            {-1.302, 0.465, 0.302, -1.453, -0.590, 1.652, -0.417, -0.547, -0.561, -0.987,
             -1.524, 0.828, 2.081, -0.179, -0.055, 1.399, 0.326, -0.279, 0.009, 0.830},
            {}, 0.0, 1.0
        });
        props.push_back({
            "Factor III volume",
            "Atchley Factor III: Molecular Size / Volume",
            {-0.733, -0.862, -3.656, 1.477, 1.891, 1.330, -1.673, 2.131, 0.533, -1.505,
              2.219, 1.299, -1.628, -3.005, 1.502, -4.760, 2.213, -0.544, 0.672, 3.097},
            {}, 0.0, 1.0
        });
        props.push_back({
            "Factor IV composition",
            "Atchley Factor IV: Amino Acid Composition",
            {1.570, -1.020, -0.259, 0.113, -0.397, 1.045, -1.474, 0.393, -0.277, 1.266,
             -1.005, -0.169, 0.421, -0.503, 0.440, 0.670, 0.908, 1.242, -2.128, -0.838},
            {}, 0.0, 1.0
        });
        props.push_back({
            "Factor V charge",
            "Atchley Factor V: Electrostatic Net Charge",
            {-0.146, -0.255, -3.242, -0.837, 0.412, 2.064, -0.078, 0.816, 1.648, -0.912,
              1.212, 0.933, -1.392, -1.853, 2.897, -2.647, 1.313, -1.262, -0.184, 1.512},
            {}, 0.0, 1.0
        });
        for (auto& p : props) p.standardize();
        return props;
    }

    static std::vector<AminoAcidProperty> from_json_file(const std::string& filepath) {
        std::ifstream f(filepath);
        if (!f.is_open()) {
            throw std::runtime_error("Cannot open amino acid properties JSON file: " + filepath);
        }
        nlohmann::json j;
        f >> j;

        std::vector<AminoAcidProperty> props;
        for (auto& [prop_name, prop_val] : j.items()) {
            AminoAcidProperty p;
            p.name = prop_name;
            p.description = prop_name;
            if (prop_val.is_object()) {
                for (int i = 0; i < 20; ++i) {
                    char aa = aa_index_to_char(i);
                    std::string key(1, aa);
                    if (prop_val.contains(key)) {
                        p.raw_values[i] = prop_val[key].get<Scalar>();
                    } else {
                        throw std::runtime_error("Property " + prop_name + " missing value for amino acid " + key);
                    }
                }
            } else if (prop_val.is_array() && prop_val.size() == 20) {
                for (int i = 0; i < 20; ++i) {
                    p.raw_values[i] = prop_val[i].get<Scalar>();
                }
            } else {
                throw std::runtime_error("Invalid property format for " + prop_name);
            }
            p.standardize();
            props.push_back(p);
        }
        std::sort(props.begin(), props.end(), [](const AminoAcidProperty& a, const AminoAcidProperty& b) {
            return a.name < b.name;
        });
        return props;
    }
};

// ============================================================================
// PRIMEMatrix: Numerically Stabilized Property-Informed Codon Substitution Matrix
// ============================================================================

// PRIMEMatrix extends MG94Matrix with the property-informed rate model. It
// inherits the whole eigensystem layout (Q, pi, eigen*, U_pi, V_pi,
// scale_factor, transition_matrix) so it can be passed to likelihood routines
// as a const MG94Matrix& — a plain upcast, not a layout-dependent
// reinterpret_cast.
class PRIMEMatrix : public MG94Matrix {
public:
    // Updates Q and performs exact self-adjoint eigendecomposition
    void update(
        Scalar alpha,
        Scalar beta,
        const std::vector<Scalar>& lambdas,
        const std::vector<Matrix>& codon_prop_diffs, // [K] each (S, S)
        const GTRParameters& gtr,
        const Eigen::Matrix<Scalar, 3, 4>& pos_nuc_freqs,
        const Vector& codon_freqs,
        const GeneticCode& code = *GeneticCode::universal(),
        Scalar max_bound = 8.0
    ) {
        int S = code.num_sense_codons;
        size_t K = lambdas.size();
        pi = codon_freqs;
        sqrt_pi.resize(S);
        inv_sqrt_pi.resize(S);
        for (int i = 0; i < S; ++i) {
            sqrt_pi(i) = std::sqrt(std::max(pi(i), 1e-12));
            inv_sqrt_pi(i) = 1.0 / sqrt_pi(i);
        }

        Scalar r_nuc[4][4] = {
            {0.0, gtr.theta_AC, 1.0, gtr.theta_AT},
            {gtr.theta_AC, 0.0, gtr.theta_CG, gtr.theta_CT},
            {1.0, gtr.theta_CG, 0.0, gtr.theta_GT},
            {gtr.theta_AT, gtr.theta_CT, gtr.theta_GT, 0.0}
        };

        Q.resize(S, S);
        Q.setZero();

        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                const auto& diff = code.diff_matrix[i][j];
                if (diff.num_diffs != 1) continue;

                Scalar nuc_rate = r_nuc[diff.diff_from[0]][diff.diff_to[0]];
                Scalar target_nuc_freq = pos_nuc_freqs(diff.diff_pos[0], diff.diff_to[0]);

                Scalar rate_modifier;
                if (diff.synonymous) {
                    rate_modifier = alpha;
                } else {
                    // Property-informed non-synonymous rate modifier:
                    // psi = exp(-M * tanh(sum_k lambda_k * Delta_k / M))
                    Scalar eta = 0.0;
                    for (size_t k = 0; k < K; ++k) {
                        eta += lambdas[k] * codon_prop_diffs[k](i, j);
                    }
                    Scalar eta_eff = max_bound * std::tanh(eta / max_bound);
                    Scalar psi = std::exp(-eta_eff);
                    rate_modifier = beta * psi;
                }

                Q(i, j) = rate_modifier * nuc_rate * target_nuc_freq;
            }
        }

        // Row diagonals
        for (int i = 0; i < S; ++i) {
            Q(i, i) = -Q.row(i).sum();
        }

        // Average substitution rate
        Scalar avg_rate = 0.0;
        for (int i = 0; i < S; ++i) {
            avg_rate -= pi(i) * Q(i, i);
        }
        scale_factor = avg_rate;

        // Symmetric transform S = diag(sqrt(pi)) * Q * diag(1/sqrt(pi))
        Matrix S_mat = sqrt_pi.asDiagonal() * Q * inv_sqrt_pi.asDiagonal();
        S_mat = 0.5 * (S_mat + S_mat.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix> solver(S_mat);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("PRIMEMatrix eigendecomposition failed");
        }

        eigenvalues = solver.eigenvalues().cwiseMin(0.0);
        eigenvectors = solver.eigenvectors();
        inv_eigenvectors = eigenvectors.transpose();
        U_pi = inv_sqrt_pi.asDiagonal() * eigenvectors;
        V_pi = inv_eigenvectors * sqrt_pi.asDiagonal();
    }
};

// ============================================================================
// Data Structures for PRIME Results
// ============================================================================

struct PropertySiteInference {
    std::string property_name;
    Scalar lambda = 0.0;
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
    Scalar adjusted_p_value = 1.0;
    Scalar log_l_null = 0.0; // LogL when lambda_k = 0
    std::string classification = "Neutral"; // "Conserved", "Changing", "Neutral"
};

struct PRIMESiteResult {
    size_t site_index = 0;
    Scalar alpha = 0.0;
    Scalar beta = 0.0;
    Scalar beta_nuisance = 0.0; // on background branches
    Scalar fel_alpha = 0.0;
    Scalar fel_beta = 0.0;
    Scalar log_l_alt = 0.0;
    Scalar log_l_fel = 0.0;
    Scalar omnibus_lrt = 0.0;
    Scalar omnibus_p_value = 1.0;
    Scalar omnibus_q_value = 1.0;
    Scalar total_branch_length = 0.0;

    std::vector<PropertySiteInference> properties;

    // Ancestral & substitution diagnostics
    size_t nonsyn_subs = 0;
    size_t unique_aa = 0;
    Scalar redundancy_ratio = 0.0;
    SiteSubstitutions subs;

    bool is_significant(Scalar threshold = 0.10) const {
        return omnibus_q_value <= threshold;
    }
};

struct PRIMEGlobalResult {
    Scalar log_likelihood = 0.0;
    Scalar aicc = 0.0;
    Scalar alpha = 1.0;
    Scalar beta = 1.0;
    std::vector<Scalar> lambdas;
    std::vector<std::string> property_names;
    Scalar baseline_tree_length = 0.0;
    Scalar gtr_log_l = 0.0;
};

struct PRIMESettings {
    Scalar p_value_threshold = 0.10;
    Scalar stage2_p_value = 0.10;
    std::string property_set = "5PROP"; // "5PROP", "Atchley", "3PROP", "4PROP", "2PROP", or JSON file
    std::string custom_property_file = "";
    Scalar shrinkage_lambda = 0.5; // L2 Ridge regularization parameter
    Scalar max_bound = 8.0;        // Soft tanh bounding ceiling
    bool full_model = true;
    bool fit_global = true;        // Run G-PRIME before S-PRIME
    bool impute_states = false;
    std::string test_branches = "All"; // "All", "Internal", "Leaves", or specific tags
    std::string test_branch_regex = "";
    std::vector<std::string> test_branch_names = {};
};

struct PRIMEResult {
    PRIMESettings settings;
    std::vector<AminoAcidProperty> properties;
    PRIMEGlobalResult global_fit;
    std::vector<PRIMESiteResult> site_results;

    size_t total_sites = 0;
    size_t invariant_sites = 0;
    size_t significant_sites = 0;
    std::map<std::string, size_t> conserved_counts;
    std::map<std::string, size_t> changing_counts;

    Scalar runtime_seconds = 0.0;

    nlohmann::json to_legacy_json(const Tree& tree, const Alignment& aln, const std::string& input_filepath = "") const;
    nlohmann::json to_modern_json(const Tree& tree, const Alignment& aln, const std::string& cli_cmd = "hyphy3 prime", const std::string& input_filepath = "") const;
};

// ============================================================================
// PRIMEAnalyzer: Main Analysis Engine
// ============================================================================

class PRIMEAnalyzer {
public:
    Tree tree;
    Alignment aln;
    MG94Parameters base_params;
    PRIMESettings settings;

    std::vector<AminoAcidProperty> properties;
    size_t K = 0; // Number of properties

    // Precomputed codon pairwise property difference matrices: [K] each (S, S)
    std::vector<Matrix> codon_prop_diffs;

    // Branch partitioning:
    std::vector<bool> is_test_branch;
    std::vector<std::string> test_branch_names;
    std::vector<std::string> bg_branch_names;
    size_t num_tested = 0;
    size_t num_background = 0;

    Scalar baseline_tree_length = 0.0;
    Scalar baseline_tested_tree_length = 0.0;
    Scalar baseline_bg_tree_length = 0.0;
    Scalar conversion_factor = 1.0;

    bool has_gtr_fit = false;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;

    PRIMEGlobalResult global_fit;

    PRIMEAnalyzer(Tree t, Alignment a, MG94Parameters params = MG94Parameters{}, PRIMESettings s = PRIMESettings{})
        : tree(std::move(t)), aln(std::move(a)), base_params(params), settings(std::move(s)) {
        init_properties();
        prepare_branch_lengths();
        precompute_codon_differences();
    }

    static PRIMEAnalyzer create_and_fit(Tree tree, Alignment aln, const PRIMESettings& settings = PRIMESettings{}) {
        PRIMEAnalyzer analyzer(std::move(tree), std::move(aln), MG94Parameters{}, settings);
        analyzer.fit_nucleotide_and_codon_baselines();
        return analyzer;
    }

    void init_properties() {
        if (!settings.custom_property_file.empty()) {
            properties = PropertyDatabase::from_json_file(settings.custom_property_file);
        } else {
            properties = PropertyDatabase::get_predefined_properties(settings.property_set);
        }
        K = properties.size();
    }

    void prepare_branch_lengths() {
        size_t num_nodes = tree.num_nodes();
        is_test_branch.assign(num_nodes, true);
        test_branch_names.clear();
        bg_branch_names.clear();

        std::unordered_set<std::string> explicit_names(
            settings.test_branch_names.begin(),
            settings.test_branch_names.end()
        );
        std::regex b_regex;
        bool use_regex = !settings.test_branch_regex.empty();
        if (use_regex) {
            b_regex = std::regex(settings.test_branch_regex);
        }

        for (size_t i = 0; i < num_nodes; ++i) {
            const auto& node = tree.nodes[i];
            if (node.id == tree.root_id) {
                is_test_branch[i] = false;
                continue;
            }

            bool is_test = true;
            if (!explicit_names.empty()) {
                is_test = explicit_names.count(node.name) > 0;
            } else if (use_regex) {
                is_test = std::regex_search(node.name, b_regex);
            } else if (settings.test_branches == "Internal") {
                is_test = !node.is_leaf;
            } else if (settings.test_branches == "Leaves") {
                is_test = node.is_leaf;
            } else if (settings.test_branches != "All") {
                is_test = (node.name.find(settings.test_branches) != std::string::npos);
            }

            is_test_branch[i] = is_test;
            if (is_test) {
                test_branch_names.push_back(node.name);
            } else {
                bg_branch_names.push_back(node.name);
            }
        }

        num_tested = test_branch_names.size();
        num_background = bg_branch_names.size();

        baseline_tree_length = 0.0;
        baseline_tested_tree_length = 0.0;
        baseline_bg_tree_length = 0.0;
        for (size_t i = 0; i < num_nodes; ++i) {
            if (tree.nodes[i].id != tree.root_id) {
                Scalar bl = tree.nodes[i].branch_length;
                baseline_tree_length += bl;
                if (is_test_branch[i]) {
                    baseline_tested_tree_length += bl;
                } else {
                    baseline_bg_tree_length += bl;
                }
            }
        }
    }

    void precompute_codon_differences() {
        const auto& code = *(aln.code ? aln.code : GeneticCode::universal());
        int S = code.num_sense_codons;
        codon_prop_diffs.resize(K);

        for (size_t k = 0; k < K; ++k) {
            codon_prop_diffs[k].resize(S, S);
            codon_prop_diffs[k].setZero();
            for (int u = 0; u < S; ++u) {
                char aa_u = code.sense_codon_aa[u];
                int idx_u = aa_char_to_index(aa_u);
                for (int v = 0; v < S; ++v) {
                    char aa_v = code.sense_codon_aa[v];
                    int idx_v = aa_char_to_index(aa_v);
                    if (idx_u >= 0 && idx_v >= 0) {
                        Scalar diff = std::abs(properties[k].z_scores[idx_u] - properties[k].z_scores[idx_v]);
                        codon_prop_diffs[k](u, v) = diff;
                    }
                }
            }
        }
    }

    void fit_nucleotide_and_codon_baselines() {
        // Step 1: Fit GTR
        GTRFitter gtr_fitter(tree, aln);
        auto gtr_res = gtr_fitter.fit();
        has_gtr_fit = true;
        gtr_rates = gtr_res.params;
        gtr_log_l = gtr_res.log_likelihood;
        gtr_aicc = gtr_res.aicc;

        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                gtr_branch_lengths[node.name] = node.branch_length;
            }
        }

        // Update tree branch lengths from GTR fit
        for (auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                auto it = gtr_branch_lengths.find(node.name);
                if (it != gtr_branch_lengths.end()) {
                    node.branch_length = it->second;
                }
            }
        }
        prepare_branch_lengths();

        // Step 2: Global MG94 Baseline Fit
        MG94Parameters base_p;
        base_p.theta_AC = gtr_rates.theta_AC;
        base_p.theta_AT = gtr_rates.theta_AT;
        base_p.theta_CG = gtr_rates.theta_CG;
        base_p.theta_CT = gtr_rates.theta_CT;
        base_p.theta_GT = gtr_rates.theta_GT;

        MG94Fitter mg_fitter(tree, aln);
        auto mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
        base_p.alpha = 1.0;
        base_p.beta = mg_res.x_opt(0);
        conversion_factor = (mg_res.x_opt.size() >= 2) ? mg_res.x_opt(1) : 1.0;
        tree = mg_res.tree;
        prepare_branch_lengths();
        global_fit.gtr_log_l = gtr_log_l;

        // Step 3: Global G-PRIME Fit (if enabled)
        if (settings.fit_global) {
            fit_global_g_prime();
        }
    }

    void fit_global_g_prime() {
        const auto& code = *(aln.code ? aln.code : GeneticCode::universal());
        int S = code.num_sense_codons;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();

        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

        // Parameter vector: [omega, lambda_0, ..., lambda_{K-1}]
        size_t dim = 1 + K;
        std::vector<Scalar> init_x(dim, 0.0);
        init_x[0] = std::clamp(base_params.beta, 0.05, 5.0);

        std::vector<Scalar> lb(dim, -5.0);
        std::vector<Scalar> ub(dim, 5.0);
        lb[0] = 0.01; ub[0] = 20.0;

        auto eval_alignment_lnl = [&](const std::vector<Scalar>& x) -> Scalar {
            Scalar omega = x[0];
            std::vector<Scalar> lams(K);
            for (size_t k = 0; k < K; ++k) lams[k] = x[1 + k];

            PRIMEMatrix pm;
            pm.update(1.0, omega, lams, codon_prop_diffs, gtr_rates, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, code, settings.max_bound);

            auto P_branches = LikelihoodEngine::build_P_branches<Matrix>(tree, pm);

            Scalar total_lnl = 0.0;
            #pragma omp parallel
            {
                auto [thread_node_L, thread_node_scale] = LikelihoodEngine::make_pattern_scratch<Vector>(num_nodes, S);

                #pragma omp for reduction(+:total_lnl) schedule(dynamic, 8)
                for (size_t p = 0; p < num_patterns; ++p) {
                    Scalar lnl = LikelihoodEngine::compute_pattern_log_likelihood(
                        tree, aln.patterns[p], leaf_to_taxon,
                        // Passing dummy MG94 with same pi
                        static_cast<const MG94Matrix&>(pm),
                        P_branches, thread_node_L, thread_node_scale
                    );
                    total_lnl += aln.patterns[p].weight * lnl;
                }
            }
            return total_lnl;
        };

        auto obj = [&](const std::vector<Scalar>& x) -> Scalar {
            return -eval_alignment_lnl(x);
        };

        // Quick optimization using Nelder-Mead Dynamic
        auto pt = NelderMeadDynamic::minimize(obj, init_x, lb, ub, 1e-3, 60);

        global_fit.beta = pt.x[0];
        global_fit.alpha = 1.0;
        global_fit.lambdas.assign(pt.x.begin() + 1, pt.x.end());
        global_fit.log_likelihood = -pt.f;
        global_fit.aicc = -2.0 * global_fit.log_likelihood + 2.0 * dim;
        global_fit.baseline_tree_length = baseline_tree_length;
        for (const auto& p : properties) {
            global_fit.property_names.push_back(p.name);
        }
    }

    using SiteProgressCallback = std::function<void(size_t site_idx, const PRIMESiteResult& res)>;
    using LiveSiteProgressCallback = std::function<void(size_t site_idx, size_t completed_sites, size_t total_sites, const PRIMESiteResult& res)>;

    // Runs Site-Level Property Inference (S-PRIME)
    PRIMEResult run(
        ProgressBar* /* progress */ = nullptr,
        bool show_progress = true,
        bool force_progress = false,
        SiteProgressCallback site_callback = nullptr,
        LiveSiteProgressCallback live_site_callback = nullptr
    ) {
        auto start_time = std::chrono::high_resolution_clock::now();
        size_t num_patterns = aln.patterns.size();
        size_t num_codons = aln.num_codons;
        size_t num_nodes = tree.num_nodes();
        const auto& code = *(aln.code ? aln.code : GeneticCode::universal());
        int S = code.num_sense_codons;

        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

        // First: Compute site-level FEL baseline for all patterns
        FELAnalyzer fel_analyzer(tree, aln, base_params);
        fel_analyzer.gtr_rates = gtr_rates;
        fel_analyzer.conversion_factor = conversion_factor;
        fel_analyzer.prepare_branch_lengths();

        std::vector<SiteResult> fel_pattern_results(num_patterns);
        for (size_t p = 0; p < num_patterns; ++p) {
            fel_pattern_results[p] = fel_analyzer.analyze_pattern(p, leaf_to_taxon);
        }

        SiteVisualizer viz(num_codons, num_patterns, SiteVisualizerConfig{
            .method_name = "S-PRIME",
            .p_value_threshold = settings.p_value_threshold,
            .pos_label = "changing",
            .neg_label = "conserved",
            .sig_label = "omnibus sig",
            .quiet = !show_progress,
            .force_tty = force_progress,
            .enable_chart = true
        });

        std::vector<PRIMESiteResult> pattern_prime_results(num_patterns);
        std::atomic<size_t> completed_sites{0};
        std::mutex callback_mutex;

        #pragma omp parallel
        {
            auto [thread_node_L, thread_node_scale] = LikelihoodEngine::make_pattern_scratch<Vector>(num_nodes, S);
            std::vector<Matrix> P_branches(num_nodes);

            #pragma omp for schedule(dynamic, 4)
            for (size_t p = 0; p < num_patterns; ++p) {
                pattern_prime_results[p] = analyze_pattern(
                    p, leaf_to_taxon, fel_pattern_results[p],
                    thread_node_L, thread_node_scale, P_branches
                );

                const auto& pat_res = pattern_prime_results[p];
                size_t w = aln.patterns[p].weight;
                size_t cur_done = completed_sites.fetch_add(w, std::memory_order_relaxed) + w;

                SiteSelectionType type = SiteSelectionType::Neutral;
                std::string detail;
                if (pat_res.omnibus_p_value <= settings.p_value_threshold) {
                    std::string cons_names, chg_names;
                    for (const auto& pr : pat_res.properties) {
                        if (pr.p_value <= settings.stage2_p_value) {
                            std::string short_name = pr.property_name.substr(0, std::min<size_t>(5, pr.property_name.size()));
                            if (pr.lambda > 0.0) {
                                cons_names += (cons_names.empty() ? "" : ",") + short_name;
                            } else if (pr.lambda < 0.0) {
                                chg_names += (chg_names.empty() ? "" : ",") + short_name;
                            }
                        }
                    }
                    if (!chg_names.empty() && cons_names.empty()) {
                        type = SiteSelectionType::Positive;
                        detail = "▲ Chg:" + chg_names + ", p=" + SiteVisualizer::format_pvalue(pat_res.omnibus_p_value);
                    } else if (!cons_names.empty() && chg_names.empty()) {
                        type = SiteSelectionType::Negative;
                        detail = "▼ Cons:" + cons_names + ", p=" + SiteVisualizer::format_pvalue(pat_res.omnibus_p_value);
                    } else {
                        type = SiteSelectionType::Significant;
                        detail = "★ p=" + SiteVisualizer::format_pvalue(pat_res.omnibus_p_value);
                        if (!cons_names.empty()) detail += " (Cons:" + cons_names + ")";
                        if (!chg_names.empty()) detail += " (Chg:" + chg_names + ")";
                    }
                }

                viz.record_pattern(p, aln.patterns[p].sites, type, pat_res.omnibus_p_value, pat_res.omnibus_lrt, detail);

                if (live_site_callback) {
                    std::lock_guard<std::mutex> lock(callback_mutex);
                    for (size_t site : aln.patterns[p].sites) {
                        PRIMESiteResult site_res = pat_res;
                        site_res.site_index = site;
                        live_site_callback(site, cur_done, num_codons, site_res);
                    }
                }
            }
        }

        viz.finish();

        // Map pattern results to site results
        std::vector<PRIMESiteResult> site_results(num_codons);
        for (size_t site = 0; site < num_codons; ++site) {
            size_t p_idx = aln.site_to_pattern[site];
            site_results[site] = pattern_prime_results[p_idx];
            site_results[site].site_index = site;

            // Ancestral substitutions for this site
            site_results[site].subs = fel_analyzer.reconstruct_ancestral_pattern(p_idx, leaf_to_taxon);
            site_results[site].subs.site_index = site;

            // Count unique amino acids and non-synonymous substitutions
            std::unordered_set<char> observed_aas;
            for (int8_t st : aln.patterns[p_idx].states) {
                if (st >= 0 && st < S) {
                    observed_aas.insert(code.sense_codon_aa[st]);
                }
            }
            site_results[site].unique_aa = observed_aas.size();
            site_results[site].nonsyn_subs = site_results[site].subs.nonsyn_test + site_results[site].subs.nonsyn_bg;
            if (site_results[site].unique_aa > 0) {
                site_results[site].redundancy_ratio = static_cast<Scalar>(site_results[site].nonsyn_subs) /
                                                     static_cast<Scalar>(site_results[site].unique_aa);
            }
        }

        // Hierarchical False Discovery Rate (HRFDR) Corrections
        // Stage 1: Benjamini-Hochberg FDR on Omnibus p-values
        std::vector<Scalar> omnibus_p_vals(num_codons);
        for (size_t site = 0; site < num_codons; ++site) {
            omnibus_p_vals[site] = site_results[site].omnibus_p_value;
        }
        std::vector<Scalar> omnibus_q_vals = hyphy::core::stats::benjamini_hochberg_fdr(omnibus_p_vals);
        for (size_t site = 0; site < num_codons; ++site) {
            site_results[site].omnibus_q_value = omnibus_q_vals[site];
        }

        // Stage 2: Holm-Bonferroni correction on individual property p-values
        size_t sig_site_count = 0;
        for (size_t site = 0; site < num_codons; ++site) {
            if (site_results[site].is_significant(settings.p_value_threshold)) {
                sig_site_count++;
            }
        }

        Scalar stage2_alpha = settings.stage2_p_value;
        if (num_codons > 0 && sig_site_count > 0) {
            // Hierarchical FDR threshold
            stage2_alpha = settings.p_value_threshold * (static_cast<Scalar>(sig_site_count) / static_cast<Scalar>(num_codons));
        }

        for (size_t site = 0; site < num_codons; ++site) {
            auto& sr = site_results[site];
            std::vector<Scalar> prop_p_vals(K);
            for (size_t k = 0; k < K; ++k) {
                prop_p_vals[k] = sr.properties[k].p_value;
            }
            std::vector<Scalar> prop_adj_p = hyphy::core::stats::holm_bonferroni(prop_p_vals);
            for (size_t k = 0; k < K; ++k) {
                sr.properties[k].adjusted_p_value = prop_adj_p[k];
            }

            // Classify each property
            for (size_t k = 0; k < K; ++k) {
                auto& p_inf = sr.properties[k];
                // Significant if omnibus passes Stage 1 AND property adjusted p-value <= stage2_alpha
                if (sr.is_significant(settings.p_value_threshold) && p_inf.adjusted_p_value <= stage2_alpha) {
                    if (p_inf.lambda > 0.0) {
                        p_inf.classification = "Conserved";
                    } else if (p_inf.lambda < 0.0) {
                        p_inf.classification = "Changing";
                    } else {
                        p_inf.classification = "Neutral";
                    }
                } else {
                    p_inf.classification = "Neutral";
                }
            }
        }

        // Invoke site_callback with finalized site results
        if (site_callback) {
            for (size_t site = 0; site < num_codons; ++site) {
                site_callback(site, site_results[site]);
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        Scalar elapsed = std::chrono::duration<Scalar>(end_time - start_time).count();

        // Compile Final Result Summary
        PRIMEResult res;
        res.settings = settings;
        res.properties = properties;
        res.global_fit = global_fit;
        res.site_results = std::move(site_results);
        res.total_sites = num_codons;
        res.runtime_seconds = elapsed;

        for (const auto& p : properties) {
            res.conserved_counts[p.name] = 0;
            res.changing_counts[p.name] = 0;
        }

        for (const auto& sr : res.site_results) {
            if (sr.alpha == 0.0 && sr.beta == 0.0 && sr.omnibus_lrt == 0.0) {
                res.invariant_sites++;
            }
            if (sr.is_significant(settings.p_value_threshold)) {
                res.significant_sites++;
                for (const auto& p_inf : sr.properties) {
                    if (p_inf.classification == "Conserved") {
                        res.conserved_counts[p_inf.property_name]++;
                    } else if (p_inf.classification == "Changing") {
                        res.changing_counts[p_inf.property_name]++;
                    }
                }
            }
        }

        return res;
    }

private:
    PRIMESiteResult analyze_pattern(
        size_t p_idx,
        const std::vector<size_t>& leaf_to_taxon,
        const SiteResult& fel_res,
        std::vector<Vector>& thread_node_L,
        std::vector<double>& thread_node_scale,
        std::vector<Matrix>& P_branches
    ) const {
        const auto& pattern = aln.patterns[p_idx];
        const auto& code = *(aln.code ? aln.code : GeneticCode::universal());
        size_t num_nodes = tree.num_nodes();
        bool has_bg = (num_background > 0);

        PRIMESiteResult res;
        res.site_index = p_idx;
        res.fel_alpha = fel_res.alpha;
        res.fel_beta = fel_res.beta;
        res.log_l_fel = fel_res.log_l_alt;
        res.properties.resize(K);
        for (size_t k = 0; k < K; ++k) {
            res.properties[k].property_name = properties[k].name;
        }

        // Fast path: Invariant site or no non-synonymous substitutions
        if ((fel_res.alpha == 0.0 && fel_res.beta == 0.0) || fel_res.beta <= 1e-4) {
            res.alpha = fel_res.alpha;
            res.beta = fel_res.beta;
            res.beta_nuisance = fel_res.beta_nuisance;
            res.log_l_alt = fel_res.log_l_alt;
            res.omnibus_lrt = 0.0;
            res.omnibus_p_value = 1.0;
            res.omnibus_q_value = 1.0;
            res.total_branch_length = 0.0;
            for (size_t k = 0; k < K; ++k) {
                res.properties[k].lambda = 0.0;
                res.properties[k].lrt = 0.0;
                res.properties[k].p_value = 1.0;
                res.properties[k].log_l_null = fel_res.log_l_alt;
            }
            return res;
        }

        // Likelihood evaluation function for pattern
        auto eval_pattern_lnl = [&](Scalar a, Scalar b, Scalar b_bg, const std::vector<Scalar>& lams) -> Scalar {
            PRIMEMatrix pm;
            pm.update(a, b, lams, codon_prop_diffs, gtr_rates, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, code, settings.max_bound);

            PRIMEMatrix pm_bg;
            if (has_bg) {
                // Background branches use standard neutral rates (lambdas = 0)
                std::vector<Scalar> zero_lams(K, 0.0);
                pm_bg.update(a, b_bg, zero_lams, codon_prop_diffs, gtr_rates, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, code, settings.max_bound);
            }

            for (size_t i = 0; i < num_nodes; ++i) {
                if (tree.nodes[i].id != tree.root_id) {
                    if (is_test_branch[i] || !has_bg) {
                        P_branches[i] = pm.transition_matrix(tree.nodes[i].branch_length);
                    } else {
                        P_branches[i] = pm_bg.transition_matrix(tree.nodes[i].branch_length);
                    }
                }
            }

            return LikelihoodEngine::compute_pattern_log_likelihood(
                tree, pattern, leaf_to_taxon,
                static_cast<const MG94Matrix&>(pm),
                P_branches, thread_node_L, thread_node_scale
            );
        };

        // L2 Ridge Shrinkage Penalty: (gamma / 2) * sum(lambda_k^2)
        auto calc_shrinkage_penalty = [&](const std::vector<Scalar>& lams) -> Scalar {
            if (settings.shrinkage_lambda <= 0.0) return 0.0;
            Scalar sum_sq = 0.0;
            for (Scalar val : lams) {
                sum_sq += val * val;
            }
            return 0.5 * settings.shrinkage_lambda * sum_sq;
        };

        // --------------------------------------------------------------------
        // 1. Fit Alternative Model (Nelder-Mead ND)
        // Parameter vector x: [alpha, beta_test, (beta_bg), lambda_0, ..., lambda_{K-1}]
        // --------------------------------------------------------------------
        size_t alt_dim = 2 + (has_bg ? 1 : 0) + K;
        std::vector<Scalar> lb_alt(alt_dim, -10.0);
        std::vector<Scalar> ub_alt(alt_dim, 10.0);
        lb_alt[0] = 1e-4; ub_alt[0] = 100.0;  // alpha
        lb_alt[1] = 1e-4; ub_alt[1] = 1000.0; // beta
        if (has_bg) {
            lb_alt[2] = 1e-4; ub_alt[2] = 1000.0; // beta_bg
        }

        auto obj_alt = [&](const std::vector<Scalar>& x) -> Scalar {
            Scalar a = x[0];
            Scalar b = x[1];
            Scalar b_bg = has_bg ? x[2] : 1.0;
            size_t lam_start = 2 + (has_bg ? 1 : 0);
            std::vector<Scalar> lams(K);
            for (size_t k = 0; k < K; ++k) lams[k] = x[lam_start + k];

            Scalar lnl = eval_pattern_lnl(a, b, b_bg, lams);
            Scalar penalty = calc_shrinkage_penalty(lams);
            return -(lnl - penalty);
        };

        // Simplex Grid Initialization
        std::vector<std::vector<Scalar>> grid_pts;

        // Point 0: The nested FEL point (lambdas = 0)
        std::vector<Scalar> pt_fel(alt_dim, 0.0);
        pt_fel[0] = std::clamp(fel_res.alpha, 1e-4, 100.0);
        pt_fel[1] = std::clamp(fel_res.beta, 1e-4, 1000.0);
        if (has_bg) pt_fel[2] = std::clamp(fel_res.beta_nuisance, 1e-4, 1000.0);
        grid_pts.push_back(pt_fel);

        // Point 1: Global G-PRIME parameters
        if (!global_fit.lambdas.empty() && global_fit.lambdas.size() == K) {
            std::vector<Scalar> pt_glob = pt_fel;
            size_t ls = 2 + (has_bg ? 1 : 0);
            for (size_t k = 0; k < K; ++k) pt_glob[ls + k] = global_fit.lambdas[k];
            grid_pts.push_back(pt_glob);
        }

        // Coordinate test vertices: +/- 1.0 on each lambda
        size_t ls = 2 + (has_bg ? 1 : 0);
        for (size_t k = 0; k < K; ++k) {
            std::vector<Scalar> p_pos = pt_fel;
            p_pos[ls + k] = 1.0;
            grid_pts.push_back(p_pos);

            std::vector<Scalar> p_neg = pt_fel;
            p_neg[ls + k] = -1.0;
            grid_pts.push_back(p_neg);
        }

        std::vector<Scalar> init_alt = pt_fel;
        Scalar best_score = obj_alt(pt_fel);
        for (const auto& g : grid_pts) {
            Scalar sc = obj_alt(g);
            if (sc < best_score) {
                best_score = sc;
                init_alt = g;
            }
        }

        auto pt_alt = NelderMeadDynamic::minimize(obj_alt, init_alt, lb_alt, ub_alt, 1e-4, 120, grid_pts);

        res.alpha = pt_alt.x[0];
        res.beta = pt_alt.x[1];
        if (has_bg) res.beta_nuisance = pt_alt.x[2];
        std::vector<Scalar> mle_lambdas(K);
        for (size_t k = 0; k < K; ++k) {
            mle_lambdas[k] = pt_alt.x[ls + k];
            res.properties[k].lambda = mle_lambdas[k];
        }

        Scalar raw_alt_lnl = eval_pattern_lnl(res.alpha, res.beta, res.beta_nuisance, mle_lambdas);
        res.log_l_alt = raw_alt_lnl;

        // Guaranteed nesting: Alt can NEVER be worse than FEL (which is lambda = 0)
        if (res.log_l_alt < res.log_l_fel) {
            res.log_l_alt = res.log_l_fel;
            res.alpha = fel_res.alpha;
            res.beta = fel_res.beta;
            for (size_t k = 0; k < K; ++k) {
                res.properties[k].lambda = 0.0;
                mle_lambdas[k] = 0.0;
            }
        }

        // Omnibus LRT & p-value vs FEL (df = K)
        bool valid_omnibus = (res.log_l_alt > -1e10 && res.log_l_fel > -1e10 &&
                              std::isfinite(res.log_l_alt) && std::isfinite(res.log_l_fel));
        Scalar omnibus_lrt = 0.0;
        if (valid_omnibus) {
            omnibus_lrt = 2.0 * (res.log_l_alt - res.log_l_fel);
            if (omnibus_lrt < 0.0) omnibus_lrt = 0.0;
        }
        res.omnibus_lrt = omnibus_lrt;

        // Chi-squared p-value for df = K
        if (!valid_omnibus) {
            res.omnibus_p_value = 1.0;
        } else if (K == 1) {
            res.omnibus_p_value = std::erfc(std::sqrt(omnibus_lrt / 2.0));
        } else if (K == 2) {
            res.omnibus_p_value = std::exp(-0.5 * omnibus_lrt);
        } else {
            // Wilson-Hilferty transformation for general df
            Scalar s_df = static_cast<Scalar>(K);
            Scalar z = std::pow(omnibus_lrt / s_df, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * s_df));
            z /= std::sqrt(2.0 / (9.0 * s_df));
            Scalar p = 0.5 * std::erfc(z / std::sqrt(2.0));
            res.omnibus_p_value = std::clamp(p, 0.0, 1.0);
        }

        // --------------------------------------------------------------------
        // 2. Fit Constrained Null Model for each Property k (lambda_k = 0)
        // --------------------------------------------------------------------
        for (size_t k = 0; k < K; ++k) {
            // Null model parameter vector: [alpha, beta, (beta_bg), lambda_{j != k}]
            // Dimension: alt_dim - 1
            size_t null_dim = alt_dim - 1;
            std::vector<Scalar> lb_null(null_dim, -10.0);
            std::vector<Scalar> ub_null(null_dim, 10.0);
            lb_null[0] = 1e-4; ub_null[0] = 100.0;
            lb_null[1] = 1e-4; ub_null[1] = 1000.0;
            if (has_bg) {
                lb_null[2] = 1e-4; ub_null[2] = 1000.0;
            }

            auto unpack_null_lams = [&](const std::vector<Scalar>& x) -> std::vector<Scalar> {
                std::vector<Scalar> lams(K, 0.0);
                size_t src_idx = 2 + (has_bg ? 1 : 0);
                for (size_t j = 0; j < K; ++j) {
                    if (j == k) {
                        lams[j] = 0.0;
                    } else {
                        lams[j] = x[src_idx++];
                    }
                }
                return lams;
            };

            auto obj_null_k = [&](const std::vector<Scalar>& x) -> Scalar {
                Scalar a = x[0];
                Scalar b = x[1];
                Scalar b_bg = has_bg ? x[2] : 1.0;
                auto lams = unpack_null_lams(x);

                Scalar lnl = eval_pattern_lnl(a, b, b_bg, lams);
                Scalar penalty = calc_shrinkage_penalty(lams);
                return -(lnl - penalty);
            };

            // Start vector 0: MLE from alternative with lambda_k clamped to 0
            std::vector<Scalar> init_null(null_dim, 0.0);
            init_null[0] = res.alpha;
            init_null[1] = res.beta;
            if (has_bg) init_null[2] = res.beta_nuisance;
            size_t dst_idx = 2 + (has_bg ? 1 : 0);
            for (size_t j = 0; j < K; ++j) {
                if (j != k) init_null[dst_idx++] = mle_lambdas[j];
            }

            // Start vector 1: All lambdas = 0 (FEL point)
            std::vector<Scalar> init_null_fel(null_dim, 0.0);
            init_null_fel[0] = pt_fel[0];
            init_null_fel[1] = pt_fel[1];
            if (has_bg) init_null_fel[2] = pt_fel[2];

            std::vector<std::vector<Scalar>> null_grid = {init_null, init_null_fel};
            auto pt_null_k = NelderMeadDynamic::minimize(obj_null_k, init_null, lb_null, ub_null, 1e-4, 80, null_grid);

            auto final_null_lams = unpack_null_lams(pt_null_k.x);
            Scalar raw_null_lnl = eval_pattern_lnl(pt_null_k.x[0], pt_null_k.x[1], has_bg ? pt_null_k.x[2] : 1.0, final_null_lams);
            res.properties[k].log_l_null = raw_null_lnl;

            // Strict nesting enforcement: Null logL cannot exceed Alt logL
            if (res.properties[k].log_l_null > res.log_l_alt) {
                res.log_l_alt = res.properties[k].log_l_null;
                res.properties[k].lambda = 0.0;
            }

            bool valid_prop = (res.log_l_alt > -1e10 && res.properties[k].log_l_null > -1e10 &&
                               std::isfinite(res.log_l_alt) && std::isfinite(res.properties[k].log_l_null));
            Scalar prop_lrt = 0.0;
            if (valid_prop) {
                prop_lrt = 2.0 * (res.log_l_alt - res.properties[k].log_l_null);
                if (prop_lrt < 0.0) prop_lrt = 0.0;
            }
            res.properties[k].lrt = prop_lrt;
            res.properties[k].p_value = valid_prop ? std::erfc(std::sqrt(prop_lrt / 2.0)) : 1.0;
        }

        res.total_branch_length = res.alpha * baseline_tested_tree_length;
        return res;
    }
};

// ============================================================================
// JSON Serialization Implementations
// ============================================================================

} // namespace hyphy::analyses
