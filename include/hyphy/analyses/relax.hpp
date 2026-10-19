#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "nlohmann/json.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <regex>
#include <memory>
#include <limits>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

struct RELAXRateDistribution {
    std::vector<Scalar> omegas;
    std::vector<Scalar> weights;
};

struct RELAXModelFit {
    Scalar log_likelihood = 0.0;
    Scalar aicc = 0.0;
    size_t parameters = 0;
    Scalar tree_scale = 1.0;
    Scalar k = 1.0;
    RELAXRateDistribution reference_distribution;
    RELAXRateDistribution test_distribution;
    std::vector<Scalar> branch_lengths;
    std::vector<Scalar> site_log_likelihoods;
};

// How branch lengths are mapped onto the per-class MG94 rate matrices.
//  SharedMixture : (default) HyPhy 2.5 BS-REL semantics. All omega classes on a branch share one
//                  synonymous rate (one conversion factor per branch partition R/T/U),
//                  chosen so that the branch length equals the mixture-expected number of
//                  substitutions: conv_part = 3 / sum_c w_c * mu_{part,c}.
//  PerClass      : legacy hyphy-3 behavior, kept for comparison. Each omega class c is
//                  normalized by its own mean rate (conv_c = 3 / mu_c), so every class implies
//                  the same expected substitutions per branch.
enum class RELAXNormalization {
    PerClass,
    SharedMixture
};

// Level at which the omega mixture is applied.
//  BranchSite : (default) every branch-site combination draws its class independently, i.e. the
//               branch transition matrix is P_b = sum_c w_c exp(conv_c * Q_c * t_b), as in HyPhy 2.5
//               BS-REL (models.generic mixture "Exp(Q1)*w1 + Exp(Q2)*w2 + ...") and hyphy-3 BUSTED.
//  Site       : legacy hyphy-3 behavior, kept for comparison. A site draws one omega class used
//               on every branch.
enum class RELAXMixtureLevel {
    Site,
    BranchSite
};

struct RELAXSettings {
    Scalar p_value_threshold = 0.05;
    RELAXNormalization normalization = RELAXNormalization::SharedMixture;
    RELAXMixtureLevel mixture_level = RELAXMixtureLevel::BranchSite;
    // Branch-set specifications. Resolution order (see RELAXAnalyzer::resolve_branch_spec):
    // explicit name list -> keyword (all/internal/leaves/unlabeled) -> tag present in the tree ->
    // comma-separated exact branch names -> regex (explicit via *_is_regex or a "regex:" prefix;
    // otherwise only as a fallback, with a warning).
    std::string test_branch_regex = "";
    std::vector<std::string> test_branch_names = {};
    bool test_is_regex = false;
    std::string reference_branch_regex = "";
    std::vector<std::string> reference_branch_names = {};
    bool reference_is_regex = false;
    bool refine_branch_lengths = true;
    // Alternating (distribution <-> branch lengths) optimization after the tree-scale fit.
    int max_joint_cycles = 5;
    Scalar joint_tolerance = 1e-3;
    bool verbose = false;
};

struct RELAXResult {
    RELAXSettings settings;

    // Phase 1: Nucleotide GTR
    Scalar gtr_log_likelihood = 0.0;
    Scalar gtr_aicc = 0.0;
    size_t gtr_parameters = 0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    // Phase 2: MG94 with separate rates for branch sets
    Scalar mg94_log_likelihood = 0.0;
    Scalar mg94_aicc = 0.0;
    size_t mg94_parameters = 0;
    Scalar mg94_omega_R = 1.0;
    Scalar mg94_omega_T = 1.0;
    std::unordered_map<std::string, Scalar> mg94_branch_lengths;

    // Phase 3 & 4: RELAX fits
    RELAXModelFit alternative_fit;
    RELAXModelFit null_fit;

    // Phase 5: Hypothesis testing
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
    Scalar k = 1.0;
    bool is_relaxed = false;
    bool is_intensified = false;
    bool is_significant = false;
    bool is_failed = false;
    std::string error_message = "";

    // Branch classification
    std::vector<std::string> test_branches;
    std::vector<std::string> reference_branches;
    std::unordered_map<std::string, std::string> branch_class; // name -> "Test" / "Reference"

    double runtime_seconds = 0.0;

nlohmann::json to_legacy_json(const Tree& tree, const Alignment& aln) const;

nlohmann::json to_modern_json(const Tree& tree, const Alignment& aln, const Provenance& prov = {}) const;
    nlohmann::json to_json(
        const Tree& tree,
        const Alignment& aln,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        if (format == JSONFormat::Legacy) {
            return to_legacy_json(tree, aln);
        }
        return to_modern_json(tree, aln, prov);
    }

    void save_json(
        const std::string& filepath,
        const Tree& tree,
        const Alignment& aln,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open file for writing: " + filepath);
        }
        nlohmann::json j = to_json(tree, aln, format, prov);
        out << j.dump(2) << "\n";
    }

};

class RELAXAnalyzer {
public:
    Tree tree;
    Alignment aln;
    RELAXSettings settings;
    std::shared_ptr<const GeneticCode> code;
    MG94Parameters base_params;

    std::vector<bool> is_test_branch;
    std::vector<int> branch_partition; // 0 = Reference, 1 = Test, 2 = Unclassified
    std::vector<size_t> leaf_to_taxon;
    std::vector<std::string> test_branch_names;
    std::vector<std::string> ref_branch_names;
    std::vector<std::string> unclassified_branch_names;
    std::unordered_map<std::string, std::string> branch_class;
    Scalar unclassified_omega = 1.0;

    RELAXAnalyzer(
        Tree t,
        Alignment a,
        RELAXSettings s = RELAXSettings{}
    ) : tree(std::move(t)), aln(std::move(a)), settings(std::move(s)) {
        code = aln.code ? aln.code : GeneticCode::universal();
        size_t num_nodes = tree.num_nodes();
        is_test_branch.assign(num_nodes, false);
        branch_partition.assign(num_nodes, 0);
        leaf_to_taxon.assign(num_nodes, static_cast<size_t>(-1));

        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        partition_branches();
    }

    static RELAXAnalyzer create(
        Tree tree,
        Alignment aln,
        RELAXSettings settings = RELAXSettings{}
    ) {
        return RELAXAnalyzer(std::move(tree), std::move(aln), std::move(settings));
    }

    // Resolves a --test / --reference specification to a per-node selection (root never selected).
    // Precedence:
    //  0. explicit name list (settings.*_branch_names): exact (case-insensitive) branch names or tags
    //  1. "regex:<re>" prefix or force_regex  -> case-insensitive regex on branch names and tags
    //  2. keyword: all | internal | leaves/terminal/tips | unlabeled/untagged
    //  3. a tag present in the tree (case-insensitive, braces optional). The RELAX aliases
    //     {T,TEST,FG,FOREGROUND} and {R,REF,REFERENCE,BG,BACKGROUND} are treated as one tag group.
    //  4. comma-separated list of exact branch names (case-sensitive; tags/keywords allowed as tokens)
    //  5. fallback: case-insensitive regex, with a warning reporting the number of matched branches.
    std::vector<bool> resolve_branch_spec(
        const std::string& spec_in,
        const std::vector<std::string>& names,
        bool force_regex,
        const char* label
    ) const {
        auto to_lower = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
            return s;
        };
        auto trim_ws = [](std::string s) {
            while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
            while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
            return s;
        };
        auto trim_tag = [&](std::string s) {
            s = trim_ws(s);
            while (!s.empty() && (s.front() == '{' || std::isspace(static_cast<unsigned char>(s.front())))) s.erase(s.begin());
            while (!s.empty() && (s.back() == '}' || std::isspace(static_cast<unsigned char>(s.back())))) s.pop_back();
            return s;
        };
        auto alias_group = [](const std::string& tl) -> int {
            if (tl == "t" || tl == "test" || tl == "fg" || tl == "foreground") return 1;
            if (tl == "r" || tl == "ref" || tl == "reference" || tl == "bg" || tl == "background") return 2;
            return 0;
        };

        std::vector<bool> sel(tree.num_nodes(), false);
        auto count = [&]() { return static_cast<size_t>(std::count(sel.begin(), sel.end(), true)); };
        auto node_tag = [&](const TreeNode& n) { return to_lower(trim_tag(n.model_tag)); };

        // 0. explicit name list
        if (!names.empty()) {
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                std::string nl = to_lower(node.name), tl = node_tag(node);
                for (const auto& nm : names) {
                    std::string q = to_lower(trim_tag(nm));
                    if (q == nl || (!tl.empty() && q == tl)) { sel[node.id] = true; break; }
                }
            }
            return sel;
        }

        std::string spec = trim_ws(spec_in);
        if (spec.empty()) return sel;

        auto apply_regex = [&](const std::string& re) -> bool {
            try {
                std::regex rx(re, std::regex::icase);
                for (const auto& node : tree.nodes) {
                    if (node.id == tree.root_id) continue;
                    if (std::regex_search(node.name, rx) ||
                        (!node.model_tag.empty() && std::regex_search(node.model_tag, rx))) {
                        sel[node.id] = true;
                    }
                }
                return true;
            } catch (const std::exception&) {
                std::cerr << "Warning: " << label << " regex '" << re << "' is invalid; no branches selected.\n";
                return false;
            }
        };

        // 1. explicit regex
        if (to_lower(spec.substr(0, 6)) == "regex:") {
            apply_regex(spec.substr(6));
            return sel;
        }
        if (force_regex) {
            apply_regex(spec);
            return sel;
        }

        // Tags present in the tree
        std::vector<std::string> tags_present;
        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            std::string tl = node_tag(node);
            if (!tl.empty() && std::find(tags_present.begin(), tags_present.end(), tl) == tags_present.end()) {
                tags_present.push_back(tl);
            }
        }
        // Selects branches carrying tag `q` (exact) or, failing that, any tag of q's alias group.
        auto select_tag = [&](const std::string& q) -> bool {
            bool exact = std::find(tags_present.begin(), tags_present.end(), q) != tags_present.end();
            int grp = alias_group(q);
            bool group_present = false;
            if (!exact && grp != 0) {
                for (const auto& tp : tags_present) if (alias_group(tp) == grp) group_present = true;
            }
            if (!exact && !group_present) return false;
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                std::string tl = node_tag(node);
                if (tl.empty()) continue;
                if (exact ? (tl == q) : (alias_group(tl) == grp)) sel[node.id] = true;
            }
            return true;
        };
        auto select_keyword = [&](const std::string& q) -> bool {
            bool is_kw = (q == "all" || q == "internal" || q == "leaves" || q == "terminal" || q == "tips" ||
                          q == "unlabeled" || q == "untagged");
            if (!is_kw) return false;
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                bool m = (q == "all") ||
                         (q == "internal" && !node.is_leaf) ||
                         ((q == "leaves" || q == "terminal" || q == "tips") && node.is_leaf) ||
                         ((q == "unlabeled" || q == "untagged") && node.model_tag.empty());
                if (m) sel[node.id] = true;
            }
            return true;
        };

        std::string spec_tag_l = to_lower(trim_tag(spec));
        // 2. keyword, 3. tag
        if (spec.find(',') == std::string::npos) {
            if (select_keyword(spec_tag_l)) return sel;
            if (select_tag(spec_tag_l)) return sel;
        }

        // 4. comma-separated exact branch names (single names included)
        {
            std::vector<std::string> unmatched;
            size_t matched_tokens = 0;
            std::stringstream ss(spec);
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                std::string t = trim_ws(tok);
                if (t.empty()) continue;
                bool hit = false;
                for (const auto& node : tree.nodes) {
                    if (node.id == tree.root_id) continue;
                    if (node.name == t) { sel[node.id] = true; hit = true; }
                }
                if (!hit) hit = select_keyword(to_lower(trim_tag(t)));
                if (!hit) hit = select_tag(to_lower(trim_tag(t)));
                if (hit) ++matched_tokens; else unmatched.push_back(t);
            }
            if (matched_tokens > 0) {
                if (!unmatched.empty()) {
                    std::cerr << "Warning: " << label << ": " << unmatched.size()
                              << " name(s) did not match any branch:";
                    for (const auto& u : unmatched) std::cerr << " '" << u << "'";
                    std::cerr << "\n";
                }
                return sel;
            }
        }

        // 5. regex fallback (with warning)
        std::fill(sel.begin(), sel.end(), false);
        apply_regex(spec);
        std::cerr << "Warning: " << label << " '" << spec << "' matched no tag or exact branch name; "
                  << "interpreted as a case-insensitive regex (" << count() << " branches matched). "
                  << "Use 'regex:<pattern>' to request regex matching explicitly.\n";
        return sel;
    }

    void partition_branches() {
        auto to_lower = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
            return s;
        };
        auto trim = [](std::string s) {
            while (!s.empty() && (std::isspace(static_cast<unsigned char>(s.front())) || s.front() == '{')) s.erase(s.begin());
            while (!s.empty() && (std::isspace(static_cast<unsigned char>(s.back())) || s.back() == '}')) s.pop_back();
            return s;
        };

        bool has_tagged_test = false;
        bool has_tagged_ref = false;

        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && !node.model_tag.empty()) {
                std::string tl = to_lower(trim(node.model_tag));
                if (tl == "t" || tl == "test" || tl == "fg" || tl == "foreground") {
                    has_tagged_test = true;
                } else if (tl == "r" || tl == "ref" || tl == "reference" || tl == "bg" || tl == "background") {
                    has_tagged_ref = true;
                }
            }
        }

        bool explicit_test = !settings.test_branch_regex.empty() || !settings.test_branch_names.empty();
        bool explicit_ref = !settings.reference_branch_regex.empty() || !settings.reference_branch_names.empty();

        std::vector<bool> test_sel, ref_sel;
        if (explicit_test) {
            test_sel = resolve_branch_spec(settings.test_branch_regex, settings.test_branch_names,
                                           settings.test_is_regex, "--test");
        }
        if (explicit_ref) {
            ref_sel = resolve_branch_spec(settings.reference_branch_regex, settings.reference_branch_names,
                                          settings.reference_is_regex, "--reference");
        }

        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;

            bool is_test = false;
            if (explicit_test) {
                is_test = test_sel[node.id];
            } else if (has_tagged_test) {
                std::string tl = to_lower(trim(node.model_tag));
                is_test = (tl == "t" || tl == "test" || tl == "fg" || tl == "foreground");
            }

            bool is_ref = false;
            if (!is_test) {
                if (explicit_ref) {
                    is_ref = ref_sel[node.id];
                } else if (has_tagged_ref) {
                    std::string tl = to_lower(trim(node.model_tag));
                    is_ref = (tl == "r" || tl == "ref" || tl == "reference" || tl == "bg" || tl == "background");
                } else {
                    is_ref = true;
                }
            }

            is_test_branch[node.id] = is_test;
            if (is_test) {
                branch_partition[node.id] = 1;
                test_branch_names.push_back(node.name);
                branch_class[node.name] = "Test";
            } else if (is_ref) {
                branch_partition[node.id] = 0;
                ref_branch_names.push_back(node.name);
                branch_class[node.name] = "Reference";
            } else {
                branch_partition[node.id] = 2;
                unclassified_branch_names.push_back(node.name);
                branch_class[node.name] = "Unclassified";
            }
        }
    }

    // Branch-length -> rate-matrix time conversion factor for each of the 3 omega classes
    // of one branch partition (Reference or Test). The unclassified partition has a single
    // omega, so both normalizations coincide there.
    //  PerClass      : conv[c] = 3 / mu_c       (exact historical expression)
    //  SharedMixture : conv[c] = 3 / mu_mix  for all c, mu_mix = sum_c w_c mu_c / sum_c w_c
    // mu_c is the mean substitution rate (per codon) of the unnormalized MG94 Q with alpha=1,
    // beta=omega_c, so the branch length is the mixture-expected substitutions per nucleotide.
    std::vector<Scalar> class_conversion_factors(
        const std::vector<MG94Matrix>& mats,
        const std::vector<Scalar>& weights
    ) const {
        std::vector<Scalar> conv(3, 1.0);
        if (settings.normalization == RELAXNormalization::SharedMixture) {
            Scalar w_sum = 0.0;
            Scalar mu_mix = 0.0;
            for (int c = 0; c < 3; ++c) {
                Scalar w = std::max(0.0, weights[c]);
                w_sum += w;
                mu_mix += w * mats[c].scale_factor;
            }
            if (w_sum > 0.0) mu_mix /= w_sum;
            Scalar shared = (mu_mix > 1e-12) ? (3.0 / mu_mix) : 1.0;
            for (int c = 0; c < 3; ++c) conv[c] = shared;
        } else {
            for (int c = 0; c < 3; ++c) {
                conv[c] = (mats[c].scale_factor > 1e-12) ? (3.0 / mats[c].scale_factor) : 1.0;
            }
        }
        return conv;
    }

    // Evaluates 3-category RELAX mixture log-likelihood
    Scalar evaluate_relax_mixture(
        const Tree& t,
        const std::vector<Scalar>& omegas_R,
        const std::vector<Scalar>& weights,
        Scalar k,
        Scalar tree_scale,
        std::vector<Scalar>* out_site_ll = nullptr,
        const std::vector<Scalar>* custom_branch_lengths = nullptr
    ) const {
        size_t num_patterns = aln.patterns.size();
        size_t num_nodes = t.num_nodes();
        size_t S = code->num_sense_codons;
        const Vector& eq_freqs = aln.codon_frequencies_f3x4;

        // Derive Test branch omegas: omega_T,i = (omega_R,i)^k
        std::vector<Scalar> omegas_T(3);
        for (int c = 0; c < 3; ++c) {
            if (omegas_R[c] <= 1e-12) {
                omegas_T[c] = (k > 0) ? 0.0 : 1.0;
            } else {
                omegas_T[c] = std::clamp(std::pow(omegas_R[c], k), 1e-6, 10000.0);
            }
        }

        // Precompute transition matrices for all 3 classes across all tree branches
        std::vector<std::vector<Matrix>> P_branches(3, std::vector<Matrix>(num_nodes));

        MG94Matrix mat_U;
        Scalar conv_U = 1.0;
        if (!unclassified_branch_names.empty()) {
            MG94Parameters p_U = base_params;
            p_U.alpha = 1.0;
            p_U.beta = unclassified_omega;
            mat_U.update(p_U, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            conv_U = (mat_U.scale_factor > 1e-12) ? (3.0 / mat_U.scale_factor) : 1.0;
        }

        std::vector<MG94Matrix> mats_R(3), mats_T(3);
        #pragma omp parallel for schedule(static)
        for (int job = 0; job < 6; ++job) {
            const int c = job % 3;
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            if (job < 3) {
                p.beta = omegas_R[c];
                mats_R[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            } else {
                p.beta = omegas_T[c];
                mats_T[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            }
        }
        const std::vector<Scalar> convs_R = class_conversion_factors(mats_R, weights);
        const std::vector<Scalar> convs_T = class_conversion_factors(mats_T, weights);

        const bool branch_site = (settings.mixture_level == RELAXMixtureLevel::BranchSite);
        std::vector<int> branch_ids;
        branch_ids.reserve(num_nodes);
        for (const auto& node : t.nodes) {
            if (node.id != t.root_id) branch_ids.push_back(node.id);
        }

        // Transition matrices for every (class, branch); classes with zero weight are never used
        // by the pruning below and are skipped. In BranchSite mode the classes are mixed inside
        // each branch matrix (stored in slot 0), so pruning runs once with a single class.
        #pragma omp parallel for schedule(dynamic)
        for (size_t bi = 0; bi < branch_ids.size(); ++bi) {
            const int nid = branch_ids[bi];
            const auto& node = t.nodes[nid];
            Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[nid] : (node.branch_length * tree_scale);
            for (int c = 0; c < 3; ++c) {
                if (weights[c] <= 1e-12) continue;
                if (branch_partition[nid] == 1) {
                    P_branches[c][nid] = mats_T[c].transition_matrix(bl * convs_T[c]);
                } else if (branch_partition[nid] == 2) {
                    P_branches[c][nid] = mat_U.transition_matrix(bl * conv_U);
                } else {
                    P_branches[c][nid] = mats_R[c].transition_matrix(bl * convs_R[c]);
                }
            }
            if (branch_site) {
                Matrix P_mix = Matrix::Zero(S, S);
                for (int c = 0; c < 3; ++c) {
                    if (weights[c] <= 1e-12) continue;
                    P_mix.noalias() += weights[c] * P_branches[c][nid];
                }
                P_branches[0][nid] = std::move(P_mix);
            }
        }

        std::vector<Scalar> class_w = weights;
        if (branch_site) {
            class_w = {1.0, 0.0, 0.0};
        }

        std::vector<Scalar> pattern_ll(num_patterns, 0.0);
        Scalar total_lnl = 0.0;

        #pragma omp parallel
        {
            auto [node_L, node_log_scale] = LikelihoodEngine::make_pattern_scratch<Vector>(num_nodes, S);

            #pragma omp for reduction(+:total_lnl) schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pll = -1e20;

                if (branch_site) {
                    LikelihoodEngine::eval_pattern_postorder(t, pattern, leaf_to_taxon, P_branches[0], static_cast<int>(S), node_L, node_log_scale);

                    Scalar class_L = eq_freqs.dot(node_L[t.root_id]);
                    if (class_L > 0.0) {
                        pll = std::log(class_L) + node_log_scale[t.root_id];
                    }
                } else {
                    std::array<Scalar, 3> terms;
                    terms.fill(-std::numeric_limits<Scalar>::infinity());
                    for (int c = 0; c < 3; ++c) {
                        if (class_w[c] <= 1e-12) continue;

                        LikelihoodEngine::eval_pattern_postorder(t, pattern, leaf_to_taxon, P_branches[c], static_cast<int>(S), node_L, node_log_scale);
                        Scalar class_L = eq_freqs.dot(node_L[t.root_id]);
                        if (class_L > 0.0) {
                            terms[c] = std::log(class_w[c]) + std::log(class_L) + node_log_scale[t.root_id];
                        }
                    }

                    Scalar lse = LikelihoodEngine::logsumexp(terms);
                    if (std::isfinite(lse)) pll = lse;
                }
                total_lnl += pattern.weight * pll;
                if (out_site_ll) {
                    pattern_ll[p] = pll;
                }
            }
        }

        if (out_site_ll) {
            out_site_ll->resize(aln.num_codons);
            for (size_t s = 0; s < aln.num_codons; ++s) {
                (*out_site_ll)[s] = pattern_ll[aln.site_to_pattern[s]];
            }
        }

        return total_lnl;
    }

    // Evaluates log-likelihood and exact analytical branch gradients d ln L / d t_b
    // under the 3-category RELAX mixture model using Inside-Outside adjoints
    std::pair<Scalar, std::vector<Scalar>> compute_relax_branch_gradients(
        const std::vector<Scalar>& branch_lengths,
        const std::vector<Scalar>& omegas_R,
        const std::vector<Scalar>& weights,
        Scalar k
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        const Vector& eq_freqs = aln.codon_frequencies_f3x4;

        // Derive Test branch omegas: omega_T,i = (omega_R,i)^k
        std::vector<Scalar> omegas_T(3);
        for (int c = 0; c < 3; ++c) {
            if (omegas_R[c] <= 1e-12) {
                omegas_T[c] = (k > 0) ? 0.0 : 1.0;
            } else {
                omegas_T[c] = std::clamp(std::pow(omegas_R[c], k), 1e-6, 10000.0);
            }
        }

        // Precompute transition matrices P_c and derivatives dP_c for all 3 classes
        std::vector<std::vector<Matrix>> P_c(3, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<Matrix>> dP_c(3, std::vector<Matrix>(num_nodes));

        MG94Matrix mat_U;
        Scalar conv_U = 1.0;
        if (!unclassified_branch_names.empty()) {
            MG94Parameters p_U = base_params;
            p_U.alpha = 1.0;
            p_U.beta = unclassified_omega;
            mat_U.update(p_U, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            conv_U = (mat_U.scale_factor > 1e-12) ? (3.0 / mat_U.scale_factor) : 1.0;
        }

        const bool shared_norm = (settings.normalization == RELAXNormalization::SharedMixture);
        const bool branch_site = (settings.mixture_level == RELAXMixtureLevel::BranchSite);
        std::vector<MG94Matrix> mats_R(3), mats_T(3);
        #pragma omp parallel for schedule(static)
        for (int job = 0; job < 6; ++job) {
            const int c = job % 3;
            // PerClass only needs matrices for classes that contribute; SharedMixture needs
            // every class's mean rate for the mixture normalization.
            if (weights[c] <= 1e-12 && !shared_norm) continue;
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            if (job < 3) {
                p.beta = omegas_R[c];
                mats_R[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            } else {
                p.beta = omegas_T[c];
                mats_T[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            }
        }
        const std::vector<Scalar> convs_R = class_conversion_factors(mats_R, weights);
        const std::vector<Scalar> convs_T = class_conversion_factors(mats_T, weights);

        std::vector<int> branch_ids;
        branch_ids.reserve(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) branch_ids.push_back(node.id);
        }

        // P and dP/dt for every contributing (class, branch), built in parallel over branches.
        // In BranchSite mode the classes are mixed into slot 0:
        //   P_mix = sum_c w_c P_c,  dP_mix = sum_c w_c conv_c Q_c P_c.
        #pragma omp parallel for schedule(dynamic)
        for (size_t bi = 0; bi < branch_ids.size(); ++bi) {
            const int nid = branch_ids[bi];
            const Scalar bl = branch_lengths[nid];
            for (int c = 0; c < 3; ++c) {
                if (weights[c] <= 1e-12) continue;
                if (branch_partition[nid] == 1) {
                    Matrix P = mats_T[c].transition_matrix(bl * convs_T[c]);
                    dP_c[c][nid] = convs_T[c] * (mats_T[c].Q * P);
                    P_c[c][nid] = std::move(P);
                } else if (branch_partition[nid] == 2) {
                    Matrix P = mat_U.transition_matrix(bl * conv_U);
                    dP_c[c][nid] = conv_U * (mat_U.Q * P);
                    P_c[c][nid] = std::move(P);
                } else {
                    Matrix P = mats_R[c].transition_matrix(bl * convs_R[c]);
                    dP_c[c][nid] = convs_R[c] * (mats_R[c].Q * P);
                    P_c[c][nid] = std::move(P);
                }
            }
            if (branch_site) {
                Matrix P_mix = Matrix::Zero(S, S);
                Matrix dP_mix = Matrix::Zero(S, S);
                for (int c = 0; c < 3; ++c) {
                    if (weights[c] <= 1e-12) continue;
                    P_mix.noalias() += weights[c] * P_c[c][nid];
                    dP_mix.noalias() += weights[c] * dP_c[c][nid];
                }
                P_c[0][nid] = std::move(P_mix);
                dP_c[0][nid] = std::move(dP_mix);
            }
        }

        std::vector<Scalar> class_w = weights;
        if (branch_site) {
            class_w = {1.0, 0.0, 0.0};
        }

        Scalar total_log_l = 0.0;
        std::vector<Scalar> grad_b(num_nodes, 0.0);

        #pragma omp parallel
        {
            std::vector<Scalar> local_grad(num_nodes, 0.0);
            Scalar local_ll = 0.0;

            #pragma omp for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                if (branch_site) {
                    auto io = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_c[0], eq_freqs, S
                    );
                    if (io.likelihood > 0.0) {
                        local_ll += pattern.weight * io.log_likelihood;
                        Scalar inv_L = pattern.weight / io.likelihood;
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int32_t vid = node.id;
                                Scalar dL_dt = io.V[vid].dot(dP_c[0][vid] * io.D[vid]);
                                local_grad[vid] += inv_L * dL_dt;
                            }
                        }
                    }
                } else {
                    std::array<LikelihoodEngine::InsideOutsideResult, 3> io_c;
                    std::array<Scalar, 3> terms;
                    terms.fill(-std::numeric_limits<Scalar>::infinity());
                    for (int c = 0; c < 3; ++c) {
                        if (class_w[c] <= 1e-12) continue;
                        io_c[c] = LikelihoodEngine::compute_inside_outside(
                            tree, pattern, leaf_to_taxon, P_c[c], eq_freqs, S
                        );
                        if (io_c[c].likelihood > 0.0) {
                            terms[c] = std::log(class_w[c]) + io_c[c].log_likelihood;
                        }
                    }

                    Scalar pat_ll = LikelihoodEngine::logsumexp(terms);
                    if (std::isfinite(pat_ll)) {
                        local_ll += pattern.weight * pat_ll;
                        std::array<Scalar, 3> post_w = {0.0, 0.0, 0.0};
                        for (int c = 0; c < 3; ++c) {
                            if (std::isfinite(terms[c])) post_w[c] = std::exp(terms[c] - pat_ll);
                        }

                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int32_t vid = node.id;
                                Scalar d_lnL_dt = 0.0;
                                for (int c = 0; c < 3; ++c) {
                                    if (post_w[c] > 1e-12 && io_c[c].likelihood > 0.0) {
                                        d_lnL_dt += post_w[c] * (io_c[c].V[vid].dot(dP_c[c][vid] * io_c[c].D[vid]) / io_c[c].likelihood);
                                    }
                                }
                                local_grad[vid] += pattern.weight * d_lnL_dt;
                            }
                        }
                    }
                }
            }

            #pragma omp critical
            {
                total_log_l += local_ll;
                for (size_t i = 0; i < num_nodes; ++i) {
                    grad_b[i] += local_grad[i];
                }
            }
        }

        return {total_log_l, grad_b};
    }

    // High-performance joint branch length refinement via exact Inside-Outside adjoint gradients and L-BFGS
    std::vector<Scalar> refine_branch_lengths(
        const std::vector<Scalar>& init_bl,
        const std::vector<Scalar>& omegas_R,
        const std::vector<Scalar>& weights,
        Scalar k,
        int max_iters = 10
    ) const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Scalar> bl = init_bl;
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id)) {
                bl[i] = std::clamp(bl[i], 1e-6, 10.0);
            }
        }

        const size_t m_history = 5;
        std::vector<std::vector<Scalar>> s_hist;
        std::vector<std::vector<Scalar>> y_hist;
        std::vector<Scalar> rho_hist;

        auto [curr_ll, curr_grad] = compute_relax_branch_gradients(bl, omegas_R, weights, k);

        Scalar best_ll = curr_ll;
        std::vector<Scalar> best_bl = bl;

        for (int iter = 0; iter < max_iters; ++iter) {
            Scalar max_g = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_g = std::max(max_g, std::abs(curr_grad[i]));
                }
            }
            if (max_g < 1e-3) break;

            // Two-loop L-BFGS recursion
            std::vector<Scalar> q = curr_grad;
            std::vector<Scalar> alpha(s_hist.size());

            for (int j = static_cast<int>(s_hist.size()) - 1; j >= 0; --j) {
                Scalar s_dot_q = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    s_dot_q += s_hist[j][i] * q[i];
                }
                alpha[j] = rho_hist[j] * s_dot_q;
                for (size_t i = 0; i < num_nodes; ++i) {
                    q[i] -= alpha[j] * y_hist[j][i];
                }
            }

            Scalar gamma = 1.0;
            if (!s_hist.empty()) {
                Scalar s_dot_y = 0.0;
                Scalar y_dot_y = 0.0;
                const auto& s_last = s_hist.back();
                const auto& y_last = y_hist.back();
                for (size_t i = 0; i < num_nodes; ++i) {
                    s_dot_y += s_last[i] * y_last[i];
                    y_dot_y += y_last[i] * y_last[i];
                }
                if (y_dot_y > 1e-12) {
                    gamma = std::clamp(s_dot_y / y_dot_y, 1e-4, 10.0);
                }
            }

            std::vector<Scalar> r(num_nodes);
            for (size_t i = 0; i < num_nodes; ++i) {
                r[i] = gamma * q[i];
            }

            for (size_t j = 0; j < s_hist.size(); ++j) {
                Scalar y_dot_r = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    y_dot_r += y_hist[j][i] * r[i];
                }
                Scalar beta = rho_hist[j] * y_dot_r;
                for (size_t i = 0; i < num_nodes; ++i) {
                    r[i] += s_hist[j][i] * (alpha[j] - beta);
                }
            }

            std::vector<Scalar> d = r;
            d[tree.root_id] = 0.0;

            Scalar dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    dir_dot_grad += d[i] * curr_grad[i];
                }
            }
            if (dir_dot_grad <= 0.0) {
                d = curr_grad;
                d[tree.root_id] = 0.0;
                dir_dot_grad = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) {
                        dir_dot_grad += d[i] * curr_grad[i];
                    }
                }
                s_hist.clear();
                y_hist.clear();
                rho_hist.clear();
            }

            // Backtracking Armijo line search
            Scalar max_step = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_step = std::max(max_step, std::abs(d[i]));
                }
            }
            Scalar step_size = s_hist.empty() ? std::min(1.0, 0.05 / std::max(1.0, max_step)) : 1.0;
            std::vector<Scalar> new_bl(num_nodes);
            Scalar new_ll = curr_ll;
            std::vector<Scalar> new_grad;

            bool line_search_ok = false;
            for (int ls = 0; ls < 8; ++ls) {
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) {
                        new_bl[i] = std::clamp(bl[i] + step_size * d[i], 1e-6, 10.0);
                    } else {
                        new_bl[i] = 0.0;
                    }
                }
                auto [test_ll, test_grad] = compute_relax_branch_gradients(new_bl, omegas_R, weights, k);
                if (test_ll > curr_ll + 1e-4 * step_size * dir_dot_grad) {
                    new_ll = test_ll;
                    new_grad = test_grad;
                    line_search_ok = true;
                    break;
                }
                step_size *= 0.5;
            }

            if (!line_search_ok) {
                break;
            }

            // Update L-BFGS history
            std::vector<Scalar> s_k(num_nodes);
            std::vector<Scalar> y_k(num_nodes);
            Scalar s_dot_y = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                s_k[i] = new_bl[i] - bl[i];
                y_k[i] = curr_grad[i] - new_grad[i];
                s_dot_y += s_k[i] * y_k[i];
            }

            if (s_dot_y > 1e-8) {
                if (s_hist.size() >= m_history) {
                    s_hist.erase(s_hist.begin());
                    y_hist.erase(y_hist.begin());
                    rho_hist.erase(rho_hist.begin());
                }
                s_hist.push_back(s_k);
                y_hist.push_back(y_k);
                rho_hist.push_back(1.0 / s_dot_y);
            }

            bl = new_bl;
            curr_ll = new_ll;
            curr_grad = new_grad;

            if (curr_ll > best_ll) {
                best_ll = curr_ll;
                best_bl = bl;
            }
        }

        return best_bl;
    }

    // ---------------------------------------------------------------------------------
    // Distribution parametrization used by the optimizers (box-bounded, no penalty walls):
    //   theta = [r0, w1, u2, a0, a1 (, lk)]
    //   w1 in [0,1], w0 = r0 * w1 (r0 in [0,1])        -> 0 <= w0 <= w1 <= 1
    //   w2 = exp(u2), u2 in [0, ln 50]                  -> 1 <= w2 <= 50
    //   p0 = a0, p1 = (1-a0) a1, p2 = (1-a0)(1-a1)      (stick-breaking, a0, a1 in [0,1])
    //   K = exp(lk), lk in [ln 1e-5, ln 50]             (alternative model only)
    // ---------------------------------------------------------------------------------
    struct RELAXDist {
        std::vector<Scalar> omegas{0.0, 0.0, 0.0};
        std::vector<Scalar> weights{0.0, 0.0, 0.0};
        Scalar k = 1.0;
    };

    static constexpr Scalar kMinK = 1e-5;
    static constexpr Scalar kMaxK = 50.0;
    static constexpr Scalar kMaxOmega = 50.0;

    static RELAXDist make_dist(Scalar w0, Scalar w1, Scalar w2, Scalar p0, Scalar p1, Scalar k = 1.0) {
        RELAXDist d;
        d.omegas = {w0, w1, w2};
        d.weights = {p0, p1, std::max(0.0, 1.0 - p0 - p1)};
        d.k = k;
        return d;
    }

    static RELAXDist decode_theta(const std::vector<Scalar>& th, bool has_k) {
        RELAXDist d;
        Scalar w1 = th[1];
        d.omegas = {th[0] * w1, w1, std::exp(th[2])};
        Scalar p0 = th[3];
        Scalar p1 = (1.0 - th[3]) * th[4];
        d.weights = {p0, p1, std::max(0.0, 1.0 - p0 - p1)};
        d.k = has_k ? std::exp(th[5]) : 1.0;
        return d;
    }

    static std::vector<Scalar> encode_theta(const RELAXDist& d, bool has_k) {
        std::vector<Scalar> th(has_k ? 6 : 5, 0.0);
        Scalar w1 = std::clamp(d.omegas[1], 0.0, 1.0);
        th[0] = (w1 > 1e-12) ? std::clamp(d.omegas[0] / w1, 0.0, 1.0) : 0.0;
        th[1] = w1;
        th[2] = std::log(std::clamp(d.omegas[2], 1.0, kMaxOmega));
        th[3] = std::clamp(d.weights[0], 0.0, 1.0);
        Scalar rest = 1.0 - th[3];
        th[4] = (rest > 1e-12) ? std::clamp(d.weights[1] / rest, 0.0, 1.0) : 0.0;
        if (has_k) th[5] = std::log(std::clamp(d.k, kMinK, kMaxK));
        return th;
    }

    static void theta_bounds(bool has_k, std::vector<Scalar>& lb, std::vector<Scalar>& ub) {
        lb = {0.0, 0.0, 0.0, 0.0, 0.0};
        ub = {1.0, 1.0, std::log(kMaxOmega), 1.0, 1.0};
        if (has_k) {
            lb.push_back(std::log(kMinK));
            ub.push_back(std::log(kMaxK));
        }
    }

    Scalar loglik_at(const RELAXDist& d, const std::vector<Scalar>& bl,
                     std::vector<Scalar>* site_ll = nullptr) const {
        return evaluate_relax_mixture(tree, d.omegas, d.weights, d.k, 1.0, site_ll, &bl);
    }

    // Nelder-Mead over theta at fixed branch lengths, restarted from the incumbent while it
    // keeps improving. `ll` must hold the log-likelihood at `th` on entry; never decreases.
    std::vector<Scalar> optimize_theta(std::vector<Scalar> th, bool has_k,
                                       const std::vector<Scalar>& bl, Scalar& ll) const {
        std::vector<Scalar> lb, ub;
        theta_bounds(has_k, lb, ub);
        auto f = [&](const std::vector<Scalar>& x) -> Scalar {
            return -loglik_at(decode_theta(x, has_k), bl);
        };
        Scalar best = -ll;
        for (int restart = 0; restart < 3; ++restart) {
            auto o = NelderMeadDynamic::minimize(f, th, lb, ub, 1e-4, 200);
            Scalar gain = best - o.f;
            if (o.f < best) {
                best = o.f;
                th = o.x;
            }
            if (gain < 1e-3) break;
        }
        ll = -best;
        return th;
    }

    // Nelder-Mead over [s, theta] with branch lengths = base_bl * s (tree-scale stage).
    std::pair<std::vector<Scalar>, Scalar> optimize_scaled(
        const std::vector<std::vector<Scalar>>& starts, bool has_k, const std::vector<Scalar>& base_bl
    ) const {
        std::vector<Scalar> lb, ub;
        theta_bounds(has_k, lb, ub);
        lb.insert(lb.begin(), 0.05);
        ub.insert(ub.begin(), 10.0);
        std::vector<Scalar> bl(base_bl.size(), 0.0);
        auto f = [&](const std::vector<Scalar>& x) -> Scalar {
            for (size_t i = 0; i < base_bl.size(); ++i) bl[i] = base_bl[i] * x[0];
            std::vector<Scalar> th(x.begin() + 1, x.end());
            return -loglik_at(decode_theta(th, has_k), bl);
        };
        std::vector<std::vector<Scalar>> grid(starts.begin() + 1, starts.end());
        auto o = NelderMeadDynamic::minimize(f, starts[0], lb, ub, 1e-4, 400, grid);
        auto x = o.x;
        Scalar fbest = o.f;
        auto o2 = NelderMeadDynamic::minimize(f, x, lb, ub, 1e-5, 250);
        if (o2.f < fbest) { x = o2.x; fbest = o2.f; }
        return {x, -fbest};
    }

    struct RELAXJointFit {
        std::vector<Scalar> theta;
        std::vector<Scalar> bl;
        Scalar ll = -std::numeric_limits<Scalar>::infinity();
        int cycles = 0;
    };

    // Alternates per-branch length refinement (L-BFGS, exact gradients) with re-optimization of
    // the distribution (and K) at the refined lengths, until the cycle gain < joint_tolerance or
    // max_joint_cycles is reached. Monotone: the log-likelihood never decreases.
    RELAXJointFit joint_optimize(std::vector<Scalar> th, std::vector<Scalar> bl, bool has_k) const {
        RELAXJointFit J;
        J.theta = std::move(th);
        J.bl = std::move(bl);
        J.ll = loglik_at(decode_theta(J.theta, has_k), J.bl);
        const int cycles = std::max(1, settings.max_joint_cycles);
        for (int cycle = 0; cycle < cycles; ++cycle) {
            Scalar start_ll = J.ll;
            if (settings.refine_branch_lengths) {
                RELAXDist d = decode_theta(J.theta, has_k);
                auto nbl = refine_branch_lengths(J.bl, d.omegas, d.weights, d.k, 30);
                Scalar nll = loglik_at(d, nbl);
                if (nll > J.ll) {
                    J.bl = std::move(nbl);
                    J.ll = nll;
                }
            }
            J.theta = optimize_theta(J.theta, has_k, J.bl, J.ll);
            J.cycles = cycle + 1;
            if (!settings.refine_branch_lengths) break;
            if (J.ll - start_ll < settings.joint_tolerance) break;
        }
        return J;
    }

    RELAXResult run(
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        bool show_progress = false,
        bool force_progress = false
    ) {
        auto start_time = std::chrono::high_resolution_clock::now();
        RELAXResult result;
        result.settings = settings;
        result.test_branches = test_branch_names;
        result.reference_branches = ref_branch_names;
        result.branch_class = branch_class;

        std::unique_ptr<ProgressBar> pbar;
        if (show_progress || force_progress) {
            pbar = std::make_unique<ProgressBar>(100, "RELAX Analysis");
        }

        // Appends wall-clock elapsed time and a log-likelihood to progress messages.
        auto with_stats = [&](const std::string& stage, Scalar ll) {
            double el = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
            std::ostringstream os;
            os << stage << " (ln L = " << std::fixed << std::setprecision(3) << ll
               << ", " << std::setprecision(1) << el << " s)";
            return os.str();
        };

        auto update_progress = [&](const std::string& stage, double frac) {
            if (pbar) {
                pbar->set_status(stage);
                pbar->update(static_cast<int>(frac * 100.0));
            } else if (progress_cb) {
                progress_cb(stage, frac);
            }
        };

        // ---------------------------------------------------------
        // Phase 1: Nucleotide GTR Model Fit
        // ---------------------------------------------------------
        update_progress("Phase 1: Fitting Nucleotide GTR model", 0.05);
        GTRFitter gtr_fitter(tree, aln);
        auto gtr_res = gtr_fitter.fit();

        // Polish the GTR fit. GTRFitter estimates the exchangeabilities once, from the initial
        // branch lengths; here rates and branch lengths are alternated to convergence. This
        // matters because HyPhy 2.5 RELAX fixes the nucleotide rates at the GTR MLE for every
        // codon fit (RELAX.bf: estimators.fixSubsetOfEstimates(relax.gtr_results, ...)).
        if (settings.refine_branch_lengths && gtr_res.tree.bridges.empty()) {
            const Alignment& nuc = gtr_fitter.nuc_aln;
            const std::vector<size_t> l2t = LikelihoodEngine::build_leaf_to_taxon_map(gtr_res.tree, nuc);
            auto eval_gtr = [&](const Tree& t, const GTRParameters& p) -> Scalar {
                GTRMatrix g;
                g.update(p, nuc.nuc_frequencies);
                return LikelihoodEngine::compute_gtr_log_likelihood(t, nuc, g, &l2t);
            };
            auto with_rates = [&](const std::array<Scalar, 5>& x) {
                GTRParameters p = gtr_res.params;
                p.theta_AC = std::exp(x[0]);
                p.theta_AT = std::exp(x[1]);
                p.theta_CG = std::exp(x[2]);
                p.theta_CT = std::exp(x[3]);
                p.theta_GT = std::exp(x[4]);
                return p;
            };
            Scalar cur = eval_gtr(gtr_res.tree, gtr_res.params);
            const Scalar lo = std::log(1e-3), hi = std::log(1e3);
            for (int cycle = 0; cycle < 10; ++cycle) {
                const Scalar start_ll = cur;
                std::array<Scalar, 5> x0 = {
                    std::log(std::max(1e-3, gtr_res.params.theta_AC)), std::log(std::max(1e-3, gtr_res.params.theta_AT)),
                    std::log(std::max(1e-3, gtr_res.params.theta_CG)), std::log(std::max(1e-3, gtr_res.params.theta_CT)),
                    std::log(std::max(1e-3, gtr_res.params.theta_GT))
                };
                auto f = [&](const std::array<Scalar, 5>& x) -> Scalar { return -eval_gtr(gtr_res.tree, with_rates(x)); };
                auto o = NelderMeadND<5>::minimize(f, x0, {lo, lo, lo, lo, lo}, {hi, hi, hi, hi, hi}, 1e-6, 400);
                if (-o.f > cur) {
                    gtr_res.params = with_rates(o.x);
                    cur = -o.f;
                }
                // GTRBranchLengthFitter exists in the committed optimizer API (unlike the newer
                // spectral-cavity sweep), keeping RELAX independent of uncommitted optimizer work.
                Tree t2 = GTRBranchLengthFitter::fit(gtr_res.tree, nuc, gtr_res.params, nuc.nuc_frequencies, 2).first;
                const Scalar c2 = eval_gtr(t2, gtr_res.params);
                if (c2 > cur) {
                    gtr_res.tree = std::move(t2);
                    cur = c2;
                }
                if (cur - start_ll < 1e-3) break;
            }
            gtr_res.log_likelihood = cur;
            const Scalar Kg = static_cast<Scalar>(5 + tree.num_nodes() - 1);
            const Scalar Ng = static_cast<Scalar>(nuc.num_nucleotides);
            gtr_res.aicc = -2.0 * cur + 2.0 * Kg * Ng / (Ng - Kg - 1.0);
        }

        result.gtr_log_likelihood = gtr_res.log_likelihood;
        result.gtr_aicc = gtr_res.aicc;
        result.gtr_parameters = 5 + (tree.num_nodes() - 1) + 2;
        result.gtr_rates = gtr_res.params;

        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                result.gtr_branch_lengths[node.name] = node.branch_length;
            }
        }

        base_params.theta_AC = gtr_res.params.theta_AC;
        base_params.theta_AT = gtr_res.params.theta_AT;
        base_params.theta_CG = gtr_res.params.theta_CG;
        base_params.theta_CT = gtr_res.params.theta_CT;
        base_params.theta_GT = gtr_res.params.theta_GT;

        update_progress("Phase 1: GTR fit complete", 0.20);

        // ---------------------------------------------------------
        // Phase 2: MG94 with Separate Rates for Branch Sets
        // ---------------------------------------------------------
        update_progress("Phase 2: Fitting MG94 with separate branch rates", 0.25);

        const size_t n_nodes = tree.num_nodes();
        std::vector<Scalar> gtr_bl(n_nodes, 0.0);
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) gtr_bl[node.id] = node.branch_length;
        }
        auto scale_bl = [](const std::vector<Scalar>& b, Scalar s) {
            std::vector<Scalar> out(b.size(), 0.0);
            for (size_t i = 0; i < b.size(); ++i) out[i] = b[i] * s;
            return out;
        };

        // Single-omega-per-set MG94 is expressed through the RELAX evaluator with one class:
        // omega_T = omega_R^k, so k = log(omega_T) / log(omega_R).
        auto mg94_k = [](Scalar w_R, Scalar w_T) -> Scalar {
            return (std::abs(std::log(std::max(1e-5, w_R))) > 1e-6)
                ? (std::log(std::max(1e-5, w_T)) / std::log(std::max(1e-5, w_R)))
                : 1.0;
        };
        const std::vector<Scalar> one_class = {1.0, 0.0, 0.0};

        auto eval_sep_mg94 = [&](const std::array<Scalar, 3>& x) -> Scalar {
            return -evaluate_relax_mixture(gtr_res.tree, {x[0], x[0], x[0]}, one_class, mg94_k(x[0], x[1]), x[2]);
        };

        auto opt_sep = NelderMeadND<3>::minimize(
            eval_sep_mg94,
            {0.1, 0.5, 1.0},
            {1e-4, 1e-4, 0.1},
            {10.0, 10.0, 10.0},
            1e-4, 100
        );

        Scalar mg_wR = opt_sep.x[0];
        Scalar mg_wT = opt_sep.x[1];
        Scalar sep_scale = opt_sep.x[2];
        Scalar mg_ll = -opt_sep.f;
        std::vector<Scalar> mg_bl = scale_bl(gtr_bl, sep_scale);

        // Per-branch refinement of the MG94 fit, alternating with (omega_R, omega_T).
        // The RELAX fits below start from these branch lengths.
        if (settings.refine_branch_lengths) {
            update_progress("Phase 2: Refining MG94 branch lengths", 0.35);
            const int cycles = std::max(1, settings.max_joint_cycles);
            for (int cycle = 0; cycle < cycles; ++cycle) {
                Scalar start_ll = mg_ll;
                Scalar kk = mg94_k(mg_wR, mg_wT);
                auto nbl = refine_branch_lengths(mg_bl, {mg_wR, mg_wR, mg_wR}, one_class, kk, 30);
                Scalar nll = evaluate_relax_mixture(tree, {mg_wR, mg_wR, mg_wR}, one_class, kk, 1.0, nullptr, &nbl);
                if (nll > mg_ll) {
                    mg_bl = std::move(nbl);
                    mg_ll = nll;
                }
                auto eval_w = [&](const std::array<Scalar, 2>& x) -> Scalar {
                    return -evaluate_relax_mixture(tree, {x[0], x[0], x[0]}, one_class, mg94_k(x[0], x[1]), 1.0, nullptr, &mg_bl);
                };
                auto o2 = NelderMeadND<2>::minimize(eval_w, {mg_wR, mg_wT}, {1e-4, 1e-4}, {10.0, 10.0}, 1e-5, 200);
                if (-o2.f > mg_ll) {
                    mg_wR = o2.x[0];
                    mg_wT = o2.x[1];
                    mg_ll = -o2.f;
                }

                // Nucleotide exchangeabilities (log scale, A<->G = 1 reference). HyPhy 2.5 RELAX
                // re-estimates them in its full MG94 refit ("Improving branch lengths, nucleotide
                // substitution biases, and global dN/dS ratios under a full codon model"); the
                // GTR estimates are starting values. The RELAX fits below inherit them.
                {
                    const MG94Parameters saved = base_params;
                    auto set_rates = [&](const std::array<Scalar, 5>& x) {
                        base_params.theta_AC = std::exp(x[0]);
                        base_params.theta_AT = std::exp(x[1]);
                        base_params.theta_CG = std::exp(x[2]);
                        base_params.theta_CT = std::exp(x[3]);
                        base_params.theta_GT = std::exp(x[4]);
                    };
                    const Scalar kk2 = mg94_k(mg_wR, mg_wT);
                    auto eval_nuc = [&](const std::array<Scalar, 5>& x) -> Scalar {
                        set_rates(x);
                        return -evaluate_relax_mixture(tree, {mg_wR, mg_wR, mg_wR}, one_class, kk2, 1.0, nullptr, &mg_bl);
                    };
                    std::array<Scalar, 5> x0 = {
                        std::log(std::max(1e-3, saved.theta_AC)), std::log(std::max(1e-3, saved.theta_AT)),
                        std::log(std::max(1e-3, saved.theta_CG)), std::log(std::max(1e-3, saved.theta_CT)),
                        std::log(std::max(1e-3, saved.theta_GT))
                    };
                    const Scalar lo = std::log(1e-3), hi = std::log(1e3);
                    auto o5 = NelderMeadND<5>::minimize(eval_nuc, x0, {lo, lo, lo, lo, lo}, {hi, hi, hi, hi, hi}, 1e-5, 300);
                    if (-o5.f > mg_ll) {
                        set_rates(o5.x);
                        mg_ll = -o5.f;
                    } else {
                        base_params = saved;
                    }
                }
                if (mg_ll - start_ll < settings.joint_tolerance) break;
            }
        }

        result.mg94_log_likelihood = mg_ll;
        result.mg94_omega_R = mg_wR;
        result.mg94_omega_T = mg_wT;

        size_t B = tree.num_nodes() - 1;
        result.mg94_parameters = 5 + B + 2 + 6; // GTR rates + branches + 2 omegas + CF3x4
        Scalar n_codons = static_cast<Scalar>(aln.num_codons);
        Scalar k_mg = static_cast<Scalar>(result.mg94_parameters);
        result.mg94_aicc = 2.0 * k_mg - 2.0 * result.mg94_log_likelihood +
            (2.0 * k_mg * (k_mg + 1.0)) / std::max(1.0, n_codons - k_mg - 1.0);

        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                result.mg94_branch_lengths[node.name] = mg_bl[node.id];
            }
        }

        update_progress(with_stats("Phase 2: MG94 separate rates fit complete", mg_ll), 0.40);

        auto with_scale = [](Scalar s, std::vector<Scalar> th) {
            th.insert(th.begin(), s);
            return th;
        };
        auto drop_k = [](const std::vector<Scalar>& th) {
            return std::vector<Scalar>(th.begin(), th.begin() + 5);
        };
        auto add_k = [](std::vector<Scalar> th, Scalar k) {
            th.resize(5);
            th.push_back(std::log(std::clamp(k, kMinK, kMaxK)));
            return th;
        };
        const Scalar w_mg = std::clamp(mg_wR, 1e-3, 0.95);

        // ---------------------------------------------------------
        // Phase 3: RELAX Null Model Fit (K = 1)
        // ---------------------------------------------------------
        update_progress("Phase 3: Fitting RELAX null model (K = 1 fixed)", 0.45);

        std::vector<std::vector<Scalar>> null_starts = {
            with_scale(1.0, encode_theta(make_dist(0.0, 0.5, 1.5, 0.6, 0.3), false)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.6, 1.2, 0.67, 0.27), false)),
            with_scale(1.0, encode_theta(make_dist(0.2 * w_mg, w_mg, 2.0, 0.3, 0.6), false)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.2, 3.0, 0.5, 0.45), false))
        };
        auto scaled_null = optimize_scaled(null_starts, false, mg_bl);
        Scalar null_s = scaled_null.first[0];
        RELAXJointFit J_null = joint_optimize(
            std::vector<Scalar>(scaled_null.first.begin() + 1, scaled_null.first.end()),
            scale_bl(mg_bl, null_s), false);

        update_progress(with_stats("Phase 3: RELAX null fit complete", J_null.ll), 0.60);

        // ---------------------------------------------------------
        // Phase 4: RELAX Alternative Model Fit (K free)
        // ---------------------------------------------------------
        update_progress("Phase 4: Fitting RELAX alternative model (K free)", 0.62);

        Scalar k_init = 1.0;
        if (std::abs(std::log(std::max(1e-5, mg_wR))) > 1e-5) {
            k_init = std::clamp(
                std::log(std::max(1e-5, mg_wT)) / std::log(std::max(1e-5, mg_wR)),
                0.01, 10.0
            );
        }

        // Start A: tree-scale grid search from generic starts
        std::vector<std::vector<Scalar>> alt_starts = {
            with_scale(1.0, encode_theta(make_dist(0.0, w_mg, 5.0, 0.35, 0.60, k_init), true)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.16, 30.0, 0.35, 0.65, 1e-4), true)),
            with_scale(1.0, encode_theta(make_dist(0.05, 0.3, 2.0, 0.40, 0.50, 0.5), true)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.5, 5.0, 0.3, 0.4, 2.0), true))
        };
        auto scaled_alt = optimize_scaled(alt_starts, true, mg_bl);
        Scalar alt_s = scaled_alt.first[0];
        RELAXJointFit J_alt = joint_optimize(
            std::vector<Scalar>(scaled_alt.first.begin() + 1, scaled_alt.first.end()),
            scale_bl(mg_bl, alt_s), true);

        // Start B: the null MLE with K = 1 (guarantees alt >= null for a monotone optimizer)
        update_progress(with_stats("Phase 4: Alternative model from null MLE; grid start reached", J_alt.ll), 0.75);
        RELAXJointFit J_alt_b = joint_optimize(add_k(J_null.theta, 1.0), J_null.bl, true);
        if (J_alt_b.ll > J_alt.ll) J_alt = std::move(J_alt_b);

        // HyPhy-style null check: re-fit the null from the alternative MLE with K := 1.
        update_progress(with_stats("Phase 4: Null model from alternative MLE; alternative at", J_alt.ll), 0.85);
        RELAXJointFit J_null_b = joint_optimize(drop_k(J_alt.theta), J_alt.bl, false);
        if (J_null_b.ll > J_null.ll) J_null = std::move(J_null_b);

        // Nesting guard: if the (improved) null beats the alternative, restart alt from it.
        if (J_null.ll > J_alt.ll) {
            RELAXJointFit J_alt_c = joint_optimize(add_k(J_null.theta, 1.0), J_null.bl, true);
            if (J_alt_c.ll > J_alt.ll) J_alt = std::move(J_alt_c);
        }

        update_progress(with_stats("Phase 4: RELAX fits complete; alternative", J_alt.ll), 0.90);

        // ---------------------------------------------------------
        // Record fits
        // ---------------------------------------------------------
        auto record_fit = [&](RELAXModelFit& fit, const RELAXJointFit& J, bool has_k, Scalar s, size_t n_params) {
            RELAXDist d = decode_theta(J.theta, has_k);
            fit.tree_scale = s;
            fit.k = d.k;
            fit.reference_distribution.omegas = d.omegas;
            fit.reference_distribution.weights = d.weights;
            std::vector<Scalar> om_T(3);
            for (int c = 0; c < 3; ++c) {
                om_T[c] = (d.omegas[c] <= 1e-12) ? ((d.k > 0) ? 0.0 : 1.0) : std::pow(d.omegas[c], d.k);
            }
            fit.test_distribution.omegas = has_k ? om_T : d.omegas;
            fit.test_distribution.weights = d.weights;
            fit.branch_lengths = J.bl;
            fit.log_likelihood = loglik_at(d, J.bl, &fit.site_log_likelihoods);
            fit.parameters = n_params;
            Scalar kp = static_cast<Scalar>(n_params);
            fit.aicc = 2.0 * kp - 2.0 * fit.log_likelihood +
                (2.0 * kp * (kp + 1.0)) / std::max(1.0, n_codons - kp - 1.0);
        };

        // GTR rates + branches + 3 omegas + 2 weights (+ k) + CF3x4
        record_fit(result.alternative_fit, J_alt, true, alt_s, 5 + B + 5 + 1 + 6);
        record_fit(result.null_fit, J_null, false, null_s, 5 + B + 5 + 6);
        Scalar alt_k = result.alternative_fit.k;

        // ---------------------------------------------------------
        // Phase 5: Hypothesis Testing & Interpretation
        // ---------------------------------------------------------
        update_progress("Phase 5: Calculating LRT & asymptotic p-value", 0.95);

        bool alt_failed = (result.alternative_fit.log_likelihood <= -1e19 || !std::isfinite(result.alternative_fit.log_likelihood));
        bool null_failed = (result.null_fit.log_likelihood <= -1e19 || !std::isfinite(result.null_fit.log_likelihood));
        if (alt_failed || null_failed) {
            result.is_failed = true;
            result.error_message = "Optimization failed: model log-likelihood underflowed or collapsed to penalty threshold (-1e20). The fit did not converge.";
            result.lrt = 0.0;
            result.p_value = 1.0;
            result.k = alt_k;
            result.is_significant = false;
            result.is_relaxed = false;
            result.is_intensified = false;
        } else {
            result.lrt = 2.0 * std::max(0.0, result.alternative_fit.log_likelihood - result.null_fit.log_likelihood);
            result.p_value = std::erfc(std::sqrt(result.lrt / 2.0));
            result.k = alt_k;

            result.is_significant = (result.p_value <= settings.p_value_threshold);
            if (result.is_significant) {
                if (result.k > 1.0) {
                    result.is_intensified = true;
                } else {
                    result.is_relaxed = true;
                }
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        result.runtime_seconds = std::chrono::duration<double>(end_time - start_time).count();

        update_progress("RELAX analysis complete", 1.0);

        return result;
    }
};

} // namespace hyphy::analyses
