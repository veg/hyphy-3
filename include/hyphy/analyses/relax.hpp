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
    int realized_cycles = 0;
    bool converged = false;
    RELAXRateDistribution reference_distribution;
    RELAXRateDistribution test_distribution;
    RELAXRateDistribution unclassified_distribution;
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
    int max_joint_cycles = 15;
    Scalar joint_tolerance = 1e-3;
    int grid_size = 0;       // 0 = default curated grid (8 null / 11 alt); > base size appends Halton points
    int starting_points = 1; // Number of top-ranked initial seeds refined via Nelder-Mead
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
    std::vector<Scalar> branch_switched_evidence_ratios;

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
    std::vector<Scalar> omegas_uncl = {0.05, 0.50, 2.0};
    std::vector<Scalar> weights_uncl = {0.50, 0.40, 0.10};

    // Neumaier compensated summation helper for deterministic across-thread reductions
    static Scalar neumaier_sum(const std::vector<Scalar>& vals) {
        Scalar sum = 0.0;
        Scalar c = 0.0;
        for (Scalar v : vals) {
            Scalar t = sum + v;
            if (std::abs(sum) >= std::abs(v)) {
                c += (sum - t) + v;
            } else {
                c += (v - t) + sum;
            }
            sum = t;
        }
        return sum + c;
    }

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
            if (node.id != tree.root_id) continue;
        }
        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            std::string tl = node_tag(node);
            if (!tl.empty() && std::find(tags_present.begin(), tags_present.end(), tl) == tags_present.end()) {
                tags_present.push_back(tl);
            }
        }
        // Selects branches carrying tag `q` (exact, case-insensitive) or, failing that, any tag of q's alias group.
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
                std::string tl_tok = to_lower(t);
                for (const auto& node : tree.nodes) {
                    if (node.id == tree.root_id) continue;
                    if (node.name == t || to_lower(node.name) == tl_tok) { sel[node.id] = true; hit = true; }
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
        auto canonical_group = [](const std::string& tl) -> std::string {
            if (tl.empty()) return "";
            if (tl == "t" || tl == "test" || tl == "fg" || tl == "foreground") return "__test__";
            if (tl == "r" || tl == "ref" || tl == "reference" || tl == "bg" || tl == "background") return "__ref__";
            return tl;
        };

        bool has_tagged_test = false;
        bool has_tagged_ref = false;
        std::vector<std::string> distinct_tag_groups;
        std::vector<std::string> non_empty_tag_groups;

        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            std::string tl = to_lower(trim(node.model_tag));
            std::string grp = canonical_group(tl);
            if (std::find(distinct_tag_groups.begin(), distinct_tag_groups.end(), grp) == distinct_tag_groups.end()) {
                distinct_tag_groups.push_back(grp);
            }
            if (!grp.empty() && std::find(non_empty_tag_groups.begin(), non_empty_tag_groups.end(), grp) == non_empty_tag_groups.end()) {
                non_empty_tag_groups.push_back(grp);
            }
            if (grp == "__test__") {
                has_tagged_test = true;
            } else if (grp == "__ref__") {
                has_tagged_ref = true;
            }
        }

        bool explicit_test = !settings.test_branch_regex.empty() || !settings.test_branch_names.empty();
        bool explicit_ref = !settings.reference_branch_regex.empty() || !settings.reference_branch_names.empty();

        // F-095: Support arbitrary single-tag Newick trees (e.g., {Marine} vs untagged) or 2-tag trees where one tag is canonical
        std::string auto_test_group;
        std::string auto_ref_group;
        if (!explicit_test && !has_tagged_test) {
            if (non_empty_tag_groups.size() == 1 && non_empty_tag_groups[0] != "__ref__") {
                auto_test_group = non_empty_tag_groups[0];
            } else if (non_empty_tag_groups.size() == 2 && has_tagged_ref) {
                auto_test_group = (non_empty_tag_groups[0] == "__ref__") ? non_empty_tag_groups[1] : non_empty_tag_groups[0];
            } else if (non_empty_tag_groups.size() >= 2) {
                throw std::runtime_error(
                    "Multiple branch tag groups exist in the tree without a designated test group. "
                    "Specify the test group explicitly via --test <group>."
                );
            }
        }
        if (!explicit_ref && !has_tagged_ref && has_tagged_test && non_empty_tag_groups.size() == 2 && distinct_tag_groups.size() == 2) {
            auto_ref_group = (non_empty_tag_groups[0] == "__test__") ? non_empty_tag_groups[1] : non_empty_tag_groups[0];
        }

        std::vector<bool> test_sel, ref_sel;
        if (explicit_test) {
            test_sel = resolve_branch_spec(settings.test_branch_regex, settings.test_branch_names,
                                           settings.test_is_regex, "--test");
        }
        if (explicit_ref) {
            ref_sel = resolve_branch_spec(settings.reference_branch_regex, settings.reference_branch_names,
                                          settings.reference_is_regex, "--reference");
        }

        const bool multi_group_tree = (distinct_tag_groups.size() > 2);

        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;

            std::string tl = to_lower(trim(node.model_tag));
            std::string grp = canonical_group(tl);

            bool is_test = false;
            if (explicit_test) {
                is_test = test_sel[node.id];
            } else if (has_tagged_test) {
                is_test = (grp == "__test__");
            } else if (!auto_test_group.empty()) {
                is_test = (grp == auto_test_group);
            }

            bool is_ref = false;
            if (!is_test) {
                if (explicit_ref) {
                    is_ref = ref_sel[node.id];
                } else if (has_tagged_ref) {
                    is_ref = (grp == "__ref__");
                } else if (!auto_ref_group.empty()) {
                    is_ref = (grp == auto_ref_group);
                } else if (!multi_group_tree) {
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
        if (ref_branch_names.empty()) {
            throw std::runtime_error(
                "No reference branches were selected for RELAX (either all branches matched the test branch set, "
                "or multiple branch tag groups exist without a designated reference group). "
                "Specify reference branches explicitly via --reference <group>."
            );
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
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        const MG94Parameters* custom_base_params = nullptr,
        const std::vector<Scalar>* custom_omegas_uncl = nullptr,
        const std::vector<Scalar>* custom_weights_uncl = nullptr
    ) const {
        size_t num_patterns = aln.patterns.size();
        size_t num_nodes = t.num_nodes();
        size_t S = code->num_sense_codons;
        const Vector& eq_freqs = aln.codon_frequencies_f3x4;
        const MG94Parameters& bp = custom_base_params ? *custom_base_params : base_params;
        const bool has_unclassified = !unclassified_branch_names.empty();
        const std::vector<Scalar>& om_U = custom_omegas_uncl ? *custom_omegas_uncl : omegas_uncl;
        const std::vector<Scalar>& wt_U = custom_weights_uncl ? *custom_weights_uncl : weights_uncl;

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

        std::vector<MG94Matrix> mats_R(3), mats_T(3), mats_U(3);
        const int n_mat_jobs = has_unclassified ? 9 : 6;
        #pragma omp parallel for schedule(static)
        for (int job = 0; job < n_mat_jobs; ++job) {
            const int c = job % 3;
            MG94Parameters p = bp;
            p.alpha = 1.0;
            if (job < 3) {
                p.beta = omegas_R[c];
                mats_R[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            } else if (job < 6) {
                p.beta = omegas_T[c];
                mats_T[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            } else {
                p.beta = om_U[c];
                mats_U[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            }
        }
        const std::vector<Scalar> convs_R = class_conversion_factors(mats_R, weights);
        const std::vector<Scalar> convs_T = class_conversion_factors(mats_T, weights);
        const std::vector<Scalar> convs_U = has_unclassified ? class_conversion_factors(mats_U, wt_U) : std::vector<Scalar>{1.0, 1.0, 1.0};

        const bool branch_site = (settings.mixture_level == RELAXMixtureLevel::BranchSite);
        std::vector<int> branch_ids;
        branch_ids.reserve(num_nodes);
        for (const auto& node : t.nodes) {
            if (node.id != t.root_id) branch_ids.push_back(node.id);
        }

        // Transition matrices for every (class, branch); classes with zero weight are never used
        // by the pruning below and are skipped. In BranchSite mode (or on unclassified branches,
        // whose 3-class BS-REL mixture (om_U, wt_U) is independent of the Test/Reference coupling)
        // the classes are mixed on the branch using br_w.
        #pragma omp parallel for schedule(static)
        for (size_t bi = 0; bi < branch_ids.size(); ++bi) {
            const int nid = branch_ids[bi];
            const auto& node = t.nodes[nid];
            Scalar bl = custom_branch_lengths ? (*custom_branch_lengths)[nid] : (node.branch_length * tree_scale);
            const int part = branch_partition[nid];
            const std::vector<Scalar>& br_w = (part == 2) ? wt_U : weights;
            if (branch_site || part == 2) {
                Matrix P_mix = Matrix::Zero(S, S);
                for (int c = 0; c < 3; ++c) {
                    if (br_w[c] <= 1e-12) continue;
                    if (part == 1) {
                        P_mix.noalias() += br_w[c] * mats_T[c].transition_matrix(bl * convs_T[c]);
                    } else if (part == 2) {
                        P_mix.noalias() += br_w[c] * mats_U[c].transition_matrix(bl * convs_U[c]);
                    } else {
                        P_mix.noalias() += br_w[c] * mats_R[c].transition_matrix(bl * convs_R[c]);
                    }
                }
                P_branches[0][nid] = P_mix;
                if (!branch_site) {
                    P_branches[1][nid] = P_mix;
                    P_branches[2][nid] = std::move(P_mix);
                }
            } else {
                for (int c = 0; c < 3; ++c) {
                    if (weights[c] <= 1e-12) continue;
                    if (part == 1) {
                        P_branches[c][nid] = mats_T[c].transition_matrix(bl * convs_T[c]);
                    } else {
                        P_branches[c][nid] = mats_R[c].transition_matrix(bl * convs_R[c]);
                    }
                }
            }
        }

        std::vector<Scalar> class_w = weights;
        if (branch_site) {
            class_w = {1.0, 0.0, 0.0};
        }

        std::vector<Scalar> pattern_ll(num_patterns, -1e20);
        std::vector<Scalar> weighted_pattern_ll(num_patterns, 0.0);
        int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            LikelihoodEngine::BatchedInsideWorkspace ws;
            ws.init(num_nodes, S, tile_size);
            Eigen::MatrixXd terms_mat(3, tile_size);

            #pragma omp for schedule(static)
            for (size_t tile = 0; tile < num_tiles; ++tile) {
                size_t p_start = tile * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);

                if (branch_site) {
                    LikelihoodEngine::compute_inside_batch_in_place(
                        t, aln, p_start, p_end, leaf_to_taxon, P_branches[0], eq_freqs, S, ws
                    );
                    for (Eigen::Index b = 0; b < B; ++b) {
                        size_t p = p_start + b;
                        Scalar pll = ws.root_log_L(b);
                        weighted_pattern_ll[p] = aln.patterns[p].weight * pll;
                        if (out_site_ll) pattern_ll[p] = pll;
                    }
                } else {
                    terms_mat.leftCols(B).setConstant(-std::numeric_limits<Scalar>::infinity());
                    for (int c = 0; c < 3; ++c) {
                        if (class_w[c] <= 1e-12) continue;
                        Scalar log_cw = std::log(class_w[c]);
                        LikelihoodEngine::compute_inside_batch_in_place(
                            t, aln, p_start, p_end, leaf_to_taxon, P_branches[c], eq_freqs, S, ws
                        );
                        for (Eigen::Index b = 0; b < B; ++b) {
                            if (ws.root_L(b) > 0.0) {
                                terms_mat(c, b) = log_cw + ws.root_log_L(b);
                            }
                        }
                    }
                    for (Eigen::Index b = 0; b < B; ++b) {
                        Scalar max_term = terms_mat.col(b).maxCoeff();
                        Scalar pll = -1e20;
                        if (max_term > -1e19) {
                            Scalar sum_exp = 0.0;
                            for (int c = 0; c < 3; ++c) {
                                if (terms_mat(c, b) > -1e19) {
                                    sum_exp += std::exp(terms_mat(c, b) - max_term);
                                }
                            }
                            if (sum_exp > 0.0) {
                                pll = max_term + std::log(sum_exp);
                            }
                        }
                        size_t p = p_start + b;
                        weighted_pattern_ll[p] = aln.patterns[p].weight * pll;
                        if (out_site_ll) pattern_ll[p] = pll;
                    }
                }
            }
        }

        Scalar total_lnl = neumaier_sum(weighted_pattern_ll);

        if (out_site_ll) {
            out_site_ll->resize(aln.num_codons);
            for (size_t s = 0; s < aln.num_codons; ++s) {
                (*out_site_ll)[s] = pattern_ll[aln.site_to_pattern[s]];
            }
        }

        return total_lnl;
    }

    // Evaluates log-likelihood and exact analytical branch gradients d ln L / d t_b
    // under the 3-category RELAX mixture model using batched Inside-Outside GEMM adjoints
    std::pair<Scalar, std::vector<Scalar>> compute_relax_branch_gradients(
        const std::vector<Scalar>& branch_lengths,
        const std::vector<Scalar>& omegas_R,
        const std::vector<Scalar>& weights,
        Scalar k,
        const MG94Parameters* custom_base_params = nullptr,
        const std::vector<Scalar>* custom_omegas_uncl = nullptr,
        const std::vector<Scalar>* custom_weights_uncl = nullptr
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        const Vector& eq_freqs = aln.codon_frequencies_f3x4;
        const MG94Parameters& bp = custom_base_params ? *custom_base_params : base_params;
        const bool has_unclassified = !unclassified_branch_names.empty();
        const std::vector<Scalar>& om_U = custom_omegas_uncl ? *custom_omegas_uncl : omegas_uncl;
        const std::vector<Scalar>& wt_U = custom_weights_uncl ? *custom_weights_uncl : weights_uncl;

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

        const bool shared_norm = (settings.normalization == RELAXNormalization::SharedMixture);
        const bool branch_site = (settings.mixture_level == RELAXMixtureLevel::BranchSite);
        std::vector<MG94Matrix> mats_R(3), mats_T(3), mats_U(3);
        const int n_mat_jobs = has_unclassified ? 9 : 6;
        #pragma omp parallel for schedule(static)
        for (int job = 0; job < n_mat_jobs; ++job) {
            const int c = job % 3;
            // PerClass only needs matrices for classes that contribute; SharedMixture needs
            // every class's mean rate for the mixture normalization.
            const std::vector<Scalar>& jw = (job >= 6) ? wt_U : weights;
            if (jw[c] <= 1e-12 && !shared_norm) continue;
            MG94Parameters p = bp;
            p.alpha = 1.0;
            if (job < 3) {
                p.beta = omegas_R[c];
                mats_R[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            } else if (job < 6) {
                p.beta = omegas_T[c];
                mats_T[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            } else {
                p.beta = om_U[c];
                mats_U[c].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            }
        }
        const std::vector<Scalar> convs_R = class_conversion_factors(mats_R, weights);
        const std::vector<Scalar> convs_T = class_conversion_factors(mats_T, weights);
        const std::vector<Scalar> convs_U = has_unclassified ? class_conversion_factors(mats_U, wt_U) : std::vector<Scalar>{1.0, 1.0, 1.0};

        std::vector<int> branch_ids;
        branch_ids.reserve(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) branch_ids.push_back(node.id);
        }

        // P and dP/dt for every contributing (class, branch), built in parallel over branches.
        // In BranchSite mode (or on unclassified branches) the classes are mixed using br_w:
        //   P_mix = sum_c br_w[c] P_c,  dP_mix = sum_c br_w[c] conv_c Q_c P_c.
        #pragma omp parallel for schedule(static)
        for (size_t bi = 0; bi < branch_ids.size(); ++bi) {
            const int nid = branch_ids[bi];
            const Scalar bl = branch_lengths[nid];
            const int part = branch_partition[nid];
            const std::vector<Scalar>& br_w = (part == 2) ? wt_U : weights;
            if (branch_site || part == 2) {
                Matrix P_mix = Matrix::Zero(S, S);
                Matrix dP_mix = Matrix::Zero(S, S);
                for (int c = 0; c < 3; ++c) {
                    if (br_w[c] <= 1e-12) continue;
                    Matrix P;
                    Matrix dP;
                    if (part == 1) {
                        P = mats_T[c].transition_matrix(bl * convs_T[c]);
                        dP = convs_T[c] * (mats_T[c].Q * P);
                    } else if (part == 2) {
                        P = mats_U[c].transition_matrix(bl * convs_U[c]);
                        dP = convs_U[c] * (mats_U[c].Q * P);
                    } else {
                        P = mats_R[c].transition_matrix(bl * convs_R[c]);
                        dP = convs_R[c] * (mats_R[c].Q * P);
                    }
                    P_mix.noalias() += br_w[c] * P;
                    dP_mix.noalias() += br_w[c] * dP;
                }
                P_c[0][nid] = P_mix;
                dP_c[0][nid] = dP_mix;
                if (!branch_site) {
                    P_c[1][nid] = P_mix;
                    P_c[2][nid] = std::move(P_mix);
                    dP_c[1][nid] = dP_mix;
                    dP_c[2][nid] = std::move(dP_mix);
                }
            } else {
                for (int c = 0; c < 3; ++c) {
                    if (weights[c] <= 1e-12) continue;
                    if (part == 1) {
                        Matrix P = mats_T[c].transition_matrix(bl * convs_T[c]);
                        dP_c[c][nid] = convs_T[c] * (mats_T[c].Q * P);
                        P_c[c][nid] = std::move(P);
                    } else {
                        Matrix P = mats_R[c].transition_matrix(bl * convs_R[c]);
                        dP_c[c][nid] = convs_R[c] * (mats_R[c].Q * P);
                        P_c[c][nid] = std::move(P);
                    }
                }
            }
        }

        std::vector<Scalar> class_w = weights;
        if (branch_site) {
            class_w = {1.0, 0.0, 0.0};
        }

        std::vector<Scalar> pattern_ll(num_patterns, 0.0);
        Matrix pattern_grad = Matrix::Zero(num_nodes, num_patterns);
        int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            std::array<LikelihoodEngine::BatchedInsideOutsideWorkspace, 3> io_c;
            int n_active_ws = branch_site ? 1 : 3;
            for (int c = 0; c < n_active_ws; ++c) io_c[c].init(num_nodes, S, tile_size);
            Eigen::MatrixXd post_w(3, tile_size);
            Matrix dP_D(S, tile_size);

            #pragma omp for schedule(static)
            for (size_t tile = 0; tile < num_tiles; ++tile) {
                size_t p_start = tile * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);

                if (branch_site) {
                    LikelihoodEngine::compute_inside_outside_batch_in_place(
                        tree, aln, p_start, p_end, leaf_to_taxon, P_c[0], eq_freqs, S, io_c[0]
                    );
                    for (Eigen::Index b = 0; b < B; ++b) {
                        size_t p = p_start + b;
                        if (io_c[0].likelihood(b) > 0.0) {
                            pattern_ll[p] = aln.patterns[p].weight * io_c[0].log_likelihood(b);
                        }
                    }
                    auto dP_D_sub = dP_D.leftCols(B);
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t vid = node.id;
                            dP_D_sub.noalias() = dP_c[0][vid] * io_c[0].D[vid].leftCols(B);
                            auto V_sub = io_c[0].V[vid].leftCols(B);
                            for (Eigen::Index b = 0; b < B; ++b) {
                                Scalar L_b = io_c[0].likelihood(b);
                                if (L_b > 0.0) {
                                    size_t p = p_start + b;
                                    pattern_grad(vid, p) = (aln.patterns[p].weight / L_b) * V_sub.col(b).dot(dP_D_sub.col(b));
                                }
                            }
                        }
                    }
                } else {
                    post_w.leftCols(B).setZero();
                    for (int c = 0; c < 3; ++c) {
                        if (class_w[c] <= 1e-12) continue;
                        LikelihoodEngine::compute_inside_outside_batch_in_place(
                            tree, aln, p_start, p_end, leaf_to_taxon, P_c[c], eq_freqs, S, io_c[c]
                        );
                    }
                    for (Eigen::Index b = 0; b < B; ++b) {
                        Scalar max_ll = -std::numeric_limits<Scalar>::infinity();
                        for (int c = 0; c < 3; ++c) {
                            if (class_w[c] > 1e-12 && io_c[c].likelihood(b) > 0.0) {
                                Scalar term = std::log(class_w[c]) + io_c[c].log_likelihood(b);
                                if (term > max_ll) max_ll = term;
                            }
                        }
                        if (max_ll > -1e19) {
                            Scalar sum_exp = 0.0;
                            for (int c = 0; c < 3; ++c) {
                                if (class_w[c] > 1e-12 && io_c[c].likelihood(b) > 0.0) {
                                    Scalar e = std::exp(std::log(class_w[c]) + io_c[c].log_likelihood(b) - max_ll);
                                    post_w(c, b) = e;
                                    sum_exp += e;
                                }
                            }
                            if (sum_exp > 0.0) {
                                size_t p = p_start + b;
                                pattern_ll[p] = aln.patterns[p].weight * (max_ll + std::log(sum_exp));
                                post_w.col(b) /= sum_exp;
                            }
                        }
                    }
                    auto dP_D_sub = dP_D.leftCols(B);
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t vid = node.id;
                            for (int c = 0; c < 3; ++c) {
                                if (class_w[c] <= 1e-12) continue;
                                dP_D_sub.noalias() = dP_c[c][vid] * io_c[c].D[vid].leftCols(B);
                                auto V_sub = io_c[c].V[vid].leftCols(B);
                                for (Eigen::Index b = 0; b < B; ++b) {
                                    Scalar pw = post_w(c, b);
                                    Scalar L_b = io_c[c].likelihood(b);
                                    if (pw > 1e-12 && L_b > 0.0) {
                                        size_t p = p_start + b;
                                        pattern_grad(vid, p) += aln.patterns[p].weight * pw * (V_sub.col(b).dot(dP_D_sub.col(b)) / L_b);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Scalar total_log_l = neumaier_sum(pattern_ll);
        std::vector<Scalar> grad_b(num_nodes, 0.0);
        std::vector<Scalar> row_vals(num_patterns, 0.0);
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i == static_cast<size_t>(tree.root_id)) continue;
            for (size_t p = 0; p < num_patterns; ++p) {
                row_vals[p] = pattern_grad(i, p);
            }
            grad_b[i] = neumaier_sum(row_vals);
        }

        return {total_log_l, grad_b};
    }

    // High-performance joint branch length refinement via exact Inside-Outside adjoint gradients and L-BFGS
    std::vector<Scalar> refine_branch_lengths(
        const std::vector<Scalar>& init_bl,
        const std::vector<Scalar>& omegas_R,
        const std::vector<Scalar>& weights,
        Scalar k,
        int max_iters = 10,
        const MG94Parameters* custom_base_params = nullptr,
        const std::vector<Scalar>* custom_omegas_uncl = nullptr,
        const std::vector<Scalar>* custom_weights_uncl = nullptr
    ) const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Scalar> bl = init_bl;
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id)) {
                bl[i] = std::clamp(bl[i], 1e-6, 1000.0);
            }
        }

        const size_t m_history = 5;
        std::vector<std::vector<Scalar>> s_hist;
        std::vector<std::vector<Scalar>> y_hist;
        std::vector<Scalar> rho_hist;

        auto [curr_ll, curr_grad] = compute_relax_branch_gradients(
            bl, omegas_R, weights, k, custom_base_params, custom_omegas_uncl, custom_weights_uncl);

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
                        new_bl[i] = std::clamp(bl[i] + step_size * d[i], 1e-6, 1000.0);
                    } else {
                        new_bl[i] = 0.0;
                    }
                }
                auto [test_ll, test_grad] = compute_relax_branch_gradients(
                    new_bl, omegas_R, weights, k, custom_base_params, custom_omegas_uncl, custom_weights_uncl);
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
    //   w2 = exp(u2), u2 in [0, ln 10000]               -> 1 <= w2 <= 10000
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
    static constexpr Scalar kMaxOmega = 10000.0;

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

    static void theta_bounds(bool has_k, std::vector<Scalar>& lb, std::vector<Scalar>& ub,
                             Scalar min_k = kMinK, Scalar max_k = kMaxK) {
        lb = {0.0, 0.0, 0.0, 0.0, 0.0};
        ub = {1.0, 1.0, std::log(kMaxOmega), 1.0, 1.0};
        if (has_k) {
            lb.push_back(std::log(min_k));
            ub.push_back(std::log(max_k));
        }
    }

    Scalar loglik_at(const RELAXDist& d, const std::vector<Scalar>& bl,
                     std::vector<Scalar>* site_ll = nullptr,
                     const MG94Parameters* custom_base_params = nullptr,
                     const std::vector<Scalar>* custom_omegas_uncl = nullptr,
                     const std::vector<Scalar>* custom_weights_uncl = nullptr) const {
        return evaluate_relax_mixture(tree, d.omegas, d.weights, d.k, 1.0, site_ll, &bl,
                                      custom_base_params, custom_omegas_uncl, custom_weights_uncl);
    }

    // Nelder-Mead over theta at fixed branch lengths, restarted from the incumbent while it
    // keeps improving. `ll` must hold the log-likelihood at `th` on entry; never decreases.
    std::vector<Scalar> optimize_theta(std::vector<Scalar> th, bool has_k,
                                       const std::vector<Scalar>& bl, Scalar& ll,
                                       const MG94Parameters* custom_base_params = nullptr,
                                       const std::vector<Scalar>* custom_omegas_uncl = nullptr,
                                       const std::vector<Scalar>* custom_weights_uncl = nullptr,
                                       Scalar min_k = kMinK, Scalar max_k = kMaxK) const {
        std::vector<Scalar> lb, ub;
        theta_bounds(has_k, lb, ub, min_k, max_k);
        if (has_k) {
            th[5] = std::clamp(th[5], lb[5], ub[5]);
        }
        auto f = [&](const std::vector<Scalar>& x) -> Scalar {
            return -loglik_at(decode_theta(x, has_k), bl, nullptr,
                              custom_base_params, custom_omegas_uncl, custom_weights_uncl);
        };
        Scalar best = f(th);
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
        size_t n_grid = starts.size();
        if (settings.grid_size > 0 && static_cast<size_t>(settings.grid_size) < n_grid) {
            n_grid = static_cast<size_t>(settings.grid_size);
        }
        const size_t n_starts = std::min<size_t>(
            n_grid, static_cast<size_t>(std::max(1, settings.starting_points))
        );
        if (n_starts <= 1) {
            std::vector<std::vector<Scalar>> grid(starts.begin() + 1, starts.begin() + n_grid);
            auto o = NelderMeadDynamic::minimize(f, starts[0], lb, ub, 1e-4, 400, grid);
            auto x = o.x;
            Scalar fbest = o.f;
            auto o2 = NelderMeadDynamic::minimize(f, x, lb, ub, 1e-5, 250);
            if (o2.f < fbest) { x = o2.x; fbest = o2.f; }
            return {x, -fbest};
        }

        // Evaluate all candidate grid points and rank them by initial objective value
        std::vector<std::pair<Scalar, size_t>> ranked;
        ranked.reserve(n_grid);
        for (size_t i = 0; i < n_grid; ++i) {
            ranked.push_back({f(starts[i]), i});
        }
        std::stable_sort(ranked.begin(), ranked.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });

        std::vector<Scalar> best_x = starts[ranked[0].second];
        Scalar best_f = ranked[0].first;
        for (size_t s_idx = 0; s_idx < n_starts; ++s_idx) {
            const auto& seed = starts[ranked[s_idx].second];
            auto o1 = NelderMeadDynamic::minimize(f, seed, lb, ub, 1e-4, 400);
            auto cand_x = o1.x;
            Scalar cand_f = o1.f;
            auto o2 = NelderMeadDynamic::minimize(f, cand_x, lb, ub, 1e-5, 250);
            if (o2.f < cand_f) { cand_x = o2.x; cand_f = o2.f; }
            if (cand_f < best_f) {
                best_f = cand_f;
                best_x = std::move(cand_x);
            }
        }
        return {best_x, -best_f};
    }

    struct RELAXJointFit {
        std::vector<Scalar> theta;
        std::vector<Scalar> bl;
        MG94Parameters fit_base;
        std::vector<Scalar> omegas_uncl = {0.05, 0.50, 2.0};
        std::vector<Scalar> weights_uncl = {0.50, 0.40, 0.10};
        Scalar ll = -std::numeric_limits<Scalar>::infinity();
        int cycles = 0;
        bool converged = false;
    };

    // Alternates per-branch length refinement (L-BFGS, exact gradients) with re-optimization of
    // the distribution (and K), unclassified distribution (if present), and GTR exchangeabilities
    // at the refined lengths, until the cycle gain < joint_tolerance (for cycle >= 2) or
    // max_joint_cycles is reached. Monotone: the log-likelihood never decreases.
    RELAXJointFit joint_optimize(
        std::vector<Scalar> th,
        std::vector<Scalar> bl,
        bool has_k,
        const MG94Parameters* init_base = nullptr,
        const std::vector<Scalar>* init_om_u = nullptr,
        const std::vector<Scalar>* init_wt_u = nullptr,
        Scalar min_k = kMinK,
        Scalar max_k = kMaxK
    ) const {
        RELAXJointFit J;
        J.theta = std::move(th);
        if (has_k) {
            J.theta[5] = std::clamp(J.theta[5], std::log(min_k), std::log(max_k));
        }
        J.bl = std::move(bl);
        J.fit_base = init_base ? *init_base : base_params;
        J.omegas_uncl = init_om_u ? *init_om_u : omegas_uncl;
        J.weights_uncl = init_wt_u ? *init_wt_u : weights_uncl;

        const bool has_unclassified = !unclassified_branch_names.empty();
        J.ll = loglik_at(decode_theta(J.theta, has_k), J.bl, nullptr,
                         &J.fit_base, &J.omegas_uncl, &J.weights_uncl);

        const int cycles = settings.refine_branch_lengths
            ? std::max(1, settings.max_joint_cycles)
            : std::max(3, settings.max_joint_cycles);
        for (int cycle = 0; cycle < cycles; ++cycle) {
            Scalar start_ll = J.ll;
            if (settings.refine_branch_lengths) {
                RELAXDist d = decode_theta(J.theta, has_k);
                auto nbl = refine_branch_lengths(J.bl, d.omegas, d.weights, d.k, 30,
                                                 &J.fit_base, &J.omegas_uncl, &J.weights_uncl);
                Scalar nll = loglik_at(d, nbl, nullptr, &J.fit_base, &J.omegas_uncl, &J.weights_uncl);
                if (nll > J.ll) {
                    J.bl = std::move(nbl);
                    J.ll = nll;
                }
            }
            J.theta = optimize_theta(J.theta, has_k, J.bl, J.ll,
                                     &J.fit_base, &J.omegas_uncl, &J.weights_uncl, min_k, max_k);

            // F-013: Optimize unclassified 3-class BS-REL mixture when unclassified branches exist
            if (has_unclassified) {
                RELAXDist d_cur = decode_theta(J.theta, has_k);
                RELAXDist d_u;
                d_u.omegas = J.omegas_uncl;
                d_u.weights = J.weights_uncl;
                std::vector<Scalar> th_u = encode_theta(d_u, false);
                std::vector<Scalar> lb_u, ub_u;
                theta_bounds(false, lb_u, ub_u);
                auto f_u = [&](const std::vector<Scalar>& x_u) -> Scalar {
                    RELAXDist du_cand = decode_theta(x_u, false);
                    return -loglik_at(d_cur, J.bl, nullptr, &J.fit_base, &du_cand.omegas, &du_cand.weights);
                };
                auto o_u = NelderMeadDynamic::minimize(f_u, th_u, lb_u, ub_u, 1e-4, 200);
                if (-o_u.f > J.ll) {
                    RELAXDist du_best = decode_theta(o_u.x, false);
                    J.omegas_uncl = du_best.omegas;
                    J.weights_uncl = du_best.weights;
                    J.ll = -o_u.f;
                }
            }

            // F-043: Re-optimize the 5 GTR nucleotide exchangeabilities inside joint_optimize
            // (active under both refine_branch_lengths and --no-branch-opt)
            {
                RELAXDist d_cur = decode_theta(J.theta, has_k);
                const MG94Parameters saved_base = J.fit_base;
                auto with_rates = [&](const std::array<Scalar, 5>& x) {
                    MG94Parameters p = saved_base;
                    p.theta_AC = std::exp(x[0]);
                    p.theta_AT = std::exp(x[1]);
                    p.theta_CG = std::exp(x[2]);
                    p.theta_CT = std::exp(x[3]);
                    p.theta_GT = std::exp(x[4]);
                    return p;
                };
                std::array<Scalar, 5> x0 = {
                    std::log(std::clamp(saved_base.theta_AC, 1e-4, 1e4)),
                    std::log(std::clamp(saved_base.theta_AT, 1e-4, 1e4)),
                    std::log(std::clamp(saved_base.theta_CG, 1e-4, 1e4)),
                    std::log(std::clamp(saved_base.theta_CT, 1e-4, 1e4)),
                    std::log(std::clamp(saved_base.theta_GT, 1e-4, 1e4))
                };
                const Scalar lo = std::log(1e-4), hi = std::log(1e4);
                auto eval_nuc = [&](const std::array<Scalar, 5>& x) -> Scalar {
                    MG94Parameters cand_p = with_rates(x);
                    return -loglik_at(d_cur, J.bl, nullptr, &cand_p, &J.omegas_uncl, &J.weights_uncl);
                };
                auto o5 = NelderMeadND<5>::minimize(eval_nuc, x0, {lo, lo, lo, lo, lo}, {hi, hi, hi, hi, hi}, 1e-4, 200);
                if (-o5.f > J.ll) {
                    J.fit_base = with_rates(o5.x);
                    J.ll = -o5.f;
                }
            }

            J.cycles = cycle + 1;
            Scalar gain = J.ll - start_ll;
            const int min_conv_cycle = settings.refine_branch_lengths ? 2 : 1;
            if (cycle >= min_conv_cycle && gain < settings.joint_tolerance) {
                J.converged = true;
                break;
            }
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
        // branch lengths; here rates and branch lengths are alternated to convergence to provide
        // high-accuracy starting values before Phase 2/3/4 refine them under the codon models.
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
            const Scalar lo = std::log(1e-4), hi = std::log(1e4);
            for (int cycle = 0; cycle < 10; ++cycle) {
                const Scalar start_ll = cur;
                std::array<Scalar, 5> x0 = {
                    std::log(std::clamp(gtr_res.params.theta_AC, 1e-4, 1e4)), std::log(std::clamp(gtr_res.params.theta_AT, 1e-4, 1e4)),
                    std::log(std::clamp(gtr_res.params.theta_CG, 1e-4, 1e4)), std::log(std::clamp(gtr_res.params.theta_CT, 1e-4, 1e4)),
                    std::log(std::clamp(gtr_res.params.theta_GT, 1e-4, 1e4))
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
        const bool has_unclassified = !unclassified_branch_names.empty();
        const std::vector<Scalar> one_class = {1.0, 0.0, 0.0};
        Scalar mg_wU = unclassified_omega;
        std::vector<Scalar> mg_om_u = {mg_wU, mg_wU, mg_wU};

        Scalar mg_wR = 0.1;
        Scalar mg_wT = 0.5;
        Scalar sep_scale = 1.0;
        Scalar mg_ll = -1e20;
        if (has_unclassified) {
            auto eval_sep_mg94_u = [&](const std::array<Scalar, 4>& x) -> Scalar {
                std::vector<Scalar> cand_u = {x[2], x[2], x[2]};
                return -evaluate_relax_mixture(gtr_res.tree, {x[0], x[0], x[0]}, one_class, mg94_k(x[0], x[1]), x[3],
                                               nullptr, nullptr, &base_params, &cand_u, &one_class);
            };
            auto opt_sep = NelderMeadND<4>::minimize(
                eval_sep_mg94_u,
                {0.1, 0.5, 0.3, 1.0},
                {1e-4, 1e-4, 1e-4, 0.1},
                {kMaxOmega, kMaxOmega, kMaxOmega, 10.0},
                1e-4, 120
            );
            mg_wR = opt_sep.x[0];
            mg_wT = opt_sep.x[1];
            mg_wU = opt_sep.x[2];
            sep_scale = opt_sep.x[3];
            mg_ll = -opt_sep.f;
            unclassified_omega = mg_wU;
            mg_om_u = {mg_wU, mg_wU, mg_wU};
        } else {
            auto eval_sep_mg94 = [&](const std::array<Scalar, 3>& x) -> Scalar {
                return -evaluate_relax_mixture(gtr_res.tree, {x[0], x[0], x[0]}, one_class, mg94_k(x[0], x[1]), x[2],
                                               nullptr, nullptr, &base_params, &mg_om_u, &one_class);
            };
            auto opt_sep = NelderMeadND<3>::minimize(
                eval_sep_mg94,
                {0.1, 0.5, 1.0},
                {1e-4, 1e-4, 0.1},
                {kMaxOmega, kMaxOmega, 10.0},
                1e-4, 100
            );
            mg_wR = opt_sep.x[0];
            mg_wT = opt_sep.x[1];
            sep_scale = opt_sep.x[2];
            mg_ll = -opt_sep.f;
        }
        std::vector<Scalar> mg_bl = scale_bl(gtr_bl, sep_scale);

        // Per-branch refinement of the MG94 fit, alternating with (omega_R, omega_T[, omega_U]).
        // The RELAX fits below start from these branch lengths.
        if (settings.refine_branch_lengths) {
            update_progress("Phase 2: Refining MG94 branch lengths", 0.35);
            const int cycles = std::max(1, settings.max_joint_cycles);
            for (int cycle = 0; cycle < cycles; ++cycle) {
                Scalar start_ll = mg_ll;
                Scalar kk = mg94_k(mg_wR, mg_wT);
                auto nbl = refine_branch_lengths(mg_bl, {mg_wR, mg_wR, mg_wR}, one_class, kk, 30,
                                                 &base_params, &mg_om_u, &one_class);
                Scalar nll = evaluate_relax_mixture(tree, {mg_wR, mg_wR, mg_wR}, one_class, kk, 1.0, nullptr, &nbl,
                                                    &base_params, &mg_om_u, &one_class);
                if (nll > mg_ll) {
                    mg_bl = std::move(nbl);
                    mg_ll = nll;
                }
                if (has_unclassified) {
                    auto eval_w3 = [&](const std::array<Scalar, 3>& x) -> Scalar {
                        std::vector<Scalar> cand_u = {x[2], x[2], x[2]};
                        return -evaluate_relax_mixture(tree, {x[0], x[0], x[0]}, one_class, mg94_k(x[0], x[1]), 1.0, nullptr, &mg_bl,
                                                       &base_params, &cand_u, &one_class);
                    };
                    auto o3 = NelderMeadND<3>::minimize(eval_w3, {mg_wR, mg_wT, mg_wU}, {1e-4, 1e-4, 1e-4}, {kMaxOmega, kMaxOmega, kMaxOmega}, 1e-5, 200);
                    if (-o3.f > mg_ll) {
                        mg_wR = o3.x[0];
                        mg_wT = o3.x[1];
                        mg_wU = o3.x[2];
                        unclassified_omega = mg_wU;
                        mg_om_u = {mg_wU, mg_wU, mg_wU};
                        mg_ll = -o3.f;
                    }
                } else {
                    auto eval_w = [&](const std::array<Scalar, 2>& x) -> Scalar {
                        return -evaluate_relax_mixture(tree, {x[0], x[0], x[0]}, one_class, mg94_k(x[0], x[1]), 1.0, nullptr, &mg_bl,
                                                       &base_params, &mg_om_u, &one_class);
                    };
                    auto o2 = NelderMeadND<2>::minimize(eval_w, {mg_wR, mg_wT}, {1e-4, 1e-4}, {kMaxOmega, kMaxOmega}, 1e-5, 200);
                    if (-o2.f > mg_ll) {
                        mg_wR = o2.x[0];
                        mg_wT = o2.x[1];
                        mg_ll = -o2.f;
                    }
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
                        return -evaluate_relax_mixture(tree, {mg_wR, mg_wR, mg_wR}, one_class, kk2, 1.0, nullptr, &mg_bl,
                                                       &base_params, &mg_om_u, &one_class);
                    };
                    std::array<Scalar, 5> x0 = {
                        std::log(std::clamp(saved.theta_AC, 1e-4, 1e4)), std::log(std::clamp(saved.theta_AT, 1e-4, 1e4)),
                        std::log(std::clamp(saved.theta_CG, 1e-4, 1e4)), std::log(std::clamp(saved.theta_CT, 1e-4, 1e4)),
                        std::log(std::clamp(saved.theta_GT, 1e-4, 1e4))
                    };
                    const Scalar lo = std::log(1e-4), hi = std::log(1e4);
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
        result.mg94_parameters = 5 + B + (has_unclassified ? 3 : 2) + 6; // GTR rates + branches + partition omegas + CF3x4
        Scalar n_codons = static_cast<Scalar>(aln.num_codons);
        Scalar k_mg = static_cast<Scalar>(result.mg94_parameters);
        result.mg94_aicc = 2.0 * k_mg - 2.0 * result.mg94_log_likelihood +
            (2.0 * k_mg * (k_mg + 1.0)) / std::max(1.0, n_codons - k_mg - 1.0);

        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                result.mg94_branch_lengths[node.name] = mg_bl[node.id];
            }
        }

        // Seed unclassified distribution from dataset Phase-2 MG94 unclassified omega
        Scalar w_u_raw = has_unclassified ? mg_wU : (0.5 * (mg_wR + mg_wT));
        Scalar w_mean = std::clamp(w_u_raw, 1e-3, 0.95);
        omegas_uncl = {std::clamp(0.1 * w_mean, 0.0, 0.9), w_mean, std::max({Scalar(1.5), Scalar(2.0) * w_mean, w_u_raw})};
        weights_uncl = {0.50, 0.40, 0.10};

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

        auto halton = [](size_t idx, size_t base) -> Scalar {
            Scalar f = 1.0;
            Scalar r = 0.0;
            size_t i = idx;
            while (i > 0) {
                f /= static_cast<Scalar>(base);
                r += f * static_cast<Scalar>(i % base);
                i /= base;
            }
            return r;
        };

        // F-014: Expanded data-scaled starting grid for null model
        std::vector<std::vector<Scalar>> null_starts = {
            with_scale(1.0, encode_theta(make_dist(0.0, 0.5, 1.5, 0.6, 0.3), false)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.6, 1.2, 0.67, 0.27), false)),
            with_scale(1.0, encode_theta(make_dist(0.2 * w_mg, w_mg, 2.0, 0.3, 0.6), false)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.2, 3.0, 0.5, 0.45), false)),
            with_scale(1.0, encode_theta(make_dist(0.05 * w_mg, 0.8 * w_mg, std::max(1.2, 3.0 * mg_wR), 0.50, 0.40), false)),
            with_scale(1.0, encode_theta(make_dist(0.1 * w_mean, w_mean, 5.0, 0.45, 0.45), false)),
            with_scale(1.0, encode_theta(make_dist(0.01, 0.15, 1.1, 0.70, 0.25), false)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.8, 10.0, 0.25, 0.65), false))
        };
        for (size_t g_idx = null_starts.size(); static_cast<int>(g_idx) < settings.grid_size; ++g_idx) {
            size_t h_id = g_idx + 1;
            Scalar w1_h = std::clamp(0.02 + 0.96 * halton(h_id, 2), 1e-3, 0.99);
            Scalar w0_h = std::clamp(w1_h * halton(h_id, 3), 0.0, w1_h);
            Scalar w2_h = std::clamp(std::exp(std::log(1.01) + std::log(kMaxOmega / 1.01) * halton(h_id, 5)), 1.01, kMaxOmega);
            Scalar p0_h = std::clamp(0.05 + 0.80 * halton(h_id, 7), 0.01, 0.90);
            Scalar p1_h = std::clamp((0.96 - p0_h) * (0.10 + 0.85 * halton(h_id, 11)), 0.01, 0.97 - p0_h);
            Scalar s_h = std::exp(std::log(0.25) + std::log(16.0) * halton(h_id, 17));
            null_starts.push_back(with_scale(s_h, encode_theta(make_dist(w0_h, w1_h, w2_h, p0_h, p1_h), false)));
        }
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

        // F-014: Expanded data-scaled starting grid covering K in {0.0, 0.05, 0.25, 0.5, 0.8, 1.0, 1.25, 2.0, 5.0, 15.0, 50.0}
        RELAXDist d_null_mle = decode_theta(J_null.theta, false);
        Scalar w0_n = std::clamp(d_null_mle.omegas[0], 0.0, 0.9);
        Scalar w1_n = std::clamp(d_null_mle.omegas[1], 1e-3, 0.99);
        Scalar w2_n = std::clamp(d_null_mle.omegas[2], 1.01, kMaxOmega);
        Scalar p0_n = std::clamp(d_null_mle.weights[0], 0.01, 0.95);
        Scalar p1_n = std::clamp(d_null_mle.weights[1], 0.01, 0.98 - p0_n);

        std::vector<std::vector<Scalar>> alt_starts = {
            with_scale(1.0, encode_theta(make_dist(0.0, w_mg, 5.0, 0.35, 0.60, k_init), true)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.16, 30.0, 0.35, 0.65, 0.0), true)),
            with_scale(1.0, encode_theta(make_dist(0.05, 0.3, 2.0, 0.40, 0.50, 0.05), true)),
            with_scale(1.0, encode_theta(make_dist(0.1 * w_mg, w_mg, 2.5, 0.45, 0.45, 0.25), true)),
            with_scale(1.0, encode_theta(make_dist(0.05, 0.3, 2.0, 0.40, 0.50, 0.5), true)),
            with_scale(1.0, encode_theta(make_dist(w0_n, w1_n, w2_n, p0_n, p1_n, 0.8), true)),
            with_scale(1.0, encode_theta(make_dist(w0_n, w1_n, w2_n, p0_n, p1_n, 1.0), true)),
            with_scale(1.0, encode_theta(make_dist(w0_n, w1_n, w2_n, p0_n, p1_n, 1.25), true)),
            with_scale(1.0, encode_theta(make_dist(0.0, 0.5, 5.0, 0.3, 0.4, 2.0), true)),
            with_scale(1.0, encode_theta(make_dist(0.1 * w_mg, w_mg, 3.0, 0.40, 0.45, 5.0), true)),
            with_scale(1.0, encode_theta(make_dist(0.05, 0.4, 2.5, 0.50, 0.40, 15.0), true)),
            with_scale(1.0, encode_theta(make_dist(w0_n, w1_n, w2_n, p0_n, p1_n, kMaxK), true))
        };
        for (size_t g_idx = alt_starts.size(); static_cast<int>(g_idx) < settings.grid_size; ++g_idx) {
            size_t h_id = g_idx + 1;
            Scalar w1_h = std::clamp(0.02 + 0.96 * halton(h_id, 2), 1e-3, 0.99);
            Scalar w0_h = std::clamp(w1_h * halton(h_id, 3), 0.0, w1_h);
            Scalar w2_h = std::clamp(std::exp(std::log(1.01) + std::log(kMaxOmega / 1.01) * halton(h_id, 5)), 1.01, kMaxOmega);
            Scalar p0_h = std::clamp(0.05 + 0.80 * halton(h_id, 7), 0.01, 0.90);
            Scalar p1_h = std::clamp((0.96 - p0_h) * (0.10 + 0.85 * halton(h_id, 11)), 0.01, 0.97 - p0_h);
            Scalar u_k = halton(h_id, 13);
            Scalar k_h = (u_k < 0.05) ? kMinK : std::clamp(std::exp(std::log(1e-3) + std::log(kMaxK / 1e-3) * ((u_k - 0.05) / 0.95)), kMinK, kMaxK);
            Scalar s_h = std::exp(std::log(0.25) + std::log(16.0) * halton(h_id, 17));
            alt_starts.push_back(with_scale(s_h, encode_theta(make_dist(w0_h, w1_h, w2_h, p0_h, p1_h, k_h), true)));
        }
        auto scaled_alt = optimize_scaled(alt_starts, true, mg_bl);
        Scalar alt_s = scaled_alt.first[0];
        RELAXJointFit J_alt = joint_optimize(
            std::vector<Scalar>(scaled_alt.first.begin() + 1, scaled_alt.first.end()),
            scale_bl(mg_bl, alt_s), true);

        // Start B: the null MLE with K = 1 (guarantees alt >= null for a monotone optimizer)
        update_progress(with_stats("Phase 4: Alternative model from null MLE; grid start reached", J_alt.ll), 0.75);
        RELAXJointFit J_alt_b = joint_optimize(
            add_k(J_null.theta, 1.0), J_null.bl, true,
            &J_null.fit_base, &J_null.omegas_uncl, &J_null.weights_uncl);
        if (J_alt_b.ll > J_alt.ll) J_alt = std::move(J_alt_b);

        // HyPhy-style null check: re-fit the null from the alternative MLE with K := 1.
        update_progress(with_stats("Phase 4: Null model from alternative MLE; alternative at", J_alt.ll), 0.82);
        RELAXJointFit J_null_b = joint_optimize(
            drop_k(J_alt.theta), J_alt.bl, false,
            &J_alt.fit_base, &J_alt.omegas_uncl, &J_alt.weights_uncl);
        if (J_null_b.ll > J_null.ll) J_null = std::move(J_null_b);

        // F-047: Opposite-side K direction stability check when |alt_ll - null_ll| < 5.0
        {
            RELAXDist d_hat = decode_theta(J_alt.theta, true);
            if (std::abs(J_alt.ll - J_null.ll) < 5.0) {
                std::vector<std::tuple<Scalar, Scalar, Scalar>> probes;
                if (d_hat.k > 1.0 + 1e-6) {
                    probes.emplace_back(kMinK, 1.0, 0.5);
                } else if (d_hat.k < 1.0 - 1e-6) {
                    probes.emplace_back(1.0, kMaxK, 2.0);
                } else {
                    probes.emplace_back(kMinK, 1.0, 0.5);
                    probes.emplace_back(1.0, kMaxK, 2.0);
                }
                for (const auto& [opp_min_k, opp_max_k, seed_k] : probes) {
                    RELAXJointFit J_opp = joint_optimize(
                        add_k(J_null.theta, seed_k), J_null.bl, true,
                        &J_null.fit_base, &J_null.omegas_uncl, &J_null.weights_uncl,
                        opp_min_k, opp_max_k);
                    if (J_opp.ll > J_alt.ll) {
                        J_alt = std::move(J_opp);
                        RELAXJointFit J_null_opp = joint_optimize(
                            drop_k(J_alt.theta), J_alt.bl, false,
                            &J_alt.fit_base, &J_alt.omegas_uncl, &J_alt.weights_uncl);
                        if (J_null_opp.ll > J_null.ll) J_null = std::move(J_null_opp);
                    }
                }
            }
        }

        // F-084: Negative-LRT retry loop (up to 5 restarts of alternative from null MLE and vice versa)
        for (int retry = 0; retry < 5 && J_alt.ll < J_null.ll; ++retry) {
            RELAXJointFit J_alt_r = joint_optimize(
                add_k(J_null.theta, 1.0), J_null.bl, true,
                &J_null.fit_base, &J_null.omegas_uncl, &J_null.weights_uncl);
            if (J_alt_r.ll > J_alt.ll) J_alt = std::move(J_alt_r);
            if (J_alt.ll > J_null.ll + 1e-3) {
                RELAXJointFit J_null_r = joint_optimize(
                    drop_k(J_alt.theta), J_alt.bl, false,
                    &J_alt.fit_base, &J_alt.omegas_uncl, &J_alt.weights_uncl);
                if (J_null_r.ll > J_null.ll) J_null = std::move(J_null_r);
            }
        }
        // Final strict nesting guard: since K=1 is in the alternative parameter space, alt_ll >= null_ll
        if (J_alt.ll < J_null.ll) {
            J_alt.theta = add_k(J_null.theta, 1.0);
            J_alt.bl = J_null.bl;
            J_alt.fit_base = J_null.fit_base;
            J_alt.omegas_uncl = J_null.omegas_uncl;
            J_alt.weights_uncl = J_null.weights_uncl;
            J_alt.ll = J_null.ll;
            J_alt.cycles = J_null.cycles;
            J_alt.converged = J_null.converged;
        }

        update_progress(with_stats("Phase 4: RELAX fits complete; alternative", J_alt.ll), 0.90);

        // ---------------------------------------------------------
        // Record fits
        // ---------------------------------------------------------
        const size_t uncl_extra_params = unclassified_branch_names.empty() ? 0 : 5;
        auto record_fit = [&](RELAXModelFit& fit, const RELAXJointFit& J, bool has_k, Scalar s, size_t n_params) {
            RELAXDist d = decode_theta(J.theta, has_k);
            fit.tree_scale = s;
            fit.k = d.k;
            fit.realized_cycles = J.cycles;
            fit.converged = J.converged;
            fit.reference_distribution.omegas = d.omegas;
            fit.reference_distribution.weights = d.weights;
            std::vector<Scalar> om_T(3);
            for (int c = 0; c < 3; ++c) {
                if (d.omegas[c] <= 1e-12) {
                    om_T[c] = (d.k > 0) ? 0.0 : 1.0;
                } else {
                    om_T[c] = std::clamp(std::pow(d.omegas[c], d.k), 1e-6, 10000.0);
                }
            }
            fit.test_distribution.omegas = has_k ? om_T : d.omegas;
            fit.test_distribution.weights = d.weights;
            if (!unclassified_branch_names.empty()) {
                fit.unclassified_distribution.omegas = J.omegas_uncl;
                fit.unclassified_distribution.weights = J.weights_uncl;
            }
            fit.branch_lengths = J.bl;
            fit.log_likelihood = loglik_at(d, J.bl, &fit.site_log_likelihoods,
                                           &J.fit_base, &J.omegas_uncl, &J.weights_uncl);
            fit.parameters = n_params;
            Scalar kp = static_cast<Scalar>(n_params);
            fit.aicc = 2.0 * kp - 2.0 * fit.log_likelihood +
                (2.0 * kp * (kp + 1.0)) / std::max(1.0, n_codons - kp - 1.0);
        };

        // GTR rates + branches + 3 omegas + 2 weights (+ 5 uncl if present) (+ k) + CF3x4
        record_fit(result.alternative_fit, J_alt, true, alt_s, 5 + B + 5 + uncl_extra_params + 1 + 6);
        record_fit(result.null_fit, J_null, false, null_s, 5 + B + 5 + uncl_extra_params + 6);
        Scalar alt_k = result.alternative_fit.k;

        // F-076: Compute branch-switched per-site log-likelihoods at the alternative MLE
        // (evaluating the alternative model with Test branches assigned K = 1, i.e., the Reference
        // distribution, while holding all other alternative-MLE parameters fixed) and store the
        // per-site evidence ratios exp(l_{alt,s} - l_{switched,s}) for legacy JSON Site Log Likelihood.constrained.
        {
            RELAXDist d_switched = decode_theta(J_alt.theta, true);
            d_switched.k = 1.0;
            std::vector<Scalar> switched_site_ll;
            loglik_at(d_switched, J_alt.bl, &switched_site_ll,
                      &J_alt.fit_base, &J_alt.omegas_uncl, &J_alt.weights_uncl);
            const auto& alt_site_ll = result.alternative_fit.site_log_likelihoods;
            if (alt_site_ll.size() == switched_site_ll.size()) {
                result.branch_switched_evidence_ratios.resize(alt_site_ll.size());
                for (size_t s = 0; s < alt_site_ll.size(); ++s) {
                    Scalar diff = std::clamp(alt_site_ll[s] - switched_site_ll[s], -700.0, 700.0);
                    result.branch_switched_evidence_ratios[s] = std::exp(diff);
                }
            }
        }

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
