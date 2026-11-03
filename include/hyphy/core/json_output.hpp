#pragma once

// Shared v3 JSON output infrastructure for all HyPhy 3 analyses.
//
// Every analysis emits the same top-level envelope:
//   $schema, schema_version, analysis, software, provenance,
//   dataset, phylogeny, model_fits, statistical_tests,
//   site_results, branch_results, timers, method_results
//
// The helpers in this header build each shared block so that a new method
// gets the full envelope for free and only contributes:
//   - an AnalysisMetadata literal (id/name/category/description/citations/settings)
//   - model_fits entries (with gtr_fit / fit_entry helpers for common shapes)
//   - its statistical_tests sub-block
//   - site_results (via ColumnarResults, which enforces equal-length columns)
//   - branch_results records and method_results payloads
//
// Adding a new method should never require re-implementing the envelope.

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/provenance.hpp"
#include "nlohmann/json.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <functional>
#include <cmath>
#include <stdexcept>

// Build-time metadata injected by CMake. Defaults keep non-CMake builds working.
#ifndef HYPHY3_VERSION
#define HYPHY3_VERSION "3.0.0"
#endif
#ifndef HYPHY3_GIT_COMMIT
#define HYPHY3_GIT_COMMIT "unknown"
#endif
#ifndef HYPHY3_BUILD_TYPE
#define HYPHY3_BUILD_TYPE "unknown"
#endif

namespace hyphy::core {

inline constexpr const char* kJSONSchemaVersion = "3.0.0";
inline constexpr const char* kSchemaBaseURL =
    "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/";

inline std::string schema_url_for(const std::string& method_id) {
    return std::string(kSchemaBaseURL) + method_id + ".v3.schema.json";
}

// ---------------------------------------------------------------------------
// Number hygiene: JSON has no NaN/Infinity. Emit null for non-finite values so
// that "not computed" is distinguishable from a real number and never leaks
// as a bare literal that a downstream parser chokes on.
// ---------------------------------------------------------------------------
inline nlohmann::json finite_or_null(double v) {
    return std::isfinite(v) ? nlohmann::json(v) : nlohmann::json(nullptr);
}

inline nlohmann::json sanitized_doubles(const std::vector<double>& vals) {
    nlohmann::json arr = nlohmann::json::array();
    arr.get_ref<nlohmann::json::array_t&>().reserve(vals.size());
    for (double v : vals) {
        arr.push_back(finite_or_null(v));
    }
    return arr;
}

// ---------------------------------------------------------------------------
// Analysis metadata block
// ---------------------------------------------------------------------------
struct AnalysisMetadata {
    std::string id;          // machine id: "fel", "meme", ... (schema enum-pinned)
    std::string name;        // display name
    std::string version;
    std::string category;    // e.g. "site_selection", "gene_selection", "branch_selection", "branch_set_selection"
    std::string description;
    std::vector<nlohmann::json> citations; // [{"citation","doi","pmid"}]
    nlohmann::json settings = nlohmann::json::object();

    nlohmann::json to_json() const {
        nlohmann::json j;
        j["id"] = id;
        j["name"] = name;
        j["version"] = version;
        j["category"] = category;
        j["description"] = description;
        if (!citations.empty()) j["citations"] = citations;
        if (!settings.empty()) j["settings"] = settings;
        return j;
    }
};

inline nlohmann::json citation(const std::string& text,
                               const std::string& doi = "",
                               const std::string& pmid = "") {
    nlohmann::json c{{"citation", text}};
    if (!doi.empty()) c["doi"] = doi;
    if (!pmid.empty()) c["pmid"] = pmid;
    return c;
}

// ---------------------------------------------------------------------------
// Software block — version/commit/build come from CMake compile definitions,
// so they can never drift from the binary that produced the output.
// ---------------------------------------------------------------------------
inline nlohmann::json software_block() {
    return {
        {"name", "hyphy"},
        {"version", HYPHY3_VERSION},
        {"git_commit", HYPHY3_GIT_COMMIT},
        {"build_type", HYPHY3_BUILD_TYPE},
        {"compiler", Provenance::detect_compiler()}
    };
}

// ---------------------------------------------------------------------------
// Provenance block.
//
// If the caller supplies a real Provenance (the CLI does), use it verbatim.
// Otherwise build an *honest* auto-provenance: real timestamp, host, cwd and
// input-file hash — but cli_command stays empty rather than fabricating a
// command line that never ran (e.g. library/Python callers). An input entry
// is only emitted when hashing actually succeeded, so a missing file can
// never produce an empty sha256 that violates the schema pattern.
// ---------------------------------------------------------------------------
inline nlohmann::json provenance_block(const Provenance& prov,
                                       const std::string& input_filepath = "",
                                       const std::string& input_format = "FASTA/NEXUS") {
    if (!prov.invocation.cli_command.empty() || !prov.inputs.empty()) {
        return prov.to_json();
    }
    Provenance auto_prov;
    auto_prov.invocation.working_directory = Provenance::get_cwd();
    if (!input_filepath.empty()) {
        size_t sz = 0;
        std::string hash = crypto::SHA256::hash_file(input_filepath, &sz);
        if (!hash.empty()) {
            auto_prov.inputs["alignment"] = {input_filepath, input_format, hash, sz};
        }
    }
    auto_prov.execution.start_time = Provenance::current_iso8601();
    auto_prov.execution.end_time = auto_prov.execution.start_time;
    auto_prov.execution.wall_time_seconds = 0.0;
#ifdef _OPENMP
    auto_prov.execution.cpu_threads = omp_get_max_threads();
#else
    auto_prov.execution.cpu_threads = 1;
#endif
    auto_prov.execution.hostname = Provenance::get_hostname();
    auto_prov.execution.os = Provenance::detect_os();
    auto_prov.execution.compiler = Provenance::detect_compiler();
    return auto_prov.to_json();
}

// ---------------------------------------------------------------------------
// Dataset block — identical shape for every codon-alignment method.
// ---------------------------------------------------------------------------
inline nlohmann::json dataset_block(const Alignment& aln) {
    nlohmann::json j;
    j["taxa_count"] = aln.num_taxa;
    j["codon_sites"] = aln.num_codons;
    j["nucleotide_sites"] = aln.num_codons * 3;
    j["unique_patterns"] = aln.patterns.size();
    j["taxa"] = aln.taxon_names;

    nlohmann::json gcode;
    gcode["id"] = aln.code ? aln.code->name : "Universal";
    gcode["name"] = aln.code ? aln.code->name : "Universal";
    gcode["sense_codons"] = aln.code ? static_cast<int>(aln.code->sense_codons.size()) : 61;
    if (aln.code) {
        gcode["stop_codons"] = aln.code->stop_codons;
    } else {
        gcode["stop_codons"] = {"TAA", "TAG", "TGA"};
    }
    j["genetic_code"] = gcode;

    j["partitions"] = nlohmann::json::array({
        {
            {"id", "default"},
            {"name", "Full Alignment"},
            {"span", nlohmann::json::array({nlohmann::json::array({1, aln.num_codons})})},
            {"sites_count", aln.num_codons},
            {"patterns_count", aln.patterns.size()}
        }
    });
    return j;
}

// ---------------------------------------------------------------------------
// Phylogeny block.
//
// branch_lengths_fn maps a node to {fit_id: length} pairs (e.g. "codon_mg94",
// "nucleotide_gtr", "unconstrained", "constrained"). tree_lengths is then the
// per-key sum over non-root nodes — computed once, identically, for all methods.
//
// Node-name discipline: nodes are keyed by display name, so a collision would
// silently overwrite an earlier entry. Rather than emit corrupt output, any
// duplicate name is a hard error. Unnamed nodes (possible only for trees built
// without the parser, which auto-names internals "Node{n}") fall back to
// "Node{id}" before the uniqueness check runs.
// ---------------------------------------------------------------------------
using BranchLengthMapFn = std::function<std::map<std::string, Scalar>(const TreeNode&)>;

// Optional per-node extras, merged into the node object (e.g. {"tested": true}).
using NodeExtrasFn = std::function<nlohmann::json(const TreeNode&)>;

inline nlohmann::json phylogeny_block(const Tree& tree,
                                      const std::string& newick_override,
                                      BranchLengthMapFn branch_lengths_fn,
                                      NodeExtrasFn node_extras_fn = nullptr) {
    nlohmann::json j;
    j["newick"] = newick_override.empty() ? tree.to_newick() : newick_override;

    nlohmann::json nodes = nlohmann::json::object();
    std::map<std::string, Scalar> totals;

    for (const auto& node : tree.nodes) {
        if (node.id != tree.root_id && node.parent_id == INVALID_INDEX) {
            continue;
        }
        const std::string key =
            node.name.empty() ? ("Node" + std::to_string(node.id)) : node.name;
        if (nodes.contains(key)) {
            throw std::runtime_error(
                "json_output: duplicate node name '" + key +
                "' in phylogeny block; node names must be unique");
        }

        nlohmann::json entry;
        entry["type"] = node.is_leaf ? "leaf" : "internal";
        if (!node.model_tag.empty()) {
            entry["model_tags"] = nlohmann::json::array({node.model_tag});
        }
        if (branch_lengths_fn) {
            auto lengths = branch_lengths_fn(node);
            if (!lengths.empty()) {
                nlohmann::json bl = nlohmann::json::object();
                for (const auto& [fit_id, len] : lengths) {
                    bl[fit_id] = finite_or_null(len);
                    if (node.id != tree.root_id && std::isfinite(len)) {
                        totals[fit_id] += len;
                    }
                }
                entry["branch_lengths"] = bl;
            }
        }
        if (node_extras_fn) {
            nlohmann::json extras = node_extras_fn(node);
            if (extras.is_object()) {
                entry.update(extras);
            }
        }
        nodes[key] = entry;
    }
    j["nodes"] = nodes;

    nlohmann::json tl = nlohmann::json::object();
    for (const auto& [fit_id, sum] : totals) {
        tl[fit_id] = finite_or_null(sum);
    }
    j["tree_lengths"] = tl;
    return j;
}

// Shared branch-length map for MG94xREV-based site methods (FEL, MEME):
// the codon-unit branch length is the nucleotide length rescaled by the
// GTR->codon conversion factor; when a nucleotide GTR fit exists its fitted
// branch lengths are reported alongside.
inline BranchLengthMapFn mg94_branch_lengths_fn(
    Scalar conversion_factor,
    bool has_gtr_fit,
    const std::unordered_map<std::string, Scalar>& gtr_branch_lengths) {
    return [conversion_factor, has_gtr_fit, &gtr_branch_lengths](const TreeNode& node) {
        std::map<std::string, Scalar> lengths;
        const Scalar mg_bl = (conversion_factor > 0.0)
                                 ? (node.branch_length / conversion_factor)
                                 : node.branch_length;
        lengths["codon_mg94"] = mg_bl;
        if (has_gtr_fit) {
            auto it = gtr_branch_lengths.find(node.name);
            lengths["nucleotide_gtr"] =
                (it != gtr_branch_lengths.end()) ? it->second : mg_bl;
        }
        return lengths;
    };
}

// ---------------------------------------------------------------------------
// Model fit entries — shared shapes so "a fit" means the same thing everywhere.
// ---------------------------------------------------------------------------
inline nlohmann::json fit_entry(const std::string& description,
                                Scalar log_likelihood,
                                Scalar aicc,
                                int parameters_count) {
    return {
        {"description", description},
        {"log_likelihood", finite_or_null(log_likelihood)},
        {"aicc", finite_or_null(aicc)},
        {"parameters_count", parameters_count}
    };
}

// The nucleotide GTR baseline fit — emitted identically by every method.
inline nlohmann::json nucleotide_gtr_fit(const Alignment& aln,
                                         const GTRParameters& rates,
                                         Scalar log_likelihood,
                                         Scalar aicc,
                                         int parameters_count) {
    nlohmann::json f = fit_entry(
        "General Time Reversible nucleotide model with empirical frequencies",
        log_likelihood, aicc, parameters_count);
    f["equilibrium_frequencies"] = {
        {"A", finite_or_null(aln.nuc_frequencies(0))},
        {"C", finite_or_null(aln.nuc_frequencies(1))},
        {"G", finite_or_null(aln.nuc_frequencies(2))},
        {"T", finite_or_null(aln.nuc_frequencies(3))}
    };
    f["substitution_rates"] = {
        {"AC", finite_or_null(rates.theta_AC)},
        {"AG", 1.0},
        {"AT", finite_or_null(rates.theta_AT)},
        {"CG", finite_or_null(rates.theta_CG)},
        {"CT", finite_or_null(rates.theta_CT)},
        {"GT", finite_or_null(rates.theta_GT)}
    };
    return f;
}

// A discrete omega distribution, shared by BUSTED/aBSREL/RELAX rate shapes:
// [{"omega": ..., "weight": ..., "annotation"?}, ...] — record-oriented, no
// {"0": {...}} index indirection.
inline nlohmann::json omega_distribution(const std::vector<Scalar>& omegas,
                                         const std::vector<Scalar>& weights,
                                         const std::vector<std::string>& annotations = {}) {
    if (omegas.size() != weights.size() ||
        (!annotations.empty() && annotations.size() != omegas.size())) {
        throw std::runtime_error(
            "json_output: omega_distribution length mismatch (omegas=" +
            std::to_string(omegas.size()) + ", weights=" +
            std::to_string(weights.size()) + ", annotations=" +
            std::to_string(annotations.size()) + ")");
    }
    nlohmann::json dist = nlohmann::json::array();
    const size_t k = omegas.size();
    for (size_t i = 0; i < k; ++i) {
        nlohmann::json entry = {
            {"omega", finite_or_null(omegas[i])},
            {"weight", finite_or_null(weights[i])}
        };
        if (i < annotations.size() && !annotations[i].empty()) {
            entry["annotation"] = annotations[i];
        }
        dist.push_back(entry);
    }
    return dist;
}

// ---------------------------------------------------------------------------
// Columnar results table — the site_results shape.
//
// All columns are added through this class, which enforces the equal-length
// invariant the schema cannot express: a mismatched column is a hard error at
// emission time, not a corrupt file.
// ---------------------------------------------------------------------------
class ColumnarResults {
public:
    struct ColumnMeta {
        std::string type;         // "integer" | "number" | "string"
        std::string description;
        std::string unit = "";    // optional

        nlohmann::json to_json() const {
            nlohmann::json m{{"type", type}, {"description", description}};
            if (!unit.empty()) m["unit"] = unit;
            return m;
        }
    };

    void add_column(const std::string& name, ColumnMeta meta, std::vector<int> data) {
        check_length(name, data.size());
        columns_[name] = meta.to_json();
        data_[name] = std::move(data);
    }
    void add_column(const std::string& name, ColumnMeta meta, std::vector<double> data) {
        check_length(name, data.size());
        columns_[name] = meta.to_json();
        data_[name] = sanitized_doubles(data);
    }
    void add_column(const std::string& name, ColumnMeta meta, std::vector<std::string> data) {
        check_length(name, data.size());
        columns_[name] = meta.to_json();
        data_[name] = std::move(data);
    }
    void add_column(const std::string& name, ColumnMeta meta, std::vector<bool> data) {
        check_length(name, data.size());
        columns_[name] = meta.to_json();
        data_[name] = std::move(data);
    }
    // Per-row JSON values (e.g. a map of branch -> EBF). data must be an array
    // whose length equals row_count.
    void add_column(const std::string& name, ColumnMeta meta, nlohmann::json data) {
        if (!data.is_array()) {
            throw std::runtime_error("json_output: column '" + name +
                                     "' data must be a JSON array");
        }
        check_length(name, data.size());
        columns_[name] = meta.to_json();
        data_[name] = std::move(data);
    }

    size_t row_count() const { return row_count_; }

    nlohmann::json to_json() const {
        return {
            {"row_count", static_cast<int>(row_count_)},
            {"columns", columns_},
            {"data", data_}
        };
    }

private:
    void check_length(const std::string& name, size_t n) {
        if (!has_columns_) {
            row_count_ = n;
            has_columns_ = true;
            return;
        }
        if (n != row_count_) {
            throw std::runtime_error(
                "json_output: column '" + name + "' has " + std::to_string(n) +
                " rows; expected " + std::to_string(row_count_));
        }
    }

    nlohmann::json columns_ = nlohmann::json::object();
    nlohmann::json data_ = nlohmann::json::object();
    size_t row_count_ = 0;
    bool has_columns_ = false;
};

// ---------------------------------------------------------------------------
// Envelope: every method starts here. Fills $schema, schema_version, analysis,
// software, provenance — the five blocks no method should hand-build.
// ---------------------------------------------------------------------------
inline nlohmann::json new_analysis_json(const AnalysisMetadata& meta,
                                        const Provenance& prov = {},
                                        const std::string& input_filepath = "") {
    nlohmann::json j;
    j["$schema"] = schema_url_for(meta.id);
    j["schema_version"] = kJSONSchemaVersion;
    j["analysis"] = meta.to_json();
    j["software"] = software_block();
    j["provenance"] = provenance_block(prov, input_filepath);
    return j;
}

} // namespace hyphy::core
