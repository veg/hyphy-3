#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/ohmic_network.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/analyses/tree_engine.hpp"
#include "hyphy/analyses/joint_align.hpp"
#if __has_include("hyphy/report/viewer_assets.hpp")
// viewer_assets.hpp is generated into the build tree by hyphy3_viewer_assets and
// is only on the include path of targets linking hyphy3_report. html_report.hpp
// itself is static, so probing it would always succeed even for consumers that
// cannot resolve its generated dependency.
#include "hyphy/report/html_report.hpp"
#define HYPHY3_HAS_HTML_REPORT 1
#endif

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::core {

/**
 * @brief Shared CLI options across HyPhy 3 analysis executables.
 */
struct CommonCLIOptions {
    std::string alignment_file;
    std::string tree_file;
    std::string output_file;
    std::string html_file;
    std::string code_name = "Universal";
    std::string branches = "All";
    std::string branch_regex;
    std::string embeddings_file;
    std::string out_poa_file;
    std::string out_stk_file;
    std::string out_enewick_file;

    bool force_unaligned = false;
    bool auto_align = false;
    bool use_ohmic = false;
    bool use_hyphaeon = false;
    bool full_model = true;
    bool all_sites = false;
    bool quiet = false;
    bool show_progress = ProgressBar::is_terminal();
    bool force_progress = false;

    int num_threads = 0;
    Scalar pvalue_threshold = 0.10;
    JSONFormat json_format = JSONFormat::ModernV3;

    // Populated by parse_common_cli_arg if an invalid --json-format was supplied
    bool has_error = false;
    std::string error_message;
};

/**
 * @brief Result status from `parse_common_cli_arg`.
 * Implicitly converts to `true` when `Handled` or `Error` so both
 * `if (parse_common_cli_arg(...))` and `status == CommonArgParseStatus::...` idioms work.
 */
struct CommonArgParseStatus {
    enum Value {
        Unhandled = 0,
        Handled = 1,
        Error = 2
    } value;

    constexpr CommonArgParseStatus(Value v) : value(v) {}
    constexpr operator bool() const { return value != Unhandled; }
    constexpr bool operator==(Value v) const { return value == v; }
    constexpr bool operator!=(Value v) const { return value != v; }
};

/**
 * @brief Try to consume a shared CLI flag at index `i`.
 *
 * Returns `CommonArgParseStatus::Handled` if `argv[i]` was recognized and valid,
 * `CommonArgParseStatus::Error` if recognized but invalid (printing an error to `std::cerr`
 * and populating `opts.has_error`), or `CommonArgParseStatus::Unhandled` otherwise.
 */
inline CommonArgParseStatus parse_common_cli_arg(int& i, int argc, char* argv[], CommonCLIOptions& opts) {
    std::string arg = argv[i];
    if ((arg == "--alignment" || arg == "-a" || arg == "--msa") && i + 1 < argc) {
        opts.alignment_file = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--unaligned" || arg == "-u" || arg == "--sequences") && i + 1 < argc) {
        opts.alignment_file = argv[++i];
        opts.force_unaligned = true;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--auto-align") {
        opts.auto_align = true;
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--out-poa" || arg == "--out-gfa") && i + 1 < argc) {
        opts.out_poa_file = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--out-stk" || arg == "--out-stockholm") && i + 1 < argc) {
        opts.out_stk_file = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--out-enewick" || arg == "--out-enwk") && i + 1 < argc) {
        opts.out_enewick_file = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--tree" || arg == "-t") && i + 1 < argc) {
        opts.tree_file = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--ohmic" || arg == "--omhic") {
        opts.use_ohmic = true;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--hyphaeon") {
        opts.use_ohmic = true;
        opts.use_hyphaeon = true;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--no-hyphaeon") {
        opts.use_hyphaeon = false;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--embeddings" && i + 1 < argc) {
        opts.use_ohmic = true;
        opts.embeddings_file = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--code" || arg == "-c") && i + 1 < argc) {
        opts.code_name = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--branches" || arg == "--test-branches") && i + 1 < argc) {
        opts.branches = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--branch-regex" || arg == "--test-branch-regex") && i + 1 < argc) {
        opts.branch_regex = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--threads" || arg == "-T") && i + 1 < argc) {
        opts.num_threads = std::stoi(argv[++i]);
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
        opts.output_file = argv[++i];
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--html") {
        if (i + 1 < argc && argv[i + 1][0] != '-') {
            opts.html_file = argv[++i];
        } else {
            opts.html_file = "__auto__";
        }
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--pvalue" || arg == "-p") && i + 1 < argc) {
        opts.pvalue_threshold = std::stod(argv[++i]);
        return CommonArgParseStatus::Handled;
    }
    if ((arg == "--json-format" || arg == "--format") && i + 1 < argc) {
        std::string fmt = argv[++i];
        if (!parse_json_format(fmt, opts.json_format)) {
            opts.has_error = true;
            opts.error_message = "Error: unknown --json-format '" + fmt + "'. Expected 'modern' or 'legacy'.";
            std::cerr << Console::danger(opts.error_message) << "\n";
            return CommonArgParseStatus::Error;
        }
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--legacy-json" || arg == "--legacy") {
        opts.json_format = JSONFormat::Legacy;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--modern-json" || arg == "--modern") {
        opts.json_format = JSONFormat::ModernV3;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--full-model") {
        opts.full_model = true;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--no-full-model" || arg == "--quick") {
        opts.full_model = false;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--all-sites") {
        opts.all_sites = true;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--progress") {
        opts.show_progress = true;
        opts.force_progress = true;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--no-progress") {
        opts.show_progress = false;
        opts.force_progress = false;
        return CommonArgParseStatus::Handled;
    }
    if (arg == "--quiet") {
        opts.quiet = true;
        opts.show_progress = false;
        return CommonArgParseStatus::Handled;
    }
    return CommonArgParseStatus::Unhandled;
}

/**
 * @brief Configure OpenMP threads and default output filenames.
 */
inline void finalize_common_cli_options(CommonCLIOptions& opts, const std::string& default_suffix) {
#ifdef _OPENMP
    if (opts.num_threads > 0) {
        omp_set_num_threads(opts.num_threads);
        LikelihoodEngine::set_user_requested_threads(opts.num_threads);
    } else if (const char* env_t = std::getenv("OMP_NUM_THREADS")) {
        int env_threads = std::atoi(env_t);
        if (env_threads > 0) {
            LikelihoodEngine::set_user_requested_threads(env_threads);
        }
    } else {
        // Default to at most 8 threads unless the user explicitly passed --threads N or OMP_NUM_THREADS
        int hw = omp_get_max_threads();
        omp_set_num_threads(std::min(hw, 8));
        LikelihoodEngine::set_user_requested_threads(0);
    }
#else
    if (opts.num_threads > 0) {
        LikelihoodEngine::set_user_requested_threads(opts.num_threads);
    }
#endif
    if (opts.output_file.empty() && !opts.alignment_file.empty()) {
        opts.output_file = opts.alignment_file + default_suffix;
    }
    if (opts.html_file == "__auto__" && !opts.output_file.empty()) {
        opts.html_file = opts.output_file + ".html";
    }
}

/**
 * @brief Bundle of loaded codon alignment, genetic code, resolved phylogeny, and Ohmic network.
 */
struct LoadedCodonDataset {
    std::shared_ptr<const GeneticCode> gcode;
    Alignment aln;
    Tree tree;
    OhmicNetwork ohmic_net;
    std::string tree_source;
};

/**
 * @brief Load GeneticCode, Codon Alignment (or co-infer MSA + Ohmic Tree from unaligned CDS when requested),
 *        and resolve Tree/OhmicNetwork with uniform error reporting.
 */
inline bool load_codon_dataset(
    const CommonCLIOptions& opts,
    LoadedCodonDataset& out,
    bool capture_ohmic_net = true,
    bool allow_fallback_code = false
) {
    try {
        out.gcode = GeneticCode::from_name(opts.code_name);
    } catch (const std::exception& e) {
        if (allow_fallback_code) {
            std::cerr << Console::warning("Warning: Unknown genetic code '" + opts.code_name + "', using Universal.\n");
            out.gcode = GeneticCode::universal();
        } else {
            std::cerr << Console::danger("Error resolving genetic code '") << opts.code_name << "': " << e.what() << "\n";
            return false;
        }
    }

    bool built_from_unaligned = false;
    if (!opts.force_unaligned) {
        try {
            out.aln = Alignment::load(opts.alignment_file, out.gcode);
        } catch (const std::exception& e) {
            std::string err = e.what();
            if (opts.auto_align &&
                (err.find("unequal lengths") != std::string::npos ||
                 err.find("length mismatch") != std::string::npos ||
                 err.find("multiple of 3") != std::string::npos)) {
                built_from_unaligned = true;
            } else {
                std::cerr << Console::danger("Error loading alignment: ") << err << "\n";
                return false;
            }
        }
    } else {
        built_from_unaligned = true;
    }

    if (built_from_unaligned) {
        try {
            auto unaligned_seqs = analyses::JointAlignEngine::load_unaligned_fasta(opts.alignment_file, out.gcode);
            if (unaligned_seqs.size() < 2) {
                throw std::runtime_error("Need at least 2 homologous sequences for de novo joint alignment.");
            }
            size_t min_c = 999999, max_c = 0;
            for (const auto& s : unaligned_seqs) {
                size_t cl = s.raw_seq.size() / 3;
                min_c = std::min(min_c, cl);
                max_c = std::max(max_c, cl);
            }
            if (!opts.quiet) {
                std::cout << Console::info("ℹ") << " "
                          << Console::bold("Unaligned homologous CDS detected (")
                          << unaligned_seqs.size() << " taxa, " << min_c << ".." << max_c
                          << " codons) — running De Novo Ohmic Co-Inference (MSA + Phylogeny + POA)...\n";
            }
            analyses::JointAlignSettings js;
            js.quiet = opts.quiet;
            js.run_marginalized_fel = false; // Only co-estimate the MAP alignment + coupled-SPR Ohmic tree + POA!
            analyses::JointAlignEngine jengine(js);
            auto [dn_aln, dn_tree] = analyses::JointAlignEngine::build_de_novo_codon_alignment_and_tree(
                unaligned_seqs, out.gcode, js
            );
            auto jres = jengine.run(dn_aln, dn_tree);

            if (!opts.out_poa_file.empty()) {
                std::ofstream pf(opts.out_poa_file);
                if (pf.is_open()) pf << jres.poa_gfa;
            }
            if (!opts.out_stk_file.empty()) {
                std::ofstream sf(opts.out_stk_file);
                if (sf.is_open()) sf << jres.stockholm_alignment;
            }
            if (!opts.out_enewick_file.empty()) {
                std::ofstream ef(opts.out_enewick_file);
                if (ef.is_open()) ef << jres.extended_newick << "\n";
            }

            out.aln = std::move(jres.map_alignment);
            if (opts.tree_file.empty()) {
                out.tree = std::move(jres.optimized_tree);
                out.ohmic_net = std::move(jres.ohmic_network);
                out.tree_source = "De Novo Ohmic Circuit + MG94 Codon Transducer + Coupled SPR (" +
                                  std::to_string(out.aln.num_codons) + " codons, " +
                                  std::to_string(out.ohmic_net.bridges.size()) + " Wheatstone bridges)";
                return true;
            }
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error building de novo alignment & phylogeny from unaligned sequences: ")
                      << e.what() << "\n";
            return false;
        }
    }

    try {
        out.tree = analyses::TreeEngine::resolve_tree(
            out.aln, opts.tree_file, opts.use_ohmic, opts.use_hyphaeon,
            opts.embeddings_file, capture_ohmic_net ? &out.ohmic_net : nullptr, &out.tree_source
        );
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error resolving phylogeny: ") << e.what() << "\n";
        return false;
    }

    return true;
}

/**
 * @brief Build a populated `Provenance` block from `CommonCLIOptions` and timing metadata.
 */
inline Provenance build_cli_provenance(
    const CommonCLIOptions& opts,
    const std::string& cli_cmd,
    const std::string& start_iso,
    double wall_time_seconds,
    const std::map<std::string, std::string>& extra_args = {}
) {
    Provenance prov;
    prov.invocation.cli_command = cli_cmd;
    prov.invocation.working_directory = Provenance::get_cwd();
    prov.invocation.arguments["alignment"] = opts.alignment_file;
    if (!opts.tree_file.empty()) prov.invocation.arguments["tree"] = opts.tree_file;
    prov.invocation.arguments["code"] = opts.code_name;
    if (!opts.branches.empty()) prov.invocation.arguments["branches"] = opts.branches;
    if (!opts.branch_regex.empty()) prov.invocation.arguments["branch_regex"] = opts.branch_regex;
    prov.invocation.arguments["pvalue_threshold"] = std::to_string(opts.pvalue_threshold);
    prov.invocation.arguments["full_model"] = opts.full_model ? "true" : "false";
    prov.invocation.arguments["json_format"] = (opts.json_format == JSONFormat::Legacy) ? "legacy" : "modern_v3";
    for (const auto& [k, v] : extra_args) {
        prov.invocation.arguments[k] = v;
    }

    if (!opts.alignment_file.empty()) {
        size_t aln_sz = 0;
        std::string aln_hash = crypto::SHA256::hash_file(opts.alignment_file, &aln_sz);
        if (!aln_hash.empty()) {
            prov.inputs["alignment"] = {opts.alignment_file, "FASTA/NEXUS", aln_hash, aln_sz};
        }
    }

    if (!opts.tree_file.empty() &&
        !analyses::TreeEngine::is_auto_tree_spec(opts.tree_file) &&
        !analyses::TreeEngine::is_ohmic_tree_spec(opts.tree_file)) {
        size_t tree_sz = 0;
        std::string tree_hash = crypto::SHA256::hash_file(opts.tree_file, &tree_sz);
        if (!tree_hash.empty()) {
            prov.inputs["tree"] = {opts.tree_file, "Newick", tree_hash, tree_sz};
        }
    }

    prov.execution.start_time = start_iso;
    prov.execution.end_time = Provenance::current_iso8601();
    prov.execution.wall_time_seconds = wall_time_seconds;
#ifdef _OPENMP
    prov.execution.cpu_threads = (opts.num_threads > 0) ? opts.num_threads : omp_get_max_threads();
#else
    prov.execution.cpu_threads = 1;
#endif
    prov.execution.hostname = Provenance::get_hostname();
    prov.execution.os = Provenance::detect_os();
    prov.execution.compiler = Provenance::detect_compiler();
    return prov;
}

/**
 * @brief Write serialized JSON string (and optional self-contained HTML report) with uniform console status output.
 */
inline bool write_cli_json_and_html(
    const CommonCLIOptions& opts,
    const std::string& json_str
) {
    try {
        std::ofstream out(opts.output_file);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot open output file for writing: " + opts.output_file);
        }
        out << json_str;
        if (json_str.empty() || json_str.back() != '\n') {
            out << "\n";
        }
        if (!opts.quiet) {
            std::string fmt_desc = (opts.json_format == JSONFormat::Legacy)
                ? "Legacy Datamonkey JSON" : "Modern HyPhy v3.0 JSON";
            std::cout << Console::success("✔") << " " << Console::bold("Saved " + fmt_desc + " report to: ")
                      << Console::brand(opts.output_file) << "\n\n";
        }
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error writing JSON: ") << e.what() << "\n";
        return false;
    }

#ifdef HYPHY3_HAS_HTML_REPORT
    if (!opts.html_file.empty()) {
        try {
            hyphy::report::write_html_report(opts.output_file, opts.html_file);
            if (!opts.quiet) {
                std::cout << Console::success("✔") << " " << Console::bold("Saved HTML report to: ")
                          << Console::brand(opts.html_file) << "\n\n";
            }
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error writing HTML report: ") << e.what() << "\n";
            return false;
        }
    }
#endif
    return true;
}

/**
 * @brief Write `nlohmann::json` report (and optional self-contained HTML report) with uniform console status output.
 */
inline bool write_cli_json_and_html(
    const CommonCLIOptions& opts,
    const nlohmann::json& j
) {
    return write_cli_json_and_html(opts, j.dump(2));
}

/**
 * @brief Format percentage `100.0 * n / tot` with 1 decimal place.
 */
inline std::string format_percentage(size_t n, size_t tot) {
    if (tot == 0) return "0.0%";
    double p = (100.0 * static_cast<double>(n)) / static_cast<double>(tot);
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(1) << p << "%";
    return ss.str();
}

} // namespace hyphy::core
