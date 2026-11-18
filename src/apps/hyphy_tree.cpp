#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/analyses/tree_engine.hpp"
#include "hyphy/core/cli_common.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;

static void print_banner() {
    Panel::print_banner(
        "TREE ENGINE",
        "Fast ML/Ohmic Phylogenetic Tree Search & Newick Tree Utilities",
        "Spectral Scaffolding • Native Polytomy Detection • Branch Labeling & Pruning (label-tree / trim-tree)",
        "Authors: Sergei L. Kosakovsky Pond & HyPhy Team • v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << "\n"
              << "  1. ML/Ohmic Tree Inference:  " << prog << " --alignment <file> [OPTIONS]\n"
              << "  2. Tree Pruning & Labeling:  " << prog << " --tree <in.nwk> [--keep-tips <file|list>] [--label <tag> --list <file|list>] -o <out.nwk>\n\n"
              << Console::bold("Tree Inference Arguments:") << "\n"
              << "  " << Console::brand("--alignment, -a") << " <file>       Path to sequence alignment (FASTA or NEXUS)\n"
              << "  " << Console::brand("--model, -m") << " <gtr|mg94>        Evolutionary model (default: gtr; use mg94 for codon branch length & omega optimization)\n"
              << "  " << Console::brand("--dnc") << " <yes|no>               Enable Ohmic Divide-and-Conquer tree search (default: yes)\n"
              << "  " << Console::brand("--dnc-min-taxa") << " <N>           Minimum taxa to trigger D&C partitioning (default: 16)\n"
              << "  " << Console::brand("--dnc-island-size") << " <N>        Target maximum taxa per Ohmic island (default: 14)\n"
              << "  " << Console::brand("--polytomies") << " <yes|no>         Allow and preserve unresolvable polytomies (default: yes)\n\n"
              << Console::bold("Tree Utility Arguments (label-tree / trim-tree / trim-label-tree):") << "\n"
              << "  " << Console::brand("--tree") << " <in.nwk>              Input Newick tree file to prune and/or label\n"
              << "  " << Console::brand("--keep-tips") << " <file|list>      Keep only specified tips (comma list, text file, or FASTA/NEXUS alignment)\n"
              << "  " << Console::brand("--trim-tips") << " <file|list>      Remove specified tips (comma list, text file, or FASTA/NEXUS alignment)\n"
              << "  " << Console::brand("--labels") << " <csv|tsv>           2-column taxon,group label map file (multi-group labeling)\n"
              << "  " << Console::brand("--label") << " <tag>                Branch label/tag to assign in Newick {tag} format (repeatable with --list; default: Foreground)\n"
              << "  " << Console::brand("--list, --branches") << " <spec>    Branches/tips to label (comma list, text file, FASTA, regex, All, Internal, Leaves)\n"
              << "  " << Console::brand("--internal-nodes") << " <strategy>  Internal branch strategy: All (clade), Parsimony, MRCA, Paraphyletic, Only, or None (default: All)\n"
              << "  " << Console::brand("--clear-tags") << "                 Strip existing {tag} annotations before labeling\n\n"
              << Console::bold("General Output & Runtime:") << "\n"
              << "  " << Console::brand("--output, -o") << " <file>          Path to output Newick tree file (use '-' for stdout)\n"
              << "  " << Console::brand("--json") << " <file>               Path to output JSON summary file\n"
              << "  " << Console::brand("--threads, -t") << " <N>            Number of OpenMP worker threads (0 = auto)\n"
              << "  " << Console::brand("--quiet, -q") << "                  Suppress console progress panels\n"
              << "  " << Console::brand("--help, -h") << "                   Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --output cd2_inferred.nwk\n"
              << "  " << prog << " --tree species.nwk --labels species_groups.csv --internal-nodes Parsimony -o species_labeled.nwk\n"
              << "  " << prog << " --tree species.nwk --keep-tips subset_taxa.txt --label Foreground --list fg_tips.txt -o subset_labeled.nwk\n\n";
}

namespace {

std::string to_lower_str(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

Tree prune_tree_tips(const Tree& in_tree, const std::unordered_set<std::string>& keep_lower) {
    if (in_tree.root_id == INVALID_INDEX || in_tree.nodes.empty()) return in_tree;

    auto serialize_pruned = [&](auto& self, int32_t nid) -> std::pair<std::string, Scalar> {
        const auto& nd = in_tree.nodes[nid];
        if (nd.is_leaf) {
            if (keep_lower.count(to_lower_str(nd.name)) > 0) {
                std::string token = nd.name;
                if (!nd.model_tag.empty()) token += "{" + nd.model_tag + "}";
                return {token, nd.branch_length};
            }
            return {"", 0.0};
        }
        std::vector<std::pair<std::string, Scalar>> kept_children;
        for (int32_t cid : nd.children) {
            auto sub = self(self, cid);
            if (!sub.first.empty()) {
                kept_children.push_back(std::move(sub));
            }
        }
        if (kept_children.empty()) {
            return {"", 0.0};
        }
        if (kept_children.size() == 1) {
            // Collapse degree-1 internal node and accumulate branch length
            Scalar merged_bl = kept_children[0].second + ((nid == in_tree.root_id) ? 0.0 : nd.branch_length);
            return {kept_children[0].first, merged_bl};
        }
        std::ostringstream oss;
        oss << "(";
        for (size_t i = 0; i < kept_children.size(); ++i) {
            if (i > 0) oss << ",";
            oss << kept_children[i].first << ":" << std::setprecision(10) << kept_children[i].second;
        }
        oss << ")" << nd.name;
        if (!nd.model_tag.empty()) oss << "{" + nd.model_tag + "}";
        return {oss.str(), (nid == in_tree.root_id) ? 0.0 : nd.branch_length};
    };

    auto root_pair = serialize_pruned(serialize_pruned, in_tree.root_id);
    if (root_pair.first.empty() || root_pair.first.find('(') == std::string::npos) {
        throw std::runtime_error("Pruning left fewer than 2 tips in the tree.");
    }
    return Tree::from_newick(root_pair.first + ";", true);
}

struct PendingLabelRule {
    std::string tag;
    std::string list_spec;
    bool has_tag = false;
    bool has_list = false;
};

} // namespace

int run_tree(int argc, char* argv[]) {
    std::string aln_path;
    std::string in_tree_path;
    std::string keep_tips_spec;
    std::string trim_tips_spec;
    std::string label_map_path;
    std::vector<PendingLabelRule> label_rules;
    std::string internal_nodes_mode = "All";
    bool clear_tags = false;
    bool label_requested = false;
    std::string output_nwk_path;
    std::string output_json_path;
    std::string model_str = "gtr";
    bool allow_polytomies = true;
    bool enable_dnc = true;
    size_t dnc_min_taxa = 16;
    size_t dnc_island_size = 14;
    bool quiet = false;
    int num_threads = 0;
    size_t spr_max_rounds = 3;
    size_t spr_max_candidates = 8;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if ((arg == "--alignment" || arg == "-a") && i + 1 < argc) {
            aln_path = argv[++i];
        } else if (arg == "--tree" && i + 1 < argc) {
            in_tree_path = argv[++i];
        } else if ((arg == "--keep-tips" || arg == "--keep") && i + 1 < argc) {
            keep_tips_spec = argv[++i];
        } else if ((arg == "--trim-tips" || arg == "--remove-tips" || arg == "--trim") && i + 1 < argc) {
            trim_tips_spec = argv[++i];
        } else if ((arg == "--labels" || arg == "--label-map" || arg == "--branch-labels" || arg == "--groups") && i + 1 < argc) {
            label_map_path = argv[++i];
            label_requested = true;
        } else if ((arg == "--label" || arg == "--label-branches" || arg == "--label-tree") && i + 1 < argc && argv[i + 1][0] != '-') {
            std::string tag_val = argv[++i];
            if (!label_rules.empty() && !label_rules.back().has_tag && label_rules.back().has_list) {
                label_rules.back().tag = tag_val;
                label_rules.back().has_tag = true;
            } else {
                label_rules.push_back({tag_val, "", true, false});
            }
            label_requested = true;
        } else if (arg == "--label-tree") {
            if (label_rules.empty()) {
                label_rules.push_back({"Foreground", "", true, false});
            }
            label_requested = true;
        } else if ((arg == "--list" || arg == "--branches" || arg == "-b") && i + 1 < argc) {
            std::string list_val = argv[++i];
            if (!label_rules.empty() && !label_rules.back().has_list) {
                label_rules.back().list_spec = list_val;
                label_rules.back().has_list = true;
            } else {
                label_rules.push_back({"Foreground", list_val, false, true});
            }
            label_requested = true;
        } else if ((arg == "--internal-nodes" || arg == "--label-mode") && i + 1 < argc) {
            internal_nodes_mode = argv[++i];
        } else if (arg == "--internal-only") {
            internal_nodes_mode = "Only";
        } else if (arg == "--leaf-only" || arg == "--tips-only") {
            internal_nodes_mode = "None";
        } else if (arg == "--clear-tags") {
            clear_tags = true;
        } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
            output_nwk_path = argv[++i];
        } else if (arg == "--json" && i + 1 < argc) {
            output_json_path = argv[++i];
        } else if ((arg == "--model" || arg == "-m") && i + 1 < argc) {
            model_str = argv[++i];
        } else if (arg == "--dnc" && i + 1 < argc) {
            std::string val = argv[++i];
            if (val == "yes" || val == "true" || val == "1") {
                enable_dnc = true;
            } else if (val == "no" || val == "false" || val == "0") {
                enable_dnc = false;
            } else {
                std::cerr << Console::danger("Error: unknown --dnc '") << val
                          << "'. Expected 'yes' or 'no'.\n";
                return 1;
            }
        } else if (arg == "--dnc-min-taxa" && i + 1 < argc) {
            dnc_min_taxa = parse_nonneg_ull("--dnc-min-taxa", argv[++i]);
        } else if (arg == "--dnc-island-size" && i + 1 < argc) {
            dnc_island_size = parse_nonneg_ull("--dnc-island-size", argv[++i]);
        } else if (arg == "--spr-rounds" && i + 1 < argc) {
            spr_max_rounds = parse_nonneg_ull("--spr-rounds", argv[++i]);
        } else if (arg == "--spr-candidates" && i + 1 < argc) {
            spr_max_candidates = parse_nonneg_ull("--spr-candidates", argv[++i]);
        } else if (arg == "--polytomies" && i + 1 < argc) {
            std::string val = argv[++i];
            if (val == "yes" || val == "true" || val == "1") {
                allow_polytomies = true;
            } else if (val == "no" || val == "false" || val == "0") {
                allow_polytomies = false;
            } else {
                std::cerr << Console::danger("Error: unknown --polytomies '") << val
                          << "'. Expected 'yes' or 'no'.\n";
                return 1;
            }
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            // 0 means "auto"; negatives and partial parses are rejected.
            num_threads = static_cast<int>(
                parse_nonneg_ull("--threads", argv[++i], std::numeric_limits<int>::max()));
        } else if (arg == "--quiet" || arg == "-q") {
            quiet = true;
        } else {
            std::cerr << "Error: Unknown argument '" << arg << "'. Use --help for usage.\n";
            return 1;
        }
    }

    // Native Tree Utility Mode (label-tree / trim-tree / trim-label-tree) when --tree <in.nwk> is provided
    if (!in_tree_path.empty()) {
        try {
            Tree tree = Tree::from_newick_file(in_tree_path, true);
            if (clear_tags) {
                for (auto& nd : tree.nodes) nd.model_tag.clear();
            }

            // If --alignment is also provided without explicit --keep-tips, allow pruning tree to alignment taxa
            if (keep_tips_spec.empty() && trim_tips_spec.empty() && !aln_path.empty()) {
                keep_tips_spec = aln_path;
            }

            if (!keep_tips_spec.empty() || !trim_tips_spec.empty()) {
                std::unordered_set<std::string> keep_lower;
                for (const auto& nd : tree.nodes) {
                    if (nd.is_leaf) keep_lower.insert(to_lower_str(nd.name));
                }
                if (!keep_tips_spec.empty()) {
                    keep_lower.clear();
                    for (const auto& item : Tree::parse_taxon_list_or_file(keep_tips_spec)) {
                        keep_lower.insert(to_lower_str(item));
                    }
                }
                if (!trim_tips_spec.empty()) {
                    for (const auto& item : Tree::parse_taxon_list_or_file(trim_tips_spec)) {
                        keep_lower.erase(to_lower_str(item));
                    }
                }
                tree = prune_tree_tips(tree, keep_lower);
            }

            if (!label_map_path.empty()) {
                auto groups = Tree::parse_label_map_file(label_map_path);
                tree.apply_label_map(groups, internal_nodes_mode, label_map_path);
            }
            for (const auto& rule : label_rules) {
                std::string tag = rule.tag.empty() ? "Foreground" : rule.tag;
                tree.apply_branch_labels(tag, rule.list_spec, internal_nodes_mode);
            }

            std::string nwk_out = tree.to_newick(true);
            if (output_nwk_path == "-") {
                std::cout << nwk_out << "\n";
            } else {
                if (output_nwk_path.empty()) {
                    output_nwk_path = in_tree_path + (label_requested ? ".labeled.nwk" : ".trimmed.nwk");
                }
                try {
                    write_file_or_throw(output_nwk_path, nwk_out + "\n");
                } catch (const std::exception& e) {
                    std::cerr << Console::danger("Error: ") << e.what() << "\n";
                    return 1;
                }
                if (!quiet) {
                    size_t n_tips = 0, n_tagged = 0;
                    for (const auto& nd : tree.nodes) {
                        if (nd.is_leaf) n_tips++;
                        if (nd.id != tree.root_id && !nd.model_tag.empty()) n_tagged++;
                    }
                    std::cout << Console::success("Saved modified Newick tree (")
                              << n_tips << " tips, " << n_tagged << " labeled branches) to: "
                              << output_nwk_path << "\n";
                }
            }
            return 0;
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error: ") << e.what() << "\n";
            return 1;
        }
    }

    if (aln_path.empty()) {
        std::cerr << "Error: Missing required argument --alignment <file> (or --tree <in.nwk> for tree utilities).\n\n";
        print_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    } else if (!std::getenv("OMP_NUM_THREADS")) {
        omp_set_num_threads(std::min(omp_get_max_threads(), 8));
    }
#endif

    try {
        Alignment aln;
        if (model_str == "mg94" || model_str == "MG94") {
            try {
                aln = Alignment::load(aln_path);
            } catch (...) {
                aln = Alignment::load_nucleotide(aln_path);
            }
        } else {
            try {
                aln = Alignment::load_nucleotide(aln_path);
            } catch (...) {
                aln = Alignment::load(aln_path);
            }
        }

        TreeEngineSettings settings;
        settings.quiet = quiet;
        settings.allow_polytomies = allow_polytomies;
        settings.enable_dnc = enable_dnc;
        settings.dnc_min_taxa = dnc_min_taxa;
        settings.dnc_target_island_size = dnc_island_size;
        settings.spr_max_rounds = spr_max_rounds;
        settings.spr_max_candidates = spr_max_candidates;

        if (model_str == "mg94" || model_str == "MG94") {
            settings.model = TreeSearchModel::MG94;
            settings.perform_model_escalation = true;
        } else {
            settings.model = TreeSearchModel::GTR;
            settings.perform_model_escalation = false;
        }

        TreeEngine engine(settings);
        TreeEngineResult res = engine.run(aln);

        if (clear_tags) {
            for (auto& nd : res.optimized_tree.nodes) nd.model_tag.clear();
        }
        if (!label_map_path.empty()) {
            auto groups = Tree::parse_label_map_file(label_map_path);
            res.optimized_tree.apply_label_map(groups, internal_nodes_mode, label_map_path);
        }
        for (const auto& rule : label_rules) {
            std::string tag = rule.tag.empty() ? "Foreground" : rule.tag;
            res.optimized_tree.apply_branch_labels(tag, rule.list_spec, internal_nodes_mode);
        }

        // Write Newick output
        if (output_nwk_path == "-") {
            std::cout << res.optimized_tree.to_newick(true) << "\n";
        } else {
            if (output_nwk_path.empty()) {
                output_nwk_path = aln_path + ".nwk";
            }
            try {
                write_file_or_throw(output_nwk_path, res.optimized_tree.to_newick(true) + "\n");
            } catch (const std::exception& e) {
                std::cerr << Console::danger("Error: ") << e.what() << "\n";
                return 1;
            }
            if (!quiet) {
                std::cout << Console::success("Saved inferred Newick tree to: ") << output_nwk_path << "\n";
            }
        }

        // Write JSON output if requested
        if (!output_json_path.empty()) {
            Provenance prov;
            prov.invocation.cli_command = "hyphy3 tree";
            for (int i = 1; i < argc; ++i) {
                prov.invocation.cli_command += " " + std::string(argv[i]);
            }
            prov.invocation.working_directory = Provenance::get_cwd();
            if (!aln_path.empty()) {
                size_t sz = 0;
                std::string hash = crypto::SHA256::hash_file(aln_path, &sz);
                if (!hash.empty()) {
                    prov.inputs["alignment"] = {aln_path, "FASTA/NEXUS", hash, sz};
                }
            }
            prov.execution.start_time = Provenance::current_iso8601();
            prov.execution.end_time = prov.execution.start_time;
            prov.execution.wall_time_seconds = res.runtime_seconds;
#ifdef _OPENMP
            prov.execution.cpu_threads = (num_threads > 0) ? num_threads : omp_get_max_threads();
#else
            prov.execution.cpu_threads = (num_threads > 0) ? num_threads : 1;
#endif
            prov.execution.hostname = Provenance::get_hostname();
            prov.execution.os = Provenance::detect_os();
            prov.execution.compiler = Provenance::detect_compiler();

            try {
                write_file_or_throw(output_json_path, engine.to_json(res, prov).dump(2) + "\n");
            } catch (const std::exception& e) {
                std::cerr << Console::danger("Error: ") << e.what() << "\n";
                return 1;
            }
            if (!quiet) {
                std::cout << Console::success("Saved analysis JSON to: ") << output_json_path << "\n";
            }
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) {
    try {
        return run_tree(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}
#endif
