#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/analyses/tree_engine.hpp"

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
              << "  " << Console::brand("--label") << " <tag>                Branch label/tag to assign in Newick {tag} format (default: Foreground)\n"
              << "  " << Console::brand("--list, --branches") << " <spec>    Branches/tips to label (comma list, text file, FASTA, regex, All, Internal, Leaves)\n"
              << "  " << Console::brand("--internal-nodes") << " <strategy>  How to label internal branches from tip sets: All (clade), MRCA, Paraphyletic, Only, or None (default: All)\n"
              << "  " << Console::brand("--clear-tags") << "                 Strip existing {tag} annotations before labeling\n\n"
              << Console::bold("General Output & Runtime:") << "\n"
              << "  " << Console::brand("--output, -o") << " <file>          Path to output Newick tree file (use '-' for stdout)\n"
              << "  " << Console::brand("--json") << " <file>               Path to output JSON summary file\n"
              << "  " << Console::brand("--threads, -t") << " <N>            Number of OpenMP worker threads\n"
              << "  " << Console::brand("--quiet, -q") << "                  Suppress console progress panels\n"
              << "  " << Console::brand("--help, -h") << "                   Show this help message\n\n"
              << Console::bold("Examples:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --output cd2_inferred.nwk\n"
              << "  " << prog << " --tree species.nwk --keep-tips subset_taxa.txt --label Foreground --list fg_tips.txt -o subset_labeled.nwk\n\n";
}

namespace {

std::string trim_ws(std::string s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::string to_lower_str(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::vector<std::string> parse_taxon_list_or_file(const std::string& spec) {
    std::vector<std::string> items;
    if (spec.empty()) return items;
    std::ifstream ifs(spec);
    if (ifs.is_open()) {
        std::string line;
        bool is_fasta = false;
        while (std::getline(ifs, line)) {
            std::string t = trim_ws(line);
            if (t.empty() || t[0] == '#') continue;
            if (t[0] == '>') {
                is_fasta = true;
                std::string hdr = trim_ws(t.substr(1));
                size_t sp = hdr.find_first_of(" \t");
                if (sp != std::string::npos) hdr = hdr.substr(0, sp);
                if (!hdr.empty()) items.push_back(hdr);
            } else if (!is_fasta) {
                std::stringstream ss(t);
                std::string tok;
                while (std::getline(ss, tok, ',')) {
                    tok = trim_ws(tok);
                    if (!tok.empty()) items.push_back(tok);
                }
            }
        }
        return items;
    }
    std::stringstream ss(spec);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = trim_ws(tok);
        if (!tok.empty()) items.push_back(tok);
    }
    return items;
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

void apply_branch_labels(
    Tree& tree,
    const std::string& label_tag,
    const std::string& list_spec,
    const std::string& internal_mode_raw
) {
    std::string imode = to_lower_str(trim_ws(internal_mode_raw));
    if (imode.empty()) imode = "all";

    std::vector<std::string> explicit_items = parse_taxon_list_or_file(list_spec);
    Tree::BranchSelection sel;
    if (explicit_items.size() > 1 || (explicit_items.size() == 1 && list_spec != explicit_items[0])) {
        sel = tree.select_branches("", "", explicit_items, false);
    } else {
        sel = tree.select_branches(list_spec.empty() ? "All" : list_spec, "", {}, false);
    }

    size_t N = tree.nodes.size();
    std::vector<bool> target(N, false);
    for (size_t i = 0; i < N; ++i) target[i] = sel.is_selected[i];

    // Compute leaf descendant counts (total vs selected) in post-order
    tree.compute_post_order();
    std::vector<int> total_leaves(N, 0);
    std::vector<int> selected_leaves(N, 0);
    int global_selected_leaves = 0;
    for (int32_t nid : tree.post_order) {
        const auto& nd = tree.nodes[nid];
        if (nd.is_leaf) {
            total_leaves[nid] = 1;
            selected_leaves[nid] = target[nid] ? 1 : 0;
            global_selected_leaves += selected_leaves[nid];
        } else {
            for (int32_t cid : nd.children) {
                total_leaves[nid] += total_leaves[cid];
                selected_leaves[nid] += selected_leaves[cid];
            }
        }
    }

    if (global_selected_leaves > 0) {
        if (imode == "all" || imode == "clade" || imode == "monophyletic") {
            for (const auto& nd : tree.nodes) {
                if (nd.id == tree.root_id || nd.is_leaf) continue;
                if (total_leaves[nd.id] > 0 && selected_leaves[nd.id] == total_leaves[nd.id]) {
                    target[nd.id] = true;
                }
            }
        } else if (imode == "only" || imode == "internal-only" || imode == "internal") {
            for (const auto& nd : tree.nodes) {
                if (nd.id == tree.root_id) continue;
                if (nd.is_leaf) {
                    target[nd.id] = false;
                } else if (total_leaves[nd.id] > 0 && selected_leaves[nd.id] == total_leaves[nd.id]) {
                    target[nd.id] = true;
                }
            }
        } else if (imode == "mrca") {
            // Find MRCA of all selected leaves and label its entire descendant clade
            int32_t mrca_id = tree.root_id;
            for (int32_t nid : tree.post_order) {
                if (selected_leaves[nid] == global_selected_leaves) {
                    mrca_id = nid;
                    break;
                }
            }
            auto mark_clade = [&](auto& self, int32_t u) -> void {
                if (u != tree.root_id) target[u] = true;
                for (int32_t c : tree.nodes[u].children) self(self, c);
            };
            mark_clade(mark_clade, mrca_id);
        } else if (imode == "paraphyletic") {
            // Label internal nodes that have >= 1 selected leaf in at least 2 children (or 1 child + ancestor with selected leaves)
            for (const auto& nd : tree.nodes) {
                if (nd.id == tree.root_id || nd.is_leaf) continue;
                int children_with_sel = 0;
                for (int32_t cid : nd.children) {
                    if (selected_leaves[cid] > 0) children_with_sel++;
                }
                bool outside_has_sel = (global_selected_leaves - selected_leaves[nd.id]) > 0;
                if ((children_with_sel >= 2) || (children_with_sel >= 1 && outside_has_sel && selected_leaves[nd.id] == total_leaves[nd.id])) {
                    target[nd.id] = true;
                }
            }
        } else if (imode == "none" || imode == "leaves" || imode == "tips") {
            for (const auto& nd : tree.nodes) {
                if (!nd.is_leaf && !sel.is_selected[nd.id]) {
                    target[nd.id] = false;
                }
            }
        }
    }

    for (auto& nd : tree.nodes) {
        if (nd.id == tree.root_id) continue;
        if (target[nd.id]) {
            nd.model_tag = label_tag;
        }
    }
}

} // namespace

int run_tree(int argc, char* argv[]) {
    std::string aln_path;
    std::string in_tree_path;
    std::string keep_tips_spec;
    std::string trim_tips_spec;
    std::string label_tag;
    std::string label_list_spec;
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
        } else if ((arg == "--label" || arg == "--label-branches" || arg == "--label-tree") && i + 1 < argc && argv[i + 1][0] != '-') {
            label_tag = argv[++i];
            label_requested = true;
        } else if (arg == "--label-tree") {
            if (label_tag.empty()) label_tag = "Foreground";
            label_requested = true;
        } else if ((arg == "--list" || arg == "--branches" || arg == "-b") && i + 1 < argc) {
            label_list_spec = argv[++i];
            label_requested = true;
            if (label_tag.empty()) label_tag = "Foreground";
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
            enable_dnc = (val == "yes" || val == "true" || val == "1");
        } else if (arg == "--dnc-min-taxa" && i + 1 < argc) {
            dnc_min_taxa = std::stoul(argv[++i]);
        } else if (arg == "--dnc-island-size" && i + 1 < argc) {
            dnc_island_size = std::stoul(argv[++i]);
        } else if (arg == "--spr-rounds" && i + 1 < argc) {
            spr_max_rounds = std::stoul(argv[++i]);
        } else if (arg == "--spr-candidates" && i + 1 < argc) {
            spr_max_candidates = std::stoul(argv[++i]);
        } else if (arg == "--polytomies" && i + 1 < argc) {
            std::string val = argv[++i];
            allow_polytomies = (val == "yes" || val == "true" || val == "1");
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
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
                    for (const auto& item : parse_taxon_list_or_file(keep_tips_spec)) {
                        keep_lower.insert(to_lower_str(item));
                    }
                }
                if (!trim_tips_spec.empty()) {
                    for (const auto& item : parse_taxon_list_or_file(trim_tips_spec)) {
                        keep_lower.erase(to_lower_str(item));
                    }
                }
                tree = prune_tree_tips(tree, keep_lower);
            }

            if (label_requested) {
                if (label_tag.empty()) label_tag = "Foreground";
                apply_branch_labels(tree, label_tag, label_list_spec, internal_nodes_mode);
            }

            std::string nwk_out = tree.to_newick(true);
            if (output_nwk_path == "-") {
                std::cout << nwk_out << "\n";
            } else {
                if (output_nwk_path.empty()) {
                    output_nwk_path = in_tree_path + (label_requested ? ".labeled.nwk" : ".trimmed.nwk");
                }
                std::ofstream ofs(output_nwk_path);
                if (!ofs.is_open()) {
                    std::cerr << Console::danger("Error: Could not open output file: ") << output_nwk_path << "\n";
                    return 1;
                }
                ofs << nwk_out << "\n";
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

        // Write Newick output
        if (output_nwk_path.empty()) {
            output_nwk_path = aln_path + ".nwk";
        }
        std::ofstream nwk_file(output_nwk_path);
        if (!nwk_file.is_open()) {
            std::cerr << Console::danger("Error: Could not open output file for Newick: ") << output_nwk_path << "\n";
            return 1;
        }
        nwk_file << res.optimized_tree.to_newick(true) << "\n";
        if (!nwk_file) {
            std::cerr << Console::danger("Error: Failed to write Newick output to: ") << output_nwk_path << "\n";
            return 1;
        }
        if (!quiet) {
            std::cout << Console::success("Saved inferred Newick tree to: ") << output_nwk_path << "\n";
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
            prov.execution.cpu_threads = (num_threads > 0) ? num_threads : 1;
            prov.execution.hostname = Provenance::get_hostname();
            prov.execution.os = Provenance::detect_os();
            prov.execution.compiler = Provenance::detect_compiler();

            std::ofstream json_file(output_json_path);
            if (!json_file.is_open()) {
                std::cerr << Console::danger("Error: Could not open output file for JSON: ") << output_json_path << "\n";
                return 1;
            }
            json_file << engine.to_json(res, prov).dump(2) << "\n";
            if (!json_file) {
                std::cerr << Console::danger("Error: Failed to write JSON output to: ") << output_json_path << "\n";
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
    return run_tree(argc, argv);
}
#endif
