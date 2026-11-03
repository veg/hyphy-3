#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/analyses/tree_engine.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace hyphy::core;
using namespace hyphy::analyses;

static void print_banner() {
    Panel::print_banner(
        "TREE ENGINE",
        "Fast ML/Ohmic Phylogenetic Tree Search",
        "Spectral Scaffolding • Native Polytomy Detection • All-Branch Quasi-Newton",
        "Authors: Sergei L. Kosakovsky Pond & HyPhy Team • v3.0.0"
    );
}

static void print_usage(const char* prog) {
    print_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " [OPTIONS]\n\n"
              << Console::bold("Required arguments:") << "\n"
              << "  " << Console::brand("--alignment, -a") << " <file>  Path to sequence alignment (FASTA or NEXUS)\n\n"
              << Console::bold("Optional arguments:") << "\n"
              << "  " << Console::brand("--model, -m") << " <gtr|mg94>   Evolutionary model (default: gtr; use mg94 for codon branch length & omega optimization)\n"
              << "  " << Console::brand("--dnc") << " <yes|no>          Enable Ohmic Divide-and-Conquer tree search (default: yes)\n"
              << "  " << Console::brand("--dnc-min-taxa") << " <N>      Minimum taxa to trigger D&C partitioning (default: 16)\n"
              << "  " << Console::brand("--dnc-island-size") << " <N>   Target maximum taxa per Ohmic island (default: 14)\n"
              << "  " << Console::brand("--polytomies") << " <yes|no>    Allow and preserve unresolvable polytomies (default: yes)\n"
              << "  " << Console::brand("--output, -o") << " <file>     Path to output Newick tree file (default: <alignment>.nwk)\n"
              << "  " << Console::brand("--json") << " <file>          Path to output JSON summary file\n"
              << "  " << Console::brand("--threads, -t") << " <N>       Number of OpenMP worker threads\n"
              << "  " << Console::brand("--quiet, -q") << "             Suppress console progress panels\n"
              << "  " << Console::brand("--help, -h") << "              Show this help message\n\n"
              << Console::bold("Example:") << "\n"
              << "  " << prog << " --alignment data/cd2.fna --output cd2_inferred.nwk\n\n";
}

int run_tree(int argc, char* argv[]) {
    std::string aln_path;
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

    if (aln_path.empty()) {
        std::cerr << "Error: Missing required argument --alignment <file>.\n\n";
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
        if (nwk_file.is_open()) {
            nwk_file << res.optimized_tree.to_newick(true) << "\n";
            if (!quiet) {
                std::cout << Console::success("Saved inferred Newick tree to: ") << output_nwk_path << "\n";
            }
        } else {
            std::cerr << Console::danger("Warning: Could not open output file for Newick: ") << output_nwk_path << "\n";
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
            if (json_file.is_open()) {
                json_file << engine.to_json(res, prov).dump(2) << "\n";
                if (!quiet) {
                    std::cout << Console::success("Saved analysis JSON to: ") << output_json_path << "\n";
                }
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
