#include "hyphy/core/console.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <cstring>

using namespace hyphy::core;

// Subcommand entry points
int run_fel(int argc, char* argv[]);
int run_contrast_fel(int argc, char* argv[]);
int run_meme(int argc, char* argv[]);
int run_busted(int argc, char* argv[]);
int run_busted_ph(int argc, char* argv[]);
int run_absrel(int argc, char* argv[]);
int run_relax(int argc, char* argv[]);
int run_fubar(int argc, char* argv[]);
int run_prime(int argc, char* argv[]);
int run_gard(int argc, char* argv[]);
int run_slac(int argc, char* argv[]);
int run_dfel(int argc, char* argv[]);
int run_ohmic(int argc, char* argv[]);
int run_tree(int argc, char* argv[]);
int run_joint_align(int argc, char* argv[]);
int run_characterize(int argc, char* argv[]);
int run_mss(int argc, char* argv[]);
int run_expcm(int argc, char* argv[]);
int run_report(int argc, char* argv[]);

void print_hyphy3_banner() {
    Panel::print_banner(
        "HYPHY 3",
        "Modern C++20 Phylogenetics Engine",
        "Next-generation molecular evolution, selection detection & hypothesis testing",
        "Authors: Sergei L. Kosakovsky Pond & HyPhy Team • https://hyphy.org • v3.0.0"
    );
}

void print_hyphy3_usage(const char* prog) {
    print_hyphy3_banner();
    std::cout << Console::bold("Usage:") << " " << prog << " " 
              << Console::accent("<analysis> | -u <unaligned.fa> | -a <alignment.nex>") << " "
              << Console::muted("[OPTIONS]") << "\n\n";

    std::cout << Console::bold("Unified Joint Alignment, Network & Selection (WALRUS):") << "\n"
              << "  " << Console::brand("walrus       ") << "  Wheatstone Alignment and Likelihood Rate Unification for Selection (aliases: characterize, auto, or pass -u/-a directly)\n\n";

    std::cout << Console::bold("Available Analyses:") << "\n"
              << "  " << Console::brand("fel          ") << "  Fixed Effects Likelihood (site-by-site selection)\n"
              << "  " << Console::brand("contrast-fel ") << "  Contrast-FEL (differential selection between branch sets)\n"
              << "  " << Console::brand("meme         ") << "  Mixed Effects Model of Evolution (episodic selection)\n"
              << "  " << Console::brand("busted       ") << "  Branch-site Unrestricted Statistical Test (gene-wide selection)\n"
              << "  " << Console::brand("busted-ph    ") << "  BUSTED-Phenotype (selection associated with a phenotype/trait)\n"
              << "  " << Console::brand("absrel       ") << "  Adaptive Branch-Site Random Effects Likelihood (lineage selection)\n"
              << "  " << Console::brand("relax        ") << "  Test for Selection Relaxation (branch set contrast)\n"
              << "  " << Console::brand("fubar        ") << "  Fast, Unconstrained Bayesian AppRoximation (Bayesian site selection)\n"
              << "  " << Console::brand("prime        ") << "  PRoperty Informed Models of Evolution (biophysical property selection)\n"
              << "  " << Console::brand("gard         ") << "  Genetic Algorithm for Recombination Detection (breakpoint screening)\n"
              << "  " << Console::brand("slac         ") << "  Single-Likelihood Ancestor Counting (counting site selection)\n"
              << "  " << Console::brand("dfel         ") << "  Directional FEL (non-stationary selection & auto-rooting)\n"
              << "  " << Console::brand("mss          ") << "  Multi-class Synonymous Substitution Models (SynREV, TSVs & fast DP selector)\n"
              << "  " << Console::brand("ohmic        ") << "  Ohmic Phylogenetic Networks & Wheatstone Pruning\n"
              << "  " << Console::brand("tree         ") << "  ML/Ohmic Tree Engine (spectral scaffold & native polytomies)\n"
              << "  " << Console::brand("joint-align  ") << "  Ohmic Co-Inference of Codon Alignment, Phylogeny & Selection\n"
              << "  " << Console::brand("expcm        ") << "  Experimentally-informed Codon Model (site-specific dN/dS analysis)\n\n";

    std::cout << Console::bold("Reports:") << "\n"
              << "  " << Console::brand("report       ") << "  Write a self-contained HTML report for an existing result JSON\n\n";

    std::cout << Console::bold("General Options:") << "\n"
              << "  " << Console::brand("--help, -h   ") << "  Show this help message\n"
              << "  " << Console::brand("--version, -v") << "  Show version information\n\n";

    std::cout << Console::bold("Examples:") << "\n"
              << "  " << prog << " -u tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa\n"
              << "  " << prog << " characterize --unaligned tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa --html ifitm3.html\n"
              << "  " << prog << " fel --alignment data/cd2.fna --tree data/cd2.nwk --html cd2.FEL.html\n"
              << "  " << prog << " dfel --alignment tests/data/adh.nex --root auto\n"
              << "  " << prog << " ohmic --alignment tests/data/cd2.fna\n"
              << "  " << prog << " joint-align --unaligned tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa\n"
              << "  " << prog << " meme --alignment tests/data/adh.nex --threads 8\n"
              << "  " << prog << " busted --alignment tests/data/adh.nex --auto-k --threads 8\n"
              << "  " << prog << " busted-ph --alignment data/cd2.fna --tree data/cd2.nwk --branches Human,Chimp,Baboon\n"
              << "  " << prog << " absrel --alignment tests/data/adh.nex --threads 8\n"
              << "  " << prog << " relax --alignment tests/data/Fig4E.nex\n"
              << "  " << prog << " fubar --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " prime --alignment data/cd2.fna --tree data/cd2.nwk --property-set 5PROP\n"
              << "  " << prog << " slac --alignment tests/data/adh.nex --samples 100\n"
              << "  " << prog << " report results.MEME.json --html results.MEME.html\n\n";
}

// Dispatch a resolved subcommand name to its run_* entry point. argc/argv are
// passed through unchanged, so argv[0] must already be the program name the
// analysis should display (e.g. "hyphy3 fel" or "hyphy3-fel").
static int dispatch_subcommand(const std::string& cmd, int argc, char* argv[]) {
    try {
        if (cmd == "walrus" || cmd == "WALRUS" || cmd == "characterize" || cmd == "auto" || cmd == "omni" || cmd == "one-button") {
            return run_characterize(argc, argv);
        } else if (cmd == "fel") {
            return run_fel(argc, argv);
        } else if (cmd == "contrast-fel" || cmd == "contrast_fel" || cmd == "cfel") {
            return run_contrast_fel(argc, argv);
        } else if (cmd == "meme") {
            return run_meme(argc, argv);
        } else if (cmd == "busted") {
            return run_busted(argc, argv);
        } else if (cmd == "busted-ph" || cmd == "busted_ph" || cmd == "BUSTED-PH") {
            return run_busted_ph(argc, argv);
        } else if (cmd == "absrel") {
            return run_absrel(argc, argv);
        } else if (cmd == "relax") {
            return run_relax(argc, argv);
        } else if (cmd == "fubar" || cmd == "FUBAR" || cmd == "fub") {
            return run_fubar(argc, argv);
        } else if (cmd == "prime" || cmd == "PRIME") {
            return run_prime(argc, argv);
        } else if (cmd == "gard" || cmd == "GARD") {
            return run_gard(argc, argv);
        } else if (cmd == "slac" || cmd == "SLAC") {
            return run_slac(argc, argv);
        } else if (cmd == "dfel" || cmd == "directional-fel" || cmd == "directional_fel") {
            return run_dfel(argc, argv);
        } else if (cmd == "mss" || cmd == "MSS" || cmd == "synrev" || cmd == "SynREV") {
            return run_mss(argc, argv);
        } else if (cmd == "ohmic" || cmd == "OHMIC" || cmd == "--ohmic") {
            return run_ohmic(argc, argv);
        } else if (cmd == "tree" || cmd == "TREE" || cmd == "tree-search" ||
                   cmd == "label-tree" || cmd == "trim-tree" || cmd == "trim-label-tree") {
            return run_tree(argc, argv);
        } else if (cmd == "joint-align" || cmd == "joint_align" || cmd == "align") {
            return run_joint_align(argc, argv);
        } else if (cmd == "expcm" || cmd == "ExpCM" || cmd == "EXPCM") {
            return run_expcm(argc, argv);
        } else if (cmd == "report") {
            return run_report(argc, argv);
        } else {
            std::cerr << "Error: Unknown analysis '" << cmd << "'.\n\n";
            print_hyphy3_usage(argv[0]);
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}

int main(int argc, char* argv[]) {
    // Alias dispatch: when invoked through a hyphenated name (hyphy3-fel,
    // hyphy-meme, ...), the basename suffix selects the subcommand so the
    // combined binary can be installed under the standalone names — e.g. as
    // symlinks — without separate executables. argv is passed through
    // unchanged, so the analysis reports "hyphy3-fel" as its program name
    // exactly as the real standalone binaries do.
    const std::string invoked = (argc > 0 && argv[0]) ? argv[0] : "";
    const std::string::size_type sep = invoked.find_last_of("/\\");
    const std::string base = (sep == std::string::npos) ? invoked : invoked.substr(sep + 1);
    for (const char* prefix : {"hyphy3-", "hyphy-"}) {
        const std::string p(prefix);
        if (base.size() > p.size() && base.compare(0, p.size(), p) == 0) {
            return dispatch_subcommand(base.substr(p.size()), argc, argv);
        }
    }

    if (argc < 2) {
        print_hyphy3_usage(argv[0]);
        return 0;
    }

    std::string cmd = argv[1];
    if (cmd == "--help" || cmd == "-h" || cmd == "help") {
        print_hyphy3_usage(argv[0]);
        return 0;
    }
    if (cmd == "--version" || cmd == "-v" || cmd == "version") {
        print_hyphy3_banner();
        std::cout << "HyPhy 3 v" HYPHY3_VERSION " (commit: " HYPHY3_GIT_COMMIT ", build: " HYPHY3_BUILD_TYPE ")\n";
        return 0;
    }

    // Direct One-Button invocation: if the user passes -u/--unaligned, -a/--alignment/--msa, or -i/--input directly to hyphy3,
    // route straight to the unified CharacterizeEngine!
    if (cmd == "-u" || cmd == "--unaligned" || cmd == "-a" || cmd == "--alignment" || cmd == "--msa" || cmd == "-i" || cmd == "--input") {
        try {
            return run_characterize(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << Console::danger("Error: ") << e.what() << "\n";
            return 1;
        }
    }

    // Shift arguments so the subcommand sees argv[0] as "hyphy3 <subcommand>"
    // and remaining arguments as its options
    std::vector<char*> sub_argv;
    std::string sub_prog = std::string(argv[0]) + " " + cmd;
    sub_argv.push_back(const_cast<char*>(sub_prog.c_str()));
    for (int i = 2; i < argc; ++i) {
        sub_argv.push_back(argv[i]);
    }
    int sub_argc = static_cast<int>(sub_argv.size());

    return dispatch_subcommand(cmd, sub_argc, sub_argv.data());
}
