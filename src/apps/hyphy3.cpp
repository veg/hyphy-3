#include "hyphy/core/console.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <cstring>

using namespace hyphy::core;

// Subcommand entry points
int run_fel(int argc, char* argv[]);
int run_meme(int argc, char* argv[]);
int run_busted(int argc, char* argv[]);
int run_absrel(int argc, char* argv[]);
int run_relax(int argc, char* argv[]);

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
              << Console::accent("<analysis>") << " "
              << Console::muted("[OPTIONS]") << "\n\n";

    std::cout << Console::bold("Available Analyses:") << "\n"
              << "  " << Console::brand("fel     ") << "  Fixed Effects Likelihood (site-by-site selection)\n"
              << "  " << Console::brand("meme    ") << "  Mixed Effects Model of Evolution (episodic selection)\n"
              << "  " << Console::brand("busted  ") << "  Branch-site Unrestricted Statistical Test (gene-wide selection)\n"
              << "  " << Console::brand("absrel  ") << "  Adaptive Branch-Site Random Effects Likelihood (lineage selection)\n"
              << "  " << Console::brand("relax   ") << "  Test for Selection Relaxation (branch set contrast)\n\n";

    std::cout << Console::bold("General Options:") << "\n"
              << "  " << Console::brand("--help, -h   ") << "  Show this help message\n"
              << "  " << Console::brand("--version, -v") << "  Show version information\n\n";

    std::cout << Console::bold("Examples:") << "\n"
              << "  " << prog << " fel --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " meme --alignment tests/data/adh.nex --threads 8\n"
              << "  " << prog << " busted --alignment tests/data/adh.nex --auto-k --threads 8\n"
              << "  " << prog << " absrel --alignment tests/data/adh.nex --threads 8\n"
              << "  " << prog << " relax --alignment tests/data/Fig4E.nex\n\n";
}

int main(int argc, char* argv[]) {
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
        return 0;
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

    if (cmd == "fel") {
        return run_fel(sub_argc, sub_argv.data());
    } else if (cmd == "meme") {
        return run_meme(sub_argc, sub_argv.data());
    } else if (cmd == "busted") {
        return run_busted(sub_argc, sub_argv.data());
    } else if (cmd == "absrel") {
        return run_absrel(sub_argc, sub_argv.data());
    } else if (cmd == "relax") {
        return run_relax(sub_argc, sub_argv.data());
    } else {
        std::cerr << "Error: Unknown analysis '" << cmd << "'.\n\n";
        print_hyphy3_usage(argv[0]);
        return 1;
    }
}
