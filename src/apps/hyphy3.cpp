#include <iostream>
#include <string>
#include <vector>
#include <cstring>

// Subcommand entry points
int run_fel(int argc, char* argv[]);
int run_meme(int argc, char* argv[]);
int run_busted(int argc, char* argv[]);
int run_absrel(int argc, char* argv[]);

void print_hyphy3_banner() {
    std::cout << "\n=======================================================\n"
              << "       HYPHY 3: Next-Generation Phylogenetics          \n"
              << "=======================================================\n"
              << " Authors:  Sergei L. Kosakovsky Pond & HyPhy Team       \n"
              << " Version:  3.0.0 (Modern C++20 Core)                   \n"
              << " Web:      https://hyphy.org | https://datamonkey.org  \n"
              << "=======================================================\n\n";
}

void print_hyphy3_usage(const char* prog) {
    print_hyphy3_banner();
    std::cout << "Usage: " << prog << " <analysis> [OPTIONS]\n\n"
              << "Available Analyses:\n"
              << "  fel         Fixed Effects Likelihood (site-by-site selection)\n"
              << "  meme        Mixed Effects Model of Evolution (episodic selection)\n"
              << "  busted      Branch-site Unrestricted Statistical Test (gene-wide selection)\n"
              << "  absrel      Adaptive Branch-Site Random Effects Likelihood (lineage selection)\n\n"
              << "General Options:\n"
              << "  --help, -h       Show this help message\n"
              << "  --version, -v    Show version information\n\n"
              << "Examples:\n"
              << "  " << prog << " fel --alignment data/cd2.fna --tree data/cd2.nwk\n"
              << "  " << prog << " meme --alignment tests/data/adh.nex --threads 8\n"
              << "  " << prog << " busted --alignment tests/data/adh.nex --auto-k --threads 8\n\n";
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
    } else {
        std::cerr << "Error: Unknown analysis '" << cmd << "'.\n\n";
        print_hyphy3_usage(argv[0]);
        return 1;
    }
}
