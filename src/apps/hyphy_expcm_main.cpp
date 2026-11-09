// Standalone entry point for hyphy3-expcm

#include "hyphy/core/console.hpp"

#include <exception>
#include <iostream>

extern int run_expcm(int argc, char* argv[]);

int main(int argc, char* argv[]) {
    try {
        return run_expcm(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << hyphy::core::Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}
