// Standalone entry point for hyphy3-expcm

extern int run_expcm(int argc, char* argv[]);

int main(int argc, char* argv[]) {
    try {
        return run_expcm(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << Console::danger("Error: ") << e.what() << "\n";
        return 1;
    }
}
