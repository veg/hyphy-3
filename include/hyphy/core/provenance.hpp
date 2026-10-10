#pragma once

#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <array>
#include <cstdint>
#include <cstring>
#include <unistd.h>
#include "nlohmann/json.hpp"

namespace hyphy::core {

namespace crypto {

// Clean, standalone FIPS 180-2 SHA-256 implementation
class SHA256 {
public:
    SHA256() { reset(); }

    void reset() {
        m_len = 0;
        m_tot_len = 0;
        m_h[0] = 0x6a09e667;
        m_h[1] = 0xbb67ae85;
        m_h[2] = 0x3c6ef372;
        m_h[3] = 0xa54ff53a;
        m_h[4] = 0x510e527f;
        m_h[5] = 0x9b05688c;
        m_h[6] = 0x1f83d9ab;
        m_h[7] = 0x5be0cd19;
    }

    void update(const uint8_t* data, size_t len) {
        for (size_t i = 0; i < len; ++i) {
            m_block[m_len++] = data[i];
            if (m_len == 64) {
                transform();
                m_tot_len += 64;
                m_len = 0;
            }
        }
    }

    void update(const std::string& str) {
        update(reinterpret_cast<const uint8_t*>(str.data()), str.size());
    }

    std::string finalize() {
        uint64_t total_bits = (m_tot_len + m_len) * 8;
        m_block[m_len++] = 0x80;
        if (m_len > 56) {
            while (m_len < 64) m_block[m_len++] = 0x00;
            transform();
            m_len = 0;
        }
        while (m_len < 56) m_block[m_len++] = 0x00;
        for (int i = 7; i >= 0; --i) {
            m_block[m_len++] = static_cast<uint8_t>((total_bits >> (i * 8)) & 0xff);
        }
        transform();

        std::ostringstream ss;
        ss << std::hex << std::setfill('0');
        for (int i = 0; i < 8; ++i) {
            ss << std::setw(8) << m_h[i];
        }
        return ss.str();
    }

    static std::string hash_string(const std::string& str) {
        SHA256 ctx;
        ctx.update(str);
        return ctx.finalize();
    }

    static std::string hash_file(const std::string& filepath, size_t* out_size = nullptr) {
        std::ifstream file(filepath, std::ios::binary);
        if (!file.is_open()) {
            return "";
        }
        SHA256 ctx;
        std::vector<uint8_t> buffer(65536);
        size_t total_bytes = 0;
        while (file.good()) {
            file.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
            std::streamsize read_count = file.gcount();
            if (read_count > 0) {
                ctx.update(buffer.data(), static_cast<size_t>(read_count));
                total_bytes += static_cast<size_t>(read_count);
            }
        }
        if (out_size) {
            *out_size = total_bytes;
        }
        return ctx.finalize();
    }

private:
    static inline uint32_t rotr(uint32_t x, uint32_t n) { return (x >> n) | (x << (32 - n)); }
    static inline uint32_t ch(uint32_t x, uint32_t y, uint32_t z) { return (x & y) ^ (~x & z); }
    static inline uint32_t maj(uint32_t x, uint32_t y, uint32_t z) { return (x & y) ^ (x & z) ^ (y & z); }
    static inline uint32_t sig0(uint32_t x) { return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22); }
    static inline uint32_t sig1(uint32_t x) { return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25); }
    static inline uint32_t theta0(uint32_t x) { return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3); }
    static inline uint32_t theta1(uint32_t x) { return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10); }

    void transform() {
        static const uint32_t K[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
        };

        uint32_t w[64];
        for (int i = 0; i < 16; ++i) {
            w[i] = (static_cast<uint32_t>(m_block[i * 4]) << 24) |
                   (static_cast<uint32_t>(m_block[i * 4 + 1]) << 16) |
                   (static_cast<uint32_t>(m_block[i * 4 + 2]) << 8) |
                   (static_cast<uint32_t>(m_block[i * 4 + 3]));
        }
        for (int i = 16; i < 64; ++i) {
            w[i] = theta1(w[i - 2]) + w[i - 7] + theta0(w[i - 15]) + w[i - 16];
        }

        uint32_t a = m_h[0], b = m_h[1], c = m_h[2], d = m_h[3];
        uint32_t e = m_h[4], f = m_h[5], g = m_h[6], h = m_h[7];

        for (int i = 0; i < 64; ++i) {
            uint32_t t1 = h + sig1(e) + ch(e, f, g) + K[i] + w[i];
            uint32_t t2 = sig0(a) + maj(a, b, c);
            h = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }

        m_h[0] += a;
        m_h[1] += b;
        m_h[2] += c;
        m_h[3] += d;
        m_h[4] += e;
        m_h[5] += f;
        m_h[6] += g;
        m_h[7] += h;
    }

    uint32_t m_h[8];
    uint8_t m_block[64];
    size_t m_len = 0;
    uint64_t m_tot_len = 0;
};

} // namespace crypto

struct InputFileProvenance {
    std::string path;
    std::string format;
    std::string sha256;
    size_t size_bytes = 0;

    nlohmann::json to_json() const {
        nlohmann::json j;
        j["path"] = path;
        j["format"] = format;
        j["sha256"] = sha256;
        j["size_bytes"] = size_bytes;
        return j;
    }
};

struct InvocationProvenance {
    std::string cli_command;
    std::string working_directory;
    std::map<std::string, std::string> arguments;

    nlohmann::json to_json() const {
        nlohmann::json j;
        j["cli_command"] = cli_command;
        j["working_directory"] = working_directory;
        j["arguments"] = arguments;
        return j;
    }
};

struct ExecutionProvenance {
    std::string start_time;
    std::string end_time;
    double wall_time_seconds = 0.0;
    int cpu_threads = 1;
    std::string hostname;
    std::string os;
    std::string compiler;

    nlohmann::json to_json() const {
        nlohmann::json j;
        j["start_time"] = start_time;
        j["end_time"] = end_time;
        j["wall_time_seconds"] = wall_time_seconds;
        j["cpu_threads"] = cpu_threads;
        j["hostname"] = hostname;
        j["os"] = os;
        j["compiler"] = compiler;
        return j;
    }
};

struct Provenance {
    InvocationProvenance invocation;
    std::map<std::string, InputFileProvenance> inputs;
    ExecutionProvenance execution;

    static std::string current_iso8601() {
        auto now = std::chrono::system_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
        auto timer = std::chrono::system_clock::to_time_t(now);
        std::tm bt;
#if defined(_WIN32)
        gmtime_s(&bt, &timer);
#else
        gmtime_r(&timer, &bt);
#endif
        std::ostringstream ss;
        ss << std::put_time(&bt, "%Y-%m-%dT%H:%M:%S")
           << '.' << std::setfill('0') << std::setw(3) << ms.count() << "Z";
        return ss.str();
    }

    static std::string detect_os() {
#if defined(__APPLE__) && defined(__MACH__)
        return "macOS";
#elif defined(__linux__)
        return "Linux";
#elif defined(_WIN32)
        return "Windows";
#else
        return "Unknown";
#endif
    }

    static std::string detect_compiler() {
#if defined(__clang__)
        return "Clang " + std::to_string(__clang_major__) + "." + std::to_string(__clang_minor__) + "." + std::to_string(__clang_patchlevel__);
#elif defined(__GNUC__)
        return "GCC " + std::to_string(__GNUC__) + "." + std::to_string(__GNUC_MINOR__) + "." + std::to_string(__GNUC_PATCHLEVEL__);
#elif defined(_MSC_VER)
        return "MSVC " + std::to_string(_MSC_VER);
#else
        return "Unknown C++ Compiler";
#endif
    }

    static std::string get_hostname() {
        char buffer[256];
        if (::gethostname(buffer, sizeof(buffer)) == 0) {
            return std::string(buffer);
        }
        return "localhost";
    }

    static std::string get_cwd() {
        char buffer[1024];
        if (::getcwd(buffer, sizeof(buffer)) != nullptr) {
            return std::string(buffer);
        }
        return "";
    }

    nlohmann::json to_json() const {
        nlohmann::json j;
        j["invocation"] = invocation.to_json();
        
        nlohmann::json in_map = nlohmann::json::object();
        for (const auto& [name, file_prov] : inputs) {
            in_map[name] = file_prov.to_json();
        }
        j["inputs"] = in_map;
        j["execution"] = execution.to_json();
        return j;
    }
};

enum class JSONFormat {
    ModernV3,
    Legacy
};

inline JSONFormat parse_json_format(const std::string& fmt) {
    if (fmt == "legacy" || fmt == "v2" || fmt == "datamonkey") {
        return JSONFormat::Legacy;
    }
    return JSONFormat::ModernV3;
}

} // namespace hyphy::core
