#pragma once
//
// HyphAeon 384D transformer embeddings bridge.
//
// hyphy-3 does not link against PyTorch. Instead, when the user explicitly opts in
// (--hyphaeon), we run a small embedded Python helper in a subprocess. The helper
// depends ONLY on the public `aeon-core` package (https://github.com/veg/HyphAeon)
// and on model weights hosted at https://huggingface.co/datamonkey/hyphaeon.
//
// There are no hard-coded filesystem paths:
//   * Python interpreter : $HYPHY3_PYTHON, otherwise `python3` on PATH
//   * Weights            : $HYPHAEON_WEIGHTS (explicit file), otherwise the Hugging Face
//                          cache used by aeon-core (~/.cache/huggingface, or $HYPHAEON_CACHE/..)
//   * Variant            : $HYPHAEON_VARIANT (default: aeon-core's default, "general")
//
// If the weights are not cached, the user is prompted (interactive terminals) before
// anything is downloaded. Non-interactive runs require HYPHY3_HYPHAEON_DOWNLOAD=1.
// Any failure is reported as an exception: HyphAeon is never silently skipped.
//

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"

#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef _WIN32
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
#endif

namespace hyphy::core::hyphaeon {

inline constexpr const char* kWeightsRepo = "https://huggingface.co/datamonkey/hyphaeon";
inline constexpr const char* kInstallHint =
    "pip install \"aeon-core @ git+https://github.com/veg/HyphAeon#subdirectory=aeon-core\"";

// Exit codes shared with the embedded helper script.
inline constexpr int kExitOk = 0;
inline constexpr int kExitMissingPackage = 3;
inline constexpr int kExitMissingWeights = 4;
inline constexpr int kExitFailure = 5;

// Embedded Python helper. The embedding computation reproduces
// stupid_trees.diffusion_batched.extract_chunked_multiscale_representations
// (taxon-representation branch only; the attention maps it also returns are unused here).
inline const char* helper_script() {
    return R"PY(
import argparse, os, struct, sys

EXIT_OK, EXIT_NO_PKG, EXIT_NO_WEIGHTS, EXIT_FAIL = 0, 3, 4, 5


def eprint(*a):
    print(*a, file=sys.stderr, flush=True)


def taxon_representations(model, c_tensor, a_tensor, tree_cache, device, chunk_size=16):
    import torch
    model.eval()
    total_sites, num_species, window_size = c_tensor.shape
    central_idx = window_size // 2
    num_nodes = num_species + 1
    static_biases = tree_cache["static_phylo_biases"]
    static_coss = tree_cache["static_rope_coss"]
    static_sins = tree_cache["static_rope_sins"]
    accum = torch.zeros((num_species, model.embed_dim), dtype=torch.float32, device=device)
    with torch.no_grad():
        for start in range(0, total_sites, chunk_size):
            end = min(start + chunk_size, total_sites)
            c_chunk = c_tensor[start:end].to(device)
            a_chunk = a_tensor[start:end].to(device)
            bsz = c_chunk.shape[0]
            x = torch.cat([model.codon_embedding(c_chunk), model.aa_embedding(a_chunk)], dim=-1)
            x = x + model.pos_embedding.unsqueeze(1) + tree_cache.get("mds_pos_static").unsqueeze(2)
            root = model.root_token.expand(bsz, 1, window_size, -1)
            x_full = torch.cat([root, x], dim=1)
            del root, x
            for i in range(len(model.col_layers)):
                col_in = (x_full.reshape(bsz * num_nodes, window_size, model.embed_dim) if window_size > 1
                          else x_full.reshape(bsz * num_nodes, model.embed_dim))
                x_full = model.col_layers[i](col_in).reshape(bsz, num_nodes, window_size, model.embed_dim)
            sub_total = bsz * window_size
            for i, layer in enumerate(model.row_layers):
                row_in = x_full.transpose(1, 2).contiguous().view(sub_total, num_nodes, model.embed_dim)
                q = layer.q_proj(row_in).view(sub_total, num_nodes, layer.num_heads, layer.head_dim).transpose(1, 2)
                k = layer.k_proj(row_in).view(sub_total, num_nodes, layer.num_heads, layer.head_dim).transpose(1, 2)
                v = layer.v_proj(row_in).view(sub_total, num_nodes, layer.num_heads, layer.head_dim).transpose(1, 2)
                cos, sin = static_coss[i], static_sins[i]
                h = layer.head_dim // 2
                q1, q2 = q[..., :h], q[..., h:]
                k1, k2 = k[..., :h], k[..., h:]
                q = torch.cat([q1 * cos - q2 * sin, q1 * sin + q2 * cos], dim=-1)
                k = torch.cat([k1 * cos - k2 * sin, k1 * sin + k2 * cos], dim=-1)
                scores = torch.matmul(q, k.transpose(-2, -1)) / (layer.head_dim ** 0.5)
                del q, k
                scores.add_(static_biases[i])
                attn = torch.softmax(scores, dim=-1)
                del scores
                out = torch.matmul(attn, v)
                del attn, v
                out = out.transpose(1, 2).contiguous().view(sub_total, num_nodes, layer.embed_dim)
                out = layer.out_proj(out) + layer.alpha_skip * row_in
                x_full = model.row_norms[i](row_in + out).reshape(bsz, window_size, num_nodes, model.embed_dim).transpose(1, 2)
            accum.add_(x_full[:, 1:, central_idx, :].sum(dim=0))
        if device.type == "mps":
            torch.mps.synchronize()
        elif device.type == "cuda":
            torch.cuda.synchronize()
    return (accum / total_sites).cpu().numpy()


def resolve_weights(W, variant, explicit, allow_download):
    """Return (weights_path, config_dict). Never touches the network unless allow_download."""
    import json
    if explicit:
        if not os.path.exists(explicit):
            eprint(f"HYPHAEON_ERROR: weights file not found: {explicit}")
            sys.exit(EXIT_FAIL)
        return explicit, {}
    from huggingface_hub import try_to_load_from_cache
    v = variant or W.DEFAULT_VARIANT
    fname = W.get_variant_filename(v)
    # aeon-core >= 0.2.0 exposes the Hugging Face hub cache as HF_HUB_CACHE; older releases had CACHE_DIR.
    hub_cache = getattr(W, "HF_HUB_CACHE", None)
    if hub_cache is not None:
        cache_dir = str(hub_cache)
    elif hasattr(W, "CACHE_DIR"):
        cache_dir = str(W.CACHE_DIR.parent / "huggingface")
    else:
        from huggingface_hub import constants as hf_constants
        cache_dir = str(hf_constants.HF_HUB_CACHE)

    def cached(name):
        try:
            p = try_to_load_from_cache(repo_id=W.HF_REPO_ID, filename=name, cache_dir=cache_dir)
            return p if isinstance(p, str) and os.path.exists(p) else None
        except Exception:
            return None

    path = cached(fname)
    if path is None:
        if not allow_download:
            eprint(f"HYPHAEON_NO_WEIGHTS: {W.HF_REPO_ID}/{fname} is not in the local cache ({cache_dir})")
            sys.exit(EXIT_NO_WEIGHTS)
        from huggingface_hub import hf_hub_download
        eprint(f"Downloading {W.HF_REPO_ID}/{fname} ...")
        path = hf_hub_download(repo_id=W.HF_REPO_ID, filename=fname, cache_dir=cache_dir)
        for cfg_name in (f"model.{v}.config.json", "config.json"):
            try:
                hf_hub_download(repo_id=W.HF_REPO_ID, filename=cfg_name, cache_dir=cache_dir)
            except Exception:
                pass
    cfg = {}
    for cfg_name in (f"model.{v}.config.json", "config.json"):
        p = cached(cfg_name)
        if p:
            with open(p) as fh:
                cfg = json.load(fh)
            break
    return path, cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--allow-download", action="store_true")
    ap.add_argument("--input")
    ap.add_argument("--output")
    a = ap.parse_args()
    try:
        import numpy as np
        import torch
        from aeon_core import weights as W
        from aeon_core.model import PhyloAxialTransformer
        from aeon_core.inference import get_device
        from aeon_core.dataset import load_alignment_and_tree
    except Exception as e:
        eprint(f"HYPHAEON_MISSING_PACKAGE: {type(e).__name__}: {e}")
        return EXIT_NO_PKG

    path, cfg = resolve_weights(W, os.environ.get("HYPHAEON_VARIANT") or None,
                                os.environ.get("HYPHAEON_WEIGHTS") or None, a.allow_download)
    if a.check:
        print(path)
        return EXIT_OK

    arch = {"embed_dim": 384, "num_layers": 6, "num_heads": 12, "window_size": 1}
    arch.update({k: cfg[k] for k in arch if k in cfg})
    device = get_device(cpu=os.environ.get("HYPHAEON_DEVICE", "").lower() == "cpu")
    model = PhyloAxialTransformer(**arch).to(device)
    model.load_state_dict(W.load_weights(weights=path, map_location=device), strict=False)
    model.eval()

    c_t, a_t, d_t, z_t, _, taxa, _ = load_alignment_and_tree(a.input, use_tn93=True, prune_duplicates=False)
    dist_mat = d_t.squeeze(0).cpu().numpy().astype(np.float64)
    mds_coords = z_t.squeeze(0).cpu().numpy()
    tree_cache = model.precompute_tree_cache(dist_mat, mds_coords)
    arr = np.asarray(taxon_representations(model, c_t, a_t, tree_cache, device), dtype=np.float64)
    with open(a.output, "wb") as fh:
        fh.write(struct.pack("<II", arr.shape[0], arr.shape[1]))
        fh.write(arr.tobytes())
    with open(a.output + ".taxa", "w") as fh:
        fh.write("\n".join(str(t) for t in taxa) + "\n")
    return EXIT_OK


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(EXIT_FAIL)
)PY";
}

namespace detail {

inline std::string env_or(const char* name, const std::string& fallback) {
    const char* v = std::getenv(name);
    return (v && *v) ? std::string(v) : fallback;
}

inline std::string tail_of_file(const std::filesystem::path& p, size_t max_lines = 12) {
    std::ifstream in(p);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) lines.push_back(line);
    std::ostringstream os;
    size_t start = lines.size() > max_lines ? lines.size() - max_lines : 0;
    for (size_t i = start; i < lines.size(); ++i) os << "    " << lines[i] << "\n";
    return os.str();
}

// Private temporary directory removed on scope exit.
struct TempDir {
    std::filesystem::path path;
    TempDir() {
#ifndef _WIN32
        std::string tmpl = (std::filesystem::temp_directory_path() / "hyphy3_hyphaeon_XXXXXX").string();
        std::vector<char> buf(tmpl.begin(), tmpl.end());
        buf.push_back('\0');
        if (!mkdtemp(buf.data())) {
            throw std::runtime_error("HyphAeon: could not create a temporary directory");
        }
        path = buf.data();
#endif
    }
    ~TempDir() {
        std::error_code ec;
        if (!path.empty()) std::filesystem::remove_all(path, ec);
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
};

// Run `python script args...` without a shell; stdout+stderr go to log_path.
inline int run_python(const std::string& python, const std::vector<std::string>& args,
                      const std::filesystem::path& log_path) {
#ifdef _WIN32
    (void)python; (void)args; (void)log_path;
    throw std::runtime_error("HyphAeon: subprocess support is not implemented on Windows");
#else
    std::vector<std::string> owned;
    owned.push_back(python);
    owned.insert(owned.end(), args.begin(), args.end());
    std::vector<char*> argv;
    for (auto& s : owned) argv.push_back(s.data());
    argv.push_back(nullptr);

    posix_spawn_file_actions_t fa;
    posix_spawn_file_actions_init(&fa);
    posix_spawn_file_actions_addopen(&fa, STDIN_FILENO, "/dev/null", O_RDONLY, 0);
    posix_spawn_file_actions_addopen(&fa, STDOUT_FILENO, log_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
    posix_spawn_file_actions_adddup2(&fa, STDOUT_FILENO, STDERR_FILENO);

    pid_t pid = 0;
    int rc = posix_spawnp(&pid, python.c_str(), &fa, nullptr, argv.data(), environ);
    posix_spawn_file_actions_destroy(&fa);
    if (rc != 0) {
        throw std::runtime_error(
            "HyphAeon: could not start the Python interpreter '" + python + "' (" + std::strerror(rc) +
            "). Set HYPHY3_PYTHON to a Python 3 interpreter that has aeon-core installed:\n    " +
            kInstallHint);
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR) throw std::runtime_error("HyphAeon: waitpid failed");
    }
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    return kExitFailure;
#endif
}

inline bool is_interactive() {
#ifdef _WIN32
    return false;
#else
    return isatty(STDIN_FILENO) && isatty(STDERR_FILENO);
#endif
}

// Decide whether a missing-weights download is allowed (prompting if interactive).
inline bool download_permitted() {
    std::string pref = env_or("HYPHY3_HYPHAEON_DOWNLOAD", "");
    for (auto& ch : pref) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (pref == "1" || pref == "yes" || pref == "true") return true;
    if (pref == "0" || pref == "no" || pref == "false") return false;
    if (!is_interactive()) return false;
    std::cerr << "\nHyphAeon model weights are not cached on this machine.\n"
              << "Download them now from " << kWeightsRepo << " (~8 MB, cached for future runs)? [y/N] "
              << std::flush;
    std::string answer;
    if (!std::getline(std::cin, answer)) return false;
    return !answer.empty() && (answer[0] == 'y' || answer[0] == 'Y');
}

}  // namespace detail

// Cheap probe: true only when aeon-core imports AND the model weights are
// already cached under HYPHY3_PYTHON (default python3), i.e. extract_embeddings
// would succeed without a download prompt. Lets callers request HyphAeon only
// when it can actually run (e.g. tests degrading gracefully in CI).
inline bool is_available() {
    try {
        detail::TempDir tmp;
        const auto script = tmp.path / "hyphaeon_check.py";
        const auto log = tmp.path / "hyphaeon_check.log";
        {
            std::ofstream s(script);
            s << helper_script();
        }
        const std::string python = detail::env_or("HYPHY3_PYTHON", "python3");
        return detail::run_python(python, {script.string(), "--check"}, log) == kExitOk;
    } catch (...) {
        return false;
    }
}

// Compute [N x 384] HyphAeon taxon embeddings for `aln`. Throws std::runtime_error with an
// actionable message on any failure (missing interpreter, package, weights, or runtime error).
inline Matrix extract_embeddings(const Alignment& aln) {
    const size_t N = aln.num_taxa;
    if (N == 0) return Matrix();

    detail::TempDir tmp;
    const auto script = tmp.path / "hyphaeon_embed.py";
    const auto fasta = tmp.path / "input.fasta";
    const auto out_bin = tmp.path / "embeddings.bin";
    const auto log = tmp.path / "hyphaeon.log";
    {
        std::ofstream s(script);
        s << helper_script();
        std::ofstream f(fasta);
        for (size_t i = 0; i < N; ++i) {
            f << ">" << aln.taxon_names[i] << "\n" << aln.sequences[i].raw_seq << "\n";
        }
    }

    const std::string python = detail::env_or("HYPHY3_PYTHON", "python3");

    // 1. Availability check (offline): packages importable + weights cached.
    int rc = detail::run_python(python, {script.string(), "--check"}, log);
    bool allow_download = false;
    if (rc == kExitMissingPackage) {
        throw std::runtime_error(
            "HyphAeon was requested (--hyphaeon) but the Python package 'aeon-core' (with torch) "
            "could not be imported by '" + python + "'.\n  Install it with:\n    " + kInstallHint +
            "\n  or point HYPHY3_PYTHON at an interpreter where it is installed. Details:\n" +
            detail::tail_of_file(log));
    }
    if (rc == kExitMissingWeights) {
        if (!detail::download_permitted()) {
            throw std::runtime_error(
                std::string("HyphAeon was requested (--hyphaeon) but the model weights are not cached, "
                            "and downloading was not approved.\n  Re-run interactively and answer 'y', set "
                            "HYPHY3_HYPHAEON_DOWNLOAD=1 to allow the download (") + kWeightsRepo +
                "), or set HYPHAEON_WEIGHTS to a local weights file.");
        }
        allow_download = true;
    } else if (rc != kExitOk) {
        throw std::runtime_error("HyphAeon availability check failed (exit " + std::to_string(rc) +
                                 "). Details:\n" + detail::tail_of_file(log));
    }

    // 2. Embedding extraction.
    std::vector<std::string> args = {script.string(), "--input", fasta.string(), "--output", out_bin.string()};
    if (allow_download) args.push_back("--allow-download");
    rc = detail::run_python(python, args, log);
    if (rc != kExitOk || !std::filesystem::exists(out_bin)) {
        throw std::runtime_error("HyphAeon embedding extraction failed (exit " + std::to_string(rc) +
                                 "). Details:\n" + detail::tail_of_file(log));
    }

    // 3. Parse [uint32 rows, uint32 cols] + float64 row-major payload.
    std::ifstream in(out_bin, std::ios::binary);
    uint32_t hdr[2] = {0, 0};
    in.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
    if (!in || hdr[0] != N || hdr[1] == 0) {
        throw std::runtime_error("HyphAeon returned " + std::to_string(hdr[0]) + " embeddings for " +
                                 std::to_string(N) + " taxa (duplicate or invalid sequences?)");
    }
    std::vector<double> buf(static_cast<size_t>(hdr[0]) * hdr[1]);
    in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size() * sizeof(double)));
    if (!in) throw std::runtime_error("HyphAeon: truncated embeddings file");
    Matrix raw(static_cast<Eigen::Index>(hdr[0]), static_cast<Eigen::Index>(hdr[1]));
    for (uint32_t r = 0; r < hdr[0]; ++r)
        for (uint32_t c = 0; c < hdr[1]; ++c)
            raw(r, c) = buf[static_cast<size_t>(r) * hdr[1] + c];

    // 4. Map rows back to the alignment's taxon order by name (exact, then first token).
    std::vector<std::string> row_names;
    {
        std::ifstream tn(out_bin.string() + ".taxa");
        std::string line;
        while (std::getline(tn, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty()) row_names.push_back(line);
        }
    }
    if (row_names.size() != N) {
        throw std::runtime_error("HyphAeon: embedding row labels are missing or incomplete");
    }
    auto first_token = [](const std::string& s) { return s.substr(0, s.find_first_of(" \t")); };
    std::unordered_map<std::string, std::vector<size_t>> by_exact, by_token;
    for (size_t r = N; r-- > 0;) {  // reverse so pop_back() yields the first occurrence
        by_exact[row_names[r]].push_back(r);
        by_token[first_token(row_names[r])].push_back(r);
    }
    Matrix M(raw.rows(), raw.cols());
    std::vector<bool> used(N, false);
    auto take = [&](std::unordered_map<std::string, std::vector<size_t>>& idx, const std::string& key) -> size_t {
        auto it = idx.find(key);
        if (it == idx.end()) return N;
        auto& v = it->second;
        while (!v.empty() && used[v.back()]) v.pop_back();
        if (v.empty()) return N;
        size_t r = v.back();
        v.pop_back();
        return r;
    };
    for (size_t i = 0; i < N; ++i) {
        const std::string& want = aln.taxon_names[i];
        size_t hit = take(by_exact, want);
        if (hit == N) hit = take(by_token, first_token(want));
        if (hit == N) {
            throw std::runtime_error("HyphAeon: no embedding returned for taxon '" + want + "'");
        }
        used[hit] = true;
        M.row(static_cast<Eigen::Index>(i)) = raw.row(static_cast<Eigen::Index>(hit));
    }
    return M;
}

}  // namespace hyphy::core::hyphaeon
