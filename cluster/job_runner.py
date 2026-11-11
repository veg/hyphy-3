#!/usr/bin/env python3
"""
Isolated Slurm compute-node runner for WALRUS web offload jobs.
Executed inside the job's scratch directory with zero shell interpolation.
"""
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

MAX_FILE_BYTES = 512 * 1024 * 1024  # 512 MB hard RLIMIT_FSIZE safety cap


def apply_limits():
    try:
        resource.setrlimit(resource.RLIMIT_FSIZE, (MAX_FILE_BYTES, MAX_FILE_BYTES))
    except Exception:
        pass
    try:
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    except Exception:
        pass


def main() -> int:
    job_dir = Path(os.environ.get("WALRUS_JOB_DIR", os.getcwd())).resolve()
    spec_path = job_dir / "spec.json"
    state_path = job_dir / "runner_state.json"

    if not spec_path.is_file():
        sys.stderr.write("Missing spec.json\n")
        return 2

    with open(spec_path, "r", encoding="utf-8") as fh:
        spec = json.load(fh)

    hyphy3_bin = spec["hyphy3_bin"]
    threads = int(spec.get("threads", 16))
    opts = spec.get("options", {})

    tmp_dir = job_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    node_name = os.uname().nodename
    started_at = time.time()
    with open(state_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "state": "running",
                "node": node_name,
                "started_at": started_at,
                "pid": os.getpid(),
            },
            fh,
        )

    seq_file = job_dir / "input.seq"
    tree_file = job_dir / "input.nwk"
    out_json = job_dir / "result.walrus.json"
    out_html = job_dir / "result.walrus.html"
    out_msa = job_dir / "result.msa.fa"
    stdout_log = job_dir / "walrus.log"

    argv = [hyphy3_bin, "walrus"]
    if opts.get("mode", "unaligned") == "aligned":
        argv.extend(["-a", str(seq_file)])
    else:
        argv.extend(["-u", str(seq_file), "--out-msa", str(out_msa)])

    if tree_file.is_file() and tree_file.stat().st_size > 0:
        argv.extend(["-t", str(tree_file)])

    argv.extend(["--code", str(opts.get("genetic_code", "Universal"))])
    argv.extend(["--threads", str(threads)])

    if opts.get("coevolution", True):
        argv.append("--coevolution")

    srv = opts.get("srv", "")
    if srv in ("Yes", "No"):
        argv.extend(["--srv", srv])

    mh = opts.get("multiple_hits", "")
    if mh in ("None", "Double", "Double+Triple"):
        argv.extend(["--multiple-hits", mh])

    argv.extend(["-o", str(out_json), "--html", str(out_html)])

    clean_env = {
        "PATH": "/usr/bin:/bin",
        "HOME": str(tmp_dir),
        "TMPDIR": str(tmp_dir),
        "OMP_NUM_THREADS": str(threads),
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "LC_ALL": "C",
    }
    if os.environ.get("LD_LIBRARY_PATH"):
        clean_env["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"]

    apply_limits()

    rc = 1
    err_msg = ""
    try:
        with open(stdout_log, "w", encoding="utf-8", buffering=1) as log_fh:
            proc = subprocess.Popen(
                argv,
                cwd=str(job_dir),
                env=clean_env,
                stdin=subprocess.DEVNULL,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                text=True,
            )
            try:
                rc = proc.wait(timeout=int(spec.get("timeout_sec", 1140)))
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
                rc = 124
                err_msg = "Job exceeded the 19-minute cluster execution time limit."
    except Exception as exc:
        rc = 1
        err_msg = f"Runner failure: {exc}"
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    completed_at = time.time()
    ok = rc == 0 and out_json.is_file() and out_json.stat().st_size > 32
    if not ok and not err_msg:
        err_msg = f"hyphy3 walrus exited with status {rc}."

    with open(state_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "state": "completed" if ok else "failed",
                "node": node_name,
                "started_at": started_at,
                "completed_at": completed_at,
                "elapsed_sec": round(completed_at - started_at, 2),
                "exit_code": rc,
                "error": err_msg,
            },
            fh,
        )

    return 0 if ok else (rc if rc != 0 else 1)


if __name__ == "__main__":
    sys.exit(main())
