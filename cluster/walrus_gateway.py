#!/usr/bin/env python3
from __future__ import annotations
"""
Zero-account, abuse-resistant HTTP gateway for offloading WALRUS jobs from the
HyPhy 3 WebAssembly workbench to the iGEM Slurm cluster (silverback).

Security & Load-Control Architecture:
1. Binds strictly to 127.0.0.1:8787 behind an outbound cloudflared HTTPS tunnel
   (zero inbound firewall ports opened on silverback).
2. Stateless HMAC-signed SHA-256 Proof-of-Work challenge with adaptive difficulty
   and single-use replay table (blocks automated bot spam without user logins).
3. Privacy-preserving subnet rate limiting (salted SHA-256 of /24 IPv4 or /48 IPv6):
   max 2 active jobs and max 12 jobs/hour per client subnet.
4. Strict in-memory biological & structural input validation before disk write:
   3..1,500 coding sequences, <= 30,000 nt/seq, <= 12 MB payload, >= 78% IUPAC
   nucleotide characters, balanced Newick parentheses, enum-only CLI flags.
5. Zero shell command interpolation (shell=False everywhere; static Slurm batch
   wrapper invoking job_runner.py).
6. Hard Slurm concurrency & priority guardrails: max 4 concurrent 16-core jobs
   (< 3.2% of silverback's 2,048 ARM64 cores), --nice=10000, --time=00:20:00,
   content-addressed result caching, and automatic 48-hour TTL storage reaper.
"""
import gzip
import hashlib
import hmac
import ipaddress
import json
import os
import re
import secrets
import shutil
import sqlite3
import subprocess
import sys
import threading
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

SERVICE_ROOT = Path(os.environ.get("WALRUS_SERVICE_ROOT", Path.home() / "walrus-service")).resolve()
JOBS_DIR = SERVICE_ROOT / "jobs"
DB_PATH = SERVICE_ROOT / "gateway" / "jobs.sqlite3"
SECRET_PATH = SERVICE_ROOT / "gateway" / "hmac_secret.bin"
HYPHY3_BIN = Path(
    os.environ.get("WALRUS_HYPHY3_BIN", SERVICE_ROOT / "hyphy-3" / "build" / "hyphy3")
).resolve()
RUNNER_PY = Path(__file__).resolve().parent / "job_runner.py"

BIND_HOST = os.environ.get("WALRUS_BIND_HOST", "127.0.0.1")
BIND_PORT = int(os.environ.get("WALRUS_BIND_PORT", "18791"))

SLURM_PARTITION = os.environ.get("WALRUS_SLURM_PARTITION", "datamonkey")
THREADS_PER_JOB = int(os.environ.get("WALRUS_THREADS_PER_JOB", "16"))
MEM_PER_JOB = os.environ.get("WALRUS_MEM_PER_JOB", "24G")
WALLTIME_SLURM = os.environ.get("WALRUS_WALLTIME", "00:20:00")
RUNNER_TIMEOUT_SEC = int(os.environ.get("WALRUS_RUNNER_TIMEOUT_SEC", "1140"))

MAX_CONCURRENT_SLURM = int(os.environ.get("WALRUS_MAX_CONCURRENT", "4"))
MAX_PENDING_QUEUE = int(os.environ.get("WALRUS_MAX_QUEUE", "25"))
MAX_ACTIVE_PER_SUBNET = int(os.environ.get("WALRUS_MAX_ACTIVE_PER_SUBNET", "2"))
MAX_HOURLY_PER_SUBNET = int(os.environ.get("WALRUS_MAX_HOURLY_PER_SUBNET", "12"))

MAX_UPLOAD_BYTES = 12 * 1024 * 1024       # 12 MB compressed/raw HTTP body
MAX_DECOMPRESSED_BYTES = 36 * 1024 * 1024 # 36 MB decompressed JSON payload
MIN_TAXA = 3
MAX_TAXA = 1500
MAX_SEQ_UNGAPPED_NT = 30000
MAX_TOTAL_NT = 15_000_000
JOB_TTL_SECONDS = 48 * 3600
MAX_JOBS_DISK_BYTES = 20 * 1024 * 1024 * 1024  # 20 GB total job cache ceiling
MAX_TELEMETRY_BYTES = 2048
MAX_TELEMETRY_PER_MIN = 45
TELEMETRY_RETENTION_DAYS = 30

ALLOWED_TELEMETRY_EVENTS = {
    "app_open",
    "sample_load",
    "custom_input",
    "stage0_done",
    "stage1_done",
    "stage2_done",
    "stage3_done",
    "ping_done",
    "render_json_done",
    "cluster_submit",
    "cluster_done",
    "export_artifact",
    "wasm_error",
    "worker_crash",
    "ui_error",
    "cluster_error",
}

ALLOWED_ORIGINS = {
    "https://veg.github.io",
    "https://spond.github.io",
    "https://hyphy.org",
    "https://www.hyphy.org",
    "https://datamonkey.org",
    "https://www.datamonkey.org",
}

ALLOWED_GENETIC_CODES = {
    "Universal",
    "Vertebrate-mtDNA",
    "Yeast-mtDNA",
    "Mold-Protozoan-mtDNA",
    "Invertebrate-mtDNA",
    "Ciliate-Nuclear",
    "Echinoderm-mtDNA",
    "Euplotid-Nuclear",
    "Alt-Yeast-Nuclear",
    "Ascidian-mtDNA",
    "Flatworm-mtDNA",
    "Blepharisma-Nuclear",
}

IUPAC_NUC_CHARS = set("ACGTRYMKSWHBVDNacgtrymkswhbvdnUu-.?")

DB_LOCK = threading.Lock()
WAKE_EVENT = threading.Event()
TELEMETRY_RATE_LOCK = threading.Lock()
TELEMETRY_RATE_BUCKETS: dict[str, tuple[int, int]] = {}


def load_or_create_secret() -> bytes:
    SECRET_PATH.parent.mkdir(parents=True, exist_ok=True)
    if SECRET_PATH.is_file():
        data = SECRET_PATH.read_bytes()
        if len(data) >= 32:
            return data
    sec = secrets.token_bytes(32)
    SECRET_PATH.write_bytes(sec)
    os.chmod(SECRET_PATH, 0o600)
    return sec


SERVER_SECRET = load_or_create_secret()


def get_db() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH), timeout=15.0, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    return conn


def init_db() -> None:
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    with DB_LOCK, get_db() as conn:
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                job_id TEXT PRIMARY KEY,
                cancel_token TEXT NOT NULL,
                cache_key TEXT NOT NULL,
                subnet_hash TEXT NOT NULL,
                label TEXT NOT NULL,
                num_taxa INTEGER NOT NULL,
                max_nt INTEGER NOT NULL,
                options_json TEXT NOT NULL,
                status TEXT NOT NULL,
                slurm_id TEXT DEFAULT '',
                node TEXT DEFAULT '',
                error TEXT DEFAULT '',
                created_at REAL NOT NULL,
                started_at REAL DEFAULT 0,
                completed_at REAL DEFAULT 0
            );
            CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status, created_at);
            CREATE INDEX IF NOT EXISTS idx_jobs_cache ON jobs(cache_key, status);
            CREATE INDEX IF NOT EXISTS idx_jobs_subnet ON jobs(subnet_hash, created_at);

            CREATE TABLE IF NOT EXISTS used_challenges (
                sig TEXT PRIMARY KEY,
                expires_at INTEGER NOT NULL
            );

            CREATE TABLE IF NOT EXISTS telemetry_events (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ts REAL NOT NULL,
                day TEXT NOT NULL,
                day_visitor TEXT NOT NULL,
                country TEXT DEFAULT '',
                event TEXT NOT NULL,
                mode TEXT DEFAULT '',
                stage TEXT DEFAULT '',
                sample_id TEXT DEFAULT '',
                n_taxa INTEGER DEFAULT 0,
                n_sites INTEGER DEFAULT 0,
                wall_ms INTEGER DEFAULT 0,
                threads INTEGER DEFAULT 0,
                err_kind TEXT DEFAULT '',
                err_msg TEXT DEFAULT ''
            );
            CREATE INDEX IF NOT EXISTS idx_telem_ts ON telemetry_events(ts);
            CREATE INDEX IF NOT EXISTS idx_telem_event ON telemetry_events(event, ts);
            CREATE INDEX IF NOT EXISTS idx_telem_day ON telemetry_events(day, day_visitor);

            CREATE TABLE IF NOT EXISTS telemetry_daily (
                day TEXT NOT NULL,
                event TEXT NOT NULL,
                total_count INTEGER NOT NULL DEFAULT 0,
                unique_visitors INTEGER NOT NULL DEFAULT 0,
                PRIMARY KEY (day, event)
            );
            """
        )


def client_subnet_hash(headers, remote_addr: str) -> str:
    raw_ip = (
        headers.get("CF-Connecting-IP")
        or headers.get("X-Forwarded-For", "").split(",")[0].strip()
        or remote_addr
        or "127.0.0.1"
    )
    try:
        ip_obj = ipaddress.ip_address(raw_ip)
        if ip_obj.version == 4:
            net = ipaddress.ip_network(f"{ip_obj}/24", strict=False)
        else:
            net = ipaddress.ip_network(f"{ip_obj}/48", strict=False)
        norm = str(net)
    except ValueError:
        norm = raw_ip
    return hmac.new(SERVER_SECRET, f"subnet:{norm}".encode("utf-8"), hashlib.sha256).hexdigest()[:24]


def daily_visitor_hash(subnet_h: str, utc_day: str) -> str:
    return hmac.new(SERVER_SECRET, f"visitor:{utc_day}:{subnet_h}".encode("utf-8"), hashlib.sha256).hexdigest()[:16]


def allow_telemetry_beacon(subnet_h: str) -> bool:
    minute_slot = int(time.time() // 60)
    with TELEMETRY_RATE_LOCK:
        if len(TELEMETRY_RATE_BUCKETS) > 4096:
            stale = [k for k, (slot, _) in TELEMETRY_RATE_BUCKETS.items() if slot < minute_slot - 1]
            for k in stale:
                TELEMETRY_RATE_BUCKETS.pop(k, None)
        slot, cnt = TELEMETRY_RATE_BUCKETS.get(subnet_h, (minute_slot, 0))
        if slot != minute_slot:
            slot, cnt = minute_slot, 0
        if cnt >= MAX_TELEMETRY_PER_MIN:
            return False
        TELEMETRY_RATE_BUCKETS[subnet_h] = (slot, cnt + 1)
        return True


def scrub_error_message(raw_msg: str) -> str:
    if not raw_msg:
        return ""
    s = re.sub(r"[\x00-\x1f\x7f]+", " ", str(raw_msg)).strip()
    # Redact quoted identifiers (e.g., user sequence or taxon names) and local file paths
    s = re.sub(r"'[^']{1,120}'", "'<id>'", s)
    s = re.sub(r'"[^"]{1,120}"', '"<id>"', s)
    s = re.sub(r"/(?:Users|home|tmp|var)/[^\s:;,)>\]]+", "<path>", s)
    return s[:240]


def leading_zero_bits(digest: bytes) -> int:
    bits = 0
    for b in digest:
        if b == 0:
            bits += 8
        else:
            bits += 8 - b.bit_length()
            break
    return bits


def current_queue_counts(conn: sqlite3.Connection):
    row = conn.execute(
        """
        SELECT
            SUM(CASE WHEN status IN ('submitted', 'running') THEN 1 ELSE 0 END) AS active,
            SUM(CASE WHEN status = 'queued' THEN 1 ELSE 0 END) AS queued
        FROM jobs
        """
    ).fetchone()
    return int(row["active"] or 0), int(row["queued"] or 0)


def compute_pow_difficulty(active: int, queued: int) -> int:
    load = active + queued
    if load < 4:
        return 18  # ~0.25s in browser WebCrypto
    if load < 10:
        return 19  # ~0.5s
    return 20      # ~1.0-1.5s under heavy queue load


def issue_pow_challenge(subnet_h: str, active: int, queued: int) -> dict:
    difficulty = compute_pow_difficulty(active, queued)
    nonce = secrets.token_hex(8)
    expires_at = int(time.time()) + 300
    body = f"{subnet_h}:{nonce}:{difficulty}:{expires_at}"
    sig = hmac.new(SERVER_SECRET, body.encode("utf-8"), hashlib.sha256).hexdigest()
    token = f"{nonce}.{difficulty}.{expires_at}.{sig}"
    return {
        "token": token,
        "difficulty": difficulty,
        "expires_at": expires_at,
    }


def verify_and_consume_pow(conn: sqlite3.Connection, subnet_h: str, pow_obj: dict) -> tuple[bool, str]:
    if not isinstance(pow_obj, dict):
        return False, "Missing Proof-of-Work object."
    token = str(pow_obj.get("token", ""))
    solution = str(pow_obj.get("solution", ""))
    if not token or not solution or len(solution) > 64:
        return False, "Invalid Proof-of-Work fields."
    parts = token.split(".")
    if len(parts) != 4:
        return False, "Malformed Proof-of-Work token."
    nonce, diff_str, exp_str, sig = parts
    try:
        difficulty = int(diff_str)
        expires_at = int(exp_str)
    except ValueError:
        return False, "Malformed Proof-of-Work numbers."

    now = int(time.time())
    if now > expires_at:
        return False, "Proof-of-Work challenge expired; please retry."
    if difficulty < 16 or difficulty > 26:
        return False, "Invalid Proof-of-Work difficulty."

    expected_body = f"{subnet_h}:{nonce}:{difficulty}:{expires_at}"
    expected_sig = hmac.new(SERVER_SECRET, expected_body.encode("utf-8"), hashlib.sha256).hexdigest()
    if not hmac.compare_digest(sig, expected_sig):
        return False, "Proof-of-Work signature mismatch (IP changed or invalid token)."

    digest = hashlib.sha256(f"{token}:{solution}".encode("utf-8")).digest()
    if leading_zero_bits(digest) < difficulty:
        return False, "Insufficient Proof-of-Work difficulty."

    conn.execute("DELETE FROM used_challenges WHERE expires_at < ?", (now,))
    exists = conn.execute("SELECT 1 FROM used_challenges WHERE sig = ?", (sig,)).fetchone()
    if exists:
        return False, "Proof-of-Work challenge was already used."
    conn.execute("INSERT INTO used_challenges (sig, expires_at) VALUES (?, ?)", (sig, expires_at))
    return True, ""


def validate_sequences_and_tree(raw_seq: str, raw_tree: str) -> tuple[bool, str, dict]:
    if not raw_seq or "\x00" in raw_seq:
        return False, "Sequence payload is empty or contains binary null bytes.", {}

    stripped = raw_seq.lstrip()
    is_nexus = stripped[:6].upper() == "#NEXUS"
    taxa_count = 0
    max_ungapped = 0
    total_ungapped = 0
    nuc_matches = 0
    alpha_total = 0

    if is_nexus:
        # Basic NEXUS sanity validation; HyPhy 3's parser handles full NEXUS grammar
        ntax_m = re.search(r"\bNTAX\s*=\s*(\d+)", raw_seq, flags=re.IGNORECASE)
        nchar_m = re.search(r"\bNCHAR\s*=\s*(\d+)", raw_seq, flags=re.IGNORECASE)
        if not ntax_m or not nchar_m:
            return False, "NEXUS input is missing NTAX or NCHAR dimensions.", {}
        taxa_count = int(ntax_m.group(1))
        max_ungapped = int(nchar_m.group(1))
        total_ungapped = taxa_count * max_ungapped
        if taxa_count < MIN_TAXA or taxa_count > MAX_TAXA:
            return (
                False,
                f"Cluster offload requires between {MIN_TAXA} and {MAX_TAXA} sequences (found NTAX={taxa_count}).",
                {},
            )
        if max_ungapped < 30 or max_ungapped > MAX_SEQ_UNGAPPED_NT:
            return (
                False,
                f"Sequence length NCHAR={max_ungapped} is outside allowed bounds (30..{MAX_SEQ_UNGAPPED_NT} nt).",
                {},
            )
    else:
        if not stripped.startswith(">"):
            return False, "Sequences must be in FASTA (>header) or #NEXUS format.", {}
        cur_len = 0
        seen_header = False
        for line in raw_seq.splitlines():
            s = line.strip()
            if not s:
                continue
            if s.startswith("("):
                break
            if s.startswith(">"):
                if seen_header:
                    if cur_len == 0:
                        return False, "FASTA contains an empty sequence record.", {}
                    max_ungapped = max(max_ungapped, cur_len)
                taxa_count += 1
                if taxa_count > MAX_TAXA:
                    return (
                        False,
                        f"Dataset exceeds the {MAX_TAXA}-sequence cluster limit.",
                        {},
                    )
                seen_header = True
                cur_len = 0
            else:
                for ch in s:
                    if ch in ("-", ".", "?", " ", "\t"):
                        continue
                    cur_len += 1
                    total_ungapped += 1
                    if ch.isalpha():
                        alpha_total += 1
                        if ch in IUPAC_NUC_CHARS:
                            nuc_matches += 1
                if cur_len > MAX_SEQ_UNGAPPED_NT:
                    return (
                        False,
                        f"A sequence exceeds the {MAX_SEQ_UNGAPPED_NT:,} nt coding-sequence limit.",
                        {},
                    )
                if total_ungapped > MAX_TOTAL_NT:
                    return (
                        False,
                        f"Total coding nucleotides exceed the {MAX_TOTAL_NT:,} nt cluster limit.",
                        {},
                    )
        if seen_header:
            if cur_len == 0:
                return False, "FASTA contains an empty trailing sequence record.", {}
            max_ungapped = max(max_ungapped, cur_len)

        if taxa_count < MIN_TAXA:
            return (
                False,
                f"At least {MIN_TAXA} homologous coding sequences are required (found {taxa_count}).",
                {},
            )
        if max_ungapped < 30:
            return False, "Coding sequences must be at least 30 nucleotides (10 codons) long.", {}
        if alpha_total > 0 and (nuc_matches / alpha_total) < 0.78:
            return (
                False,
                "Uploaded sequences look like amino-acid protein sequences rather than nucleotide coding sequences (CDS). WALRUS requires nucleotide CDS (A/C/G/T).",
                {},
            )

    clean_tree = (raw_tree or "").strip()
    if clean_tree:
        if "\x00" in clean_tree or len(clean_tree) > 512 * 1024:
            return False, "Tree string exceeds 512 KB or contains invalid characters.", {}
        if clean_tree.count("(") != clean_tree.count(")"):
            return False, "Provided Newick tree has unbalanced parentheses.", {}
        if not clean_tree.endswith(";"):
            clean_tree += ";"

    return True, "", {
        "num_taxa": taxa_count,
        "max_nt": max_ungapped,
        "clean_tree": clean_tree,
    }


def sanitize_options(raw_opts: dict) -> tuple[bool, str, dict]:
    if not isinstance(raw_opts, dict):
        raw_opts = {}
    mode = str(raw_opts.get("mode", "unaligned"))
    if mode not in ("unaligned", "aligned"):
        return False, "options.mode must be 'unaligned' or 'aligned'.", {}

    gcode = str(raw_opts.get("genetic_code", "Universal"))
    if gcode not in ALLOWED_GENETIC_CODES:
        return False, f"Unsupported genetic code '{gcode}'.", {}

    coev = bool(raw_opts.get("coevolution", True))
    srv = str(raw_opts.get("srv", ""))
    if srv not in ("", "Yes", "No"):
        return False, "options.srv must be 'Yes', 'No', or empty.", {}

    mh = str(raw_opts.get("multiple_hits", ""))
    if mh not in ("", "None", "Double", "Double+Triple"):
        return False, "Invalid options.multiple_hits value.", {}

    return True, "", {
        "mode": mode,
        "genetic_code": gcode,
        "coevolution": coev,
        "srv": srv,
        "multiple_hits": mh,
    }


def parse_stage_progress(log_text: str, status: str) -> tuple[str, int, str]:
    if status == "completed":
        return "Completed — Full Dossier Ready", 100, ""
    if status == "queued":
        return "Queued for 16-core ARM64 compute node", 4, ""
    if status == "submitted":
        return "Allocated in Slurm — starting compute node", 8, ""
    if status == "canceled":
        return "Canceled", 0, ""
    if status == "failed":
        return "Failed", 100, ""

    lines = [ln.strip() for ln in log_text.splitlines() if ln.strip()]
    tail = "\n".join(lines[-10:]) if lines else ""
    joined = "\n".join(lines[-40:]).lower()

    stage_label = "Stage 0: 6-Frame QA & Homology Admission"
    pct = 15
    if "epistasis" in joined or "co-evol" in joined or "simes" in joined:
        stage_label = "Stage 3: Epistasis & Omnibus Forensics"
        pct = 90
    elif "absrel" in joined or "relax" in joined or "lof" in joined:
        stage_label = "Stage 3: Lineage Regimes (RELAX / aBSREL / LoF)"
        pct = 78
    elif "meme" in joined or "fel" in joined or "busted" in joined or "stage 3" in joined:
        stage_label = "Stage 3: Unified Codon Selection Pass"
        pct = 62
    elif "recomb" in joined or "spr" in joined or "tree" in joined or "gc3" in joined or "stage 2" in joined:
        stage_label = "Stage 2: MG94×REV Phylogeny & Recombination"
        pct = 42
    elif "poa" in joined or "pair-hmm" in joined or "align" in joined or "stage 1" in joined:
        stage_label = "Stage 1: Codon Pair-HMM + POA Alignment"
        pct = 26

    return stage_label, pct, tail


def submit_job_to_slurm(job_id: str) -> tuple[bool, str]:
    job_dir = JOBS_DIR / job_id
    slurm_script = job_dir / "run_walrus.slurm"
    slurm_out = job_dir / "slurm.out"

    # Static Slurm batch script: zero user-supplied strings interpolated
    script_body = (
        "#!/bin/bash\n"
        "#SBATCH --job-name=walrus_web\n"
        f"#SBATCH --partition={SLURM_PARTITION}\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task={THREADS_PER_JOB}\n"
        f"#SBATCH --mem={MEM_PER_JOB}\n"
        f"#SBATCH --time={WALLTIME_SLURM}\n"
        "#SBATCH --nice=10000\n"
        "#SBATCH --no-requeue\n"
        f"#SBATCH --output={slurm_out}\n"
        f"#SBATCH --error={slurm_out}\n"
        f"export WALRUS_JOB_DIR={job_dir}\n"
        f"exec /usr/bin/python3 {RUNNER_PY}\n"
    )
    slurm_script.write_text(script_body, encoding="utf-8")
    os.chmod(slurm_script, 0o700)

    try:
        res = subprocess.run(
            ["/usr/bin/sbatch", "--parsable", str(slurm_script)],
            cwd=str(job_dir),
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except Exception as exc:
        return False, f"sbatch invocation error: {exc}"

    if res.returncode != 0:
        return False, (res.stderr or res.stdout or "sbatch failed").strip()[:240]

    slurm_id = res.stdout.strip().split(";")[0].strip()
    if not slurm_id.isdigit():
        return False, f"Unexpected sbatch output: {slurm_id[:120]}"
    return True, slurm_id


def query_active_slurm_ids(slurm_ids: list[str]) -> dict[str, tuple[str, str]]:
    """Return mapping slurm_id -> (slurm_state_code, nodelist) for currently queued/running jobs."""
    if not slurm_ids:
        return {}
    try:
        res = subprocess.run(
            [
                "/usr/bin/squeue",
                "--noheader",
                "-j",
                ",".join(slurm_ids),
                "-o",
                "%i|%t|%N",
            ],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
        )
    except Exception:
        return {}
    out = {}
    if res.returncode == 0:
        for line in res.stdout.splitlines():
            parts = line.strip().split("|")
            if len(parts) >= 3:
                out[parts[0].strip()] = (parts[1].strip(), parts[2].strip())
    return out


def scheduler_loop() -> None:
    last_reap = 0.0
    while True:
        WAKE_EVENT.wait(timeout=2.5)
        WAKE_EVENT.clear()
        try:
            reconcile_and_dispatch()
            now = time.time()
            if now - last_reap > 300.0:
                reap_expired_jobs()
                last_reap = now
        except Exception as exc:
            sys.stderr.write(f"[scheduler] error: {exc}\n")


def reconcile_and_dispatch() -> None:
    with DB_LOCK, get_db() as conn:
        active_rows = conn.execute(
            "SELECT job_id, slurm_id, status, started_at FROM jobs WHERE status IN ('submitted', 'running')"
        ).fetchall()

    if active_rows:
        sids = [r["slurm_id"] for r in active_rows if r["slurm_id"]]
        live_map = query_active_slurm_ids(sids)
        now = time.time()
        with DB_LOCK, get_db() as conn:
            for r in active_rows:
                jid = r["job_id"]
                sid = r["slurm_id"]
                job_dir = JOBS_DIR / jid
                rstate_file = job_dir / "runner_state.json"
                rstate = {}
                if rstate_file.is_file():
                    try:
                        rstate = json.loads(rstate_file.read_text(encoding="utf-8"))
                    except Exception:
                        rstate = {}

                if rstate.get("state") in ("completed", "failed"):
                    conn.execute(
                        """
                        UPDATE jobs
                        SET status = ?, node = ?, started_at = ?, completed_at = ?, error = ?
                        WHERE job_id = ? AND status != 'canceled'
                        """,
                        (
                            rstate["state"],
                            rstate.get("node", ""),
                            float(rstate.get("started_at") or r["started_at"] or now),
                            float(rstate.get("completed_at") or now),
                            str(rstate.get("error", "")),
                            jid,
                        ),
                    )
                    continue

                if sid in live_map:
                    st_code, node = live_map[sid]
                    if st_code in ("R", "CG") and r["status"] != "running":
                        conn.execute(
                            "UPDATE jobs SET status = 'running', node = ?, started_at = ? WHERE job_id = ?",
                            (node, float(rstate.get("started_at") or now), jid),
                        )
                    elif node:
                        conn.execute("UPDATE jobs SET node = ? WHERE job_id = ?", (node, jid))
                else:
                    # Not in squeue and no completed runner_state yet: wait briefly for NFS flush, else inspect files
                    out_json = job_dir / "result.walrus.json"
                    if out_json.is_file() and out_json.stat().st_size > 32:
                        conn.execute(
                            "UPDATE jobs SET status = 'completed', completed_at = ? WHERE job_id = ?",
                            (now, jid),
                        )
                    elif now - float(r["started_at"] or now) > 15.0:
                        slurm_err = ""
                        slurm_out = job_dir / "slurm.out"
                        if slurm_out.is_file():
                            slurm_err = slurm_out.read_text(encoding="utf-8", errors="replace")[-300:].strip()
                        conn.execute(
                            "UPDATE jobs SET status = 'failed', completed_at = ?, error = ? WHERE job_id = ?",
                            (
                                now,
                                slurm_err or "Slurm job terminated before producing output (time or memory limit).",
                                jid,
                            ),
                        )

    # Dispatch queued jobs up to MAX_CONCURRENT_SLURM
    with DB_LOCK, get_db() as conn:
        active_cnt, _ = current_queue_counts(conn)
        slots = max(0, MAX_CONCURRENT_SLURM - active_cnt)
        if slots <= 0:
            return
        next_jobs = conn.execute(
            "SELECT job_id FROM jobs WHERE status = 'queued' ORDER BY created_at ASC LIMIT ?",
            (slots,),
        ).fetchall()

    for row in next_jobs:
        jid = row["job_id"]
        ok, info = submit_job_to_slurm(jid)
        now = time.time()
        with DB_LOCK, get_db() as conn:
            if ok:
                conn.execute(
                    "UPDATE jobs SET status = 'submitted', slurm_id = ?, started_at = ? WHERE job_id = ? AND status = 'queued'",
                    (info, now, jid),
                )
            else:
                conn.execute(
                    "UPDATE jobs SET status = 'failed', completed_at = ?, error = ? WHERE job_id = ?",
                    (now, info, jid),
                )


def reap_expired_jobs() -> None:
    now = time.time()
    cutoff = now - JOB_TTL_SECONDS
    telem_cutoff = now - (TELEMETRY_RETENTION_DAYS * 86400)
    with DB_LOCK, get_db() as conn:
        old_rows = conn.execute(
            "SELECT job_id FROM jobs WHERE created_at < ? AND status IN ('completed', 'failed', 'canceled')",
            (cutoff,),
        ).fetchall()
        for r in old_rows:
            shutil.rmtree(JOBS_DIR / r["job_id"], ignore_errors=True)
            conn.execute("DELETE FROM jobs WHERE job_id = ?", (r["job_id"],))

        # Roll up telemetry events older than TELEMETRY_RETENTION_DAYS into compact daily totals
        conn.execute(
            """
            INSERT INTO telemetry_daily (day, event, total_count, unique_visitors)
            SELECT day, event, COUNT(*), COUNT(DISTINCT day_visitor)
            FROM telemetry_events
            WHERE ts < ?
            GROUP BY day, event
            ON CONFLICT(day, event) DO UPDATE SET
                total_count = telemetry_daily.total_count + excluded.total_count,
                unique_visitors = MAX(telemetry_daily.unique_visitors, excluded.unique_visitors)
            """,
            (telem_cutoff,),
        )
        conn.execute("DELETE FROM telemetry_events WHERE ts < ?", (telem_cutoff,))

    # Enforce total storage ceiling across jobs/
    entries = []
    total_bytes = 0
    for child in JOBS_DIR.iterdir():
        if not child.is_dir():
            continue
        sz = sum(f.stat().st_size for f in child.rglob("*") if f.is_file())
        mtime = child.stat().st_mtime
        entries.append((mtime, sz, child.name, child))
        total_bytes += sz

    if total_bytes > MAX_JOBS_DISK_BYTES:
        entries.sort(key=lambda x: x[0])  # oldest first
        with DB_LOCK, get_db() as conn:
            for _, sz, jid, path in entries:
                if total_bytes <= int(MAX_JOBS_DISK_BYTES * 0.8):
                    break
                st = conn.execute("SELECT status FROM jobs WHERE job_id = ?", (jid,)).fetchone()
                if st and st["status"] in ("queued", "submitted", "running"):
                    continue
                shutil.rmtree(path, ignore_errors=True)
                conn.execute("DELETE FROM jobs WHERE job_id = ?", (jid,))
                total_bytes -= sz


def build_telemetry_stats(conn: sqlite3.Connection) -> dict:
    now = time.time()
    utc_today = time.strftime("%Y-%m-%d", time.gmtime(now))
    t_24h = now - 86400.0
    t_7d = now - (7 * 86400.0)

    today_uv = conn.execute(
        "SELECT COUNT(DISTINCT day_visitor) AS c FROM telemetry_events WHERE day = ?",
        (utc_today,),
    ).fetchone()["c"]

    daily_uv_rows = conn.execute(
        """
        SELECT day, COUNT(DISTINCT day_visitor) AS uv, COUNT(*) AS events
        FROM telemetry_events
        WHERE ts >= ?
        GROUP BY day
        ORDER BY day DESC
        """,
        (t_7d,),
    ).fetchall()

    countries_7d = conn.execute(
        """
        SELECT country, COUNT(DISTINCT day_visitor) AS uv, COUNT(*) AS events
        FROM telemetry_events
        WHERE ts >= ? AND country != ''
        GROUP BY country
        ORDER BY uv DESC, events DESC
        LIMIT 12
        """,
        (t_7d,),
    ).fetchall()

    events_24h = {
        r["event"]: int(r["c"])
        for r in conn.execute(
            "SELECT event, COUNT(*) AS c FROM telemetry_events WHERE ts >= ? GROUP BY event ORDER BY c DESC",
            (t_24h,),
        ).fetchall()
    }

    stage_perf_7d = [
        {
            "event": r["event"],
            "runs": int(r["runs"]),
            "mean_wall_ms": int(round(r["avg_ms"] or 0)),
            "max_wall_ms": int(r["max_ms"] or 0),
            "mean_taxa": round(float(r["avg_taxa"] or 0.0), 1),
            "mean_sites": round(float(r["avg_sites"] or 0.0), 1),
        }
        for r in conn.execute(
            """
            SELECT event, COUNT(*) AS runs, AVG(wall_ms) AS avg_ms, MAX(wall_ms) AS max_ms,
                   AVG(NULLIF(n_taxa, 0)) AS avg_taxa, AVG(NULLIF(n_sites, 0)) AS avg_sites
            FROM telemetry_events
            WHERE ts >= ? AND event IN ('stage0_done', 'stage1_done', 'stage2_done', 'stage3_done', 'ping_done', 'cluster_done')
            GROUP BY event
            ORDER BY runs DESC
            """,
            (t_7d,),
        ).fetchall()
    ]

    top_samples_7d = [
        {"sample_id": r["sample_id"], "loads": int(r["c"])}
        for r in conn.execute(
            """
            SELECT sample_id, COUNT(*) AS c
            FROM telemetry_events
            WHERE ts >= ? AND event = 'sample_load' AND sample_id != ''
            GROUP BY sample_id
            ORDER BY c DESC
            LIMIT 10
            """,
            (t_7d,),
        ).fetchall()
    ]

    top_errors_7d = [
        {
            "event": r["event"],
            "stage": r["stage"],
            "err_kind": r["err_kind"],
            "err_msg": r["err_msg"],
            "count": int(r["c"]),
            "last_seen_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(float(r["last_ts"]))),
        }
        for r in conn.execute(
            """
            SELECT event, stage, err_kind, err_msg, COUNT(*) AS c, MAX(ts) AS last_ts
            FROM telemetry_events
            WHERE ts >= ? AND event IN ('wasm_error', 'worker_crash', 'ui_error', 'cluster_error')
            GROUP BY event, stage, err_kind, err_msg
            ORDER BY c DESC, last_ts DESC
            LIMIT 20
            """,
            (t_7d,),
        ).fetchall()
    ]

    cluster_jobs = {
        r["status"]: int(r["c"])
        for r in conn.execute("SELECT status, COUNT(*) AS c FROM jobs GROUP BY status").fetchall()
    }

    return {
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(now)),
        "privacy": "cookieless daily-rotating HMAC subnet visitor IDs; zero sequence or filename retention",
        "today_utc": utc_today,
        "today_unique_visitors": int(today_uv or 0),
        "daily_activity_7d": [
            {"day": r["day"], "unique_visitors": int(r["uv"]), "events": int(r["events"])}
            for r in daily_uv_rows
        ],
        "top_countries_7d": [
            {"country": r["country"], "unique_visitors": int(r["uv"]), "events": int(r["events"])}
            for r in countries_7d
        ],
        "events_last_24h": events_24h,
        "stage_performance_7d": stage_perf_7d,
        "top_samples_7d": top_samples_7d,
        "cluster_jobs_48h_window": cluster_jobs,
        "top_errors_7d": top_errors_7d,
    }


class WalrusGatewayHandler(BaseHTTPRequestHandler):
    server_version = "WalrusClusterGateway/1.1"

    def log_message(self, fmt: str, *args) -> None:
        sys.stdout.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {fmt % args}\n")
        sys.stdout.flush()

    def _cors_origin(self) -> str:
        origin = self.headers.get("Origin", "")
        if origin in ALLOWED_ORIGINS:
            return origin
        if origin.startswith("http://localhost:") or origin.startswith("http://127.0.0.1:"):
            return origin
        return "https://veg.github.io"

    def _send_cors_headers(self) -> None:
        self.send_header("Access-Control-Allow-Origin", self._cors_origin())
        self.send_header("Vary", "Origin")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header(
            "Access-Control-Allow-Headers",
            "Content-Type, Content-Encoding, Accept-Encoding",
        )
        self.send_header("Access-Control-Max-Age", "86400")

    def _send_json(self, status: int, payload: dict) -> None:
        raw = json.dumps(payload, separators=(",", ":")).encode("utf-8")
        accept_gz = "gzip" in self.headers.get("Accept-Encoding", "") and len(raw) > 1024
        if accept_gz:
            raw = gzip.compress(raw, compresslevel=5)
        self.send_response(status)
        self._send_cors_headers()
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        if accept_gz:
            self.send_header("Content-Encoding", "gzip")
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def _send_file(self, status: int, filepath: Path, content_type: str) -> None:
        if not filepath.is_file():
            self._send_json(HTTPStatus.NOT_FOUND, {"error": "Requested output artifact not found."})
            return
        data = filepath.read_bytes()
        accept_gz = "gzip" in self.headers.get("Accept-Encoding", "") and len(data) > 1024
        if accept_gz:
            data = gzip.compress(data, compresslevel=5)
        self.send_response(status)
        self._send_cors_headers()
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "public, max-age=3600")
        if accept_gz:
            self.send_header("Content-Encoding", "gzip")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _read_json_body(self, max_bytes: int = MAX_UPLOAD_BYTES) -> tuple[bool, str, dict]:
        cl_str = self.headers.get("Content-Length", "")
        if not cl_str.isdigit():
            return False, "Missing Content-Length header.", {}
        cl = int(cl_str)
        if cl <= 0 or cl > max_bytes:
            return False, f"Payload exceeds {max_bytes} byte limit.", {}
        body = self.rfile.read(cl)
        if self.headers.get("Content-Encoding", "").lower() == "gzip":
            try:
                with gzip.GzipFile(fileobj=__import__("io").BytesIO(body)) as gz:
                    body = gz.read(MAX_DECOMPRESSED_BYTES + 1)
                if len(body) > MAX_DECOMPRESSED_BYTES:
                    return False, "Decompressed payload exceeds 36 MB limit.", {}
            except Exception:
                return False, "Invalid gzip payload.", {}
        try:
            obj = json.loads(body.decode("utf-8"))
        except Exception:
            return False, "Malformed UTF-8 JSON body.", {}
        if not isinstance(obj, dict):
            return False, "JSON root must be an object.", {}
        return True, "", obj

    def do_OPTIONS(self) -> None:
        self.send_response(HTTPStatus.NO_CONTENT)
        self._send_cors_headers()
        self.end_headers()

    def do_HEAD(self) -> None:
        self.send_response(HTTPStatus.OK)
        self._send_cors_headers()
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.end_headers()

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        path = parsed.path.rstrip("/") or "/"

        if path in ("/", "/v1/status"):
            with DB_LOCK, get_db() as conn:
                active, queued = current_queue_counts(conn)
            self._send_json(
                HTTPStatus.OK,
                {
                    "service": "WALRUS iGEM Cluster Offload",
                    "online": HYPHY3_BIN.is_file(),
                    "cluster": "silverback (ARM64 Neoverse Slurm)",
                    "threads_per_job": THREADS_PER_JOB,
                    "max_concurrent_jobs": MAX_CONCURRENT_SLURM,
                    "active_jobs": active,
                    "queued_jobs": queued,
                    "telemetry": True,
                    "limits": {
                        "max_taxa": MAX_TAXA,
                        "max_seq_nt": MAX_SEQ_UNGAPPED_NT,
                        "max_upload_mb": MAX_UPLOAD_BYTES // (1024 * 1024),
                        "walltime_minutes": 20,
                        "retention_hours": 48,
                    },
                },
            )
            return

        if path == "/v1/stats":
            with DB_LOCK, get_db() as conn:
                summary = build_telemetry_stats(conn)
            self._send_json(HTTPStatus.OK, summary)
            return

        if path == "/v1/challenge":
            subnet_h = client_subnet_hash(self.headers, self.client_address[0])
            with DB_LOCK, get_db() as conn:
                active, queued = current_queue_counts(conn)
            ch = issue_pow_challenge(subnet_h, active, queued)
            ch["active_jobs"] = active
            ch["queued_jobs"] = queued
            self._send_json(HTTPStatus.OK, ch)
            return

        m_job = re.fullmatch(r"/v1/jobs/(wlr_[A-Za-z0-9_-]{16,32})(?:/(result\.json|result\.html|result\.msa\.fa))?", path)
        if m_job:
            jid, artifact = m_job.group(1), m_job.group(2)
            with DB_LOCK, get_db() as conn:
                row = conn.execute("SELECT * FROM jobs WHERE job_id = ?", (jid,)).fetchone()
            if not row:
                self._send_json(HTTPStatus.NOT_FOUND, {"error": "Job ID not found or expired (48h retention)."})
                return

            job_dir = JOBS_DIR / jid
            if artifact == "result.json":
                self._send_file(HTTPStatus.OK, job_dir / "result.walrus.json", "application/json; charset=utf-8")
                return
            if artifact == "result.html":
                self._send_file(HTTPStatus.OK, job_dir / "result.walrus.html", "text/html; charset=utf-8")
                return
            if artifact == "result.msa.fa":
                self._send_file(HTTPStatus.OK, job_dir / "result.msa.fa", "text/plain; charset=utf-8")
                return

            # Status query
            status = row["status"]
            log_path = job_dir / "walrus.log"
            log_txt = ""
            if log_path.is_file():
                try:
                    with open(log_path, "rb") as fh:
                        fh.seek(0, os.SEEK_END)
                        sz = fh.tell()
                        fh.seek(max(0, sz - 8192))
                        log_txt = fh.read().decode("utf-8", errors="replace")
                except Exception:
                    log_txt = ""

            stage_label, pct, log_tail = parse_stage_progress(log_txt, status)
            now = time.time()
            started = float(row["started_at"] or 0.0)
            completed = float(row["completed_at"] or 0.0)
            elapsed = (
                round(completed - started, 1)
                if completed > 0 and started > 0
                else (round(now - started, 1) if started > 0 and status in ("running", "submitted") else 0.0)
            )

            qpos = 0
            if status == "queued":
                with DB_LOCK, get_db() as conn:
                    qrow = conn.execute(
                        "SELECT COUNT(*) AS c FROM jobs WHERE status = 'queued' AND created_at <= ?",
                        (row["created_at"],),
                    ).fetchone()
                    qpos = int(qrow["c"] or 1)

            self._send_json(
                HTTPStatus.OK,
                {
                    "job_id": jid,
                    "label": row["label"],
                    "num_taxa": row["num_taxa"],
                    "max_nt": row["max_nt"],
                    "options": json.loads(row["options_json"]),
                    "status": status,
                    "stage": stage_label,
                    "progress_pct": pct,
                    "queue_position": qpos,
                    "node": row["node"],
                    "threads": THREADS_PER_JOB,
                    "created_at": row["created_at"],
                    "elapsed_sec": elapsed,
                    "log_tail": log_tail,
                    "error": row["error"],
                    "has_json": (job_dir / "result.walrus.json").is_file(),
                    "has_html": (job_dir / "result.walrus.html").is_file(),
                    "has_msa": (job_dir / "result.msa.fa").is_file(),
                },
            )
            return

        self._send_json(HTTPStatus.NOT_FOUND, {"error": "Unknown endpoint."})

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        path = parsed.path.rstrip("/")

        if path == "/v1/telemetry":
            subnet_h = client_subnet_hash(self.headers, self.client_address[0])
            if not allow_telemetry_beacon(subnet_h):
                self.send_response(HTTPStatus.NO_CONTENT)
                self._send_cors_headers()
                self.end_headers()
                return

            ok, _, body = self._read_json_body(max_bytes=MAX_TELEMETRY_BYTES)
            if not ok:
                self.send_response(HTTPStatus.NO_CONTENT)
                self._send_cors_headers()
                self.end_headers()
                return

            ev = str(body.get("event", "")).strip()
            if ev not in ALLOWED_TELEMETRY_EVENTS:
                self.send_response(HTTPStatus.NO_CONTENT)
                self._send_cors_headers()
                self.end_headers()
                return

            now = time.time()
            utc_day = time.strftime("%Y-%m-%d", time.gmtime(now))
            visitor = daily_visitor_hash(subnet_h, utc_day)
            raw_cc = str(self.headers.get("CF-IPCountry", "")).strip().upper()
            country = raw_cc if (len(raw_cc) == 2 and raw_cc.isalpha() and raw_cc != "XX") else ""

            mode = re.sub(r"[^a-z0-9_-]+", "", str(body.get("mode", "")).lower())[:16]
            stage = re.sub(r"[^a-z0-9_-]+", "", str(body.get("stage", "")).lower())[:24]
            sample_id = (
                re.sub(r"[^A-Za-z0-9._-]+", "", str(body.get("sample_id", "")))[:64]
                if ev == "sample_load"
                else ""
            )

            def _clamp_int(val, hi: int) -> int:
                try:
                    return max(0, min(hi, int(round(float(val)))))
                except Exception:
                    return 0

            n_taxa = _clamp_int(body.get("n_taxa", 0), 100_000)
            n_sites = _clamp_int(body.get("n_sites", 0), 1_000_000)
            wall_ms = _clamp_int(body.get("wall_ms", 0), 86_400_000)
            threads = _clamp_int(body.get("threads", 0), 512)
            err_kind = re.sub(r"[^A-Za-z0-9_.:-]+", "", str(body.get("err_kind", "")))[:48]
            err_msg = scrub_error_message(str(body.get("err_msg", "")))

            try:
                with DB_LOCK, get_db() as conn:
                    conn.execute(
                        """
                        INSERT INTO telemetry_events (
                            ts, day, day_visitor, country, event, mode, stage,
                            sample_id, n_taxa, n_sites, wall_ms, threads, err_kind, err_msg
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            now,
                            utc_day,
                            visitor,
                            country,
                            ev,
                            mode,
                            stage,
                            sample_id,
                            n_taxa,
                            n_sites,
                            wall_ms,
                            threads,
                            err_kind,
                            err_msg,
                        ),
                    )
            except Exception:
                pass

            self.send_response(HTTPStatus.NO_CONTENT)
            self._send_cors_headers()
            self.end_headers()
            return

        m_cancel = re.fullmatch(r"/v1/jobs/(wlr_[A-Za-z0-9_-]{16,32})/cancel", path)
        if m_cancel:
            jid = m_cancel.group(1)
            ok, err, body = self._read_json_body()
            if not ok:
                self._send_json(HTTPStatus.BAD_REQUEST, {"error": err})
                return
            token = str(body.get("cancel_token", ""))
            with DB_LOCK, get_db() as conn:
                row = conn.execute("SELECT * FROM jobs WHERE job_id = ?", (jid,)).fetchone()
                if not row:
                    self._send_json(HTTPStatus.NOT_FOUND, {"error": "Job not found."})
                    return
                if not token or not hmac.compare_digest(token, row["cancel_token"]):
                    self._send_json(HTTPStatus.FORBIDDEN, {"error": "Invalid cancel token."})
                    return
                sid = row["slurm_id"]
                conn.execute(
                    "UPDATE jobs SET status = 'canceled', completed_at = ? WHERE job_id = ? AND status IN ('queued', 'submitted', 'running')",
                    (time.time(), jid),
                )
            if sid:
                subprocess.run(["/usr/bin/scancel", sid], capture_output=True, timeout=6, check=False)
            WAKE_EVENT.set()
            self._send_json(HTTPStatus.OK, {"job_id": jid, "status": "canceled"})
            return

        if path == "/v1/jobs":
            if not HYPHY3_BIN.is_file():
                self._send_json(HTTPStatus.SERVICE_UNAVAILABLE, {"error": "Cluster hyphy3 binary is not ready."})
                return

            ok, err, body = self._read_json_body()
            if not ok:
                self._send_json(HTTPStatus.BAD_REQUEST, {"error": err})
                return

            subnet_h = client_subnet_hash(self.headers, self.client_address[0])
            now = time.time()

            with DB_LOCK, get_db() as conn:
                _, queued_cnt = current_queue_counts(conn)
                if queued_cnt >= MAX_PENDING_QUEUE:
                    self._send_json(
                        HTTPStatus.TOO_MANY_REQUESTS,
                        {"error": f"Cluster queue is currently full ({MAX_PENDING_QUEUE} pending jobs). Please retry in a few minutes."},
                    )
                    return

                sub_active = conn.execute(
                    "SELECT COUNT(*) AS c FROM jobs WHERE subnet_hash = ? AND status IN ('queued', 'submitted', 'running')",
                    (subnet_h,),
                ).fetchone()["c"]
                if int(sub_active or 0) >= MAX_ACTIVE_PER_SUBNET:
                    self._send_json(
                        HTTPStatus.TOO_MANY_REQUESTS,
                        {"error": f"Limit of {MAX_ACTIVE_PER_SUBNET} concurrent active cluster jobs per network reached. Wait for your running job to finish or cancel it."},
                    )
                    return

                sub_hour = conn.execute(
                    "SELECT COUNT(*) AS c FROM jobs WHERE subnet_hash = ? AND created_at > ?",
                    (subnet_h, now - 3600.0),
                ).fetchone()["c"]
                if int(sub_hour or 0) >= MAX_HOURLY_PER_SUBNET:
                    self._send_json(
                        HTTPStatus.TOO_MANY_REQUESTS,
                        {"error": f"Hourly rate limit ({MAX_HOURLY_PER_SUBNET} cluster jobs/hour) reached."},
                    )
                    return

                pow_ok, pow_err = verify_and_consume_pow(conn, subnet_h, body.get("pow"))
                if not pow_ok:
                    self._send_json(HTTPStatus.FORBIDDEN, {"error": pow_err})
                    return

            raw_seq = str(body.get("sequences", ""))
            raw_tree = str(body.get("tree", ""))
            raw_label = re.sub(r"[^A-Za-z0-9._ -]+", "_", str(body.get("label", "dataset.fa")))[:96] or "dataset.fa"

            opt_ok, opt_err, clean_opts = sanitize_options(body.get("options", {}))
            if not opt_ok:
                self._send_json(HTTPStatus.BAD_REQUEST, {"error": opt_err})
                return

            seq_ok, seq_err, meta = validate_sequences_and_tree(raw_seq, raw_tree)
            if not seq_ok:
                self._send_json(HTTPStatus.BAD_REQUEST, {"error": seq_err})
                return

            # Content-addressed cache lookup
            h = hashlib.sha256()
            h.update(raw_seq.strip().encode("utf-8"))
            h.update(b"\x00")
            h.update(meta["clean_tree"].encode("utf-8"))
            h.update(b"\x00")
            h.update(json.dumps(clean_opts, sort_keys=True).encode("utf-8"))
            cache_key = h.hexdigest()

            with DB_LOCK, get_db() as conn:
                cached = conn.execute(
                    "SELECT job_id FROM jobs WHERE cache_key = ? AND status = 'completed' ORDER BY completed_at DESC LIMIT 1",
                    (cache_key,),
                ).fetchone()
                if cached:
                    cjid = cached["job_id"]
                    cdir = JOBS_DIR / cjid
                    if (cdir / "result.walrus.json").is_file():
                        self._send_json(
                            HTTPStatus.OK,
                            {
                                "job_id": cjid,
                                "cancel_token": "",
                                "status": "completed",
                                "cached": True,
                            },
                        )
                        return

            job_id = "wlr_" + secrets.token_urlsafe(15).replace("-", "A").replace("_", "B")
            cancel_token = secrets.token_urlsafe(18)
            job_dir = JOBS_DIR / job_id
            job_dir.mkdir(parents=True, exist_ok=False)
            os.chmod(job_dir, 0o700)

            (job_dir / "input.seq").write_text(raw_seq, encoding="utf-8")
            if meta["clean_tree"]:
                (job_dir / "input.nwk").write_text(meta["clean_tree"] + "\n", encoding="utf-8")

            spec = {
                "job_id": job_id,
                "hyphy3_bin": str(HYPHY3_BIN),
                "threads": THREADS_PER_JOB,
                "timeout_sec": RUNNER_TIMEOUT_SEC,
                "options": clean_opts,
            }
            (job_dir / "spec.json").write_text(json.dumps(spec, indent=2), encoding="utf-8")

            with DB_LOCK, get_db() as conn:
                conn.execute(
                    """
                    INSERT INTO jobs (
                        job_id, cancel_token, cache_key, subnet_hash, label,
                        num_taxa, max_nt, options_json, status, created_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, 'queued', ?)
                    """,
                    (
                        job_id,
                        cancel_token,
                        cache_key,
                        subnet_h,
                        raw_label,
                        meta["num_taxa"],
                        meta["max_nt"],
                        json.dumps(clean_opts),
                        now,
                    ),
                )

            WAKE_EVENT.set()
            self._send_json(
                HTTPStatus.ACCEPTED,
                {
                    "job_id": job_id,
                    "cancel_token": cancel_token,
                    "status": "queued",
                    "cached": False,
                    "threads": THREADS_PER_JOB,
                },
            )
            return

        self._send_json(HTTPStatus.NOT_FOUND, {"error": "Unknown endpoint."})


class ReusableThreadingHTTPServer(ThreadingHTTPServer):
    allow_reuse_address = True


def main() -> None:
    init_db()
    sched = threading.Thread(target=scheduler_loop, name="walrus-scheduler", daemon=True)
    sched.start()
    WAKE_EVENT.set()

    httpd = ReusableThreadingHTTPServer((BIND_HOST, BIND_PORT), WalrusGatewayHandler)
    print(f"WALRUS Cluster Gateway listening on http://{BIND_HOST}:{BIND_PORT} (hyphy3={HYPHY3_BIN})", flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
