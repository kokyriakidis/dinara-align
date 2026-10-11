# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""What every benchmark shares: the rivals' pins and builds, the CPU every tool is built for, the runners'
`seq` protocol, and the rivals' kept results.

Every rival is cloned at a pinned commit into `benchmarks/.cache/` and built there; nothing of theirs is
vendored. Every tool is built for the same CPU (see `set_cpu`): by default the host's own, each compiler
told so, Rust by `-C target-cpu`, C and C++ by `-march` or `-mcpu`, and Mojo by `--target-cpu`, so no tool
runs a narrower instruction set than another. Each table names the CPU it was built for and the processor
it ran on.

A pinned rival's times do not change, so each is timed once on a machine and its results kept: by
`measure` for the runners that speak the `seq` protocol, by `kept_rows` for the in-process drivers. A
later run replays them and times dinara-align alone, until a pin, a runner's source, a workload or the
CPU changes, or a fresh run is asked for.
"""

import hashlib
import json
import os
import platform
import random
import resource
import shutil
import signal
import subprocess
import sys
import threading
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CACHE = HERE / ".cache"
DATA = CACHE / "data"
RESULTS = CACHE / "results"

RIVALS = {
    "pa-bench": ("https://github.com/pairwise-alignment/pa-bench", "af7a50d0c3aa9518a141c6e1511c2126ab416848"),
    "astar-pairwise-aligner": (
        "https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner",
        "bf2e14e0cbc3a9a03600dcda0641d7f89e401e63",
    ),
    # v2.3.6: the modes' and batches' drivers, and its regression set (see `scripts/test_wfa.py`).
    "WFA2-lib": ("https://github.com/smarco/WFA2-lib", "bcf473a6561fa297934a80eaa7c04b4ee412360c"),
    "SSW": (
        "https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library",
        "a66636b79ef36ac178122437053be0d8ef345271",
    ),
    "parasail": ("https://github.com/jeffdaily/parasail", "fb985ee4f2302c72c87d7d56443712b2d920b106"),
    # Edlib v1.2.7.
    "edlib": ("https://github.com/Martinsos/edlib", "ec2310eda1841ab48c14cd3d866778a4f5eb1491"),
    "ksw2": ("https://github.com/lh3/ksw2", "289609bd9e5381a13b16239d0a7703f1ff03f9ca"),
    # Bound-and-align (Doblas et al., Bioinformatics 2025), x86-64 only.
    "QuickEd": ("https://github.com/maxdoblas/QuickEd", "2b1d44b8204a17fa0752fbe0b69b081c86a56e10"),
    # Striped bit-vector and difference-recurrence kernels (Shao and Ruan, Bioinformatics 2024), x86-64 only.
    "bsalign": ("https://github.com/ruanjue/bsalign", "1707aabea62f03eca77bdc059c5a02ccb7f9c0fa"),
}
"""Each rival's repository and the commit its numbers were taken at."""

DNA = "ACGT"

# region Inputs


def mutate(sequence: str, rate: float, alphabet: str, rng: random.Random) -> str:
    """Substitutions, deletions and insertions in equal thirds, at `rate` edits per letter."""
    out = []
    for letter in sequence:
        roll = rng.random()
        if roll < rate / 3:
            out.append(rng.choice([other for other in alphabet if other != letter]))
        elif roll < 2 * rate / 3:
            continue
        elif roll < rate:
            out.extend((letter, rng.choice(alphabet)))
        else:
            out.append(letter)
    return "".join(out)


# endregion Inputs

# region Builds


def fetch(name: str) -> Path:
    """Clones a rival at its pinned commit, or moves an existing clone there."""
    url, commit = RIVALS[name]
    destination = CACHE / name
    if not destination.exists():
        subprocess.run(["git", "clone", "--quiet", "--filter=blob:none", url, str(destination)], check=True)
    if subprocess.run(["git", "-C", str(destination), "cat-file", "-e", commit], capture_output=True).returncode:
        subprocess.run(["git", "-C", str(destination), "fetch", "--quiet", "origin"], check=True)
    subprocess.run(["git", "-C", str(destination), "checkout", "--quiet", commit], check=True)
    return destination


def nightly_environment(install: bool) -> dict | None:
    """An environment whose `cargo` honours A*PA's pinned nightly, or `None` when there is none.

    A system `rustup` serves as is. Without one, `--install-rust` puts a private copy under the
    cache, leaving the shell profile and any system Rust untouched.
    """
    if shutil.which("rustup"):
        return dict(os.environ)
    private = CACHE / "rust"
    environment = dict(os.environ, RUSTUP_HOME=str(private / "rustup"), CARGO_HOME=str(private / "cargo"))
    environment["PATH"] = f"{private / 'cargo' / 'bin'}{os.pathsep}{environment['PATH']}"
    if (private / "cargo" / "bin" / "rustup").exists():
        return environment
    if not install:
        return None
    installer = subprocess.run(
        ["curl", "--proto", "=https", "--tlsv1.2", "-sSf", "https://sh.rustup.rs"], check=True, capture_output=True
    )
    subprocess.run(
        ["sh", "-s", "--", "-y", "--no-modify-path", "--default-toolchain", "none", "--profile", "minimal"],
        input=installer.stdout,
        env=environment,
        check=True,
        capture_output=True,
    )
    return environment


CPU = "native"
"""The CPU every tool is built for: `native`, the host's own, or a name all three compilers know, such
as `x86-64-v3` or `generic`; see `set_cpu`."""


def set_cpu(name: str) -> None:
    """Builds every tool after this for the CPU `name`."""
    global CPU
    CPU = name


def arm() -> bool:
    """Whether the host is an ARM machine."""
    return platform.machine().lower() in ("arm64", "aarch64")


def x86() -> bool:
    """Whether the host is an x86-64 machine, where the rivals with x86 vector code alone build."""
    return platform.machine().lower() in ("x86_64", "amd64")


def c_cpu_flag() -> str:
    """The C and C++ compilers' flag for `CPU`: `-mcpu` on ARM, `-march` elsewhere."""
    return f"-mcpu={CPU}" if arm() else f"-march={CPU}"


def build_environment(environment: dict) -> dict:
    """`environment` with every compiler a Rust build runs told to build for `CPU`: rustc through
    `RUSTFLAGS`, a build script's C and C++ through `CFLAGS` and `CXXFLAGS`, each added to what is set."""
    out = dict(environment)
    out["RUSTFLAGS"] = f"{out.get('RUSTFLAGS', '')} -C target-cpu={CPU}".strip()
    for name in ("CFLAGS", "CXXFLAGS"):
        out[name] = f"{out.get(name, '')} {c_cpu_flag()}".strip()
    return out


def mojo_cpu_args() -> list:
    """Mojo's flags for `CPU`: none for `native`, its default, else `--target-cpu`."""
    return [] if CPU == "native" else ["--target-cpu", CPU]


def processor() -> str:
    """The host's processor, as the operating system names it."""
    try:
        if sys.platform == "darwin":
            return subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, check=True
            ).stdout.strip()
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except (OSError, subprocess.CalledProcessError):
        pass
    return platform.processor() or platform.machine()


def build_note() -> str:
    """The line each table opens with: what every tool was built for and what it ran on."""
    return (
        f"Every tool built for `{CPU}` (Rust `-C target-cpu={CPU}`, C and C++ `{c_cpu_flag()}`, Mojo "
        f"{'its default, the host' if CPU == 'native' else '--target-cpu ' + CPU}), run on {processor()}.\n\n"
    )


def cargo_runner(crate: str, environment: dict) -> Path:
    """Builds one Rust runner in release mode and returns its binary."""
    target = CACHE / f"target-{crate}"
    subprocess.run(
        ["cargo", "build", "--release", "--quiet"],
        cwd=HERE / crate,
        env=dict(build_environment(environment), CARGO_TARGET_DIR=str(target)),
        check=True,
    )
    return target / "release" / f"{crate}-runner"


def mojo_binary(source: Path, name: str) -> Path:
    """Builds a dinara-align runner, `source`, against this checkout, unless it already is.

    A build takes seconds, so the binary is kept while it is newer than every source and the toolchain's
    lock, and was built with the same command, which a stamp beside it records.
    """
    binary = CACHE / "bin" / name
    binary.parent.mkdir(parents=True, exist_ok=True)
    accelerator = os.environ.get("MOJO_ACCELERATOR", "").split()
    command = ["mojo", "build", "-I", str(ROOT), str(source), "-o", str(binary), *mojo_cpu_args(), *accelerator]
    stamp = binary.with_suffix(".command")
    sources = [source, ROOT / "pixi.lock", *(ROOT / "dinara_align").rglob("*.mojo")]
    if (
        binary.exists()
        and stamp.exists()
        and stamp.read_text() == " ".join(command)
        and all(item.stat().st_mtime < binary.stat().st_mtime for item in sources)
    ):
        return binary
    subprocess.run(command, check=True)
    stamp.write_text(" ".join(command))
    return binary


def mojo_runner() -> Path:
    """The `seq` runner, `ours.mojo` (see `run_tool`)."""
    return mojo_binary(HERE / "ours.mojo", "dinara-align-runner")


def c_driver(driver: Path, name: str, extra: list = []) -> Path:
    """A C++ driver, `driver`, over Edlib, WFA2-lib, parasail and SSW, and on x86-64 KSW2, built as `name` for
    `CPU` with `extra` flags besides: the modes' rivals and the batches'."""
    edlib = fetch("edlib")
    wfa = fetch("WFA2-lib")
    parasail = fetch("parasail")
    ssw = fetch("SSW")
    # Each CPU its own parasail build: it dispatches its kernels at run time either way, but its own C follows.
    parasail_build = parasail / f"build-{CPU}"
    if not (parasail_build / "libparasail.a").exists():
        parasail_build.mkdir(exist_ok=True)
        subprocess.run(
            ["cmake", "-DCMAKE_POLICY_VERSION_MINIMUM=3.5", "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_SHARED_LIBS=OFF",
             f"-DCMAKE_C_FLAGS={c_cpu_flag()}", ".."],
            cwd=parasail_build, check=True, capture_output=True,
        )
        subprocess.run(["make", "-j4", "parasail"], cwd=parasail_build, check=True, capture_output=True)
    # WFA2-lib's Makefile builds for the host's own CPU; given the chosen one instead, kept per CPU.
    wfa_library = wfa / "lib" / f"libwfa-{CPU}.a"
    if not wfa_library.exists():
        subprocess.run(["make", "clean"], cwd=wfa, check=True, capture_output=True)
        subprocess.run(["make", "setup", "lib_wfa", f"CC_FLAGS=-Wall -fPIE -O3 {c_cpu_flag()}"], cwd=wfa, check=True, capture_output=True)
        (wfa / "lib" / "libwfa.a").rename(wfa_library)
    build = CACHE / "drivers"
    build.mkdir(parents=True, exist_ok=True)
    binary = build / f"{name}-{CPU}"
    flags = ["-O3", c_cpu_flag(), *extra]
    sources = [str(ssw / "src" / "ssw.c")]
    defines = []
    if x86():
        ksw2 = fetch("ksw2")
        sources += [str(ksw2 / "ksw2_extz2_sse.c"), str(ksw2 / "ksw2_extd2_sse.c"), str(ksw2 / "kalloc.c")]
        defines += ["-DWITH_KSW2", f"-I{ksw2}"]
    objects = []
    for source in sources:
        target = build / (Path(source).stem + f"-{CPU}.o")
        subprocess.run(["cc", "-c", *flags, *defines, f"-I{ssw / 'src'}", source, "-o", str(target)], check=True)
        objects.append(str(target))
    subprocess.run(
        [
            "c++", "-std=c++17", *flags, *defines,
            str(driver), str(edlib / "edlib" / "src" / "edlib.cpp"), *objects,
            f"-I{edlib / 'edlib' / 'include'}", f"-I{wfa}", f"-I{ssw / 'src'}", f"-I{parasail}", f"-I{parasail_build}",
            str(wfa_library), str(parasail_build / "libparasail.a"),
            "-lz", "-lm", "-lpthread", "-o", str(binary),
        ],
        check=True,
    )
    return binary


# endregion Builds

# region Kept results


KEPT_ROWS = CACHE / "rival-results.json"
"""The in-process drivers' rows from earlier runs (see `kept_rows`)."""
REMEASURE = False
"""Whether the rivals run again even where results for the same inputs are kept."""


def fingerprint(*paths: Path) -> str:
    """A hash of every file at or under `paths`, by name and content: a changed workload, driver or lock
    is a different fingerprint, and a rebuilt binary or a touched file is not."""
    digest = hashlib.sha1()
    for path in paths:
        files = [path] if path.is_file() else sorted(item for item in path.rglob("*") if item.is_file())
        for file in files:
            digest.update(file.name.encode())
            digest.update(file.read_bytes())
    return digest.hexdigest()[:16]


def kept_rows(tool: str, pins: list, inputs: list, measured):
    """`measured()`'s rows for a rival's in-process driver, measured once and kept.

    A rival is pinned, so the same commits and driver built for the same CPU give the same times on the
    same workload files: its rows are kept under all of those, `pins` the rival repositories it builds
    from and `inputs` the files its driver and its workloads are, and replayed until any of them
    changes, or `REMEASURE` asks for a fresh run. dinara-align, whose code changes, runs every time.
    """
    commits = "-".join(RIVALS[name][1][:12] for name in pins)
    key = "\t".join([tool, commits, fingerprint(*inputs), CPU])
    kept = json.loads(KEPT_ROWS.read_text()) if KEPT_ROWS.exists() else {}
    if key in kept and not REMEASURE:
        print(f"{tool}: kept from an earlier run", file=sys.stderr, flush=True)
        return kept[key]
    rows = measured()
    kept[key] = rows
    KEPT_ROWS.parent.mkdir(parents=True, exist_ok=True)
    KEPT_ROWS.write_text(json.dumps(kept))
    return rows


# endregion Kept results

# region The seq protocol


GRACE = 1.0
"""The seconds past its budget a runner may take, starting up included, before it is stopped."""

MEMORY_CAP = 8 << 30
"""Bytes a capped tool may take before it is stopped as its budget would stop it: WFA's keep-every-front
mode grows with the square of the distance, BSAlign's and KSW2's tracebacks with the matrix, past any
machine's memory on the longest divergent pairs, and faster than a budget of seconds stops them."""

CAPPED = ("wfa", "bsalign", "ksw2")
"""The tools run under `MEMORY_CAP`, by their names' first word."""

KEPT = CACHE / "pa-bench-rivals.json"
"""The `seq` runners' results from earlier runs (see `measure`)."""

TIMER = Path("/usr/bin/time")
"""GNU time, which on Linux reports a runner's own peak memory (see `run_tool`)."""


DONE = "done"
"""A run that aligned every pair."""
BUDGET = "budget"
"""A run its time budget stopped: QuickEd's "t.o." where it aligned nothing."""
MEMORY = "memory"
"""A run `MEMORY_CAP` stopped: QuickEd's "o.o.m."."""
FAILED = "failed"
"""A run that ended on an error of its own."""


def run_tool(binary: Path, tool: str, path: Path, budget: float) -> tuple[list[tuple[float, int, int]], int, str]:
    """One tool over one `.seq` file, `<runner> seq <tool> <budget> <file>`: each pair's time, cost and
    growth of peak memory it finished, its peak resident bytes, and how the run ended (`DONE`, `BUDGET`,
    `MEMORY` or `FAILED`).

    The runner prints a row, `tool, file, seconds, cost[, growth]`, as each pair finishes and checks its
    budget only between pairs, so one still busy a grace past the budget is stopped, keeping the rows it
    printed. A pair's growth is what aligning it added to the most the process had held before it, A*PA2's
    memory measure, so the input and the runtime count for nothing.

    The peak memory is the runner's own. On macOS `wait4` reports it, even when stopped, where
    `subprocess` would discard it. On Linux a forked process starts from its parent's peak, so a runner
    forked from this harness, which holds whole datasets, would report the harness's: there GNU time, a
    small process, forks the runner and reports its peak, and a runner stopped mid-pair has none to
    report, -1. The runner starts a session of its own, so stopping it stops GNU time and the runner
    together. A tool in `CAPPED` stops at `MEMORY_CAP`, its pairs so far kept: on Linux the kernel
    refuses it the memory, so it fails before its budget runs out, and on macOS, which does not enforce
    that limit, its resident size is watched.
    """
    timed = sys.platform.startswith("linux") and TIMER.exists()
    command = [str(binary), "seq", tool, str(budget), str(path)]
    if timed:
        command = [str(TIMER), "-f", "peak %M", *command]
    capped = tool.replace(":", "-").split("-")[0] in CAPPED

    def limit() -> None:
        """Caps the runner's address space at `MEMORY_CAP`, run in the child before it starts."""
        resource.setrlimit(resource.RLIMIT_AS, (MEMORY_CAP, MEMORY_CAP))

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE if timed else subprocess.DEVNULL,
        text=True,
        start_new_session=True,
        preexec_fn=limit if capped and sys.platform.startswith("linux") else None,
    )
    stopped = threading.Event()
    exhausted = threading.Event()

    def stop() -> None:
        """Kills the runner's whole session, GNU time with it, and marks the run stopped."""
        stopped.set()
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass

    timer = threading.Timer(budget + GRACE, stop)
    timer.start()
    if capped and sys.platform == "darwin":

        def watch() -> None:
            """Stops the runner once its resident size passes `MEMORY_CAP`, polling every 50 ms."""
            while process.poll() is None and not stopped.is_set():
                resident = subprocess.run(["ps", "-o", "rss=", "-p", str(process.pid)], capture_output=True, text=True)
                if resident.stdout.strip() and int(resident.stdout) * 1024 > MEMORY_CAP:
                    exhausted.set()
                    stop()
                    return
                stopped.wait(0.05)

        threading.Thread(target=watch, daemon=True).start()
    output = process.stdout.read()
    process.stdout.close()
    errors = process.stderr.read() if timed else ""
    if timed:
        process.stderr.close()
    _, status, usage = os.wait4(process.pid, 0)
    timer.cancel()
    process.returncode = os.waitstatus_to_exitcode(status)
    if timed:
        # GNU time counts kilobytes.
        reported = [line for line in errors.splitlines() if line.startswith("peak ")]
        peak = int(reported[-1].split()[1]) * 1024 if reported else -1
    else:
        # macOS counts bytes.
        peak = usage.ru_maxrss
    rows = []
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) in (4, 5) and fields[0] == tool:
            rows.append((float(fields[2]), int(fields[3]), int(fields[4]) if len(fields) == 5 else -1))
    if exhausted.is_set():
        status = MEMORY
    elif stopped.is_set():
        status = BUDGET
    elif process.returncode != 0:
        # Under the cap, a runner the kernel refused memory dies of it; anything else is its own error.
        status = MEMORY if capped and sys.platform.startswith("linux") else FAILED
    elif len(rows) < path.read_bytes().count(b"\n") // 2:
        # A runner stops starting pairs once its budget is spent, and ends as it would at the last.
        status = BUDGET
    else:
        status = DONE
    return rows, peak, status


RUNNER_SOURCES = {
    "target-quicked": ("QuickEd", ["quicked/runner.c"]),
    "target-bsalign": ("bsalign", ["bsalign/runner.c"]),
    "target-astarpa": ("astar-pairwise-aligner", ["astarpa/src/main.rs", "astarpa/Cargo.lock"]),
    "target-pa-wrapper": ("pa-bench", ["pa-wrapper/src/main.rs", "pa-wrapper/Cargo.lock"]),
}
"""Each rival runner's build directory, the rival it builds from and our runner's sources."""


def identity(binary: Path) -> str:
    """What a rival's results depend on: its pinned commit, our runner's source and the CPU it was built
    for, not the binary's age. A rebuild of unchanged sources keeps every result; moving a pin, editing the
    runner, or building for another CPU starts afresh."""
    build = next(part for part in binary.parents if part.name in RUNNER_SOURCES)
    rival, sources = RUNNER_SOURCES[build.name]
    source = b"".join((HERE / item).read_bytes() for item in sources)
    return f"{RIVALS[rival][1][:12]}-{hashlib.sha1(source).hexdigest()[:12]}-{CPU}"


def replay(rows: list, status: str, budget: float) -> tuple[list[tuple], str]:
    """What a run with a smaller `budget` would have kept of a longer run's pairs: the runner starts a
    pair while what it has spent is under the budget, and a pair still running a grace past the budget
    is stopped."""
    kept, spent = [], 0.0
    for row in rows:
        seconds = row[0]
        if spent >= budget:
            return kept, BUDGET
        if spent + seconds > budget + GRACE:
            return kept, BUDGET
        kept.append(tuple(row))
        spent += seconds
    return kept, status


def measure(binary: Path, tool: str, path: Path, budget: float, kept: dict | None) -> tuple[list[tuple], int, str]:
    """`run_tool`, or for a rival what an earlier run of it on the same file already says.

    The rivals are pinned, so each one is run on a file once, at the largest budget anyone has asked for,
    and kept for good under the file's content and the rival's build (see `identity`): a smaller budget
    replays that run (see `replay`), and a larger one reuses it if it aligned every pair, running again
    only if it did not. A rival's peak memory is that of
    the run kept. dinara-align, `kept` None, runs every time.
    """
    if kept is None:
        return run_tool(binary, tool, path, budget)
    # By the sample's content, not its name alone: a dataset built again under the same name, as the GIAB
    # pairs were from a later assembly, would otherwise replay a rival's rows from the pairs it replaced.
    key = "\t".join([tool, path.name, fingerprint(path), identity(binary)])
    entry = kept.get(key)
    pairs = path.read_bytes().count(b"\n") // 2
    if entry is not None and not REMEASURE:
        # Results kept before a run's ending was told apart know only whether it was stopped.
        status = entry.get("status", BUDGET if entry.get("stopped") else DONE)
        whole = len(entry["rows"]) == pairs and status == DONE
        if budget <= entry["budget"] or whole or status == MEMORY:
            rows, status = replay([tuple(row) for row in entry["rows"]], status, budget)
            return rows, entry["peak"], status
    rows, peak, status = run_tool(binary, tool, path, budget)
    kept[key] = {"budget": budget, "rows": rows, "peak": peak, "status": status, "stopped": status != DONE}
    return rows, peak, status


def load_kept() -> dict:
    """The `seq` runners' kept results."""
    return json.loads(KEPT.read_text()) if KEPT.exists() else {}


# endregion The seq protocol

# region Cells


def duration(seconds: float) -> str:
    """A time in the unit that keeps three significant figures readable."""
    if seconds < 1e-3:
        return f"{seconds * 1e6:.0f} µs"
    if seconds < 0.9995:
        return f"{seconds * 1e3:.3g} ms"
    return f"{seconds:.3g} s"


def megabytes(peak: int, status: str) -> str:
    """A peak's table cell in megabytes, a lower bound for a run stopped, QuickEd's "o.o.m." for one the
    memory cap stopped, a dash for none reported."""
    if status == MEMORY:
        return "o.o.m."
    if peak < 0:
        return "—"
    return f"{'≥ ' if status != DONE else ''}{peak / 2**20:.0f} MB"


def growth(rows: list) -> str:
    """The median and largest growth of peak memory over a tool's pairs, A*PA2's memory table: what
    aligning a pair added to the most the process had held before it, so the runtime and the input count
    for nothing, and a pair within an earlier pair's peak counts zero."""
    grown = sorted(row[2] for row in rows if len(row) > 2 and row[2] >= 0)
    if not grown:
        return "—"

    def size(value: int) -> str:
        """`value` bytes in megabytes, one decimal below ten."""
        return f"{value / 2**20:.1f}" if value < 10 * 2**20 else f"{value / 2**20:.0f}"

    return f"{size(grown[len(grown) // 2])} / {size(grown[-1])} MB"


# endregion Cells
