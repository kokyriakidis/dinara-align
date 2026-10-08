#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Times dinara-align against other exact global DNA aligners, and checks that every answer agrees.

    pixi run bench                    # 1k and 10k workloads, about a minute
    pixi run bench --full             # adds the 100k DNA pairs, about five minutes
    pixi run bench --install-rust     # also installs the nightly A*PA needs, privately

Every rival is cloned at a pinned commit into `benchmarks/.cache/` and built there; nothing of
theirs is vendored. Inputs are generated from a fixed seed, so two runs time the same pairs. Each
runner prints `tool, workload, task, device, seconds, answer` rows, and the answers within a
workload must all be equal: a speed table over disagreeing aligners would compare different
questions, so a disagreement fails the run.

Every runner times each measurement the same way, warm and in-process: a call shorter than a tenth
of a second is repeated in twenty batches and the fastest batch's average kept, a longer one is
timed once. So one run of each runner is enough.

Every tool is built for the same CPU, `--cpu` (see `set_cpu`): by default the host's own, each
compiler told so, Rust by `-C target-cpu`, C and C++ by `-march` or `-mcpu`, and Mojo by
`--target-cpu`, so no tool runs a narrower instruction set than another. A*PA's kernels take their
SIMD width at compile time, so a build for the baseline CPU would run its 256-bit vectors as
pairs of 128-bit ones. Each table names the CPU it was built for and the processor it ran on.
"""

import argparse
import os
import platform
import random
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CACHE = HERE / ".cache"
DATA = CACHE / "data"
RESULTS = CACHE / "results"

RIVALS = {
    "hyalite": ("https://github.com/Psy-Fer/hyalite", "0189bcbfaf9e2fa57c7fb07ad7356e02c1d259f1"),
    "pa-bench": ("https://github.com/pairwise-alignment/pa-bench", "af7a50d0c3aa9518a141c6e1511c2126ab416848"),
    "astar-pairwise-aligner": (
        "https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner",
        "bf2e14e0cbc3a9a03600dcda0641d7f89e401e63",
    ),
    # v2.3.6, for its regression set alone (see `scripts/test_wfa.py`).
    "WFA2-lib": ("https://github.com/smarco/WFA2-lib", "bcf473a6561fa297934a80eaa7c04b4ee412360c"),
    # The local and overlap aligners of `local_bench.py`.
    "SSW": (
        "https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library",
        "a66636b79ef36ac178122437053be0d8ef345271",
    ),
    "parasail": ("https://github.com/jeffdaily/parasail", "fb985ee4f2302c72c87d7d56443712b2d920b106"),
    "abPOA": ("https://github.com/yangao07/abPOA", "2e095ba7f8de6bc62aaaa9a684cda781b8098ee0"),
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


def generate(full: bool) -> None:
    """Writes the four workload files, one `name, first, second` row per pair."""
    DATA.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)

    def random_sequence(length: int) -> str:
        """`length` random bases."""
        return "".join(rng.choice(DNA) for _ in range(length))

    def homologous_pairs(path: str, name: str, count: int, length: int, rate: float) -> None:
        """Writes `count` pairs of `length` bases to `path`, each second its first mutated at `rate`."""
        with open(DATA / path, "w") as out:
            for _ in range(count):
                first = random_sequence(length)
                out.write(f"{name}\t{first}\t{mutate(first, rate, DNA, rng)}\n")

    # Batches scored with the DNA default, sized like short-read and long-read work.
    homologous_pairs("dna_reads.tsv", "reads-150bp", 10_000, 150, 0.02)
    homologous_pairs("dna_kilobase.tsv", "reads-1kbp", 1_000, 1_000, 0.10)

    lengths = [(1_000, "1k"), (10_000, "10k")] + ([(100_000, "100k")] if full else [])
    # One divergence suffices for the affine pairs: a full sweep's time does not depend on it.
    with open(DATA / "dna_affine.tsv", "w") as affine:
        for length, label in lengths:
            first = random_sequence(length)
            affine.write(f"affine-{label}\t{first}\t{mutate(first, 0.05, DNA, rng)}\n")
    # A*PA's time does depend on it, so edit distance sweeps three.
    with open(DATA / "dna_edit.tsv", "w") as edit:
        for length, label in lengths:
            for rate in (0.01, 0.05, 0.15):
                first = random_sequence(length)
                edit.write(f"edit-{label}-{round(rate * 100)}%\t{first}\t{mutate(first, rate, DNA, rng)}\n")


# endregion Inputs

# region Rivals


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
        f"{'its default, the host' if CPU == 'native' else '--target-cpu ' + CPU}), run on {processor()}, "
        "inputs in memory, timed warm in-process.\n\n"
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


def mojo_runner() -> Path:
    """Builds the dinara-align runner against this checkout's sources, unless it already is.

    A build takes seconds, so the binary is kept while it is newer than every source and the toolchain's
    lock, and was built with the same command, which a stamp beside it records.
    """
    binary = CACHE / "bin" / "dinara-align-runner"
    binary.parent.mkdir(parents=True, exist_ok=True)
    accelerator = os.environ.get("MOJO_ACCELERATOR", "").split()
    command = ["mojo", "build", "-I", str(ROOT), str(HERE / "ours.mojo"), "-o", str(binary), *mojo_cpu_args(), *accelerator]
    stamp = binary.with_suffix(".command")
    sources = [HERE / "ours.mojo", ROOT / "pixi.lock", *(ROOT / "dinara_align").rglob("*.mojo")]
    if (
        binary.exists()
        and stamp.exists()
        and stamp.read_text() == " ".join(command)
        and all(source.stat().st_mtime < binary.stat().st_mtime for source in sources)
    ):
        return binary
    subprocess.run(command, check=True)
    stamp.write_text(" ".join(command))
    return binary


# endregion Rivals

# region Report


def run(binary: Path, label: str, repeat: int) -> list[list[str]]:
    """Runs one runner `repeat` times and keeps each measurement's fastest time.

    The minimum rather than the mean, because everything slowing a run down is noise from outside
    the aligner; an answer that changes between repetitions is kept too, so `report` sees it.
    """
    best: dict[tuple, list[str]] = {}
    for attempt in range(repeat):
        print(f"running {label} ({attempt + 1}/{repeat}) ...", file=sys.stderr, flush=True)
        outcome = subprocess.run([str(binary), str(DATA)], check=True, capture_output=True, text=True)
        for line in outcome.stdout.splitlines():
            row = line.split("\t")
            if len(row) != 6:
                continue
            key = (row[0], row[1], row[2], row[3], row[5])
            if key not in best or float(row[4]) < float(best[key][4]):
                best[key] = row
    return list(best.values())


def duration(seconds: float) -> str:
    """A time in the unit that keeps three significant figures readable."""
    if seconds < 1e-3:
        return f"{seconds * 1e6:.0f} µs"
    if seconds < 0.9995:
        return f"{seconds * 1e3:.3g} ms"
    return f"{seconds:.3g} s"


def report(rows: list[list[str]]) -> bool:
    """Writes the raw rows and a pivoted table, and answers whether every workload agreed."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(RESULTS / "results.tsv", "w") as raw:
        raw.write("tool\tworkload\ttask\tdevice\tseconds\tanswer\n")
        raw.writelines("\t".join(row) + "\n" for row in rows)

    answers = defaultdict(set)
    for tool, workload, _task, _device, _seconds, answer in rows:
        answers[workload].add(answer)
    disagreeing = sorted(workload for workload, seen in answers.items() if len(seen) > 1)

    columns = []
    for tool, _workload, _task, device, _seconds, _answer in rows:
        column = f"{tool} ({device})" if tool == "dinara-align" else tool
        if column not in columns:
            columns.append(column)
    cells = {}
    order = []
    for tool, workload, task, device, seconds, answer in rows:
        key = (workload, task)
        if key not in order:
            order.append(key)
        cells[key, f"{tool} ({device})" if tool == "dinara-align" else tool] = duration(float(seconds))

    lines = [
        "| workload | task | " + " | ".join(columns) + " | agree |",
        "| :-- | :-- | " + " | ".join("--:" for _ in columns) + " | :-: |",
    ]
    for workload, task in order:
        timings = " | ".join(cells.get(((workload, task), column), "—") for column in columns)
        lines.append(f"| {workload} | {task} | {timings} | {'✗' if workload in disagreeing else '✓'} |")
    table = build_note() + "\n".join(lines) + "\n"
    (RESULTS / "results.md").write_text(table)
    print(table)
    for workload in disagreeing:
        print(f"DISAGREEMENT on {workload}: {sorted(answers[workload])}", file=sys.stderr)
    return not disagreeing


# endregion Report


def main() -> None:
    """Generates the workloads, runs every tool it can build on them, and writes the report; exits nonzero
    when the tools' answers disagree."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full", action="store_true", help="include the 100k DNA pairs")
    parser.add_argument("--install-rust", action="store_true", help="install A*PA's nightly under the cache")
    parser.add_argument("--repeat", type=int, default=1, help="runs per tool, keeping the fastest (default 1)")
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    options = parser.parse_args()
    set_cpu(options.cpu)

    generate(options.full)
    # dinara-align runs first: it writes the DNA scoring every other runner reads.
    rows = run(mojo_runner(), "dinara-align", options.repeat)

    if shutil.which("cargo"):
        fetch("hyalite")
        rows += run(cargo_runner("hyalite", dict(os.environ)), "hyalite", options.repeat)
    else:
        print("skipping hyalite: no `cargo` on PATH", file=sys.stderr)

    nightly = nightly_environment(options.install_rust)
    if nightly is None:
        print("skipping A*PA: it needs `rustup` for its pinned nightly; pass --install-rust", file=sys.stderr)
    else:
        fetch("astar-pairwise-aligner")
        rows += run(cargo_runner("astarpa", nightly), "A*PA", options.repeat)

    if shutil.which("cargo"):
        fetch("pa-bench")
        rows += run(cargo_runner("pa-wrapper", dict(os.environ)), "Edlib and WFA2-lib", options.repeat)
    else:
        print("skipping Edlib and WFA2-lib: no `cargo` on PATH", file=sys.stderr)

    sys.exit(0 if report(rows) else 1)


if __name__ == "__main__":
    main()
