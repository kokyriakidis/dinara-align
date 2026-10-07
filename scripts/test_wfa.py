#!/usr/bin/env python3
"""Holds dinara-align to WFA2-lib's regression set: its pairs and its own results for each mode.

    pixi run test-wfa

Clones WFA2-lib at the commit pinned in `benchmarks/run.py` into the benchmarks' cache, the first time
only, builds `tests/wfa_utest.mojo` and runs it on the clone's `tests` directory. Nothing of WFA2-lib
is built: only its data files are read.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "benchmarks"))

from run import fetch  # noqa: E402


def main() -> None:
    wfa = fetch("WFA2-lib")
    binary = ROOT / "build" / "wfa_utest"
    binary.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["mojo", "build", "-I", str(ROOT), str(ROOT / "tests" / "wfa_utest.mojo"), "-o", str(binary)], check=True
    )
    sys.exit(subprocess.run([str(binary), str(wfa / "tests")]).returncode)


if __name__ == "__main__":
    main()
