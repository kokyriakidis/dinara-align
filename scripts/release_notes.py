#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Prints the section of CHANGELOG.md for one release, the notes of its GitHub release.

    python3 scripts/release_notes.py v0.2.0
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    """Prints the notes under the version's `##` heading, the first argument with or without its `v`."""
    version = sys.argv[1].removeprefix("v")
    lines = (ROOT / "CHANGELOG.md").read_text().splitlines()
    out, inside = [], False
    for line in lines:
        if line.startswith("## "):
            if inside:
                break
            inside = line.removeprefix("## ").split()[0].strip("[]") == version
            continue
        if inside:
            out.append(line)
    if not out:
        sys.exit(f"no section for {version} in CHANGELOG.md")
    print("\n".join(out).strip())


if __name__ == "__main__":
    main()
