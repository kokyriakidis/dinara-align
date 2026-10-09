#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Packs build/python/dinara_align (see `pixi run build-python`) into a wheel in build/dist, for `pip
install`: the package's Python, its extension and the Mojo runtime libraries beside it.

    pixi run build-wheel

The extension finds the interpreter's C API in the process that loads it rather than linking one
Python's, so one wheel serves every CPython 3: tagged `py3-none` and the platform, macOS 11 on Apple
silicon, or on Linux the glibc the runtime libraries need (see `scripts/bundle.sh`).
"""

import base64
import hashlib
import platform
import sys
import tomllib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "build" / "python" / "dinara_align"
DIST = ROOT / "build" / "dist"


def platform_tag() -> str:
    """The wheel's platform tag for this machine; exits for a platform with none known."""
    machine = platform.machine().lower()
    if sys.platform == "darwin":
        return "macosx_11_0_arm64"
    if machine in ("x86_64", "amd64"):
        return "manylinux_2_35_x86_64"
    if machine in ("aarch64", "arm64"):
        return "manylinux_2_35_aarch64"
    sys.exit(f"build_wheel: no wheel tag known for {sys.platform} {machine}")


def digest(data: bytes) -> str:
    """A `RECORD` entry's hash of `data`: SHA-256, URL-safe base64 without padding."""
    return "sha256=" + base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()


def main() -> None:
    """Writes the wheel, its metadata and `RECORD` with it, and prints its path."""
    if not (PACKAGE / "_dinara.so").exists():
        sys.exit("build_wheel: run `pixi run build-python` first")
    project = tomllib.loads((ROOT / "pixi.toml").read_text())["workspace"]
    version = project["version"]
    tag = f"py3-none-{platform_tag()}"
    info = f"dinara_align-{version}.dist-info"
    metadata = "\n".join(
        [
            "Metadata-Version: 2.1",
            "Name: dinara-align",
            f"Version: {version}",
            f"Summary: {project['description']}",
            "Home-page: https://github.com/kokyriakidis/dinara-align",
            f"Author: {project['authors'][0]}",
            f"License: {project['license']}",
            "Requires-Python: >=3.9",
            "",
        ]
    )
    wheel = f"Wheel-Version: 1.0\nGenerator: dinara-align build_wheel\nRoot-Is-Purelib: false\nTag: {tag}\n"
    DIST.mkdir(parents=True, exist_ok=True)
    target = DIST / f"dinara_align-{version}-{tag}.whl"
    records = []
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(PACKAGE.iterdir()):
            # The licenses go in the metadata below, where installers look for them.
            if path.name in ("__pycache__", "LICENSE", "NOTICE"):
                continue
            data = path.read_bytes()
            name = f"dinara_align/{path.name}"
            info_entry = zipfile.ZipInfo(name)
            # Libraries keep their execute bit, which the loader needs on some systems.
            info_entry.external_attr = (0o755 if path.suffix in (".so", ".dylib") else 0o644) << 16
            info_entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info_entry, data)
            records.append(f"{name},{digest(data)},{len(data)}")
        # The license and its notice, where installers look for them.
        for source in ("LICENSE", "NOTICE"):
            data = (ROOT / source).read_bytes()
            name = f"{info}/licenses/{source}"
            archive.writestr(name, data)
            records.append(f"{name},{digest(data)},{len(data)}")
        for name, text in ((f"{info}/METADATA", metadata), (f"{info}/WHEEL", wheel)):
            archive.writestr(name, text)
            records.append(f"{name},{digest(text.encode())},{len(text.encode())}")
        records.append(f"{info}/RECORD,,")
        archive.writestr(f"{info}/RECORD", "\n".join(records) + "\n")
    print(target)


if __name__ == "__main__":
    main()
