#!/usr/bin/env bash
# Builds dinara-align's C API (see c/dinara.h) into build/c: libdinara for the platform's oldest CPU, or
# `target-cpu`, beside the Mojo runtime libraries it loads and the header.
#
#     pixi run build-c [target-cpu]
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
bash "$root/scripts/bundle.sh" --library "$root/c/dinara.mojo" "$root/build/c" "$@"
cp "$root/c/dinara.h" "$root/build/c/"
