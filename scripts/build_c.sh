#!/usr/bin/env bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Builds dinara-align's C API (see c/dinara.h) into build/c: libdinara for the platform's oldest CPU, or
# `target-cpu`, beside the Mojo runtime libraries it loads and the header.
#
#     pixi run build-c [target-cpu]
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
bash "$root/scripts/bundle.sh" --library "$root/c/dinara.mojo" "$root/build/c" "$@"
cp "$root/c/dinara.h" "$root/build/c/"
