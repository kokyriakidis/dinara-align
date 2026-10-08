#!/usr/bin/env bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Builds the Python package (see python/dinara_align) into build/python: the package's Python, its
# extension `_dinara` for the platform's oldest CPU or `target-cpu`, and the Mojo runtime libraries it
# loads, so `PYTHONPATH=build/python python3 -c "import dinara_align"` works without pixi or Mojo.
#
#     pixi run build-python [target-cpu]
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
out=$root/build/python/dinara_align
mkdir -p "$out"
bash "$root/scripts/bundle.sh" --library "$root/python/_dinara.mojo" "$out" "$@"
# Python imports an extension by its module's name, the bundle names a library `lib<name>`.
for library in "$out"/lib_dinara.so "$out"/lib_dinara.dylib; do
    [ -e "$library" ] && mv "$library" "$out/_dinara.so"
done
cp "$root/python/dinara_align/"*.py "$out/"
echo "$out: import dinara_align with $(dirname "$out") on PYTHONPATH"
