#!/usr/bin/env bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Builds a Mojo program for a baseline CPU and puts beside it the runtime libraries it loads from the pixi
# environment, so the folder runs on a machine with neither pixi nor Mojo:
#
#     pixi run bundle <program.mojo> [out-dir] [target-cpu]
#     bash scripts/bundle.sh --library <exports.mojo> [out-dir] [target-cpu]
#
# With `--library` the Mojo file's `@export` functions become a shared library, `lib<name>.so` or
# `lib<name>.dylib`, beside the same runtime libraries; a program linking it needs only its folder on
# its own search path (see `pixi run build-c`).
#
# `mojo build` has no static mode and writes the environment's absolute library path into the binary ahead
# of any `$ORIGIN` passed through `-Xlinker`, so the bundled copies would be ignored wherever that path
# exists. The binary's search path is rewritten instead: `$ORIGIN` on Linux, `@executable_path` on macOS.
# Linux targets still need glibc 2.35 or later (Ubuntu 22.04, Debian 12, RHEL 10).
set -euo pipefail

library=false
if [ "${1:-}" = --library ]; then
    library=true
    shift
fi
program=${1:?usage: pixi run bundle <program.mojo> [out-dir] [target-cpu]}
name=$(basename "$program" .mojo)
out=${2:-build/bundle/$name}
root=$(cd "$(dirname "$0")/.." && pwd)
env_lib=${CONDA_PREFIX:?run through pixi: pixi run bundle ...}/lib

# The oldest CPU each platform's builds run on; the host's own would bake in its newest instructions.
case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) cpu=x86-64 ;;
    Linux-aarch64) cpu=generic ;;
    Darwin-arm64) cpu=apple-m1 ;;
    *) echo "bundle: no baseline CPU known for $(uname -sm)" >&2; exit 1 ;;
esac
cpu=${3:-$cpu}

mkdir -p "$out"
if $library; then
    [ "$(uname -s)" = Darwin ] && file=lib$name.dylib || file=lib$name.so
    # shellcheck disable=SC2086 # `$MOJO_ACCELERATOR` is a list of flags, empty for the host's own GPU.
    mojo build --emit shared-lib --target-cpu "$cpu" ${MOJO_ACCELERATOR:-} -I "$root" "$program" -o "$out/$file"
else
    file=$name
    # shellcheck disable=SC2086
    mojo build --target-cpu "$cpu" ${MOJO_ACCELERATOR:-} -I "$root" "$program" -o "$out/$file"
fi

if [ "$(uname -s)" = Linux ]; then
    # `ldd` lists every library the binary loads, its libraries' libraries too.
    for path in $(ldd "$out/$file" | awk -v lib="$env_lib/" 'index($3, lib) == 1 { print $3 }'); do
        copy=$out/$(basename "$path")
        cp -L "$path" "$copy"
        chmod u+w "$copy"
        # The environment's libstdc++ carries debug information: 24 MB, 2.8 MB without.
        strip --strip-unneeded "$copy"
        patchelf --set-rpath '$ORIGIN' "$copy"
    done
    patchelf --set-rpath '$ORIGIN' "$out/$file"
    left=$(env -u LD_LIBRARY_PATH ldd "$out/$file" | grep -F "$env_lib/" || true)
else
    # Mojo's dylibs name each other through `@rpath`, which a dylib resolves through the search paths of the
    # images that load it, so only the binary's own needs rewriting; the closure is followed by hand.
    pending=$out/$file
    while [ -n "$pending" ]; do
        current=${pending%%$'\n'*}
        [ "$current" = "$pending" ] && pending= || pending=${pending#*$'\n'}
        for dep in $(otool -L "$current" | awk 'NR > 1 && $1 ~ /^@rpath\// { sub("@rpath/", "", $1); print $1 }'); do
            [ -e "$out/$dep" ] && continue
            cp "$env_lib/$dep" "$out/$dep"
            pending=${pending:+$pending$'\n'}$out/$dep
        done
    done
    for path in $(otool -l "$out/$file" | awk '$1 == "cmd" && $2 == "LC_RPATH" { getline; getline; print $2 }'); do
        install_name_tool -delete_rpath "$path" "$out/$file"
    done
    if $library; then
        # Found through the program's search path, its runtime libraries through its own folder.
        install_name_tool -id "@rpath/$file" "$out/$file"
        install_name_tool -add_rpath @loader_path "$out/$file"
    else
        install_name_tool -add_rpath @executable_path "$out/$file"
    fi
    # Editing load commands voids the signature, and Apple silicon runs nothing unsigned.
    codesign --force --sign - "$out/$file"
    left=$(otool -l "$out/$file" | grep -F "$env_lib" || true)
fi

if [ -n "$left" ]; then
    echo "bundle: $out/$file still loads from the pixi environment:" >&2
    echo "$left" >&2
    exit 1
fi
# The licenses travel with every bundle: MPL-2.0, and Apache-2.0 for the parts taken from AffineGaps (see NOTICE).
cp "$root/LICENSE" "$root/NOTICE" "$out/"
mkdir -p "$out/LICENSES"
cp "$root/LICENSES/Apache-2.0.txt" "$out/LICENSES/"
runtime=$(ls "$out" | grep -c "^lib" || true)
$library && runtime=$((runtime - 1))
echo "$out: $file for $cpu with $runtime runtime libraries, $(du -sh "$out" | cut -f1)"
