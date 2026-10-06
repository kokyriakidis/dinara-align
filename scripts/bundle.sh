#!/usr/bin/env bash
# Builds a Mojo program for a baseline CPU and puts beside it the runtime libraries it loads from the pixi
# environment, so the folder runs on a machine with neither pixi nor Mojo:
#
#     pixi run bundle <program.mojo> [out-dir] [target-cpu]
#
# `mojo build` has no static mode and writes the environment's absolute library path into the binary ahead
# of any `$ORIGIN` passed through `-Xlinker`, so the bundled copies would be ignored wherever that path
# exists. The binary's search path is rewritten instead: `$ORIGIN` on Linux, `@executable_path` on macOS.
# Linux targets still need glibc 2.35 or later (Ubuntu 22.04, Debian 12, RHEL 10).
set -euo pipefail

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
mojo build --target-cpu "$cpu" -I "$root" "$program" -o "$out/$name"

if [ "$(uname -s)" = Linux ]; then
    # `ldd` lists every library the binary loads, its libraries' libraries too.
    for path in $(ldd "$out/$name" | awk -v lib="$env_lib/" 'index($3, lib) == 1 { print $3 }'); do
        copy=$out/$(basename "$path")
        cp -L "$path" "$copy"
        chmod u+w "$copy"
        # The environment's libstdc++ carries debug information: 24 MB, 2.8 MB without.
        strip --strip-unneeded "$copy"
        patchelf --set-rpath '$ORIGIN' "$copy"
    done
    patchelf --set-rpath '$ORIGIN' "$out/$name"
    left=$(env -u LD_LIBRARY_PATH ldd "$out/$name" | grep -F "$env_lib/" || true)
else
    # Mojo's dylibs name each other through `@rpath`, which a dylib resolves through the executable's search
    # path, so only the executable's needs rewriting; the closure is followed by hand.
    pending=$out/$name
    while [ -n "$pending" ]; do
        file=${pending%%$'\n'*}
        [ "$file" = "$pending" ] && pending= || pending=${pending#*$'\n'}
        for dep in $(otool -L "$file" | awk 'NR > 1 && $1 ~ /^@rpath\// { sub("@rpath/", "", $1); print $1 }'); do
            [ -e "$out/$dep" ] && continue
            cp "$env_lib/$dep" "$out/$dep"
            pending=${pending:+$pending$'\n'}$out/$dep
        done
    done
    for path in $(otool -l "$out/$name" | awk '$1 == "cmd" && $2 == "LC_RPATH" { getline; getline; print $2 }'); do
        install_name_tool -delete_rpath "$path" "$out/$name"
    done
    install_name_tool -add_rpath @executable_path "$out/$name"
    # Editing load commands voids the signature, and Apple silicon runs nothing unsigned.
    codesign --force --sign - "$out/$name"
    left=$(otool -l "$out/$name" | grep -F "$env_lib" || true)
fi

if [ -n "$left" ]; then
    echo "bundle: $out/$name still loads from the pixi environment:" >&2
    echo "$left" >&2
    exit 1
fi
echo "$out: $name for $cpu with $(($(ls "$out" | wc -l) - 1)) libraries, $(du -sh "$out" | cut -f1)"
