#!/usr/bin/env bash
# Builds the command-line aligner (see cli/dinara_align_cli.mojo) into build/cli as `dinara-align`, for
# the platform's oldest CPU or `target-cpu`, beside the runtime libraries it loads, so the folder runs
# without pixi or Mojo.
#
#     pixi run build-cli [target-cpu]
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
bash "$root/scripts/bundle.sh" "$root/cli/dinara_align_cli.mojo" "$root/build/cli" "$@"
mv "$root/build/cli/dinara_align_cli" "$root/build/cli/dinara-align"
echo "$root/build/cli/dinara-align"
