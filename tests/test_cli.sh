#!/usr/bin/env bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# The command-line aligner (see cli/dinara_align_cli.mojo) on a small reference and reads: each output
# format, both strands, a cap, and a refusal.
#
#     pixi run test-cli
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cli=$root/build/cli/dinara-align
work=$(mktemp -d)
trap 'rm -r "$work"' EXIT
printf ">chr\nTTTTTTTTTTACGTACGTTTGCAGGATCCAGTTTTTTTTTT\n" >"$work/ref.fa"
printf "@r1\nACGTACGTTTGCAGGATCCAG\n+\nIIIIIIIIIIIIIIIIIIIII\n@r2\nCTGGATCCTGCAAACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIII\n" >"$work/reads.fq"
failures=0

expect() {
    if [ "$1" != "$2" ]; then
        printf 'failed: %s\n  got:  %s\n  want: %s\n' "$3" "$1" "$2" >&2
        failures=$((failures + 1))
    fi
}

expect "$("$cli" -r ACGTACGTTTGCA -q ACGTCGTTTTGCA | tail -1)" \
    "$(printf 'query\treference\t+\t2\t-2\t0\t13\t0\t13\t4=1D2=1I6=')" "a pair given whole"
expect "$("$cli" --costs affine:4,6,2 --mode infix "$work/ref.fa" "$work/reads.fq" | sed -n 2p)" \
    "$(printf 'r1\tchr\t+\t0\t0\t10\t31\t0\t21\t21=')" "a read placed in a reference"
expect "$("$cli" --costs affine:4,6,2 --mode infix --both-strands --format sam "$work/ref.fa" "$work/reads.fq" | tail -1 | cut -f1-6,12-13)" \
    "$(printf 'r2\t16\tchr\t11\t255\t21=\tNM:i:0\tMD:Z:21')" "SAM, the read on the reverse strand"
expect "$("$cli" --costs affine:4,6,2 --mode local:2 --format paf "$work/ref.fa" "$work/reads.fq" | head -1 | cut -f1-12)" \
    "$(printf 'r1\t21\t0\t21\t+\tchr\t41\t10\t31\t21\t21\t255')" "PAF, a local alignment"
expect "$("$cli" --distance --mode infix --max-cost 3 "$work/ref.fa" "$work/reads.fq" | tail -1)" \
    "$(printf 'r2\tchr\t*')" "a cap"
if "$cli" --costs affine:0,6,2 -r A -q A 2>/dev/null; then
    expect "accepted" "refused" "costs a search cannot run by"
fi

# Options and values the command line would drop without a word, and bytes that are not a sequence's.
for refused in "--mode local:2 --match-score 5" "--mode global:7" "--costs affine:4,6,2:9" "--mode infix --zdrop 3" \
    "--mode infix --match-score -3"; do
    # shellcheck disable=SC2086
    if "$cli" -r ACGTACGTTTGCA -q ACGTCGTTTTGCA $refused >/dev/null 2>&1; then
        expect "accepted" "refused" "$refused"
    fi
done
if "$cli" -r GT -q "$(printf '\370\216\216\216')" --mode reference-in-query >/dev/null 2>&1; then
    expect "accepted" "refused" "a sequence that is not ASCII"
fi

if [ "$failures" -gt 0 ]; then
    echo "$failures checks failed" >&2
    exit 1
fi
echo "command line: every check passed"
