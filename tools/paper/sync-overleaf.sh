#!/usr/bin/env bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Pushes paper/, as committed, to the Overleaf project behind the `overleaf` remote (Overleaf's Git
# integration), so that the project holds exactly the paper's files.
#
#     git remote add overleaf https://git.overleaf.com/<project id>    # once; Overleaf: Menu > Git
#     pixi run paper-sync            # after committing a change to paper/
#     pixi run paper-sync --force    # the first time, to replace the new project's placeholder files
#
# The sync is one way: paper/ is generated from tools/paper/, so an edit made in Overleaf would be lost at
# the next `pixi run paper`. Git refuses a plain push over such an edit, and the script stops there.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

if ! git remote get-url overleaf >/dev/null 2>&1; then
    echo "No overleaf remote: git remote add overleaf https://git.overleaf.com/<project id> (Overleaf: Menu > Git)" >&2
    exit 1
fi
if ! git diff --quiet HEAD -- paper || [ -n "$(git ls-files --others --exclude-standard paper)" ]; then
    echo "paper/ has changes that are not committed; commit them, then sync" >&2
    exit 1
fi

# paper/'s own history, its files at the top level, as Overleaf expects them.
split=$(git subtree split --prefix=paper HEAD 2>/dev/null)
if ! git push overleaf "$split:refs/heads/master" "$@"; then
    echo "Overleaf refused the push: the project has changes made there. Look at them (git fetch overleaf)," >&2
    echo "and if they can go, push again with --force." >&2
    exit 1
fi
echo "Overleaf now has paper/ at $(git rev-parse --short HEAD)."
