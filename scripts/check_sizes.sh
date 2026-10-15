#!/usr/bin/env bash
# Refuse to commit files larger than 1 MB (data belong in data/ or in the Zenodo deposit, never in git).
set -euo pipefail
limit=$((1024 * 1024)); bad=0
while IFS= read -r f; do
    [[ -f $f ]] || continue
    s=$(stat -c %s "$f")
    if (( s > limit )); then echo "too large for git ($s bytes): $f"; bad=1; fi
done < <(git diff --cached --name-only)
exit $bad
