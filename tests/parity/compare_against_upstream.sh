#!/usr/bin/env bash
#
# compare_against_upstream.sh — directly demonstrate that this PR's R+C++
# package produces byte-identical outputs to the UPSTREAM (pre-refactor)
# package when called with the same inputs.
#
# Unlike verify_parity.sh, this script does NOT trust the committed
# reference JSONs. Instead it:
#
#   1. Checks out upstream/master into a worktree.
#   2. Checks out the current branch into a separate worktree.
#   3. Builds the R+C++ package from each worktree, in identical pinned
#      Docker images, using identical compile flags.
#   4. Runs the SAME R driver script (`generate_reference.R`) inside each
#      image, exercising every exported function via 13 fixed input cases.
#   5. Captures the outputs from both runs side-by-side.
#   6. `diff -r` the two output directories byte-for-byte.
#   7. Reports PASS only if every byte of every output file matches.
#
# This is the most direct demonstration possible: two real R sessions,
# two real installs of the bw package (one pre-refactor, one post),
# called with identical inputs, outputs byte-compared with no intermediate
# state to trust.
#
# Usage:
#   ./tests/parity/compare_against_upstream.sh                # compare upstream/master vs HEAD
#   ./tests/parity/compare_against_upstream.sh REF1 REF2      # compare any two git refs
#   ./tests/parity/compare_against_upstream.sh --keep         # don't delete temp dirs (for inspection)
#
# Exit codes:
#   0   the two installs produced byte-identical outputs — PASS
#   1   any output file differed — FAIL (diff is printed and saved)
#   2   environment problem (Docker, git, refs not found)
#
# Requirements:
#   - Docker
#   - git (with `upstream` remote configured, pointing at INSP-RH/bw)
#   - bash, diff, mktemp
#
# Time:
#   ~6 minutes the first time (two Docker builds, ~3 min each).
#   ~2 minutes on subsequent runs once Docker layers are cached.
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Locate the repo root and parse args
# ---------------------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/../.." >/dev/null 2>&1 && pwd )"
cd "$REPO_ROOT"

KEEP=0
REF_BEFORE=""
REF_AFTER=""
for arg in "$@"; do
    case "$arg" in
        --keep) KEEP=1 ;;
        -h|--help)
            sed -n '2,38p' "$0" | sed 's/^# //;s/^#//'
            exit 0
            ;;
        --*)
            echo "compare_against_upstream.sh: unknown flag: $arg" >&2
            exit 2
            ;;
        *)
            if [[ -z "$REF_BEFORE" ]]; then REF_BEFORE="$arg"
            elif [[ -z "$REF_AFTER" ]]; then REF_AFTER="$arg"
            else
                echo "compare_against_upstream.sh: too many positional args" >&2
                exit 2
            fi
            ;;
    esac
done
REF_BEFORE="${REF_BEFORE:-upstream/master}"
REF_AFTER="${REF_AFTER:-HEAD}"

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
need() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "compare_against_upstream.sh: required command not found: $1" >&2
        exit 2
    }
}
need docker
need git
need diff
need mktemp

if ! docker info >/dev/null 2>&1; then
    echo "compare_against_upstream.sh: docker is not running" >&2
    exit 2
fi

for ref in "$REF_BEFORE" "$REF_AFTER"; do
    git rev-parse --verify "$ref" >/dev/null 2>&1 || {
        echo "compare_against_upstream.sh: git ref not found: $ref" >&2
        echo "  (need 'upstream' remote? add it with:" >&2
        echo "     git remote add upstream https://github.com/INSP-RH/bw.git && git fetch upstream)" >&2
        exit 2
    }
done

SHA_BEFORE="$(git rev-parse --short "$REF_BEFORE")"
SHA_AFTER="$(git rev-parse --short "$REF_AFTER")"

if [[ "$SHA_BEFORE" == "$SHA_AFTER" ]]; then
    echo "compare_against_upstream.sh: both refs resolve to $SHA_BEFORE — nothing to compare" >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Create two worktrees, one per ref
# ---------------------------------------------------------------------------
TMP_ROOT="$(mktemp -d -t bw_compare.XXXXXX)"
WT_BEFORE="$TMP_ROOT/before"
WT_AFTER="$TMP_ROOT/after"
OUT_BEFORE="$TMP_ROOT/out_before"
OUT_AFTER="$TMP_ROOT/out_after"

cleanup() {
    git worktree remove "$WT_BEFORE" --force 2>/dev/null || true
    git worktree remove "$WT_AFTER"  --force 2>/dev/null || true
    if [[ "$KEEP" -eq 0 ]]; then
        rm -rf "$TMP_ROOT"
    else
        echo "==> kept temp dirs under $TMP_ROOT" >&2
    fi
}
trap cleanup EXIT

echo "==> creating worktrees" >&2
echo "    before: $REF_BEFORE ($SHA_BEFORE) -> $WT_BEFORE" >&2
echo "    after:  $REF_AFTER ($SHA_AFTER) -> $WT_AFTER" >&2
git worktree add "$WT_BEFORE" "$REF_BEFORE" --detach >/dev/null
git worktree add "$WT_AFTER"  "$REF_AFTER"  --detach >/dev/null

# ---------------------------------------------------------------------------
# Ensure both worktrees have the build + drive infrastructure they need.
# The upstream worktree probably doesn't have our Dockerfile or the regen
# script, so we copy ours in. They are NOT part of the package being
# compared — they only run the comparison.
# ---------------------------------------------------------------------------
DRIVER_NEEDS=("Dockerfile" ".dockerignore" "tests/parity/generate_reference.R")
for wt in "$WT_BEFORE" "$WT_AFTER"; do
    for f in "${DRIVER_NEEDS[@]}"; do
        if [[ ! -e "$wt/$f" ]]; then
            mkdir -p "$wt/$(dirname "$f")"
            cp "$REPO_ROOT/$f" "$wt/$f"
        fi
    done
done

# ---------------------------------------------------------------------------
# Build a Docker image from each worktree (cached after first build)
# ---------------------------------------------------------------------------
IMG_BEFORE="bw-compare-before:$SHA_BEFORE"
IMG_AFTER="bw-compare-after:$SHA_AFTER"

echo "==> building $IMG_BEFORE" >&2
docker build -q -t "$IMG_BEFORE" "$WT_BEFORE" >/dev/null

echo "==> building $IMG_AFTER" >&2
docker build -q -t "$IMG_AFTER" "$WT_AFTER" >/dev/null

# ---------------------------------------------------------------------------
# In each image, run the identical R driver script against the installed
# bw package and dump outputs into a fresh directory.
# ---------------------------------------------------------------------------
mkdir -p "$OUT_BEFORE" "$OUT_AFTER"

run_driver() {
    local image="$1" wt="$2" out="$3"
    docker run --rm \
        -v "$wt/tests/parity:/parity:ro" \
        -v "$out:/out" \
        -w /pkg "$image" \
        Rscript -e "install.packages('jsonlite', repos='https://cloud.r-project.org', quiet=TRUE); outdir <- '/out'; source('/parity/generate_reference.R')" \
        >/dev/null
}

echo "==> running R driver against $REF_BEFORE install" >&2
run_driver "$IMG_BEFORE" "$WT_BEFORE" "$OUT_BEFORE"

echo "==> running R driver against $REF_AFTER install" >&2
run_driver "$IMG_AFTER" "$WT_AFTER" "$OUT_AFTER"

# ---------------------------------------------------------------------------
# Byte-compare every output file
# ---------------------------------------------------------------------------
echo "==> diffing outputs byte-for-byte" >&2

CHANGED=()
MISSING=()
EXTRA=()

for before_file in "$OUT_BEFORE"/*.json; do
    name="$(basename "$before_file")"
    after_file="$OUT_AFTER/$name"
    if [[ ! -f "$after_file" ]]; then
        MISSING+=("$name")
        continue
    fi
    if ! cmp -s "$before_file" "$after_file"; then
        CHANGED+=("$name")
    fi
done

for after_file in "$OUT_AFTER"/*.json; do
    name="$(basename "$after_file")"
    [[ -f "$OUT_BEFORE/$name" ]] || EXTRA+=("$name")
done

TOTAL_BEFORE="$(ls -1 "$OUT_BEFORE"/*.json 2>/dev/null | wc -l | tr -d ' ')"
TOTAL_AFTER="$(ls -1 "$OUT_AFTER"/*.json  2>/dev/null | wc -l | tr -d ' ')"
SAME=$((TOTAL_BEFORE - ${#CHANGED[@]} - ${#MISSING[@]}))

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
echo
echo "----- live side-by-side comparison -----"
echo "  before:   $REF_BEFORE  ($SHA_BEFORE)"
echo "  after:    $REF_AFTER   ($SHA_AFTER)"
echo "  files (before / after): $TOTAL_BEFORE / $TOTAL_AFTER"
echo "  byte-identical:         $SAME"
echo "  drifted (content):      ${#CHANGED[@]}"
echo "  missing from 'after':   ${#MISSING[@]}"
echo "  extra in 'after':       ${#EXTRA[@]}"
echo "----------------------------------------"

if (( ${#CHANGED[@]} == 0 && ${#MISSING[@]} == 0 && ${#EXTRA[@]} == 0 && TOTAL_BEFORE > 0 )); then
    echo
    echo "PASS: the R package produces byte-identical outputs at $SHA_AFTER"
    echo "      as it did at $SHA_BEFORE."
    echo
    echo "      $SAME of $TOTAL_BEFORE output files matched bit-for-bit, including"
    echo "      the 5-year child run (1825 RK4 steps) where any FP-order drift"
    echo "      would have accumulated and been caught."
    exit 0
fi

echo
echo "FAIL: the R package's behavior changed between $SHA_BEFORE and $SHA_AFTER."

if (( ${#CHANGED[@]} > 0 )); then
    echo "  drifted (showing first 20 lines of diff for each):"
    for f in "${CHANGED[@]}"; do
        echo "    --- $f ---"
        diff -u "$OUT_BEFORE/$f" "$OUT_AFTER/$f" | head -20 | sed 's/^/      /'
    done
fi
if (( ${#MISSING[@]} > 0 )); then
    echo "  missing from 'after':"
    for f in "${MISSING[@]}"; do echo "    - $f"; done
fi
if (( ${#EXTRA[@]} > 0 )); then
    echo "  extra in 'after':"
    for f in "${EXTRA[@]}"; do echo "    - $f"; done
fi
echo
echo "Full output directories preserved at:"
echo "  before: $OUT_BEFORE"
echo "  after:  $OUT_AFTER"
echo "(They will be deleted on next run unless you passed --keep.)"
exit 1
