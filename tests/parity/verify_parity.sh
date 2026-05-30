#!/usr/bin/env bash
#
# verify_parity.sh — prove the R+C++ library is unchanged in TWO dimensions:
#
#   1. Numerical outputs: the 13 reference JSON files under
#      tests/parity/reference_values/ are byte-identical to what the current
#      source produces.
#
#   2. R API surface: the exported function list, signatures, and default
#      values captured in tests/parity/api_surface/exports.txt are
#      byte-identical to what the current package exports.
#
# Together these are the bit-level contract for "what an R caller can observe
# from this package today." Any code change that alters either dimension
# trips the gate; the change is held to that lens in review.
#
# This script is the executable backing of docs/QA_R_BEHAVIOR_UNCHANGED.md.
#
# Usage:
#   ./tests/parity/verify_parity.sh                 # verify both, return exit code
#   ./tests/parity/verify_parity.sh --keep          # don't delete the temp dir
#   ./tests/parity/verify_parity.sh --update        # overwrite committed snapshots
#                                                     (only after deliberate change)
#   ./tests/parity/verify_parity.sh --save-artifacts
#                                                   # on PASS, write the verification
#                                                   # log + (empty) diffs to
#                                                   # tests/parity/last_verification/
#                                                   # as committed proof artifacts
#
# Exit codes:
#   0   both contracts hold — PASS
#   1   one or both contracts violated — FAIL (full diff saved on failure)
#   2   environment problem (Docker missing, repo layout unexpected, ...)
#
# Requirements:
#   - Docker (regeneration runs inside the pinned r-base:4.5.3 image)
#   - bash, diff, mktemp (standard on macOS and Linux)
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Locate the repo root (this script is invoked from anywhere)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/../.." >/dev/null 2>&1 && pwd )"

REF_DIR="$REPO_ROOT/tests/parity/reference_values"
API_DIR="$REPO_ROOT/tests/parity/api_surface"
GEN_REFS="$REPO_ROOT/tests/parity/generate_reference.R"
GEN_API="$REPO_ROOT/tests/parity/snapshot_api.R"
DOCKERFILE="$REPO_ROOT/Dockerfile"

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
KEEP=0
UPDATE=0
SAVE_ARTIFACTS=0
for arg in "$@"; do
    case "$arg" in
        --keep)            KEEP=1 ;;
        --update)          UPDATE=1 ;;
        --save-artifacts)  SAVE_ARTIFACTS=1 ;;
        -h|--help)
            sed -n '2,30p' "$0" | sed 's/^# //;s/^#//'
            exit 0
            ;;
        *)
            echo "verify_parity.sh: unknown arg: $arg" >&2
            echo "Try --help" >&2
            exit 2
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
need() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "verify_parity.sh: required command not found: $1" >&2
        exit 2
    }
}
need docker
need diff
need mktemp

for path in "$REF_DIR" "$API_DIR" "$GEN_REFS" "$GEN_API" "$DOCKERFILE"; do
    [[ -e "$path" ]] || {
        echo "verify_parity.sh: expected file/dir missing: $path" >&2
        echo "  (is this script being run from a checkout of the bw repo?)" >&2
        exit 2
    }
done

if ! docker info >/dev/null 2>&1; then
    echo "verify_parity.sh: docker is not running" >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Build the pinned Docker image (cached after first run)
# ---------------------------------------------------------------------------
IMAGE_TAG="bw-parity-verify"
echo "==> building $IMAGE_TAG (cached after first build)" >&2
docker build -q -t "$IMAGE_TAG" "$REPO_ROOT" >/dev/null

# ---------------------------------------------------------------------------
# Regenerate BOTH snapshots into a temp directory.
# Layout: $TMP_OUT/references/*.json    +    $TMP_OUT/api_surface/exports.txt
# (NEVER into the committed dirs, unless --update was explicitly passed)
# ---------------------------------------------------------------------------
TMP_OUT="$(mktemp -d -t bw_parity.XXXXXX)"
mkdir -p "$TMP_OUT/references" "$TMP_OUT/api_surface"

cleanup() {
    if [[ "$KEEP" -eq 0 ]]; then
        rm -rf "$TMP_OUT"
    else
        echo "==> kept temp dir: $TMP_OUT" >&2
    fi
}
trap cleanup EXIT

echo "==> regenerating reference values + API surface snapshot" >&2
docker run --rm \
    -v "$REPO_ROOT/tests/parity:/parity:ro" \
    -v "$TMP_OUT:/out" \
    -w /pkg "$IMAGE_TAG" \
    Rscript -e "install.packages('jsonlite', repos='https://cloud.r-project.org', quiet=TRUE); outdir <- '/out/references'; source('/parity/generate_reference.R'); outdir <- '/out/api_surface'; source('/parity/snapshot_api.R')" \
    >/dev/null

# ---------------------------------------------------------------------------
# --update mode: overwrite committed snapshots with the regenerated ones.
# Use after a *deliberate* behavior change. Caller is expected to commit
# the resulting diff.
# ---------------------------------------------------------------------------
if [[ "$UPDATE" -eq 1 ]]; then
    echo "==> --update: overwriting committed snapshots with regenerated outputs" >&2
    for f in "$TMP_OUT/references"/*.json; do
        cp "$f" "$REF_DIR/$(basename "$f")"
    done
    cp "$TMP_OUT/api_surface/exports.txt" "$API_DIR/exports.txt"
    echo "    done. inspect 'git diff tests/parity/' and commit if intentional." >&2
    exit 0
fi

# ---------------------------------------------------------------------------
# Diff regenerated snapshots against committed ones, byte-for-byte
# ---------------------------------------------------------------------------
echo "==> diffing regenerated snapshots vs committed contracts" >&2

DIFF_LOG="$(mktemp -t bw_parity_diff.XXXXXX)"
trap 'cleanup; rm -f "$DIFF_LOG"' EXIT

# --- 1. numerical outputs ----------------------------------------------
REF_CHANGED=()
REF_MISSING=()
REF_EXTRA=()

for committed in "$REF_DIR"/*.json; do
    name="$(basename "$committed")"
    regen="$TMP_OUT/references/$name"
    if [[ ! -f "$regen" ]]; then
        REF_MISSING+=("$name")
        continue
    fi
    if ! cmp -s "$committed" "$regen"; then
        REF_CHANGED+=("$name")
        {
            echo "===== references/$name (- committed, + regenerated) ====="
            diff -u "$committed" "$regen" || true
            echo
        } >> "$DIFF_LOG"
    fi
done

for regen in "$TMP_OUT/references"/*.json; do
    name="$(basename "$regen")"
    [[ -f "$REF_DIR/$name" ]] || REF_EXTRA+=("$name")
done

REF_TOTAL="$(ls -1 "$REF_DIR"/*.json | wc -l | tr -d ' ')"
REF_SAME=$((REF_TOTAL - ${#REF_CHANGED[@]} - ${#REF_MISSING[@]}))
REF_OK=0
if (( ${#REF_CHANGED[@]} == 0 && ${#REF_MISSING[@]} == 0 && ${#REF_EXTRA[@]} == 0 )); then
    REF_OK=1
fi

# --- 2. R API surface --------------------------------------------------
API_OK=0
if cmp -s "$API_DIR/exports.txt" "$TMP_OUT/api_surface/exports.txt"; then
    API_OK=1
else
    {
        echo "===== api_surface/exports.txt (- committed, + regenerated) ====="
        diff -u "$API_DIR/exports.txt" "$TMP_OUT/api_surface/exports.txt" || true
        echo
    } >> "$DIFF_LOG"
fi

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
echo
echo "----- R parity verification report -----"
echo "  numerical outputs:"
echo "    committed reference files:   $REF_TOTAL"
echo "    byte-identical to current:   $REF_SAME"
echo "    drifted (content changed):   ${#REF_CHANGED[@]}"
echo "    missing from regeneration:   ${#REF_MISSING[@]}"
echo "    extra in regeneration:       ${#REF_EXTRA[@]}"
echo "  R API surface:"
if (( API_OK )); then
    echo "    exports.txt:                 byte-identical"
else
    echo "    exports.txt:                 DRIFTED"
fi
echo "----------------------------------------"

if (( REF_OK && API_OK )); then
    echo
    echo "PASS: both contracts hold."
    echo "      - R+C++ outputs are byte-identical to the committed reference."
    echo "      - R API surface is byte-identical to the committed snapshot."
    echo "      No R-visible behavior change in this checkout."

    # --- on --save-artifacts, write proof artifacts into the repo --------
    if (( SAVE_ARTIFACTS )); then
        ART_DIR="$REPO_ROOT/tests/parity/last_verification"
        mkdir -p "$ART_DIR"

        # 1. Per-contract diffs (empty on success) — proof of no difference
        : > "$ART_DIR/outputs.diff"
        for committed in "$REF_DIR"/*.json; do
            name="$(basename "$committed")"
            regen="$TMP_OUT/references/$name"
            diff -u "$committed" "$regen" >> "$ART_DIR/outputs.diff" || true
        done

        : > "$ART_DIR/api_surface.diff"
        diff -u "$API_DIR/exports.txt" "$TMP_OUT/api_surface/exports.txt" \
            >> "$ART_DIR/api_surface.diff" || true

        # 2. Per-case file inventory (sha256 of each committed snapshot —
        #    a quick "what was proven" manifest, stable across runs)
        {
            echo "# Verification manifest"
            echo "# Files proven byte-identical to current source via verify_parity.sh"
            echo "# (Re)generate via: ./tests/parity/verify_parity.sh --save-artifacts"
            echo
            echo "numerical outputs:"
            (cd "$REF_DIR" && shasum -a 256 *.json | sort | sed 's/^/  /')
            echo
            echo "R API surface:"
            (cd "$API_DIR" && shasum -a 256 exports.txt | sed 's/^/  /')
        } > "$ART_DIR/manifest.txt"

        # 3. Human-readable summary report
        {
            echo "Parity verification — PASS"
            echo "=========================="
            echo
            echo "All R-visible snapshots are byte-identical to the committed contracts."
            echo
            echo "Numerical outputs:"
            echo "  $REF_TOTAL committed reference files, $REF_SAME byte-identical, 0 drifted."
            echo "R API surface:"
            echo "  exports.txt byte-identical."
            echo
            echo "outputs.diff      empty (no per-byte difference in any of the $REF_TOTAL JSONs)"
            echo "api_surface.diff  empty (no difference in exports.txt)"
            echo "manifest.txt      sha256 of every contract file at the moment of proof"
            echo
            echo "Reproduce: ./tests/parity/verify_parity.sh"
        } > "$ART_DIR/verification.log"

        echo
        echo "==> wrote proof artifacts to tests/parity/last_verification/"
        echo "    - verification.log   (human-readable summary)"
        echo "    - outputs.diff       (empty — no per-byte difference in 13 JSONs)"
        echo "    - api_surface.diff   (empty — no difference in exports.txt)"
        echo "    - manifest.txt       (sha256 of every contract file)"
    fi

    exit 0
fi

echo
echo "FAIL: R-visible behavior or API has changed."
if (( ${#REF_CHANGED[@]} > 0 )); then
    echo "  numerical outputs changed:"
    for f in "${REF_CHANGED[@]}"; do echo "    - $f"; done
fi
if (( ${#REF_MISSING[@]} > 0 )); then
    echo "  numerical outputs missing from regeneration:"
    for f in "${REF_MISSING[@]}"; do echo "    - $f"; done
fi
if (( ${#REF_EXTRA[@]} > 0 )); then
    echo "  numerical outputs extra in regeneration:"
    for f in "${REF_EXTRA[@]}"; do echo "    - $f"; done
fi
if (( ! API_OK )); then
    echo "  R API surface: signatures, exports, or defaults changed."
fi

PERM_DIFF="$REPO_ROOT/r_parity_diff.log"
cp "$DIFF_LOG" "$PERM_DIFF"
echo
echo "Full diff saved to: $PERM_DIFF"
echo "If the change was *intentional*, run:"
echo "    $0 --update"
echo "and commit the resulting snapshot updates with a message that explains"
echo "what changed and why."
exit 1
