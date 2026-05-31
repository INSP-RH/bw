# QA — Proof That R Behavior Is Unchanged

**Status:** baseline established. As of the commit that introduces this
document, running `./tests/parity/verify_parity.sh` produces:

```
----- R parity verification report -----
  numerical outputs:
    committed reference files:   13
    byte-identical to current:   13
    drifted (content changed):   0
    missing from regeneration:   0
    extra in regeneration:       0
  R API surface:
    exports.txt:                 byte-identical
----------------------------------------

PASS: both contracts hold.
      - R+C++ outputs are byte-identical to the committed reference.
      - R API surface is byte-identical to the committed snapshot.
      No R-visible behavior change in this checkout.
```

The gate enforces **two** byte-level contracts on what an R caller can
observe from this package:

1. **Numerical outputs** — the 13 reference JSON files under
   `tests/parity/reference_values/` are the same bits the current source
   produces for 13 fixed input cases.
2. **R API surface** — the exported function list, signatures, and default
   values captured in `tests/parity/api_surface/exports.txt` are the same
   bytes the current package exports.

Any code change that touches the R/C++ math layer, the `R/` wrappers, or
`NAMESPACE` must produce the same exit-code-zero report, or the change is
considered to have altered R-visible behavior and must be reviewed under
that lens.

---

## Why this exists

The `bw` library is a published clinical model with an existing user base.
Future refactors of `src/` — for example, isolating the pure-C++ math from
the Rcpp wrappers so it becomes reachable from other languages — are
intended to be structurally invisible: same algorithm, same constants, same
operation order. But **"structurally invisible" is not the same as "provably
unchanged"**. Without a mechanical guarantee, "we believe this is
behavior-preserving" is a claim that can quietly degrade across review.

This document, the catalog it points at, and the script it references are
that mechanical guarantee. They make "no R behavior change" a CI check, not
an assertion.

The scope of this PR is **only the guarantee mechanism**. No `R/`, `src/`,
`man/`, or `DESCRIPTION` files are touched. The proof script returns exit
code 0 against this commit, establishing the baseline; subsequent PRs that
touch the math layer are measured against it.

---

## What "byte-identical" means here and what it does not mean

The verification asserts that the snapshot files written by the R+C++
pipeline today are *the same bytes* as the files committed to the
repository — for both the numerical output snapshots and the API surface
snapshot. That entails identity at all of:

- the bit-level value of every floating-point number produced (the JSONs
  are serialised at 17 significant digits — full IEEE-754 `double` precision),
- the order of fields and array elements,
- the names, argument order, and default values of every exported R function,
- the file size, line endings, and trailing whitespace.

It does **not** assert:

- that the library is mathematically correct relative to the published model
  (that is the job of the existing `tests/testthat/` suite),
- that the library will behave identically on different hardware or with a
  different `libm` (the Docker image pins R, the compiler, glibc, and the
  toolchain — change the image and you change the contract).

The strict reading is the right one: the contract is "this specific source
tree, compiled inside this specific image, produces these specific bytes."
Any drift means *something* in that chain changed, and a human has to look
at the diff and decide whether the change was intentional.

---

## The proof methodology

```
                                                       ┌─ committed contracts ───────────────────────┐
                                                       │ tests/parity/reference_values/*.json (13)   │
                                                       │   17-sig-digit numerical outputs            │
                                                       │ tests/parity/api_surface/exports.txt        │
                                                       │   R exports, signatures, defaults           │
                                                       └────────────────┬────────────────────────────┘
                                                                        │
              ┌────────────────────────────────────────┐                │
              │ tests/parity/generate_reference.R      │                │ byte-compare
              │   13 fixed input cases, deterministic  │                │ both files
              │ tests/parity/snapshot_api.R            │                │
              │   walks the installed package namespace│                │
              └────────────────┬───────────────────────┘                │
                               │                                        │
                               │ runs inside r-base:4.5.3               │
                               ▼                                        │
                  ┌────────────────────────────┐                        │
                  │ R + C++ source (current)   │ ── produces ── temp/refs/*.json ──┐
                  │   src/, R/, NAMESPACE, ... │ ── produces ── temp/api/exports.txt ─┘
                  └────────────────────────────┘
```

The script's job is to make the right-hand path mechanical: regenerate
both snapshots, diff each against its committed twin, report. It does not
interpret the result.

The chain is reproducible because every link is pinned:

| link | pinned by |
|---|---|
| OS, glibc | `r-base:4.5.3` Docker image (defined by `Dockerfile`) |
| R version | `r-base:4.5.3` |
| C++ compiler, optimization flags | R's `Makevars` defaults inherited from `r-base:4.5.3` |
| C++ source | this repo's `src/` |
| R source, package metadata | this repo's `R/`, `DESCRIPTION`, `NAMESPACE` |
| Algorithm inputs | `tests/parity/generate_reference.R` (13 cases) |
| Output precision | `jsonlite::write_json(..., digits = 17)` |
| Output layout | `generate_reference.R`'s `serialize_model()` (rows = individuals) |
| API snapshot format | `tests/parity/snapshot_api.R` (deterministic via `base::args()`) |

Change any link and the contract changes — by design.

---

## How to run the proof

From the repo root:

```bash
./tests/parity/verify_parity.sh
```

The script:

1. Builds the pinned Docker image (cached after first run; subsequent
   invocations finish in a few seconds plus regeneration time).
2. Inside the image, regenerates **both** snapshots from the current
   source into a fresh temp directory: the 13 reference JSONs from
   `generate_reference.R`, plus `exports.txt` from `snapshot_api.R`.
3. Byte-compares each regenerated file against its committed twin.
4. Prints a per-contract summary report and exits 0 (both contracts hold)
   or 1 (either contract drifted).

The committed snapshots are **never touched** unless `--update` is passed
explicitly. Use `--update` only when a drift is the intentional result of
a change a human has reviewed and accepted, and commit the resulting diff
with a message that explains what changed.

For full usage:

```bash
./tests/parity/verify_parity.sh --help
```

To leave a **committed proof artifact** in the repo after a successful
run (used by the refactor commits to ship visible evidence that the gate
held), pass `--save-artifacts`:

```bash
./tests/parity/verify_parity.sh --save-artifacts
```

This writes, on a passing run only:

- `tests/parity/last_verification/verification.log` — human-readable summary
- `tests/parity/last_verification/outputs.diff` — empty file (proves zero
  per-byte difference across all 13 reference JSONs)
- `tests/parity/last_verification/api_surface.diff` — empty file (proves
  zero difference in `exports.txt`)
- `tests/parity/last_verification/manifest.txt` — sha256 of every contract
  file at the moment of proof

These are committed alongside the code change they prove. A reviewer who
doesn't want to run Docker themselves can read these directly and verify
the hashes against `tests/parity/reference_values/` and
`tests/parity/api_surface/`.

See [`REFERENCE_VALUES_CATALOG.md`](REFERENCE_VALUES_CATALOG.md) for what
each of the 13 cases exercises.

---

## How this gets used in code review

For any PR that changes files under `src/`, `R/`, `tests/parity/generate_reference.R`,
or the `Dockerfile`:

1. The PR author runs `./tests/parity/verify_parity.sh` locally before
   pushing.
2. If the script returns 0, the PR is by mechanical proof a
   no-behavior-change PR (or at least, no change to any of the 13 contracted
   outputs). Reviewers focus on code structure, not numerical correctness.
3. If the script returns 1, the PR is by mechanical proof a behavior change.
   The PR description must explain *why*, and the regenerated references
   must be committed in the same PR (via `--update`).

CI integration (a GitHub Actions workflow that runs the script in a clean
environment on every push) is a natural follow-up; the script provides
everything CI needs.

---

## Current pass status

The script has been run end-to-end against this branch as of the commit that
introduces this document. Verbatim output:

```
==> building bw-parity-verify (cached after first build)
==> regenerating reference values + API surface snapshot
==> diffing regenerated snapshots vs committed contracts

----- R parity verification report -----
  numerical outputs:
    committed reference files:   13
    byte-identical to current:   13
    drifted (content changed):   0
    missing from regeneration:   0
    extra in regeneration:       0
  R API surface:
    exports.txt:                 byte-identical
----------------------------------------

PASS: both contracts hold.
      - R+C++ outputs are byte-identical to the committed reference.
      - R API surface is byte-identical to the committed snapshot.
      No R-visible behavior change in this checkout.
```

Exit code: 0.

Per-file confirmation (13 of 13):

| file | size | status |
|---|---:|---|
| `adult_female_baseline.json` | 43,145 B | identical |
| `adult_male_baseline.json` | 48,254 B | identical |
| `adult_male_deficit.json` | 66,261 B | identical |
| `adult_male_ei2400.json` | 48,841 B | identical |
| `adult_multi.json` | 89,825 B | identical |
| `child_female_730d.json` | 58,061 B | identical |
| `child_male_365d.json` | 29,038 B | identical |
| `child_male_5yr.json` | 146,044 B | identical (1,825 RK4 steps over 5 years) |
| `energy_exponential.json` | 233 B | identical |
| `energy_linear.json` | 88 B | identical |
| `energy_logarithmic.json` | 204 B | identical |
| `energy_stepwise_l.json` | 92 B | identical |
| `energy_stepwise_r.json` | 92 B | identical |

This is the baseline. Every future change is measured against it.

---

## What invalidates each contract

**Numerical-output contract** breaks — and the gate (correctly) fails — when:

- The C++ source in `src/` is changed in a way that alters the numerical
  output of any of the 13 cases. This is the common case the proof is
  designed to catch.
- The R source in `R/` is changed in a way that alters argument defaults,
  argument handling, or output structure for any function in any of the
  13 cases.
- The Docker image's R version, compiler, or `libm` changes.
- The output-serialisation code in `tests/parity/generate_reference.R`
  changes (e.g. someone reduces `digits = 17` to `digits = 15`, or changes
  array layout).
- The input cases in `tests/parity/generate_reference.R` change.

**API-surface contract** breaks — and the gate (correctly) fails — when:

- `NAMESPACE` adds, removes, or renames an exported function.
- Any exported function's argument list changes (added, removed, renamed,
  reordered, or default value altered).
- The package version in `DESCRIPTION` changes (intentional — version bumps
  are deliberate contract changes).
- The S3 / S4 method dispatch surface changes.
- The snapshot format in `tests/parity/snapshot_api.R` changes.

The first two bullets of each contract are the cases the proof exists to
discipline. The first bullet of the API contract is particularly important
because R has no compile-time signature checks — silent argument-list
changes break downstream callers at runtime.

The third, fourth, and fifth are *intentional changes to the contract itself*.
When they happen, the right response is:

1. Run `./tests/parity/verify_parity.sh --update` to refresh the committed
   references.
2. In the same commit that changes the toolchain / serialisation / inputs,
   also commit the regenerated reference files.
3. The commit message must explain that the contract changed and why.

This keeps "the contract changed" a deliberate, reviewable action rather
than something that quietly happens.

---

## Why this is the right shape for this library

Three properties make `bw` a particularly clean fit for byte-level golden
testing:

1. **Deterministic numerics.** The library uses RK4 with a fixed step size
   on deterministic ODEs. No randomness, no parallel reductions, no
   non-deterministic BLAS. The same inputs produce the same bits every
   run, given the same toolchain.
2. **Small, stable input space.** The 13 cases cover every code path. The
   library does not have a combinatorial input space that would force
   tolerance-based testing.
3. **High-stakes correctness.** It's a published clinical model.
   Tolerance-based tests would silently allow drift that flips a
   clinically meaningful bit; byte-level testing makes any such drift
   impossible to land accidentally.

For libraries that fail any of these properties — non-deterministic BLAS,
floating-point reductions across threads, or huge input spaces — the
appropriate gate is a tolerance-based comparison with carefully chosen
`rtol`/`atol`. That is not this library, and this document deliberately
takes the stricter contract that fits.

---

## For skeptics: hand-verify without trusting any automation

If you don't want to trust `verify_parity.sh`, follow the step-by-step
walkthrough in [`MANUAL_VERIFICATION.md`](MANUAL_VERIFICATION.md). It
breaks the proof into primitive commands (`git diff`, `cmp`, `diff -u`,
`shasum`, `R CMD INSTALL`) you can read and run by hand, with the expected
output for each step. Includes a step where you deliberately break a
constant and confirm the gate catches it.

## References

- [`MANUAL_VERIFICATION.md`](MANUAL_VERIFICATION.md) — step-by-step
  walkthrough for hand-checking the parity claim without trusting any
  automation.
- [`REFERENCE_VALUES_CATALOG.md`](REFERENCE_VALUES_CATALOG.md) — what each
  of the 13 cases exercises.
- [`../tests/parity/verify_parity.sh`](../tests/parity/verify_parity.sh) —
  the executable that produces the report above.
- [`../tests/parity/generate_reference.R`](../tests/parity/generate_reference.R) —
  the 13 numerical-output input cases.
- [`../tests/parity/snapshot_api.R`](../tests/parity/snapshot_api.R) —
  the API surface dumper.
- [`../tests/parity/reference_values/`](../tests/parity/reference_values) —
  the 13 frozen numerical-output JSON files.
- [`../tests/parity/api_surface/exports.txt`](../tests/parity/api_surface/exports.txt) —
  the frozen R API surface.
- [`../tests/parity/last_verification/`](../tests/parity/last_verification) —
  the most recent successful proof artifacts (log, empty diffs, sha256
  manifest). Updated by `--save-artifacts` on each refactor commit that
  needs to ship visible evidence.
- [`../Dockerfile`](../Dockerfile) — the pinned R/C++ toolchain.
