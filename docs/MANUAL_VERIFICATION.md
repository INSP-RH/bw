# Manual Verification — Hand-Check the "No R Behavior Change" Claim

This is the **skeptic's walkthrough**. The PR claims R-visible behavior is
byte-identical before and after the kernel-extraction refactor. The
`tests/parity/verify_parity.sh` script automates the proof, but if you don't
want to take the automation's word for it, this document walks you through
every check by hand, command by command, with the output you should see at
each step.

Estimated time: ~15 minutes once Docker has built the image (which it does
itself, once, on first call — about 3 minutes the first time).

Prerequisites:
- Docker installed and running.
- A checkout of this PR's branch (`refactor/extract-cpp-kernel`).
- Optional: R 4.5+ installed locally if you want to verify outside Docker
  too. Not required for any of the checks below.

> **What you're proving.** Two contracts: (1) every exported R function
> produces the same numerical output for the same inputs, bit-for-bit; (2)
> every exported R function has the same name, argument list, and default
> values. If both hold, no R caller can observe the refactor.

---

## Step 1 — The PR touches no R-facing files

Before running any code, prove the surface change is limited. From the repo root:

```bash
git diff --stat upstream/master..HEAD -- 'R/**' 'man/**' DESCRIPTION NAMESPACE
```

**Expected output:** empty (no diff). Confirm:

```bash
git diff --stat upstream/master..HEAD -- 'R/**' 'man/**' DESCRIPTION NAMESPACE | wc -l
```

**Expected:** `0`.

This means every exported function lives in the same R source file, with the same signature and the same roxygen documentation as before. The Rcpp dispatch table (`src/RcppExports.cpp`) should also be unchanged:

```bash
git diff --stat upstream/master..HEAD -- src/RcppExports.cpp
```

**Expected:** empty.

Now look at the files this PR DOES change:

```bash
git diff --stat upstream/master..HEAD
```

**Expected:** ~24 files. All additions or deletions inside `src/` (kernel restructure), plus new `tests/parity/` content and `docs/`. Nothing else touched.

---

## Step 2 — Inspect the new src/ layout

Confirm the math has been moved into a pure-C++ subtree with no Rcpp dependency. Look at one kernel file:

```bash
head -30 src/kernel/energy.cpp
```

**Expected:** the file's `#include` block contains only standard-library headers (`<cstddef>`, `<cmath>`, etc.) and `bw/energy.hpp`. **It does not include `<Rcpp.h>`.** Confirm:

```bash
grep -c 'Rcpp.h' src/kernel/*.cpp src/include/bw/*.hpp
```

**Expected:** `0` for every file. Each file should report `0`.

Now confirm the Rcpp shim files DO include Rcpp:

```bash
grep -l 'Rcpp.h' src/*_rcpp.cpp
```

**Expected:** all three (`src/adult_rcpp.cpp`, `src/child_rcpp.cpp`, `src/energy_rcpp.cpp`).

Read one of them to see how thin it is — for example:

```bash
wc -l src/energy_rcpp.cpp src/kernel/energy.cpp src/include/bw/energy.hpp
```

**Expected:** the `_rcpp.cpp` shim is short (~75 lines); the kernel is the bulk of the code; the header declares the public surface. Together they replace the original ~112-line `energy_build.cpp`.

---

## Step 3 — Build the R package the normal way

Confirm `R CMD INSTALL` still works. Inside the pinned Docker environment so you don't need R locally:

```bash
docker build -t bw-manual-verify .
docker run --rm -v "$(pwd):/pkg:ro" -w /tmp bw-manual-verify \
    sh -c "cp -r /pkg /tmp/bw && R CMD INSTALL --no-docs /tmp/bw"
```

**Expected:** the install completes with `* DONE (bw)` at the end. No errors. The `kernel/*.cpp` files get compiled alongside the top-level `*_rcpp.cpp` shims and linked into `bw.so`.

If you have R locally and prefer to install on the host:

```bash
R CMD INSTALL .
```

Same expected result.

---

## Step 4 — Run the existing testthat suite

Confirm the R package's own test suite still passes (this is upstream's own
verification of correctness; it should be unchanged by this PR):

```bash
docker run --rm bw-manual-verify Rscript -e "testthat::test_local()"
```

**Expected:** all tests pass with no failures.

This proves the refactor doesn't break R-side behavior the upstream
maintainer already cared about enough to write tests for.

---

## Step 5 — Hand-check a single numerical output

The bit-level claim. Pick the simplest reference case — `adult_male_baseline`
— compute it yourself and compare against the committed JSON.

In a fresh Docker R session:

```bash
docker run --rm bw-manual-verify Rscript -e '
  library(bw)
  library(jsonlite)
  result <- adult_weight(76, 1.73, 36, "male")
  cat(toJSON(result$Body_Weight[1, 1:5], digits = 17), "\n")
'
```

**Expected:** prints something like:

```
[76,76.0001234567,76.0002...]   (first 5 values of Body_Weight over time)
```

The exact 17-digit values aren't worth eyeballing alone. Instead, compare
against what's in the committed reference:

```bash
docker run --rm bw-manual-verify Rscript -e '
  ref <- jsonlite::fromJSON("/pkg/tests/parity/reference_values/adult_male_baseline.json")
  cat(jsonlite::toJSON(ref$Body_Weight[1, 1:5], digits = 17), "\n")
' -e ":" --args _   # workaround: mount the pkg
```

Or more straightforwardly, use `cmp` to compare on-disk:

```bash
# Generate a fresh adult_male_baseline.json from current source, byte-diff it
docker run --rm -v "$(pwd):/pkg:ro" -v "/tmp:/out" -w /pkg bw-manual-verify sh -c '
  R CMD INSTALL --no-docs . >/dev/null 2>&1 &&
  Rscript -e "
    library(bw); library(jsonlite)
    r <- adult_weight(76, 1.73, 36, \"male\")
    s <- lapply(r, function(x) if (is.matrix(x)) lapply(seq_len(nrow(x)), function(i) as.numeric(x[i,])) else if (is.character(x) && length(x)==1) x else if (is.logical(x) && length(x)==1) x else as.numeric(x))
    write_json(s, \"/out/fresh.json\", auto_unbox=TRUE, digits=17)
  "
'
cmp tests/parity/reference_values/adult_male_baseline.json /tmp/fresh.json && echo "BYTE-IDENTICAL" || echo "DIFFER"
```

**Expected:** `BYTE-IDENTICAL`.

The committed reference was generated by this same procedure from the
pre-refactor source. Byte-identity here means the refactor produced exactly
the same numerical output the original C++ did, down to the last bit of
every double-precision value.

If you want to repeat for the other 12 cases, the recipes are in
[`tests/parity/generate_reference.R`](../tests/parity/generate_reference.R).

---

## Step 6 — Hand-check the R API surface

Confirm that no exported function's name, argument list, or defaults
changed. Get the current API:

```bash
docker run --rm bw-manual-verify Rscript -e '
  ns <- asNamespace("bw")
  for (name in sort(getNamespaceExports("bw"))) {
    fn <- get(name, envir = ns)
    if (is.function(fn)) {
      a <- args(fn)
      cat(name, "  ", deparse(a, width.cutoff = 500L)[1], "\n", sep = "")
    }
  }
'
```

Open `tests/parity/api_surface/exports.txt` in your editor and compare
visually against the output above. Each function and its full argument list
(with defaults) should match exactly.

For a stricter check, regenerate the snapshot fresh and diff:

```bash
docker run --rm -v "$(pwd):/pkg:ro" -v "/tmp:/out" -w /pkg bw-manual-verify sh -c '
  R CMD INSTALL --no-docs . >/dev/null 2>&1 &&
  Rscript -e "outdir <- \"/out\"; source(\"/pkg/tests/parity/snapshot_api.R\")"
'
diff -u tests/parity/api_surface/exports.txt /tmp/exports.txt && echo "API SURFACE IDENTICAL"
```

**Expected:** `API SURFACE IDENTICAL`, no diff output.

This proves the package has the same 8 exports, every function has the same
signature, every default value is identical, and no S3 / S4 methods were
added or removed.

---

## Step 7 — Run the full automated gate (now that you trust it)

Having done steps 1–6 by hand, you can now run the automated gate as a
sanity convenience:

```bash
./tests/parity/verify_parity.sh
```

**Expected last lines:**

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
```

Exit code: `0`.

If you do not trust this script, you can read it in full at
[`tests/parity/verify_parity.sh`](../tests/parity/verify_parity.sh) — it's
~250 lines of bash, mostly `cmp` and `diff` invocations against the same
files you compared by hand in steps 5 and 6.

---

## Step 8 — Verify the committed proof artifacts match reality

The PR also commits artifacts under `tests/parity/last_verification/` as
human-readable proof. Confirm these are not lying.

The two diff files should be 0 bytes:

```bash
ls -la tests/parity/last_verification/{outputs,api_surface}.diff
```

**Expected:** both files exist and are `0 B`. An empty diff IS the proof of
no difference — `diff -u a b` produces zero output when `a` and `b` are
identical.

The manifest contains sha256 hashes of every contract file as of the moment
the artifacts were generated. Verify them now:

```bash
cd tests/parity/reference_values && shasum -a 256 *.json | sort
```

Compare against `numerical outputs:` section of
`tests/parity/last_verification/manifest.txt`. Every hash should match.

Same for the API surface:

```bash
shasum -a 256 tests/parity/api_surface/exports.txt
```

Compare against the `R API surface:` section of `manifest.txt`.

If all hashes match, the committed proof artifacts accurately describe the
current state of the contract files. Anyone reviewing the PR on GitHub can
see these artifacts without running anything.

---

## Step 9 — Prove the gate actually fails when it should

You should also confirm the gate is not just rubber-stamping. Deliberately
break byte-identity and verify the gate fails:

```bash
# Make a real numerical change to the kernel: bump a constant slightly
sed -i.bak 's/roF = 9440.727/roF = 9440.728/' src/include/bw/adult.hpp 2>/dev/null \
  || sed -i.bak 's/roF = 9440.727/roF = 9440.728/' src/kernel/adult.cpp
./tests/parity/verify_parity.sh
echo "Exit code: $?"
```

**Expected:** the gate returns exit 1, lists which reference cases drifted
(adult cases will fail because `roF` is in the fat-mass formula), and saves
a full diff to `r_parity_diff.log`.

Now restore:

```bash
mv src/include/bw/adult.hpp.bak src/include/bw/adult.hpp 2>/dev/null \
  || mv src/kernel/adult.cpp.bak src/kernel/adult.cpp
./tests/parity/verify_parity.sh
```

**Expected:** PASS, exit 0.

This proves the gate is sensitive — a one-millionth-of-a-percent change to
a physical constant trips it. The gate is real, not theater.

---

## Step 10 — Spot-check the kernel itself

If you want to read the math directly, the kernel is concentrated in three
files:

```bash
wc -l src/kernel/*.cpp
```

**Expected:**

```
 131 src/kernel/energy.cpp
 490 src/kernel/child.cpp
 658 src/kernel/adult.cpp
1279 total
```

Compare with the pre-refactor (run `git show upstream/master:src/`):

```bash
git show upstream/master -- 'src/*.cpp' 'src/*.h' | grep -c '^+++ b/src/'
```

Each kernel file mirrors the original `*.cpp` / `*.h` of the same module —
same constants, same loops, same operation order. Types changed (Rcpp
vectors → `const double*` + `std::size_t`); arithmetic did not. You can
diff a function-by-function (e.g. `Adult::dL` from the original vs the
kernel's `dL`) — they should be the same expression with different types
on the parameters.

---

## What success looks like, all together

If steps 1–9 all returned the expected outputs, you have hand-verified:

- (Step 1) The PR touches no R-package files: no `R/`, no `man/`, no
  `DESCRIPTION`, no `NAMESPACE`, no `RcppExports.cpp`.
- (Step 2) The new C++ kernel under `src/kernel/` contains no Rcpp
  dependency.
- (Step 3) The R package builds and installs the same as before.
- (Step 4) The existing testthat suite passes.
- (Step 5) The numerical output of at least one reference case is bit-for-bit
  identical to the committed reference. (Repeat for any other cases if you
  want — the script does all 13.)
- (Step 6) The R API surface is character-for-character identical to the
  committed snapshot.
- (Step 7) The automated gate agrees: PASS, exit 0.
- (Step 8) The committed proof artifacts (empty diffs + sha256 manifest)
  accurately reflect the current state.
- (Step 9) The gate is not rubber-stamping — it fails when you deliberately
  break a constant.

At this point, you've personally confirmed the PR's correctness claim
without relying on a single piece of automation you didn't read.

---

## If anything diverges from the expected output

Stop and inspect. The most likely causes, in decreasing order:

1. **You're not on the branch `refactor/extract-cpp-kernel`.** Run
   `git branch --show-current` to confirm.
2. **Docker isn't running** — `verify_parity.sh` exits with code 2 in that case.
3. **You modified something between steps and forgot.** Run
   `git status --short`; if anything's modified, reset and try again.
4. **You're on a different platform with a different `libm`.** The contract
   pins the `r-base:4.5.3` Docker image specifically because `libm`'s
   handling of `exp`/`pow`/`log` is not bit-stable across implementations.
   If you're running the R commands OUTSIDE Docker on a different OS, you
   may legitimately see drift in the last few bits of long-horizon
   integrations (child 5-year run especially). Stay inside Docker for
   bit-exact comparison.

Report any other divergence as a PR comment — it's either a real bug in the
refactor (in which case the gate caught it) or a documentation gap (in
which case this file needs updating).
