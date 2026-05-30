# Reference Values Catalog

This document describes the 13 frozen reference cases stored under
`tests/parity/reference_values/`. They are the bit-level contract for "what
R+C++ produces today." They are produced by
[`tests/parity/generate_reference.R`](../tests/parity/generate_reference.R)
running inside the pinned `r-base:4.5.3` Docker image, and serialized to JSON
at **17 significant digits** (full IEEE-754 double precision) via the `jsonlite`
package.

[`tests/parity/verify_parity.sh`](../tests/parity/verify_parity.sh) regenerates
these files from the current source and asserts byte-for-byte equality with
the committed copies. The acceptance criterion for any code change that
touches the math layer is: **this script returns exit code 0**.

See [`QA_R_BEHAVIOR_UNCHANGED.md`](QA_R_BEHAVIOR_UNCHANGED.md) for the proof
methodology and current pass/fail status.

---

## File listing

| # | file | size | exercises |
|---:|---|---:|---|
| 1 | `adult_male_baseline.json` | 48,254 B | adult RK4, default `EI`, steady state, male |
| 2 | `adult_male_ei2400.json` | 48,841 B | adult RK4, explicit `EI=2400`, male |
| 3 | `adult_female_baseline.json` | 43,145 B | adult RK4, default `EI`, steady state, female (sex parameter path) |
| 4 | `adult_male_deficit.json` | 66,261 B | adult RK4, `EIchange=-500 kcal/day` over 365 days, weight-loss trajectory, male |
| 5 | `adult_multi.json` | 89,825 B | adult RK4, two individuals in one call (multi-individual vectorization, mixed sex) |
| 6 | `child_male_365d.json` | 29,038 B | child RK4, age 6, male, default intake, 1-year horizon |
| 7 | `child_female_730d.json` | 58,061 B | child RK4, age 10, female, default intake, 2-year horizon |
| 8 | `child_male_5yr.json` | 146,044 B | child RK4, age 6, male, 5-year horizon (1825 integration steps — longest stress test) |
| 9 | `energy_linear.json` | 88 B | `EnergyBuilder` linear interpolation |
| 10 | `energy_exponential.json` | 233 B | `EnergyBuilder` exponential interpolation |
| 11 | `energy_logarithmic.json` | 204 B | `EnergyBuilder` logarithmic interpolation |
| 12 | `energy_stepwise_l.json` | 92 B | `EnergyBuilder` step-left interpolation |
| 13 | `energy_stepwise_r.json` | 92 B | `EnergyBuilder` step-right interpolation |

Total: 13 files, ~530 KB. The exact inputs that produced each file are
preserved in `tests/parity/generate_reference.R`.

---

## What each case exercises (in detail)

### Adult model (`adult_weight()`)

The adult model integrates four coupled ODEs (adaptive thermogenesis, ECF,
glycogen, lean mass) via RK4 with `dt = 1 day`. The choice of cases below
exercises every code path in the adult model and every meaningful input
permutation of the public API.

**Case 1 — `adult_male_baseline.json`** — `adult_weight(76, 1.73, 36, "male")`.
Default everything: no explicit `EI`, no `EIchange`, no `NAchange`, baseline
`PAL = 1.5`. Should produce a steady-state-ish trajectory. Tests the
`EI`-from-RMR×PAL default path and the all-zero `EIchange`/`NAchange` default
matrices.

**Case 2 — `adult_male_ei2400.json`** — `adult_weight(76, 1.73, 36, "male", EI = 2400)`.
Same subject, but with explicit baseline energy intake. Tests the explicit-`EI`
branch of the Adult constructor.

**Case 3 — `adult_female_baseline.json`** — `adult_weight(65, 1.65, 30, "female")`.
Female sex path through `estimate_rmr`, `estimate_fat`, `estimate_ecf` (each
of which has separate male/female formulas), and through the
sex-weighted blend in the RK4.

**Case 4 — `adult_male_deficit.json`** — `adult_weight(90, 1.80, 40, "male", EIchange = matrix(-500, 1, 365), NAchange = matrix(0, 1, 365))`.
Time-varying `EIchange` (constant −500 kcal/day over a year). Exercises the
day-indexed `EIchange[idx, :]` lookup at every RK4 step. Produces a
weight-loss trajectory whose endpoint is sensitive to any drift in the
lean-mass / fat-mass partitioning.

**Case 5 — `adult_multi.json`** — `adult_weight(c(76, 65), c(1.73, 1.65), c(36, 30), c("male", "female"))`.
Two individuals (one male, one female) in a single call. Exercises the cohort
broadcast paths and verifies that per-individual results don't cross-contaminate.

### Child model (`child_weight()`)

The child model integrates two coupled ODEs (FFM, FM) via RK4 with `dt = 1 day`,
with sex- and age-dependent reference tables. Each case tests a different
horizon to surface any drift that grows with step count.

**Case 6 — `child_male_365d.json`** — `child_weight(6, "male", days = 365)`.
Age 6 male over one year. Exercises male sex path through 7 sex-blended
parameter sets (`GROWTH_*`, `EB_*`, `K_CHILD`, `DELTA_MAX`, `FFM_REF_MALE`,
`FM_REF_MALE`) and the reference-EI precomputation.

**Case 7 — `child_female_730d.json`** — `child_weight(10, "female", days = 730)`.
Age 10 female over two years. Female sex path; ages span an FFM/FM reference
table boundary (age 10 → 12) and exercises the piecewise-linear interpolation
of those tables at every RK4 stage.

**Case 8 — `child_male_5yr.json`** — `child_weight(6, "male", days = 1825)`.
Age 6 male over five years. **1,825 integration steps** — the longest stress
test. Any algorithmic drift that compounds with step count surfaces here.
Sensitive to the exact floating-point ordering inside the RK4 update.

### Energy interpolation (`EnergyBuilder`)

Tests the C++ `EnergyBuilder` directly (via `bw:::EnergyBuilder`, bypassing
the R wrapper which drops the first column). Inputs:
`energy = matrix(c(2000, 2200, 1800), nrow = 1)`, `time = c(0, 5, 10)`.

**Cases 9–13** — one file per interpolation method: `Linear`, `Exponential`,
`Logarithmic`, `Stepwise_L`, `Stepwise_R`. Each produces an 11-element vector
(times 0..10). These are closed-form expressions — no iteration, no integration
— so they are the cleanest possible regression check on the C++ math layer.

Example content (`energy_linear.json`, full file):

```json
{"energy":[[2000,2040,2080,2120,2160,2200,2120,2040,1960,1880,1800]],"method":"Linear"}
```

---

## Why these specific cases and not others

The 13 cases were chosen to be **minimal but sufficient**: every code path
through `src/adult_weight.cpp`, `src/child_weight.cpp`, and `src/energy_build.cpp`
is exercised by at least one case, and the long-horizon child run gives
sub-ULP sensitivity to any change in operation order inside the RK4 step.

Adding more cases would make the suite slower without improving discriminating
power. Removing cases would create blind spots — for example, deleting Case 4
(`adult_male_deficit`) would silently allow any change to the `EIchange`
day-indexing logic to pass the gate.

The cases are not a random sample; they are a **proof set**. Each one earns
its place by exercising a path that no other case covers.

---

## How to add a new case

If a new code path is introduced (e.g. a new public function in `R/`), add a
case here:

1. Edit `tests/parity/generate_reference.R` to add a block that calls the new
   function with fixed inputs and writes its output to a new JSON file under
   `outdir`. Use `digits = 17` for `write_json`.
2. Run `./tests/parity/verify_parity.sh --update` to materialise the new JSON.
3. Inspect the new file. Spot-check that the numbers are sensible.
4. Commit the new JSON plus the script change to the same PR.
5. Add a row to the file listing above and a description to the appropriate
   section.

From that point on, the new case is part of the bit-level contract, and any
future change to the path it exercises is gated by `verify_parity.sh`.

---

## How to retire a case

Don't, unless the public API surface that the case exercises is itself
removed. The cases are cheap (~530 KB total, runs in well under a minute);
preserving them costs nothing, removing them creates blind spots.

If a public function is genuinely removed, delete the corresponding case in
the same commit as the function removal, and update the listing above.
