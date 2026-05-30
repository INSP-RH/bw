#!/usr/bin/env Rscript
# Generate R reference values for the parity gate.
#
# Run via:  ./tests/parity/verify_parity.sh   (which calls this inside Docker)
#
# Produces 13 JSON files under tests/parity/reference_values/ containing the
# numerical output of adult_weight(), child_weight(), and the C++
# EnergyBuilder at 17 significant digits (full IEEE-754 double precision).
# These files are the byte-level contract for "what R+C++ produces today";
# see docs/QA_R_BEHAVIOR_UNCHANGED.md.
#
# Requires: the bw package installed, plus jsonlite.

library(bw)
library(jsonlite)

# Allow the caller (verify_parity.sh) to override the output directory.
if (!exists("outdir")) outdir <- "tests/parity/reference_values"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Serialize a model result to a JSON-friendly structure.
# Matrices become lists of row vectors (one row per individual); scalars and
# vectors pass through unchanged.
serialize_model <- function(model) {
  lapply(model, function(x) {
    if (is.matrix(x)) {
      lapply(seq_len(nrow(x)), function(i) as.numeric(x[i, ]))
    } else if (is.character(x) && length(x) == 1) {
      x
    } else if (is.logical(x) && length(x) == 1) {
      x
    } else {
      as.numeric(x)
    }
  })
}

cat("Generating adult model references...\n")

# Case 1 — male, baseline steady state
r1 <- adult_weight(76, 1.73, 36, "male")
write_json(serialize_model(r1),
           file.path(outdir, "adult_male_baseline.json"),
           auto_unbox = TRUE, digits = 17)

# Case 2 — male with explicit EI
r2 <- adult_weight(76, 1.73, 36, "male", EI = rep(2400, 1))
write_json(serialize_model(r2),
           file.path(outdir, "adult_male_ei2400.json"),
           auto_unbox = TRUE, digits = 17)

# Case 3 — female, baseline
r3 <- adult_weight(65, 1.65, 30, "female")
write_json(serialize_model(r3),
           file.path(outdir, "adult_female_baseline.json"),
           auto_unbox = TRUE, digits = 17)

# Case 4 — male with negative energy change (weight loss)
ei_change <- matrix(rep(-500, 365), nrow = 1, ncol = 365)
na_change <- matrix(rep(0, 365),   nrow = 1, ncol = 365)
r4 <- adult_weight(90, 1.80, 40, "male", EIchange = ei_change, NAchange = na_change)
write_json(serialize_model(r4),
           file.path(outdir, "adult_male_deficit.json"),
           auto_unbox = TRUE, digits = 17)

# Case 5 — multiple individuals, mixed sex
r5 <- adult_weight(c(76, 65), c(1.73, 1.65), c(36, 30), c("male", "female"))
write_json(serialize_model(r5),
           file.path(outdir, "adult_multi.json"),
           auto_unbox = TRUE, digits = 17)

cat("Generating child model references...\n")

# Case 6 — male child, age 6, default intake, 1 year
r6 <- child_weight(6, "male", days = 365)
write_json(serialize_model(r6),
           file.path(outdir, "child_male_365d.json"),
           auto_unbox = TRUE, digits = 17)

# Case 7 — female child, age 10, default intake, 2 years
r7 <- child_weight(10, "female", days = 730)
write_json(serialize_model(r7),
           file.path(outdir, "child_female_730d.json"),
           auto_unbox = TRUE, digits = 17)

# Case 8 — male child, 5 years (1825 RK4 steps, longest stress test)
r8 <- child_weight(6, "male", days = 1825)
write_json(serialize_model(r8),
           file.path(outdir, "child_male_5yr.json"),
           auto_unbox = TRUE, digits = 17)

cat("Generating energy interpolation references...\n")

# Cases 9-13 — all 5 deterministic methods.
# Use bw:::EnergyBuilder directly (not the R wrapper, which drops the first
# column), so the references exercise the C++ kernel as-is.
energy   <- matrix(c(2000, 2200, 1800), nrow = 1)
time_vec <- c(0, 5, 10)
for (method in c("Linear", "Exponential", "Logarithmic", "Stepwise_L", "Stepwise_R")) {
  result <- bw:::EnergyBuilder(energy, time_vec, method)
  write_json(list(energy = list(as.numeric(result[1, ])),
                  method = method),
             file.path(outdir, paste0("energy_", tolower(method), ".json")),
             auto_unbox = TRUE, digits = 17)
}

cat(sprintf("Done. %d reference files written to %s/\n",
            length(list.files(outdir, pattern = "\\.json$")), outdir))
