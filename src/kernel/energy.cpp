//
//  kernel/energy.cpp
//
//  Pure C++ implementation of the deterministic interpolation methods
//  originally in src/energy_build.cpp. Math is preserved verbatim, only
//  the types/containers were translated from Rcpp sugar to raw pointers.
//
//  Layout convention: matrices are column-major (matching R / Rcpp), so
//  element (row r, col c) of an (nrow × ncols) matrix lives at
//  buffer[r + c * nrow].
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
//----------------------------------------------------------------------------------------

#include "bw/energy.hpp"

#include <cmath>
#include <cstddef>
#include <cstring>

namespace bw {
namespace energy {

bool parse_method(const char* name, Method* out) {
  if (name == nullptr || out == nullptr) return false;
  if (std::strcmp(name, "Linear") == 0)      { *out = Method::Linear;      return true; }
  if (std::strcmp(name, "Exponential") == 0) { *out = Method::Exponential; return true; }
  if (std::strcmp(name, "Logarithmic") == 0) { *out = Method::Logarithmic; return true; }
  if (std::strcmp(name, "Stepwise_L") == 0)  { *out = Method::Stepwise_L;  return true; }
  if (std::strcmp(name, "Stepwise_R") == 0)  { *out = Method::Stepwise_R;  return true; }
  return false;
}

std::size_t output_ncols(const double* Time, std::size_t ntimes) {
  // Matches original: int days = floor(Time(Time.size()-1));
  // Output has (days + 1) columns.
  const int days = static_cast<int>(std::floor(Time[ntimes - 1]));
  return static_cast<std::size_t>(days + 1);
}

void build_deterministic(const double* Energy,
                         std::size_t nrow,
                         std::size_t ntimes,
                         const double* Time,
                         Method method,
                         double* Evalues) {

  // Number of times to calculate (matches original `int days`).
  const int days = static_cast<int>(std::floor(Time[ntimes - 1]));
  int j = 0; // indicator of time value we are taking

  // To avoid logarithm starting at 0 we displace the exponential to let for
  // a maximum y2 - y1 of 1000.
  const double K = 5000.0;
  const double log_K = std::log(K);

  // Case; exponential; logarithmic or stepwise
  for (int i = 0; i < days; i++) {

    // Per-iteration scalars used by the deterministic methods. These are
    // computed once outside the per-row loop, identically to how Rcpp sugar
    // evaluates the scalar subexpressions in the original code.
    const double dt = Time[j + 1] - Time[j];
    const double di = static_cast<double>(i) - Time[j];

    const std::size_t col_i   = static_cast<std::size_t>(i)     * nrow;
    const std::size_t col_j   = static_cast<std::size_t>(j)     * nrow;
    const std::size_t col_jp1 = static_cast<std::size_t>(j + 1) * nrow;

    switch (method) {
      case Method::Linear:
        for (std::size_t r = 0; r < nrow; ++r) {
          Evalues[r + col_i] =
              (Energy[r + col_jp1] - Energy[r + col_j]) / dt * di
              + Energy[r + col_j];
        }
        break;

      case Method::Stepwise_L:
        for (std::size_t r = 0; r < nrow; ++r) {
          Evalues[r + col_i] = Energy[r + col_j];
        }
        break;

      case Method::Stepwise_R:
        for (std::size_t r = 0; r < nrow; ++r) {
          Evalues[r + col_i] = Energy[r + col_jp1];
        }
        break;

      case Method::Exponential:
        for (std::size_t r = 0; r < nrow; ++r) {
          Evalues[r + col_i] =
              std::exp(
                (std::log(Energy[r + col_jp1] - Energy[r + col_j] + K) - log_K)
                / dt * di
                + log_K)
              - K + Energy[r + col_j];
        }
        break;

      case Method::Logarithmic:
        for (std::size_t r = 0; r < nrow; ++r) {
          Evalues[r + col_i] =
              1000.0 * std::log(
                (std::exp((Energy[r + col_jp1] - Energy[r + col_j]) / 1000.0) - 1.0)
                / dt * di
                + 1.0)
              + Energy[r + col_j];
        }
        break;
    }

    // Update to next time
    if (i + 1 >= Time[j + 1]) {
      j = j + 1;
    }
  }

  // Last day: Evalues(_, ncol-1) = Energy(_, Energy.ncol() - 1)
  const std::size_t ncols_out = output_ncols(Time, ntimes);
  const std::size_t last_out_col = (ncols_out - 1) * nrow;
  const std::size_t last_in_col  = (ntimes     - 1) * nrow;
  for (std::size_t r = 0; r < nrow; ++r) {
    Evalues[r + last_out_col] = Energy[r + last_in_col];
  }
}

} // namespace energy
} // namespace bw
