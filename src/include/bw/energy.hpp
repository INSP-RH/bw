//
//  bw/energy.hpp
//
//  Pure C++ public interface for the energy interpolation kernel.
//  Contains NO Rcpp dependency; safe to consume from any C++17 caller.
//
//  Math is preserved verbatim from the original Rcpp implementation; do not
//  reorder floating-point operations or substitute equivalent identities.
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
//----------------------------------------------------------------------------------------

#ifndef BW_ENERGY_HPP
#define BW_ENERGY_HPP

#include <cstddef>

namespace bw {
namespace energy {

enum class Method {
  Linear,
  Exponential,
  Logarithmic,
  Stepwise_L,
  Stepwise_R
};

// Parse one of "Linear", "Exponential", "Logarithmic", "Stepwise_L",
// "Stepwise_R". Returns true on match and writes the enum to *out.
// "Brownian" is intentionally NOT recognized here — it requires R RNG and
// is handled in the Rcpp wrapper.
bool parse_method(const char* name, Method* out);

// Compute the number of output columns for a given Time vector:
//   ncols_out = floor(Time[ntimes - 1]) + 1
// Caller uses this to size the destination buffer.
std::size_t output_ncols(const double* Time, std::size_t ntimes);

// Deterministic interpolation methods (Linear / Exponential / Logarithmic /
// Stepwise_L / Stepwise_R).
//
// Layout: Energy is column-major with shape (nrow, ntimes); element (r, c) is
// at Energy[r + c * nrow]. Evalues must be pre-allocated with shape
// (nrow, output_ncols(Time, ntimes)) in the same column-major layout.
//
// Math is byte-identical (per IEEE-754) to the original EnergyBuilder()
// implementation for these five methods.
void build_deterministic(const double* Energy,
                         std::size_t nrow,
                         std::size_t ntimes,
                         const double* Time,
                         Method method,
                         double* Evalues);

} // namespace energy
} // namespace bw

#endif // BW_ENERGY_HPP
