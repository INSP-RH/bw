//
//  energy_rcpp.cpp
//
//  Thin Rcpp marshalling layer around bw::energy::build_deterministic().
//  The Brownian-bridge path is kept here because it consumes R's RNG via
//  Rcpp::rnorm(); the pure C++ kernel handles only the deterministic
//  methods (Linear / Exponential / Logarithmic / Stepwise_L / Stepwise_R).
//
//  IMPORTANT: the exported function name and signature MUST stay byte-
//  identical to the original so RcppExports.cpp regenerates unchanged.
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
//----------------------------------------------------------------------------------------

#include <Rcpp.h>
#include <string>

#include "bw/energy.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix EnergyBuilder(NumericMatrix Energy, NumericVector Time,
                            std::string interpol) {

  //Number of times to calculate
  int days = floor(Time(Time.size()-1));

  //Numeric matrix to return
  NumericMatrix Evalues(Energy.nrow(), days + 1);

  //Brownian bridge — keeps R's RNG, so it stays in the Rcpp layer.
  if (interpol.compare("Brownian") == 0) {

    for (int j = 0; j < (Time.size()-1); j++){

      //Get times
      double T = Time(j+1);
      double t = Time(j);

      //Simulate W brownian path
      NumericMatrix W(Energy.nrow(), (T - t) + 1); //By default W(_, 0) = 0;
      for (int i = 1; i < (T - t + 1); i++){
        W(_, i) = W(_,i-1) + rnorm(Energy.nrow());
      }

      //Get brownian bridge
      for (int i = 0 ; i < (T - t + 1); i++){
        Evalues(_,i + t) = Energy(_,j)*( (T - t) - i )/(T - t) + Energy(_,j+1)*i/(T-t) +
          W(_, i) -  (i/(T-t))*W(_, (T-t));
      }
    }

    return Evalues;
  }

  //Deterministic methods — dispatch into the pure C++ kernel.
  bw::energy::Method method;
  if (!bw::energy::parse_method(interpol.c_str(), &method)) {
    // Unknown interpolation: preserve original behavior (no branch matched →
    // Evalues filled with zeros for the first `days` cols, last col set below).
    Evalues(_, Evalues.ncol() - 1) = Energy(_, Energy.ncol() - 1);
    return Evalues;
  }

  bw::energy::build_deterministic(
      Energy.begin(),
      static_cast<std::size_t>(Energy.nrow()),
      static_cast<std::size_t>(Time.size()),
      Time.begin(),
      method,
      Evalues.begin());

  return Evalues;
}
