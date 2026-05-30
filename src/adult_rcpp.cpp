// Thin Rcpp marshalling layer for the adult body-weight model.
//
// All math lives in src/kernel/adult.cpp; this file only:
//   * converts NumericVector / NumericMatrix inputs into std::vector<double>
//     and bw::Matrix (with the transposed [day x individual] layout the
//     kernel expects — same convention the original code used),
//   * runs the model,
//   * marshals the AdultResult back into a Named Rcpp::List with the
//     EXACT shape (NumericMatrix sized nind x (nsims+1), StringMatrix for
//     BMI_Category, scalar Correct_Values and Model_Type) the R caller
//     expects.
//
// The [[Rcpp::export]] signatures here MUST byte-match the originals so
// RcppExports.cpp does not regenerate.

#include <Rcpp.h>

#include "bw/adult.hpp"

#include <cstddef>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

inline std::vector<double> to_std(const NumericVector& v) {
    return std::vector<double>(v.begin(), v.end());
}

// Convert an Rcpp NumericMatrix into a row-major bw::Matrix preserving
// (nrow, ncol) and value layout — so M_rcpp(i, j) == M_bw(i, j).
inline bw::Matrix to_std(const NumericMatrix& m) {
    bw::Matrix out(static_cast<std::size_t>(m.nrow()),
                   static_cast<std::size_t>(m.ncol()));
    for (int i = 0; i < m.nrow(); ++i) {
        for (int j = 0; j < m.ncol(); ++j) {
            out(static_cast<std::size_t>(i),
                static_cast<std::size_t>(j)) = m(i, j);
        }
    }
    return out;
}

inline NumericMatrix to_rcpp(const bw::Matrix& m) {
    NumericMatrix out(static_cast<int>(m.nrow), static_cast<int>(m.ncol));
    for (std::size_t i = 0; i < m.nrow; ++i) {
        for (std::size_t j = 0; j < m.ncol; ++j) {
            out(static_cast<int>(i), static_cast<int>(j)) = m(i, j);
        }
    }
    return out;
}

inline StringMatrix to_rcpp(const bw::StringMatrix& m) {
    StringMatrix out(static_cast<int>(m.nrow), static_cast<int>(m.ncol));
    for (std::size_t i = 0; i < m.nrow; ++i) {
        for (std::size_t j = 0; j < m.ncol; ++j) {
            out(static_cast<int>(i), static_cast<int>(j)) = m(i, j);
        }
    }
    return out;
}

inline List to_list(const bw::AdultResult& r) {
    return List::create(
        Named("Time")                   = NumericVector(r.Time.begin(), r.Time.end()),
        Named("Age")                    = to_rcpp(r.Age),
        Named("Adaptive_Thermogenesis") = to_rcpp(r.Adaptive_Thermogenesis),
        Named("Extracellular_Fluid")    = to_rcpp(r.Extracellular_Fluid),
        Named("Glycogen")               = to_rcpp(r.Glycogen),
        Named("Fat_Mass")               = to_rcpp(r.Fat_Mass),
        Named("Lean_Mass")              = to_rcpp(r.Lean_Mass),
        Named("Body_Weight")            = to_rcpp(r.Body_Weight),
        Named("Body_Mass_Index")        = to_rcpp(r.Body_Mass_Index),
        Named("BMI_Category")           = to_rcpp(r.BMI_Category),
        Named("Energy_Intake")          = to_rcpp(r.Energy_Intake),
        Named("Correct_Values")         = r.Correct_Values,
        Named("Model_Type")             = r.Model_Type);
}

} // anonymous

// [[Rcpp::export]]
List adult_weight_wrapper(NumericVector bw, NumericVector ht, NumericVector age,
                          NumericVector sex, NumericMatrix EIchange,
                          NumericMatrix NAchange, NumericVector PAL,
                          NumericVector pcarb_base, NumericVector pcarb, double dt,
                          double days, bool checkValues) {
    bw::Adult Person(to_std(bw), to_std(ht), to_std(age), to_std(sex),
                     to_std(EIchange), to_std(NAchange),
                     to_std(PAL), to_std(pcarb), to_std(pcarb_base),
                     dt, checkValues);
    return to_list(Person.rk4(days));
}

// [[Rcpp::export]]
List adult_weight_wrapper_EI(NumericVector bw, NumericVector ht, NumericVector age,
                             NumericVector sex, NumericMatrix EIchange,
                             NumericMatrix NAchange, NumericVector PAL,
                             NumericVector pcarb_base, NumericVector pcarb, double dt,
                             NumericVector extradata, double days,
                             bool checkValues, bool isEnergy) {
    bw::Adult Person(to_std(bw), to_std(ht), to_std(age), to_std(sex),
                     to_std(EIchange), to_std(NAchange),
                     to_std(PAL), to_std(pcarb), to_std(pcarb_base),
                     dt, to_std(extradata), checkValues, isEnergy);
    return to_list(Person.rk4(days));
}

// [[Rcpp::export]]
List adult_weight_wrapper_EI_fat(NumericVector bw, NumericVector ht, NumericVector age,
                                 NumericVector sex, NumericMatrix EIchange,
                                 NumericMatrix NAchange, NumericVector PAL,
                                 NumericVector pcarb_base, NumericVector pcarb, double dt,
                                 NumericVector input_EI, NumericVector input_fat,
                                 double days, bool checkValues) {
    bw::Adult Person(to_std(bw), to_std(ht), to_std(age), to_std(sex),
                     to_std(EIchange), to_std(NAchange),
                     to_std(PAL), to_std(pcarb), to_std(pcarb_base),
                     dt, to_std(input_EI), to_std(input_fat), checkValues);
    return to_list(Person.rk4(days));
}
