//
//  child_rcpp.cpp — thin Rcpp wrapper around the pure-C++ bw::Child kernel.
//
//  Preserves the exact [[Rcpp::export]] signatures from the original
//  child_weight_wrapper.cpp so RcppExports.cpp regenerates byte-identical.
//
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
//----------------------------------------------------------------------------------------

#include <Rcpp.h>
#include "bw/child.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

using namespace Rcpp;

// Convert a NumericVector to std::vector<double> (independent of Rcpp memory).
static std::vector<double> to_std(const NumericVector& v) {
    return std::vector<double>(v.begin(), v.end());
}

// Convert NumericMatrix (column-major) to row-major std::vector<double>.
// out[r*ncol + c] = mat(r, c)
static std::vector<double> matrix_to_row_major(const NumericMatrix& mat) {
    std::size_t nr = static_cast<std::size_t>(mat.nrow());
    std::size_t nc = static_cast<std::size_t>(mat.ncol());
    std::vector<double> out(nr * nc);
    for (std::size_t c = 0; c < nc; ++c) {
        for (std::size_t r = 0; r < nr; ++r) {
            out[r*nc + c] = mat(r, c);
        }
    }
    return out;
}

// Convert a kernel result (row-major-by-column-major: data[i + nind*j] for
// matrix(i,j)) back into an Rcpp NumericMatrix with the same shape.
// Layout in the kernel matches NumericMatrix's column-major (i + nind*j),
// so we can just copy byte-equivalent values into the NumericMatrix.
static NumericMatrix make_matrix(const std::vector<double>& data,
                                 std::size_t nind, std::size_t nsteps) {
    NumericMatrix mat(nind, nsteps);
    for (std::size_t j = 0; j < nsteps; ++j) {
        for (std::size_t i = 0; i < nind; ++i) {
            mat(i, j) = data[i + nind*j];
        }
    }
    return mat;
}

// [[Rcpp::export]]
List child_weight_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, NumericMatrix input_EIntake, double days, double dt, bool checkValues){

    std::vector<double> age_v = to_std(age);
    std::vector<double> sex_v = to_std(sex);
    std::vector<double> ffm_v = to_std(FFM);
    std::vector<double> fm_v  = to_std(FM);
    std::vector<double> ei_v  = matrix_to_row_major(input_EIntake);

    bw::Child Person(age_v, sex_v, ffm_v, fm_v,
                     ei_v,
                     static_cast<std::size_t>(input_EIntake.nrow()),
                     static_cast<std::size_t>(input_EIntake.ncol()),
                     dt, checkValues);

    bw::ChildRK4Result R = Person.rk4(days - 1); // days - 1: see original wrapper

    NumericVector TIME(R.nsteps);
    for (std::size_t j = 0; j < R.nsteps; ++j) TIME[j] = R.Time[j];

    return List::create(
        Named("Time")           = TIME,
        Named("Age")            = make_matrix(R.Age,           R.nind, R.nsteps),
        Named("Fat_Free_Mass")  = make_matrix(R.Fat_Free_Mass, R.nind, R.nsteps),
        Named("Fat_Mass")       = make_matrix(R.Fat_Mass,      R.nind, R.nsteps),
        Named("Body_Weight")    = make_matrix(R.Body_Weight,   R.nind, R.nsteps),
        Named("Correct_Values") = R.Correct_Values,
        Named("Model_Type")     = "Children"
    );
}

// [[Rcpp::export]]
List child_weight_wrapper_richardson(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, double K, double Q, double A, double B, double nu, double C, double days, double dt, bool checkValues){

    std::vector<double> age_v = to_std(age);
    std::vector<double> sex_v = to_std(sex);
    std::vector<double> ffm_v = to_std(FFM);
    std::vector<double> fm_v  = to_std(FM);

    bw::Child Person(age_v, sex_v, ffm_v, fm_v,
                     K, Q, A, B, nu, C, dt, checkValues);

    bw::ChildRK4Result R = Person.rk4(days - 1);

    NumericVector TIME(R.nsteps);
    for (std::size_t j = 0; j < R.nsteps; ++j) TIME[j] = R.Time[j];

    return List::create(
        Named("Time")           = TIME,
        Named("Age")            = make_matrix(R.Age,           R.nind, R.nsteps),
        Named("Fat_Free_Mass")  = make_matrix(R.Fat_Free_Mass, R.nind, R.nsteps),
        Named("Fat_Mass")       = make_matrix(R.Fat_Mass,      R.nind, R.nsteps),
        Named("Body_Weight")    = make_matrix(R.Body_Weight,   R.nind, R.nsteps),
        Named("Correct_Values") = R.Correct_Values,
        Named("Model_Type")     = "Children"
    );
}

// [[Rcpp::export]]
NumericMatrix intake_reference_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, double days,  double dt){

    // Empty energy intake matrix
    NumericMatrix EI(1,1);

    std::vector<double> age_v = to_std(age);
    std::vector<double> sex_v = to_std(sex);
    std::vector<double> ffm_v = to_std(FFM);
    std::vector<double> fm_v  = to_std(FM);
    std::vector<double> ei_v  = matrix_to_row_major(EI);

    bw::Child Person(age_v, sex_v, ffm_v, fm_v,
                     ei_v,
                     static_cast<std::size_t>(EI.nrow()),
                     static_cast<std::size_t>(EI.ncol()),
                     dt, false);

    std::size_t nind = static_cast<std::size_t>(age.size());
    int ncols = static_cast<int>(std::floor(days/dt) + 1);
    NumericMatrix EnergyIntake(static_cast<int>(nind), ncols);

    // Get energy matrix:
    //   for (double i = 0; i < floor(days/dt) + 1; i++)
    //       EnergyIntake(_, i) = Person.IntakeReference(age + dt*i/365.0);
    // Iterate exactly with double `i` to mirror original loop bounds.
    for (double i = 0; i < std::floor(days/dt) + 1; i += 1.0) {
        std::vector<double> tvec(nind);
        for (std::size_t k = 0; k < nind; ++k) {
            tvec[k] = age[k] + dt*i/365.0;
        }
        std::vector<double> col = Person.IntakeReference(tvec);
        int col_idx = static_cast<int>(i);
        for (std::size_t k = 0; k < nind; ++k) {
            EnergyIntake(static_cast<int>(k), col_idx) = col[k];
        }
    }

    return EnergyIntake;
}

// [[Rcpp::export]]
List mass_reference_wrapper(NumericVector age, NumericVector sex){

    NumericMatrix EI(1,1);
    NumericMatrix inputFM(1,1);
    NumericMatrix inputFFM(1,1);

    std::vector<double> age_v = to_std(age);
    std::vector<double> sex_v = to_std(sex);
    // Mirror reference: inputFFM and inputFM are 1x1 matrices in the original.
    // They are NOT used (FFMReference / FMReference don't read FFM/FM), so
    // their content is moot; we still construct the kernel object with empty
    // FFM/FM vectors that won't be touched by the called methods.
    std::vector<double> ffm_v;
    std::vector<double> fm_v;
    // size them to age.size() so getParameters' nind matches and any incidental
    // access is safe. Reference passed NumericMatrix(1,1) but the Rcpp class
    // stored them in NumericVector slots — implicit conversion of NumericMatrix
    // to NumericVector flattens column-major.
    // However the only data members the called functions use here are sex,
    // age, and the per-individual constant tables.
    ffm_v.assign(age.size(), 0.0);
    fm_v.assign(age.size(), 0.0);
    std::vector<double> ei_v  = matrix_to_row_major(EI);

    bw::Child Person(age_v, sex_v, ffm_v, fm_v,
                     ei_v,
                     static_cast<std::size_t>(EI.nrow()),
                     static_cast<std::size_t>(EI.ncol()),
                     1.0, false);

    std::vector<double> FMref  = Person.FMReference(age_v);
    std::vector<double> FFMref = Person.FFMReference(age_v);

    NumericVector FM(FMref.size());
    NumericVector FFM(FFMref.size());
    for (std::size_t k = 0; k < FMref.size(); ++k)  FM[k]  = FMref[k];
    for (std::size_t k = 0; k < FFMref.size(); ++k) FFM[k] = FFMref[k];

    return List::create(Named("FM")  = FM,
                        Named("FFM") = FFM);
}
