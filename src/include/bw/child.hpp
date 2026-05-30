//
//  child.hpp — pure C++ public interface for the Child model.
//
//  Translated from src/child_weight.h (Rcpp version). NO Rcpp here.
//  Arithmetic order is preserved verbatim from the reference impl.
//
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
//----------------------------------------------------------------------------------------

#ifndef BW_CHILD_HPP
#define BW_CHILD_HPP

#include <cstddef>
#include <vector>

namespace bw {

// Output of rk4(): per-individual time series.
// All matrices are stored row-major: row i (individual), col j (step).
struct ChildRK4Result {
    std::size_t nind;
    std::size_t nsteps;        // == nsims + 1
    std::vector<double> Time;  // length nsteps
    std::vector<double> Age;          // nind * nsteps
    std::vector<double> Fat_Free_Mass;// nind * nsteps
    std::vector<double> Fat_Mass;     // nind * nsteps
    std::vector<double> Body_Weight;  // nind * nsteps
    bool Correct_Values;
};

class Child {
public:
    // Default (classic) constructor — EIntake matrix path.
    // EIntake is row-major: rows are time slots, cols are individuals; nrows = ei_nrows.
    Child(const std::vector<double>& input_age,
          const std::vector<double>& input_sex,
          const std::vector<double>& input_FFM,
          const std::vector<double>& input_FM,
          const std::vector<double>& input_EIntake,
          std::size_t ei_nrows,
          std::size_t ei_ncols,
          double input_dt,
          bool checkValues);

    // Richardson (generalized logistic) constructor.
    Child(const std::vector<double>& input_age,
          const std::vector<double>& input_sex,
          const std::vector<double>& input_FFM,
          const std::vector<double>& input_FM,
          double input_K, double input_Q, double input_A,
          double input_B, double input_nu, double input_C,
          double input_dt,
          bool checkValues);

    ~Child();

    // RK4 over `days`. Result laid out row-major: matrix[i, j] -> data[i + nind*j]
    // to match Rcpp's column-major NumericMatrix(nind, nsims+1).
    ChildRK4Result rk4(double days);

    // Reference helpers (callable directly).
    std::vector<double> IntakeReference(const std::vector<double>& t);
    std::vector<double> FFMReference(const std::vector<double>& t);
    std::vector<double> FMReference(const std::vector<double>& t);

    std::size_t n_individuals() const { return nind; }

private:
    // Inputs / configuration
    std::vector<double> age;
    std::vector<double> sex;
    std::vector<double> FFM;
    std::vector<double> FM;
    std::vector<double> EIntake; // (ei_nrows x ei_ncols), row-major: r*ei_ncols + c
    std::size_t ei_nrows;
    std::size_t ei_ncols;
    bool check;

    // Private unchanging constants
    double rhoFM;
    double deltamin;
    double P;
    double h;
    double dt;
    bool generalized_logistic;

    std::size_t nind;

    // Per-individual constants
    std::vector<double> K;
    std::vector<double> deltamax;

    // Robinson's curve constants
    double K_logistic;
    double Q_logistic;
    double A_logistic;
    double B_logistic;
    double nu_logistic;
    double C_logistic;

    // g, DYNAMICS PAPER
    std::vector<double> A, tA, tauA;
    std::vector<double> B, tB, tauB;
    std::vector<double> D, tD, tauD;

    // g, IMPACT PAPER
    std::vector<double> A1, tA1, tauA1;
    std::vector<double> B1, tB1, tauB1;
    std::vector<double> D1, tD1, tauD1;

    // EB, IMPACT PAPER
    std::vector<double> A_EB, tA_EB, tauA_EB;
    std::vector<double> B_EB, tB_EB, tauB_EB;
    std::vector<double> D_EB, tD_EB, tauD_EB;

    // Linear regression coefficients (unused in math, but kept for parity)
    std::vector<double> ffm_beta0, ffm_beta1;
    std::vector<double> fm_beta0,  fm_beta1;

    // Private functions
    void build();
    void getParameters();

    std::vector<double> general_ode(const std::vector<double>& t,
                                    const std::vector<double>& input_A,
                                    const std::vector<double>& input_B,
                                    const std::vector<double>& input_D,
                                    const std::vector<double>& input_tA,
                                    const std::vector<double>& input_tB,
                                    const std::vector<double>& input_tD,
                                    const std::vector<double>& input_tauA,
                                    const std::vector<double>& input_tauB,
                                    const std::vector<double>& input_tauD);

    std::vector<double> Growth_dynamic(const std::vector<double>& t);
    std::vector<double> Growth_impact(const std::vector<double>& t);
    std::vector<double> EB_impact(const std::vector<double>& t);
    std::vector<double> cRhoFFM(const std::vector<double>& input_FFM);
    std::vector<double> cP(const std::vector<double>& FFMv, const std::vector<double>& FMv);
    std::vector<double> Delta(const std::vector<double>& t);
    std::vector<double> Expenditure(const std::vector<double>& t,
                                    const std::vector<double>& FFMv,
                                    const std::vector<double>& FMv);
    std::vector<double> Intake(const std::vector<double>& t);

    // Returns matrix [2 x nind] row-major: out[0*nind+i] = dFFM, out[1*nind+i] = dFM
    void dMass(const std::vector<double>& t,
               const std::vector<double>& FFMv,
               const std::vector<double>& FMv,
               std::vector<double>& dFFM_out,
               std::vector<double>& dFM_out);
};

} // namespace bw

#endif // BW_CHILD_HPP
