// Pure C++ public header for the adult body-weight model.
//
// This header MUST NOT include <Rcpp.h>. It exposes the adult RK4 kernel
// using only standard C++ types. The Rcpp marshaling layer lives in
// src/adult_rcpp.cpp.

#ifndef BW_ADULT_HPP
#define BW_ADULT_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace bw {

// Row-major dense matrix of doubles. Stored as nrow rows of ncol entries.
// All matrices used by the adult kernel are double-precision.
struct Matrix {
    std::size_t nrow;
    std::size_t ncol;
    std::vector<double> data; // row-major: data[i*ncol + j]

    Matrix() : nrow(0), ncol(0) {}
    Matrix(std::size_t r, std::size_t c)
        : nrow(r), ncol(c), data(r * c, 0.0) {}
    Matrix(std::size_t r, std::size_t c, double fill)
        : nrow(r), ncol(c), data(r * c, fill) {}

    double& operator()(std::size_t i, std::size_t j) { return data[i * ncol + j]; }
    double  operator()(std::size_t i, std::size_t j) const { return data[i * ncol + j]; }
};

// Row-major dense matrix of std::string for BMI category output.
struct StringMatrix {
    std::size_t nrow;
    std::size_t ncol;
    std::vector<std::string> data;

    StringMatrix() : nrow(0), ncol(0) {}
    StringMatrix(std::size_t r, std::size_t c)
        : nrow(r), ncol(c), data(r * c) {}

    std::string& operator()(std::size_t i, std::size_t j) { return data[i * ncol + j]; }
    const std::string& operator()(std::size_t i, std::size_t j) const { return data[i * ncol + j]; }
};

// Output of Adult::rk4. Each numeric matrix is shaped (nind, nsims+1) — one
// row per individual, one column per timestep. This matches the original
// Rcpp NumericMatrix layout the wrapper used to return.
struct AdultResult {
    std::vector<double> Time;            // length nsims+1
    Matrix              Age;             // (nind, nsims+1)
    Matrix              Adaptive_Thermogenesis;
    Matrix              Extracellular_Fluid;
    Matrix              Glycogen;
    Matrix              Fat_Mass;
    Matrix              Lean_Mass;
    Matrix              Body_Weight;
    Matrix              Body_Mass_Index;
    StringMatrix        BMI_Category;
    Matrix              Energy_Intake;
    bool                Correct_Values;
    std::string         Model_Type;
};

// Pure-C++ adult body-weight model.
//
// EIchange and NAchange are passed in already-transposed form (rows = day
// index, cols = individual) — the same layout used by the original Rcpp
// kernel. The R wrapper transposes its input before calling.
class Adult {
public:
    // Constructor when initial energy intake is estimated from Mifflin-St Jeor.
    Adult(const std::vector<double>& weight,
          const std::vector<double>& height,
          const std::vector<double>& age_yrs,
          const std::vector<double>& sexstring,
          const Matrix& input_EIchange,
          const Matrix& input_NAchange,
          const std::vector<double>& physicalactivity,
          const std::vector<double>& percentc,
          const std::vector<double>& percentb,
          double dt, bool checkValues);

    // Constructor when initial energy intake OR initial fat is user-supplied.
    Adult(const std::vector<double>& weight,
          const std::vector<double>& height,
          const std::vector<double>& age_yrs,
          const std::vector<double>& sexstring,
          const Matrix& input_EIchange,
          const Matrix& input_NAchange,
          const std::vector<double>& physicalactivity,
          const std::vector<double>& percentc,
          const std::vector<double>& percentb,
          double dt,
          const std::vector<double>& extradata,
          bool checkValues, bool isEnergy);

    // Constructor when both initial energy intake AND initial fat are supplied.
    Adult(const std::vector<double>& weight,
          const std::vector<double>& height,
          const std::vector<double>& age_yrs,
          const std::vector<double>& sexstring,
          const Matrix& input_EIchange,
          const Matrix& input_NAchange,
          const std::vector<double>& physicalactivity,
          const std::vector<double>& percentc,
          const std::vector<double>& percentb,
          double dt,
          const std::vector<double>& input_EI,
          const std::vector<double>& input_fat,
          bool checkValues);

    ~Adult();

    // Run the RK4 integrator for `days` days; returns full trajectories.
    AdultResult rk4(double days);

    // Public state mirrors the original NumericVector public members.
    std::vector<double> bw;
    std::vector<double> ht;
    std::vector<double> age;
    std::vector<double> sex;
    std::vector<double> EI;
    std::vector<double> PAL;
    std::vector<double> fat;
    std::vector<double> lean;
    std::vector<double> steadystate;
    std::vector<double> G_base;
    std::vector<double> ecfinit;
    std::vector<double> CIb;
    std::vector<double> pcarb;
    std::vector<double> pcarb_base;
    Matrix              EIchange;
    Matrix              NAchange;

private:
    // Private state.
    std::vector<double> kG;
    std::vector<double> K;
    std::vector<double> rmr;
    std::vector<double> delta;
    std::vector<double> atinit;

    double roG;
    double Na;
    double zetaNa;
    double zetaCI;
    double roF;
    double roL;
    double gammaF;
    double gammaL;
    double etaF;
    double etaL;
    double betaTEF;
    double betaAT;
    double tauAT;
    double C;
    double alfa1;
    double alfa2;
    double rmrbw;
    double rmrage;
    double rmrht;
    double rmr_m;
    double rmr_f;
    int    nind;
    double dt;
    bool   check;

    // Auxiliary builders.
    void getRMR();
    void getParameters();
    void getBaselineMass();
    void getCaloricSteadyState();
    void getEnergy();
    void getDelta();
    void getK();
    void getCarbConstants();
    void getATinit();
    void getECFinit();

    void build(const std::vector<double>& weight,
               const std::vector<double>& height,
               const std::vector<double>& age_yrs,
               const std::vector<double>& sexstring,
               const Matrix& input_EIchange,
               const Matrix& input_NAchange,
               const std::vector<double>& physicalactivity,
               const std::vector<double>& percentc,
               const std::vector<double>& percentb,
               double input_dt, bool checkValues);

    void build(const std::vector<double>& weight,
               const std::vector<double>& height,
               const std::vector<double>& age_yrs,
               const std::vector<double>& sexstring,
               const Matrix& input_EIchange,
               const Matrix& input_NAchange,
               const std::vector<double>& physicalactivity,
               const std::vector<double>& percentc,
               const std::vector<double>& percentb,
               double input_dt,
               const std::vector<double>& extradata,
               bool checkValues, bool isEnergy);

    void build(const std::vector<double>& weight,
               const std::vector<double>& height,
               const std::vector<double>& age_yrs,
               const std::vector<double>& sexstring,
               const Matrix& input_EIchange,
               const Matrix& input_NAchange,
               const std::vector<double>& physicalactivity,
               const std::vector<double>& percentc,
               const std::vector<double>& percentb,
               double input_dt,
               const std::vector<double>& input_EI,
               const std::vector<double>& input_fat,
               bool checkValues);

    // Per-step derivative helpers. Each returns a length-nind vector and
    // preserves the original arithmetic verbatim.
    std::vector<double> TotalIntake(double t);
    std::vector<std::string> BMIClassifier(const std::vector<double>& BMI);
    std::vector<double> CI(double t);
    std::vector<double> R(double t,
                          const std::vector<double>& L,
                          const std::vector<double>& G,
                          const std::vector<double>& AT,
                          const std::vector<double>& ECF);
    std::vector<double> fatMass(const std::vector<double>& L);
    std::vector<double> deltaEI(double t);
    std::vector<double> deltaNA(double t);
    std::vector<double> TEF(double t);
    std::vector<double> dAT(double t, const std::vector<double>& AT);
    std::vector<double> dECF(double t, const std::vector<double>& ECF);
    std::vector<double> dG(double t, const std::vector<double>& G);
    std::vector<double> dL(double t,
                           const std::vector<double>& L,
                           const std::vector<double>& G,
                           const std::vector<double>& AT,
                           const std::vector<double>& ECF);
};

} // namespace bw

#endif // BW_ADULT_HPP
