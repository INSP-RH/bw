// Pure C++ implementation of the adult body-weight RK4 model.
//
// This file MUST NOT include <Rcpp.h>. The math is translated VERBATIM
// from the historical src/adult_weight.cpp — only the container types
// changed (Rcpp::NumericVector -> std::vector<double>, Rcpp::NumericMatrix
// -> bw::Matrix). Operation order, parenthesization, constants, and the
// sequence of RK4 sub-integrations (AT, ECF, Glycogen, Lean) are byte-for-
// byte preserved to keep the parity gate green.

#include "bw/adult.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace bw {

// ---- small elementwise helpers --------------------------------------------
namespace {

inline std::vector<double> vsub(const std::vector<double>& a,
                                const std::vector<double>& b) {
    std::vector<double> r(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) r[i] = a[i] - b[i];
    return r;
}

inline std::vector<double> vadd(const std::vector<double>& a,
                                const std::vector<double>& b) {
    std::vector<double> r(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) r[i] = a[i] + b[i];
    return r;
}

inline std::vector<double> vmul(const std::vector<double>& a,
                                const std::vector<double>& b) {
    std::vector<double> r(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) r[i] = a[i] * b[i];
    return r;
}

inline std::vector<double> vscale(const std::vector<double>& a, double s) {
    std::vector<double> r(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) r[i] = a[i] * s;
    return r;
}

inline std::vector<double> vaxpy(const std::vector<double>& y,
                                 double a,
                                 const std::vector<double>& x) {
    // returns y + a*x elementwise
    std::vector<double> r(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) r[i] = y[i] + a * x[i];
    return r;
}

// Average of two equal-length vectors: 0.5*(a + b).
inline std::vector<double> vavg(const std::vector<double>& a,
                                const std::vector<double>& b) {
    std::vector<double> r(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) r[i] = 0.5 * (a[i] + b[i]);
    return r;
}

// Extract column j of an (nind x N) matrix as a length-nind vector.
inline std::vector<double> col(const Matrix& M, std::size_t j) {
    std::vector<double> r(M.nrow);
    for (std::size_t i = 0; i < M.nrow; ++i) r[i] = M(i, j);
    return r;
}

// Extract row i of an (R x C) matrix as a length-C vector.
inline std::vector<double> row(const Matrix& M, std::size_t i) {
    std::vector<double> r(M.ncol);
    for (std::size_t j = 0; j < M.ncol; ++j) r[j] = M(i, j);
    return r;
}

// Write a length-nind vector into column j of an (nind x N) matrix.
inline void set_col(Matrix& M, std::size_t j, const std::vector<double>& v) {
    for (std::size_t i = 0; i < M.nrow; ++i) M(i, j) = v[i];
}

inline void set_col(StringMatrix& M, std::size_t j,
                    const std::vector<std::string>& v) {
    for (std::size_t i = 0; i < M.nrow; ++i) M(i, j) = v[i];
}

} // anonymous

// ---- constructors ---------------------------------------------------------

Adult::Adult(const std::vector<double>& weight, const std::vector<double>& height,
             const std::vector<double>& age_yrs, const std::vector<double>& sexstring,
             const Matrix& input_EIchange, const Matrix& input_NAchange,
             const std::vector<double>& physicalactivity,
             const std::vector<double>& percentc, const std::vector<double>& percentb,
             double input_dt, bool checkValues) {
    build(weight, height, age_yrs, sexstring, input_EIchange, input_NAchange,
          physicalactivity, percentc, percentb, input_dt, checkValues);
}

Adult::Adult(const std::vector<double>& weight, const std::vector<double>& height,
             const std::vector<double>& age_yrs, const std::vector<double>& sexstring,
             const Matrix& input_EIchange, const Matrix& input_NAchange,
             const std::vector<double>& physicalactivity,
             const std::vector<double>& percentc, const std::vector<double>& percentb,
             double input_dt, const std::vector<double>& extradata,
             bool checkValues, bool isEnergy) {
    build(weight, height, age_yrs, sexstring, input_EIchange, input_NAchange,
          physicalactivity, percentc, percentb, input_dt, extradata, checkValues, isEnergy);
}

Adult::Adult(const std::vector<double>& weight, const std::vector<double>& height,
             const std::vector<double>& age_yrs, const std::vector<double>& sexstring,
             const Matrix& input_EIchange, const Matrix& input_NAchange,
             const std::vector<double>& physicalactivity,
             const std::vector<double>& percentc, const std::vector<double>& percentb,
             double input_dt, const std::vector<double>& input_EI,
             const std::vector<double>& input_fat, bool checkValues) {
    build(weight, height, age_yrs, sexstring, input_EIchange, input_NAchange,
          physicalactivity, percentc, percentb, input_dt, input_EI, input_fat, checkValues);
}

Adult::~Adult() {}

// ---- build variants -------------------------------------------------------

void Adult::build(const std::vector<double>& weight, const std::vector<double>& height,
                  const std::vector<double>& age_yrs, const std::vector<double>& sexstring,
                  const Matrix& input_EIchange, const Matrix& input_NAchange,
                  const std::vector<double>& physicalactivity,
                  const std::vector<double>& percentc, const std::vector<double>& percentb,
                  double input_dt, bool checkValues) {
    dt         = input_dt;
    bw         = weight;
    ht         = height;
    age        = age_yrs;
    sex        = sexstring;
    EIchange   = input_EIchange;
    NAchange   = input_NAchange;
    PAL        = physicalactivity;
    pcarb      = percentc;
    pcarb_base = percentb;
    check      = checkValues;

    getParameters();
    getRMR();
    getATinit();
    getECFinit();
    getBaselineMass();
    getCaloricSteadyState();
    getEnergy();
    getDelta();
    getK();
    getCarbConstants();
}

void Adult::build(const std::vector<double>& weight, const std::vector<double>& height,
                  const std::vector<double>& age_yrs, const std::vector<double>& sexstring,
                  const Matrix& input_EIchange, const Matrix& input_NAchange,
                  const std::vector<double>& physicalactivity,
                  const std::vector<double>& percentc, const std::vector<double>& percentb,
                  double input_dt, const std::vector<double>& extradata,
                  bool checkValues, bool isEnergy) {
    dt         = input_dt;
    bw         = weight;
    ht         = height;
    age        = age_yrs;
    sex        = sexstring;
    EIchange   = input_EIchange;
    NAchange   = input_NAchange;
    PAL        = physicalactivity;
    pcarb      = percentc;
    pcarb_base = percentb;
    check      = checkValues;

    getParameters();
    getRMR();
    getATinit();
    getECFinit();

    if (isEnergy) {
        EI = extradata;
        getBaselineMass();
    } else {
        getCaloricSteadyState();
        getEnergy();
        fat = extradata;
        // lean = bw - (ecfinit + fat + 3.7*G_base)
        lean.assign(bw.size(), 0.0);
        for (std::size_t i = 0; i < bw.size(); ++i) {
            lean[i] = bw[i] - (ecfinit[i] + fat[i] + 3.7 * G_base[i]);
        }
    }

    getDelta();
    getK();
    getCarbConstants();
}

void Adult::build(const std::vector<double>& weight, const std::vector<double>& height,
                  const std::vector<double>& age_yrs, const std::vector<double>& sexstring,
                  const Matrix& input_EIchange, const Matrix& input_NAchange,
                  const std::vector<double>& physicalactivity,
                  const std::vector<double>& percentc, const std::vector<double>& percentb,
                  double input_dt, const std::vector<double>& input_EI,
                  const std::vector<double>& input_fat, bool checkValues) {
    dt         = input_dt;
    bw         = weight;
    ht         = height;
    age        = age_yrs;
    sex        = sexstring;
    EIchange   = input_EIchange;
    NAchange   = input_NAchange;
    PAL        = physicalactivity;
    pcarb      = percentc;
    pcarb_base = percentb;
    check      = checkValues;

    getParameters();
    getRMR();
    getATinit();
    getECFinit();

    EI  = input_EI;
    fat = input_fat;
    // lean = bw - (ecfinit + fat + 3.7*G_base)
    lean.assign(bw.size(), 0.0);
    for (std::size_t i = 0; i < bw.size(); ++i) {
        lean[i] = bw[i] - (ecfinit[i] + fat[i] + 3.7 * G_base[i]);
    }

    getDelta();
    getK();
    getCarbConstants();
}

// ---- parameters -----------------------------------------------------------

void Adult::getParameters() {
    nind    = static_cast<int>(bw.size());

    roG     = 4206.501;
    Na      = 3220;
    zetaNa  = 3000;
    zetaCI  = 4000;
    roF     = 9440.727;
    roL     = 1816.444;
    gammaF  = 3.107075;
    gammaL  = 21.98853;
    etaF    = 179.2543;
    etaL    = 229.4455;
    betaTEF = 0.1;
    betaAT  = 0.14;
    tauAT   = 14.0;
    C       = 10.4 * (roL / roF);
    alfa1   = -(1 + etaL / roL) * C;
    alfa2   = -(1 + etaF / roF);
    rmrbw   = 9.99;
    rmrage  = 4.92;
    rmrht   = 625.0;
    rmr_m   = 5.0;
    rmr_f   = 161.0;
    G_base.assign(static_cast<std::size_t>(nind), 0.5);
}

void Adult::getRMR() {
    // rmr = (rmrbw*bw + rmrht*ht - rmrage*age + rmr_m)*(1-sex)
    //     + (rmrbw*bw + rmrht*ht - rmrage*age - rmr_f)*sex;
    rmr.assign(bw.size(), 0.0);
    for (std::size_t i = 0; i < bw.size(); ++i) {
        rmr[i] = (rmrbw * bw[i] + rmrht * ht[i] - rmrage * age[i] + rmr_m) * (1 - sex[i])
               + (rmrbw * bw[i] + rmrht * ht[i] - rmrage * age[i] - rmr_f) * sex[i];
    }
}

void Adult::getCaloricSteadyState() {
    // steadystate = rmr * PAL
    steadystate.assign(rmr.size(), 0.0);
    for (std::size_t i = 0; i < rmr.size(); ++i) {
        steadystate[i] = rmr[i] * PAL[i];
    }
}

void Adult::getATinit() {
    atinit.assign(static_cast<std::size_t>(nind), 0.0);
}

void Adult::getEnergy() {
    EI = steadystate;
}

void Adult::getDelta() {
    // delta = ((1.0 - betaTEF)*PAL - 1.0) * rmr / bw
    delta.assign(bw.size(), 0.0);
    for (std::size_t i = 0; i < bw.size(); ++i) {
        delta[i] = ((1.0 - betaTEF) * PAL[i] - 1.0) * rmr[i] / bw[i];
    }
}

void Adult::getECFinit() {
    // ecfinit = (0.025*age + 9.57*ht + 0.191*bw - 12.4)*(1.0-sex)
    //         + (-4.0 + 5.98*ht + 0.167*bw)*sex
    ecfinit.assign(bw.size(), 0.0);
    for (std::size_t i = 0; i < bw.size(); ++i) {
        ecfinit[i] = (0.025 * age[i] + 9.57 * ht[i] + 0.191 * bw[i] - 12.4) * (1.0 - sex[i])
                   + (-4.0 + 5.98 * ht[i] + 0.167 * bw[i]) * sex[i];
    }
}

void Adult::getBaselineMass() {
    // fat = (bw * (0.14*age + 37.31*log(bw/(ht^2)) - 103.94)/100.0)*(1-sex)
    //     + (bw * (0.14*age + 39.96*log(bw/(ht^2)) - 102.01)/100.0)*sex
    fat.assign(bw.size(), 0.0);
    for (std::size_t i = 0; i < bw.size(); ++i) {
        const double bmi_arg = bw[i] / std::pow(ht[i], 2.0);
        const double l       = std::log(bmi_arg);
        fat[i] = (bw[i] * (0.14 * age[i] + 37.31 * l - 103.94) / 100.0) * (1 - sex[i])
               + (bw[i] * (0.14 * age[i] + 39.96 * l - 102.01) / 100.0) * sex[i];
    }
    // lean = bw - (ecfinit + fat + 3.7*G_base)
    lean.assign(bw.size(), 0.0);
    for (std::size_t i = 0; i < bw.size(); ++i) {
        lean[i] = bw[i] - (ecfinit[i] + fat[i] + 3.7 * G_base[i]);
    }
}

std::vector<double> Adult::TEF(double t) {
    // betaTEF * deltaEI(t)
    std::vector<double> de = deltaEI(t);
    for (auto& v : de) v *= betaTEF;
    return de;
}

std::vector<double> Adult::dG(double t, const std::vector<double>& G) {
    // (CI(t) - kG*pow(G, 2.0)) / roG
    std::vector<double> ci = CI(t);
    std::vector<double> r(G.size());
    for (std::size_t i = 0; i < G.size(); ++i) {
        r[i] = (ci[i] - kG[i] * std::pow(G[i], 2.0)) / roG;
    }
    return r;
}

std::vector<double> Adult::dAT(double t, const std::vector<double>& AT) {
    // (betaAT * deltaEI(t) - AT) * (1.0 / tauAT)
    std::vector<double> de = deltaEI(t);
    std::vector<double> r(AT.size());
    for (std::size_t i = 0; i < AT.size(); ++i) {
        r[i] = (betaAT * de[i] - AT[i]) * (1.0 / tauAT);
    }
    return r;
}

std::vector<double> Adult::dECF(double t, const std::vector<double>& ECF) {
    // (deltaNA(t) - zetaNa*(ECF - ecfinit) - zetaCI*(1.0 - CI(t)/CIb)) / Na
    std::vector<double> dn = deltaNA(t);
    std::vector<double> ci = CI(t);
    std::vector<double> r(ECF.size());
    for (std::size_t i = 0; i < ECF.size(); ++i) {
        r[i] = (dn[i] - zetaNa * (ECF[i] - ecfinit[i]) - zetaCI * (1.0 - ci[i] / CIb[i])) / Na;
    }
    return r;
}

void Adult::getCarbConstants() {
    // CIb = pcarb_base * EI
    // kG  = CIb / (G_base^2)
    CIb.assign(EI.size(), 0.0);
    kG.assign(EI.size(), 0.0);
    for (std::size_t i = 0; i < EI.size(); ++i) {
        CIb[i] = pcarb_base[i] * EI[i];
        kG[i]  = CIb[i] / std::pow(G_base[i], 2.0);
    }
}

std::vector<double> Adult::CI(double t) {
    // pcarb * TotalIntake(t)
    std::vector<double> ti = TotalIntake(t);
    for (std::size_t i = 0; i < ti.size(); ++i) ti[i] *= pcarb[i];
    return ti;
}

std::vector<double> Adult::TotalIntake(double t) {
    // EI + deltaEI(t)
    std::vector<double> de = deltaEI(t);
    for (std::size_t i = 0; i < de.size(); ++i) de[i] += EI[i];
    return de;
}

void Adult::getK() {
    // K = (rmr * PAL) - gammaL*lean - gammaF*fat - delta*bw
    K.assign(rmr.size(), 0.0);
    for (std::size_t i = 0; i < rmr.size(); ++i) {
        K[i] = (rmr[i] * PAL[i]) - gammaL * lean[i] - gammaF * fat[i] - delta[i] * bw[i];
    }
}

std::vector<double> Adult::fatMass(const std::vector<double>& L) {
    // fat * exp(roL * (L - lean) / (roF * C))
    std::vector<double> r(L.size());
    for (std::size_t i = 0; i < L.size(); ++i) {
        r[i] = fat[i] * std::exp(roL * (L[i] - lean[i]) / (roF * C));
    }
    return r;
}

std::vector<double> Adult::dL(double t,
                              const std::vector<double>& L,
                              const std::vector<double>& G,
                              const std::vector<double>& AT,
                              const std::vector<double>& ECF) {
    // R(t, L, G, AT, ECF) * (C/roL)
    std::vector<double> r = R(t, L, G, AT, ECF);
    const double scale = C / roL;
    for (auto& v : r) v *= scale;
    return r;
}

std::vector<double> Adult::R(double t,
                             const std::vector<double>& L,
                             const std::vector<double>& G,
                             const std::vector<double>& AT,
                             const std::vector<double>& ECF) {
    // F      = fatMass(L)
    // weight = L + F + ECF + 3.7*G
    // R3     = K + delta*weight + TEF(t) + AT - TotalIntake(t) + roG*dG(t, G)
    // return (R3 + gammaL*L + gammaF*F)/(alfa1 + alfa2*F)
    std::vector<double> F  = fatMass(L);
    std::vector<double> ti = TotalIntake(t);
    std::vector<double> dg = dG(t, G);
    std::vector<double> tef = TEF(t);
    std::vector<double> r(L.size());
    for (std::size_t i = 0; i < L.size(); ++i) {
        double weight = L[i] + F[i] + ECF[i] + 3.7 * G[i];
        double R3     = K[i] + delta[i] * weight + tef[i] + AT[i] - ti[i] + roG * dg[i];
        r[i] = (R3 + gammaL * L[i] + gammaF * F[i]) / (alfa1 + alfa2 * F[i]);
    }
    return r;
}

std::vector<std::string> Adult::BMIClassifier(const std::vector<double>& BMI) {
    std::vector<std::string> classification(BMI.size());
    for (std::size_t i = 0; i < BMI.size(); ++i) {
        classification[i] = "Unknown";
        if (BMI[i] < 18.5) {
            classification[i] = "Underweight";
        } else if (BMI[i] >= 18.5 && BMI[i] < 25) {
            classification[i] = "Normal";
        } else if (BMI[i] >= 25 && BMI[i] < 30) {
            classification[i] = "Pre-Obese";
        } else if (BMI[i] >= 30) {
            classification[i] = "Obese";
        }
    }
    return classification;
}

// ---- RK4 integrator -------------------------------------------------------

AdultResult Adult::rk4(double days) {
    // Match the original's loop-bound math: int(min(ceil(days/dt), nrow(EIchange)-1)).
    const int nsims = static_cast<int>(
        std::min(std::ceil(days / dt),
                 static_cast<double>(EIchange.nrow) - 1.0));

    const std::size_t N = static_cast<std::size_t>(nind);
    const std::size_t T = static_cast<std::size_t>(nsims + 1);

    Matrix AT(N, T);
    Matrix ECF(N, T);
    Matrix GLY(N, T);
    Matrix L(N, T);
    Matrix F(N, T);
    Matrix BW(N, T);
    Matrix BMI(N, T);
    Matrix TEI(N, T);
    Matrix AGE(N, T);
    StringMatrix CAT(N, T);
    std::vector<double> TIME(T, 0.0);

    // Initial states.
    set_col(AT,  0, atinit);
    set_col(ECF, 0, ecfinit);
    set_col(GLY, 0, G_base);
    set_col(L,   0, lean);
    {
        std::vector<double> F0 = fatMass(lean);
        set_col(F, 0, F0);
    }
    set_col(BW, 0, bw);
    {
        std::vector<double> bmi0(N);
        for (std::size_t i = 0; i < N; ++i) bmi0[i] = bw[i] / std::pow(ht[i], 2.0);
        set_col(BMI, 0, bmi0);
        set_col(CAT, 0, BMIClassifier(bmi0));
    }
    set_col(TEI, 0, EI);
    TIME[0] = 0.0;
    set_col(AGE, 0, age);

    bool correctVals = true;

    std::vector<double> k1, k2, k3, k4;

    for (int i = 1; i <= nsims; ++i) {
        // (check block is a no-op in the original; preserved here as such.)
        if (!correctVals) break;

        const double t_prev = TIME[static_cast<std::size_t>(i - 1)];
        const std::size_t ip = static_cast<std::size_t>(i);
        const std::size_t im = static_cast<std::size_t>(i - 1);

        // ----- Adaptive thermogenesis ---------------------------------
        {
            std::vector<double> ATprev = col(AT, im);
            k1 = dAT(t_prev,            ATprev);
            k2 = dAT(t_prev + 0.5 * dt, vaxpy(ATprev, 0.5 * dt, k1));
            k3 = dAT(t_prev + 0.5 * dt, vaxpy(ATprev, 0.5 * dt, k2));
            k4 = dAT(t_prev + dt,       vaxpy(ATprev, dt,       k3));
            std::vector<double> next(N);
            for (std::size_t j = 0; j < N; ++j) {
                next[j] = ATprev[j] + dt * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
            }
            set_col(AT, ip, next);
        }

        // ----- Extracellular fluid ------------------------------------
        {
            std::vector<double> ECFprev = col(ECF, im);
            k1 = dECF(t_prev,            ECFprev);
            k2 = dECF(t_prev + 0.5 * dt, vaxpy(ECFprev, 0.5 * dt, k1));
            k3 = dECF(t_prev + 0.5 * dt, vaxpy(ECFprev, 0.5 * dt, k2));
            k4 = dECF(t_prev + dt,       vaxpy(ECFprev, dt,       k3));
            std::vector<double> next(N);
            for (std::size_t j = 0; j < N; ++j) {
                next[j] = ECFprev[j] + dt * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
            }
            set_col(ECF, ip, next);
        }

        // ----- Glycogen -----------------------------------------------
        {
            std::vector<double> GLYprev = col(GLY, im);
            k1 = dG(t_prev,            GLYprev);
            k2 = dG(t_prev + 0.5 * dt, vaxpy(GLYprev, 0.5 * dt, k1));
            k3 = dG(t_prev + 0.5 * dt, vaxpy(GLYprev, 0.5 * dt, k2));
            k4 = dG(t_prev + dt,       vaxpy(GLYprev, dt,       k3));
            std::vector<double> next(N);
            for (std::size_t j = 0; j < N; ++j) {
                next[j] = GLYprev[j] + dt * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
            }
            set_col(GLY, ip, next);
        }

        // ----- Lean mass (uses already-updated AT/ECF/GLY for k2/k3/k4)
        {
            std::vector<double> Lprev   = col(L,   im);
            std::vector<double> GLYprev = col(GLY, im);
            std::vector<double> ATprev  = col(AT,  im);
            std::vector<double> ECFprev = col(ECF, im);
            std::vector<double> GLYnew  = col(GLY, ip);
            std::vector<double> ATnew   = col(AT,  ip);
            std::vector<double> ECFnew  = col(ECF, ip);

            // k1 at TIME(i-1) with previous values
            k1 = dL(t_prev, Lprev, GLYprev, ATprev, ECFprev);
            // k2 at TIME(i-1) + dt/2 using 0.5*(new + old) averages
            k2 = dL(t_prev + 0.5 * dt,
                    vaxpy(Lprev, 0.5 * dt, k1),
                    vavg(GLYnew, GLYprev),
                    vavg(ATnew,  ATprev),
                    vavg(ECFnew, ECFprev));
            // k3 same time/averages but using k2
            k3 = dL(t_prev + 0.5 * dt,
                    vaxpy(Lprev, 0.5 * dt, k2),
                    vavg(GLYnew, GLYprev),
                    vavg(ATnew,  ATprev),
                    vavg(ECFnew, ECFprev));
            // k4 at TIME(i-1)+dt using the new values directly
            k4 = dL(t_prev + dt,
                    vaxpy(Lprev, dt, k3),
                    GLYnew, ATnew, ECFnew);

            std::vector<double> next(N);
            for (std::size_t j = 0; j < N; ++j) {
                next[j] = Lprev[j] + dt * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
            }
            set_col(L, ip, next);
        }

        // ----- F, BW, BMI, CAT, TIME, AGE, TEI -----------------------
        {
            std::vector<double> Lcur = col(L, ip);
            std::vector<double> Fcur = fatMass(Lcur);
            set_col(F, ip, Fcur);

            std::vector<double> ECFcur = col(ECF, ip);
            std::vector<double> GLYcur = col(GLY, ip);

            std::vector<double> BWcur(N);
            for (std::size_t j = 0; j < N; ++j) {
                BWcur[j] = Fcur[j] + Lcur[j] + ECFcur[j] + 3.7 * GLYcur[j];
            }
            set_col(BW, ip, BWcur);

            std::vector<double> BMIcur(N);
            for (std::size_t j = 0; j < N; ++j) {
                BMIcur[j] = BWcur[j] / std::pow(ht[j], 2.0);
            }
            set_col(BMI, ip, BMIcur);
            set_col(CAT, ip, BMIClassifier(BMIcur));

            TIME[ip] = t_prev + dt;

            std::vector<double> AGEprev = col(AGE, im);
            std::vector<double> AGEcur(N);
            for (std::size_t j = 0; j < N; ++j) {
                AGEcur[j] = AGEprev[j] + dt / 365.0;
            }
            set_col(AGE, ip, AGEcur);

            std::vector<double> TI = TotalIntake(TIME[ip]);
            set_col(TEI, ip, TI);
        }
    }

    AdultResult out;
    out.Time                   = std::move(TIME);
    out.Age                    = std::move(AGE);
    out.Adaptive_Thermogenesis = std::move(AT);
    out.Extracellular_Fluid    = std::move(ECF);
    out.Glycogen               = std::move(GLY);
    out.Fat_Mass               = std::move(F);
    out.Lean_Mass              = std::move(L);
    out.Body_Weight            = std::move(BW);
    out.Body_Mass_Index        = std::move(BMI);
    out.BMI_Category           = std::move(CAT);
    out.Energy_Intake          = std::move(TEI);
    out.Correct_Values         = correctVals;
    out.Model_Type             = "Adult";
    return out;
}

std::vector<double> Adult::deltaEI(double t) {
    // EIchange(floor(t/dt), _) — row across all individuals at the
    // (floor) day index. Matches the original Rcpp behavior verbatim.
    const std::size_t idx = static_cast<std::size_t>(std::floor(t / dt));
    return row(EIchange, idx);
}

std::vector<double> Adult::deltaNA(double t) {
    const std::size_t idx = static_cast<std::size_t>(std::floor(t / dt));
    return row(NAchange, idx);
}

} // namespace bw
