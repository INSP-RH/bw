//
//  kernel/child.cpp — pure C++ implementation of the Child weight-dynamics model.
//
//  Translated verbatim from src/child_weight.cpp. Arithmetic order is preserved.
//
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
//----------------------------------------------------------------------------------------

#include "bw/child.hpp"

#include <cmath>
#include <cstddef>
#include <vector>
#include <algorithm>

namespace bw {

// ----- Constructors --------------------------------------------------------

Child::Child(const std::vector<double>& input_age,
             const std::vector<double>& input_sex,
             const std::vector<double>& input_FFM,
             const std::vector<double>& input_FM,
             const std::vector<double>& input_EIntake,
             std::size_t ei_nrows_,
             std::size_t ei_ncols_,
             double input_dt,
             bool checkValues)
    : age(input_age),
      sex(input_sex),
      FFM(input_FFM),
      FM(input_FM),
      EIntake(input_EIntake),
      ei_nrows(ei_nrows_),
      ei_ncols(ei_ncols_),
      check(checkValues),
      generalized_logistic(false)
{
    dt = input_dt;
    build();
}

Child::Child(const std::vector<double>& input_age,
             const std::vector<double>& input_sex,
             const std::vector<double>& input_FFM,
             const std::vector<double>& input_FM,
             double input_K, double input_Q, double input_A,
             double input_B, double input_nu, double input_C,
             double input_dt,
             bool checkValues)
    : age(input_age),
      sex(input_sex),
      FFM(input_FFM),
      FM(input_FM),
      EIntake(),
      ei_nrows(0),
      ei_ncols(0),
      check(checkValues),
      generalized_logistic(true)
{
    dt = input_dt;
    K_logistic = input_K;
    A_logistic = input_A;
    Q_logistic = input_Q;
    B_logistic = input_B;
    nu_logistic = input_nu;
    C_logistic = input_C;
    build();
}

Child::~Child() {}

void Child::build() {
    getParameters();
}

// ----- general_ode ---------------------------------------------------------
//
// Reference:
//   return input_A*exp(-(t-input_tA)/input_tauA) +
//          input_B*exp(-0.5*pow((t-input_tB)/input_tauB,2)) +
//          input_D*exp(-0.5*pow((t-input_tD)/input_tauD,2));
//
// Element-wise; per-element order preserved.
std::vector<double> Child::general_ode(const std::vector<double>& t,
                                       const std::vector<double>& input_A,
                                       const std::vector<double>& input_B,
                                       const std::vector<double>& input_D,
                                       const std::vector<double>& input_tA,
                                       const std::vector<double>& input_tB,
                                       const std::vector<double>& input_tD,
                                       const std::vector<double>& input_tauA,
                                       const std::vector<double>& input_tauB,
                                       const std::vector<double>& input_tauD)
{
    std::vector<double> out(nind);
    for (std::size_t i = 0; i < nind; ++i) {
        out[i] = input_A[i]*std::exp(-(t[i]-input_tA[i])/input_tauA[i]) +
                 input_B[i]*std::exp(-0.5*std::pow((t[i]-input_tB[i])/input_tauB[i], 2)) +
                 input_D[i]*std::exp(-0.5*std::pow((t[i]-input_tD[i])/input_tauD[i], 2));
    }
    return out;
}

std::vector<double> Child::Growth_dynamic(const std::vector<double>& t) {
    return general_ode(t, A, B, D, tA, tB, tD, tauA, tauB, tauD);
}

std::vector<double> Child::Growth_impact(const std::vector<double>& t) {
    return general_ode(t, A1, B1, D1, tA1, tB1, tD1, tauA1, tauB1, tauD1);
}

std::vector<double> Child::EB_impact(const std::vector<double>& t) {
    return general_ode(t, A_EB, B_EB, D_EB, tA_EB, tB_EB, tD_EB, tauA_EB, tauB_EB, tauD_EB);
}

// ----- cRhoFFM -------------------------------------------------------------
// return 4.3*input_FFM + 837.0;
std::vector<double> Child::cRhoFFM(const std::vector<double>& input_FFM) {
    std::vector<double> out(nind);
    for (std::size_t i = 0; i < nind; ++i) {
        out[i] = 4.3*input_FFM[i] + 837.0;
    }
    return out;
}

// ----- cP ------------------------------------------------------------------
// NumericVector rhoFFM = cRhoFFM(FFM);
// NumericVector C      = 10.4 * rhoFFM / rhoFM;
// return C/(C + FM);
std::vector<double> Child::cP(const std::vector<double>& FFMv,
                              const std::vector<double>& FMv)
{
    std::vector<double> rhoFFM = cRhoFFM(FFMv);
    std::vector<double> C(nind);
    for (std::size_t i = 0; i < nind; ++i) {
        C[i] = 10.4 * rhoFFM[i] / rhoFM;
    }
    std::vector<double> out(nind);
    for (std::size_t i = 0; i < nind; ++i) {
        out[i] = C[i]/(C[i] + FMv[i]);
    }
    return out;
}

// ----- Delta ---------------------------------------------------------------
// return deltamin + (deltamax - deltamin)*(1.0 / (1.0 + pow((t / P),h)));
std::vector<double> Child::Delta(const std::vector<double>& t) {
    std::vector<double> out(nind);
    for (std::size_t i = 0; i < nind; ++i) {
        out[i] = deltamin + (deltamax[i] - deltamin)*(1.0 / (1.0 + std::pow((t[i] / P), h)));
    }
    return out;
}

// ----- FFMReference --------------------------------------------------------
std::vector<double> Child::FFMReference(const std::vector<double>& t) {
    // ffm_ref[row][i], 17 rows, nind cols
    std::vector<std::vector<double>> ffm_ref(17, std::vector<double>(nind));
    for (std::size_t i = 0; i < nind; ++i) {
        ffm_ref[0][i]  = 10.134*(1-sex[i])+9.477*sex[i];
        ffm_ref[1][i]  = 12.099*(1 - sex[i]) + 11.494*sex[i];
        ffm_ref[2][i]  = 14.0*(1 - sex[i]) + 13.2*sex[i];
        ffm_ref[3][i]  = 16.0*(1 - sex[i]) + 14.7*sex[i];
        ffm_ref[4][i]  = 17.9*(1 - sex[i]) + 16.3*sex[i];
        ffm_ref[5][i]  = 19.9*(1 - sex[i]) + 18.2*sex[i];
        ffm_ref[6][i]  = 22.0*(1 - sex[i]) + 20.5*sex[i];
        ffm_ref[7][i]  = 24.4*(1 - sex[i]) + 23.3*sex[i];
        ffm_ref[8][i]  = 27.5*(1 - sex[i]) + 26.4*sex[i];
        ffm_ref[9][i]  = 29.5*(1 - sex[i]) + 28.5*sex[i];
        ffm_ref[10][i] = 33.2*(1 - sex[i]) + 32.4*sex[i];
        ffm_ref[11][i] = 38.1*(1 - sex[i]) + 36.1*sex[i];
        ffm_ref[12][i] = 43.6*(1 - sex[i]) + 38.9*sex[i];
        ffm_ref[13][i] = 49.1*(1 - sex[i]) + 40.7*sex[i];
        ffm_ref[14][i] = 54.0*(1 - sex[i]) + 41.7*sex[i];
        ffm_ref[15][i] = 57.7*(1 - sex[i]) + 42.3*sex[i];
        ffm_ref[16][i] = 60.0*(1 - sex[i]) + 42.6*sex[i];
    }

    std::vector<double> ffm_ref_t(nind);
    int jmin;
    int jmax;
    double diff;
    for (std::size_t i = 0; i < nind; ++i) {
        if (t[i] >= 18.0) {
            ffm_ref_t[i] = ffm_ref[16][i];
        } else {
            jmin = static_cast<int>(std::floor(t[i]));
            jmin = std::max(jmin, 2);
            jmin = jmin - 2;
            jmax = std::min(jmin + 1, 17);
            diff = t[i] - std::floor(t[i]);
            ffm_ref_t[i] = ffm_ref[jmin][i] + diff*(ffm_ref[jmax][i] - ffm_ref[jmin][i]);
        }
    }
    return ffm_ref_t;
}

// ----- FMReference ---------------------------------------------------------
std::vector<double> Child::FMReference(const std::vector<double>& t) {
    std::vector<std::vector<double>> fm_ref(17, std::vector<double>(nind));
    for (std::size_t i = 0; i < nind; ++i) {
        fm_ref[0][i]  = 2.456*(1-sex[i])+ 2.433*sex[i];
        fm_ref[1][i]  = 2.576*(1 - sex[i]) + 2.606*sex[i];
        fm_ref[2][i]  = 2.7*(1 - sex[i]) + 2.8*sex[i];
        fm_ref[3][i]  = 2.7*(1 - sex[i]) + 2.9*sex[i];
        fm_ref[4][i]  = 2.8*(1 - sex[i]) + 3.2*sex[i];
        fm_ref[5][i]  = 2.9*(1 - sex[i]) + 3.7*sex[i];
        fm_ref[6][i]  = 3.3*(1 - sex[i]) + 4.3*sex[i];
        fm_ref[7][i]  = 3.7*(1 - sex[i]) + 5.2*sex[i];
        fm_ref[8][i]  = 4.8*(1 - sex[i]) + 7.2*sex[i];
        fm_ref[9][i]  = 5.9*(1 - sex[i]) + 8.5*sex[i];
        fm_ref[10][i] = 6.7*(1 - sex[i]) + 9.2*sex[i];
        fm_ref[11][i] = 7.0*(1 - sex[i]) + 10.0*sex[i];
        fm_ref[12][i] = 7.2*(1 - sex[i]) + 11.3*sex[i];
        fm_ref[13][i] = 7.5*(1 - sex[i]) + 12.8*sex[i];
        fm_ref[14][i] = 8.0*(1 - sex[i]) + 14.0*sex[i];
        fm_ref[15][i] = 8.4*(1 - sex[i]) + 14.3*sex[i];
        fm_ref[16][i] = 8.8*(1 - sex[i]) + 14.3*sex[i];
    }

    std::vector<double> fm_ref_t(nind);
    int jmin;
    int jmax;
    double diff;
    for (std::size_t i = 0; i < nind; ++i) {
        if (t[i] >= 18.0) {
            fm_ref_t[i] = fm_ref[16][i];
        } else {
            jmin = static_cast<int>(std::floor(t[i]));
            jmin = std::max(jmin, 2);
            jmin = jmin - 2;
            jmax = std::min(jmin + 1, 17);
            diff = t[i] - std::floor(t[i]);
            fm_ref_t[i] = fm_ref[jmin][i] + diff*(fm_ref[jmax][i] - fm_ref[jmin][i]);
        }
    }
    return fm_ref_t;
}

// ----- IntakeReference -----------------------------------------------------
// return EB + K + (22.4 + delta)*FFMref + (4.5 + delta)*FMref +
//        230.0/rhoFFM*(p*EB + growth) + 180.0/rhoFM*((1-p)*EB - growth);
std::vector<double> Child::IntakeReference(const std::vector<double>& t) {
    std::vector<double> EB     = EB_impact(t);
    std::vector<double> FFMref = FFMReference(t);
    std::vector<double> FMref  = FMReference(t);
    std::vector<double> delta  = Delta(t);
    std::vector<double> growth = Growth_dynamic(t);
    std::vector<double> p      = cP(FFMref, FMref);
    std::vector<double> rhoFFM = cRhoFFM(FFMref);
    std::vector<double> out(nind);
    for (std::size_t i = 0; i < nind; ++i) {
        out[i] = EB[i] + K[i] + (22.4 + delta[i])*FFMref[i] + (4.5 + delta[i])*FMref[i] +
                 230.0/rhoFFM[i]*(p[i]*EB[i] + growth[i]) +
                 180.0/rhoFM*((1-p[i])*EB[i] - growth[i]);
    }
    return out;
}

// ----- Expenditure ---------------------------------------------------------
// NumericVector Expend = K + (22.4 + delta)*FFM + (4.5 + delta)*FM +
//                        0.24*DeltaI + (230.0/rhoFFM *p + 180.0/rhoFM*(1.0-p))*Intakeval +
//                        growth*(230.0/rhoFFM - 180.0/rhoFM);
// return Expend/(1.0 + 230.0/rhoFFM *p + 180.0/rhoFM*(1.0-p));
std::vector<double> Child::Expenditure(const std::vector<double>& t,
                                       const std::vector<double>& FFMv,
                                       const std::vector<double>& FMv)
{
    std::vector<double> delta     = Delta(t);
    std::vector<double> Iref      = IntakeReference(t);
    std::vector<double> Intakeval = Intake(t);
    std::vector<double> DeltaI(nind);
    for (std::size_t i = 0; i < nind; ++i) {
        DeltaI[i] = Intakeval[i] - Iref[i];
    }
    std::vector<double> p      = cP(FFMv, FMv);
    std::vector<double> rhoFFM = cRhoFFM(FFMv);
    std::vector<double> growth = Growth_dynamic(t);

    std::vector<double> out(nind);
    for (std::size_t i = 0; i < nind; ++i) {
        double Expend = K[i] + (22.4 + delta[i])*FFMv[i] + (4.5 + delta[i])*FMv[i] +
                        0.24*DeltaI[i] + (230.0/rhoFFM[i] *p[i] + 180.0/rhoFM*(1.0-p[i]))*Intakeval[i] +
                        growth[i]*(230.0/rhoFFM[i] - 180.0/rhoFM);
        out[i] = Expend/(1.0 + 230.0/rhoFFM[i] *p[i] + 180.0/rhoFM*(1.0-p[i]));
    }
    return out;
}

// ----- dMass ---------------------------------------------------------------
// Mass(0,_) = (1.0*p*(Intake(t) - expend) + growth)/rhoFFM;
// Mass(1,_) = ((1.0 - p)*(Intake(t) - expend) - growth)/rhoFM;
//
// Note: Intake(t) is invoked TWICE in the reference. Calls are pure; we
// preserve the call count to keep semantics identical (no observable effect,
// but defensive against future side-effects).
void Child::dMass(const std::vector<double>& t,
                  const std::vector<double>& FFMv,
                  const std::vector<double>& FMv,
                  std::vector<double>& dFFM_out,
                  std::vector<double>& dFM_out)
{
    std::vector<double> rhoFFM = cRhoFFM(FFMv);
    std::vector<double> p      = cP(FFMv, FMv);
    std::vector<double> growth = Growth_dynamic(t);
    std::vector<double> expend = Expenditure(t, FFMv, FMv);

    std::vector<double> intake_ffm = Intake(t);
    for (std::size_t i = 0; i < nind; ++i) {
        dFFM_out[i] = (1.0*p[i]*(intake_ffm[i] - expend[i]) + growth[i])/rhoFFM[i];
    }
    std::vector<double> intake_fm = Intake(t);
    for (std::size_t i = 0; i < nind; ++i) {
        dFM_out[i] = ((1.0 - p[i])*(intake_fm[i] - expend[i]) - growth[i])/rhoFM;
    }
}

// ----- rk4 -----------------------------------------------------------------
ChildRK4Result Child::rk4(double days) {

    int nsims = static_cast<int>(std::floor(days/dt));
    std::size_t nsteps = static_cast<std::size_t>(nsims) + 1;

    ChildRK4Result R;
    R.nind = nind;
    R.nsteps = nsteps;
    R.Time.assign(nsteps, 0.0);
    R.Age.assign(nind*nsteps, 0.0);
    R.Fat_Free_Mass.assign(nind*nsteps, 0.0);
    R.Fat_Mass.assign(nind*nsteps, 0.0);
    R.Body_Weight.assign(nind*nsteps, 0.0);
    R.Correct_Values = true;

    // Initial states — column 0 (matrix is nind rows, nsteps cols, column-major
    // layout: element (i, j) -> data[i + nind*j]).
    for (std::size_t i = 0; i < nind; ++i) {
        R.Fat_Free_Mass[i + nind*0] = FFM[i];
        R.Fat_Mass    [i + nind*0] = FM[i];
        R.Body_Weight [i + nind*0] = FFM[i] + FM[i];
        R.Age         [i + nind*0] = age[i];
    }
    R.Time[0] = 0.0;

    // Temporary state vectors
    std::vector<double> prevFFM(nind), prevFM(nind), prevAge(nind);
    std::vector<double> tmpFFM(nind), tmpFM(nind), tmpAge(nind);
    std::vector<double> k1_ffm(nind), k1_fm(nind);
    std::vector<double> k2_ffm(nind), k2_fm(nind);
    std::vector<double> k3_ffm(nind), k3_fm(nind);
    std::vector<double> k4_ffm(nind), k4_fm(nind);

    for (int istep = 1; istep <= nsims; ++istep) {
        std::size_t j = static_cast<std::size_t>(istep);

        // Pull previous column
        for (std::size_t i = 0; i < nind; ++i) {
            prevFFM[i] = R.Fat_Free_Mass[i + nind*(j-1)];
            prevFM[i]  = R.Fat_Mass    [i + nind*(j-1)];
            prevAge[i] = R.Age         [i + nind*(j-1)];
        }

        // k1 = dMass(AGE[i-1], FFM[i-1], FM[i-1])
        dMass(prevAge, prevFFM, prevFM, k1_ffm, k1_fm);

        // k2 = dMass(AGE[i-1] + 0.5 * dt/365.0, FFM[i-1] + 0.5 * k1_ffm, FM[i-1] + 0.5 * k1_fm)
        for (std::size_t i = 0; i < nind; ++i) {
            tmpAge[i] = prevAge[i] + 0.5 * dt/365.0;
            tmpFFM[i] = prevFFM[i] + 0.5 * k1_ffm[i];
            tmpFM[i]  = prevFM[i]  + 0.5 * k1_fm[i];
        }
        dMass(tmpAge, tmpFFM, tmpFM, k2_ffm, k2_fm);

        // k3 = dMass(AGE[i-1] + 0.5 * dt/365.0, FFM[i-1] + 0.5 * k2_ffm, FM[i-1] + 0.5 * k2_fm)
        for (std::size_t i = 0; i < nind; ++i) {
            tmpAge[i] = prevAge[i] + 0.5 * dt/365.0;
            tmpFFM[i] = prevFFM[i] + 0.5 * k2_ffm[i];
            tmpFM[i]  = prevFM[i]  + 0.5 * k2_fm[i];
        }
        dMass(tmpAge, tmpFFM, tmpFM, k3_ffm, k3_fm);

        // k4 = dMass(AGE[i-1] + dt/365.0, FFM[i-1] + k3_ffm, FM[i-1] + k3_fm)
        for (std::size_t i = 0; i < nind; ++i) {
            tmpAge[i] = prevAge[i] + dt/365.0;
            tmpFFM[i] = prevFFM[i] + k3_ffm[i];
            tmpFM[i]  = prevFM[i]  + k3_fm[i];
        }
        dMass(tmpAge, tmpFFM, tmpFM, k4_ffm, k4_fm);

        // Update — DO NOT factor or rearrange.
        for (std::size_t i = 0; i < nind; ++i) {
            double new_ffm = prevFFM[i] + dt*(k1_ffm[i] + 2.0*k2_ffm[i] + 2.0*k3_ffm[i] + k4_ffm[i])/6.0;
            double new_fm  = prevFM[i]  + dt*(k1_fm[i]  + 2.0*k2_fm[i]  + 2.0*k3_fm[i]  + k4_fm[i])/6.0;
            R.Fat_Free_Mass[i + nind*j] = new_ffm;
            R.Fat_Mass    [i + nind*j] = new_fm;
            R.Body_Weight [i + nind*j] = new_ffm + new_fm;
        }

        R.Time[j] = R.Time[j-1] + dt;

        // Age — accumulated, NOT recomputed.
        for (std::size_t i = 0; i < nind; ++i) {
            R.Age[i + nind*j] = R.Age[i + nind*(j-1)] + dt/365.0;
        }
    }

    return R;
}

// ----- getParameters -------------------------------------------------------
void Child::getParameters() {
    rhoFM    = 9.4*1000.0;
    deltamin = 10.0;
    P        = 12.0;
    h        = 10.0;

    nind = age.size();

    auto two = [this](double male_v, double female_v) {
        std::vector<double> out(nind);
        for (std::size_t i = 0; i < nind; ++i) {
            out[i] = male_v*(1 - sex[i]) + female_v*sex[i];
        }
        return out;
    };

    // NOTE: keep the exact (1 - sex) vs (1-sex) literals — they're arithmetically
    // identical but for parity readability we mirror the reference verbatim above.
    // The values are pure constants times sex, so order is fine either way.

    ffm_beta0 = two(2.9, 3.8);
    ffm_beta1 = two(2.9, 2.3);
    fm_beta0  = two(1.2, 0.56);
    fm_beta1  = two(0.41, 0.74);
    K         = two(800, 700);
    deltamax  = two(19, 17);
    A         = two(3.2, 2.3);
    B         = two(9.6, 8.4);
    D         = two(10.1, 1.1);
    tA        = two(4.7, 4.5);
    tB        = two(12.5, 11.7);
    tD        = two(15.0, 16.2);
    tauA      = two(2.5, 1.0);
    tauB      = two(1.0, 0.9);
    tauD      = two(1.5, 0.7);
    A_EB      = two(7.2, 16.5);
    B_EB      = two(30, 47.0);
    D_EB      = two(21, 41.0);
    tA_EB     = two(5.6, 4.8);
    tB_EB     = two(9.8, 9.1);
    tD_EB     = two(15.0, 13.5);
    tauA_EB   = two(15, 7.0);
    tauB_EB   = two(1.5, 1.0);
    tauD_EB   = two(2.0, 1.5);
    A1        = two(3.2, 2.3);
    B1        = two(9.6, 8.4);
    D1        = two(10.0, 1.1);
    tA1       = two(4.7, 4.5);
    tB1       = two(12.5, 11.7);
    tD1       = two(15.0, 16.0);
    tauA1     = two(1.0, 1.0);
    tauB1     = two(0.94, 0.94);
    tauD1     = two(0.69, 0.69);
}

// ----- Intake --------------------------------------------------------------
// Intake in calories
std::vector<double> Child::Intake(const std::vector<double>& t) {
    std::vector<double> out(nind);
    if (generalized_logistic) {
        for (std::size_t i = 0; i < nind; ++i) {
            out[i] = A_logistic + (K_logistic - A_logistic)/std::pow(C_logistic + Q_logistic*std::exp(-B_logistic*t[i]), 1/nu_logistic);
        }
    } else {
        // EIntake row lookup: timeval = floor(365.0 * (t[0] - age[0]) / dt).
        // CRITICAL: the IEEE-754 quirk (365.0 * (1.0/365.0) == 0.9999999999999999)
        // is part of the reference; preserve this expression verbatim.
        int timeval = static_cast<int>(std::floor(365.0*(t[0] - age[0])/dt));
        // EIntake is row-major in our storage: row r, col c -> EIntake[r*ei_ncols + c]
        // The reference uses NumericMatrix(timeval, _) — full row of width nind == ei_ncols.
        for (std::size_t i = 0; i < nind; ++i) {
            out[i] = EIntake[static_cast<std::size_t>(timeval)*ei_ncols + i];
        }
    }
    return out;
}

} // namespace bw
