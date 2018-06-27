//
//  child_weight.cpp
//
//  This is a function that calculates
//  weight change for children using the dynamic
//  weight model by Kevin D. Hall et al. and
//  Runge Kutta method to solve the ODE system.
//
//  Input:
//  age             .-  Years since individual first arrived to Earth
//  sex             .-  Either 1 = "female" or 0 = "male"
//  FFM             .-  Fat Free Mass (kg) of the individual
//  FM              .-  Fat Mass (kg) of the individual
//  input_EIntake   .-  Energy intake (kcal) of individual per day
//  days            .-  Days to model (integer)
//  dt              .-  Time step used to solve the ODE system numerically
//  K               .-  Richardson parameter
//  Q               .-  Richardson parameter
//  A               .-  Richardson parameter
//  B               .-  Richardson parameter
//  nu              .-  Richardson parameter
//  C               .-  Richardson parameter
//  Note:
//  Weight = FFM + FM. No extracellular fluid or glycogen is considered
//  Please see child_weight.hpp for additional information
//
//  Authors:
//  Dalia Camacho-García-Formentí
//  Rodrigo Zepeda-Tello
//
// References:
//
//  Deurenberg, Paul, Jan A Weststrate, and Jaap C Seidell. 1991. “Body Mass Index as a Measure of Body Fatness:
//      Age-and Sex-Specific Prediction Formulas.” British Journal of Nutrition 65 (2). Cambridge University Press: 105–14.
//
//  Ellis, Kenneth J, Roman J Shypailo, Steven A Abrams, and William W Wong. 2000. “The Reference Child and Adolescent Models of
//      Body Composition: A Contemporary Comparison.” Annals of the New York Academy of Sciences 904 (1). Wiley Online Library: 374–82.
//
//  Fomon, Samuel J, Ferdinand Haschke, Ekhard E Ziegler, and Steven E Nelson. 1982.
//      “Body Composition of Reference Children from Birth to Age 10 Years.” The American Journal of
//      Clinical Nutrition 35 (5). Am Soc Nutrition: 1169–75.
//
//  Hall, Kevin D, Nancy F Butte, Boyd A Swinburn, and Carson C Chow. 2013. “Dynamics of Childhood Growth
//      and Obesity: Development and Validation of a Quantitative Mathematical Model.” The Lancet Diabetes & Endocrinology 1 (2). Elsevier: 97–105.
//
//  Haschke, F. 1989. “Body Composition During Adolescence.” Body Composition Measurements in Infants and Children.
//      Ross Laboratories Columbus, OH, 76–83.
//
//  Katan, Martijn B, Janne C De Ruyter, Lothar DJ Kuijper, Carson C Chow, Kevin D Hall, and Margreet R Olthof. 2016.
//      “Impact of Masked Replacement of Sugar-Sweetened with Sugar-Free Beverages on Body Weight Increases with Initial Bmi:
//      Secondary Analysis of Data from an 18 Month Double–Blind Trial in Children.” PloS One 11 (7). Public Library of Science: e0159771.
//
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute,
// sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
// is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies
// or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
// BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//----------------------------------------------------------------------------------------


#include "child_weight.h"

//Default (classic) constructor for energy matrix
Child::Child(NumericVector input_age, NumericVector input_sex, NumericVector input_FFM, NumericVector input_FM, NumericMatrix input_EIntake,
             double input_dt, bool checkValues){
    age   = input_age;
    sex   = input_sex;
    FM    = input_FM;
    FFM   = input_FFM;
    dt    = input_dt;
    EIntake = input_EIntake;
    check = checkValues;
    generalized_logistic = false;
    build();
}

//Constructor which uses Richard's curve with the parameters of https://en.wikipedia.org/wiki/Generalised_logistic_function
Child::Child(NumericVector input_age, NumericVector input_sex, NumericVector input_FFM, NumericVector input_FM, double input_K,
             double input_Q, double input_A, double input_B, double input_nu, double input_C, 
             double input_dt, bool checkValues){
    age   = input_age;
    sex   = input_sex;
    FM    = input_FM;
    FFM   = input_FFM;
    dt    = input_dt;
    K_logistic = input_K;
    A_logistic = input_A;
    Q_logistic = input_Q;
    B_logistic = input_B;
    nu_logistic = input_nu;
    C_logistic = input_C;
    check = checkValues;
    generalized_logistic = true;
    build();
}

Child::~Child(void){
    
}

void Child::build(){
    getParameters();
}

//General function for expressing growth and eb terms
NumericVector Child::general_ode(NumericVector t, NumericVector input_A, NumericVector input_B,
                                 NumericVector input_D, NumericVector input_tA,
                                 NumericVector input_tB, NumericVector input_tD,
                                 NumericVector input_tauA, NumericVector input_tauB,
                                 NumericVector input_tauD){
    
    return input_A*exp(-(t-input_tA)/input_tauA ) +
            input_B*exp(-0.5*pow((t-input_tB)/input_tauB,2)) +
            input_D*exp(-0.5*pow((t-input_tD)/input_tauD,2));
}

NumericVector Child::Growth_dynamic(NumericVector t){
    return general_ode(t, A, B, D, tA, tB, tD, tauA, tauB, tauD);
}

NumericVector Child::Growth_impact(NumericVector t){
    return general_ode(t, A1, B1, D1, tA1, tB1, tD1, tauA1, tauB1, tauD1);
}

NumericVector Child::EB_impact(NumericVector t){
    return general_ode(t, A_EB, B_EB, D_EB, tA_EB, tB_EB, tD_EB, tauA_EB, tauB_EB, tauD_EB);
}

NumericVector Child::cRhoFFM(NumericVector input_FFM){
    return 4.3*input_FFM + 837.0;
}

NumericVector Child::cP(NumericVector FFM, NumericVector FM){
    NumericVector rhoFFM = cRhoFFM(FFM);
    NumericVector C      = 10.4 * rhoFFM / rhoFM;
    return C/(C + FM);
}

NumericVector Child::Delta(NumericVector t){
    return deltamin + (deltamax - deltamin)*(1.0 / (1.0 + pow((t / P),h)));
}

NumericVector Child::FFMReference(NumericVector t){ 
  /*  return ffm_beta0 + ffm_beta1*t; */
   NumericMatrix ffm_ref(17,nind);
    ffm_ref(0,_)   = 10.134*(1-sex)+9.477*sex;
    ffm_ref(1,_)   = 12.099*(1 - sex) + 11.494*sex;
    ffm_ref(2,_)   = 14.0*(1 - sex) + 13.2*sex;
    ffm_ref(3,_)   = 16.0*(1 - sex) + 14.7*sex;
    ffm_ref(4,_)   = 17.9*(1 - sex) + 16.3*sex;
    ffm_ref(5,_)   = 19.9*(1 - sex) + 18.2*sex;
    ffm_ref(6,_)   = 22.0*(1 - sex) + 20.5*sex;
    ffm_ref(7,_)   = 24.4*(1 - sex) + 23.3*sex;
    ffm_ref(8,_)   = 27.5*(1 - sex) + 26.4*sex;
    ffm_ref(9,_)   = 29.5*(1 - sex) + 28.5*sex;
    ffm_ref(10,_)  = 33.2*(1 - sex) + 32.4*sex;
    ffm_ref(11,_)  = 38.1*(1 - sex) + 36.1*sex;
    ffm_ref(12,_)  = 43.6*(1 - sex) + 38.9*sex;
    ffm_ref(13,_)  = 49.1*(1 - sex) + 40.7*sex;
    ffm_ref(14,_)  = 54.0*(1 - sex) + 41.7*sex;
    ffm_ref(15,_)  = 57.7*(1 - sex) + 42.3*sex;
    ffm_ref(16,_)  = 60.0*(1 - sex) + 42.6*sex;
 
 NumericVector ffm_ref_t(nind);
 int jmin;
 int jmax;
 double diff;
 for(int i=0;i<nind;i++){
  if(t(i)>=18.0){
     ffm_ref_t(i)=ffm_ref(16,i);
  }else{
   jmin=floor(t(i));
   jmin=std::max(jmin,2);
   jmin=jmin-2;
   jmax= std::min(jmin+1,17);
   diff= t(i)-floor(t(i));
   ffm_ref_t(i)=ffm_ref(jmin,i)+diff*(ffm_ref(jmax,i)-ffm_ref(jmin,i));
    }
    }
  return ffm_ref_t;
}

NumericVector Child::FMReference(NumericVector t){
   /* return fm_beta0 + fm_beta1*t;*/
    NumericMatrix fm_ref(17,nind);
    fm_ref(0,_)   = 2.456*(1-sex)+ 2.433*sex;
    fm_ref(1,_)   = 2.576*(1 - sex) + 2.606*sex;
    fm_ref(2,_)   = 2.7*(1 - sex) + 2.8*sex;
    fm_ref(3,_)   = 2.7*(1 - sex) + 2.9*sex;
    fm_ref(4,_)   = 2.8*(1 - sex) + 3.2*sex;
    fm_ref(5,_)   = 2.9*(1 - sex) + 3.7*sex;
    fm_ref(6,_)   = 3.3*(1 - sex) + 4.3*sex;
    fm_ref(7,_)   = 3.7*(1 - sex) + 5.2*sex;
    fm_ref(8,_)   = 4.8*(1 - sex) + 7.2*sex;
    fm_ref(9,_)   = 5.9*(1 - sex) + 8.5*sex;
    fm_ref(10,_)  = 6.7*(1 - sex) + 9.2*sex;
    fm_ref(11,_)  = 7.0*(1 - sex) + 10.0*sex;
    fm_ref(12,_)  = 7.2*(1 - sex) + 11.3*sex;
    fm_ref(13,_)  = 7.5*(1 - sex) + 12.8*sex;
    fm_ref(14,_)  = 8.0*(1 - sex) + 14.0*sex;
    fm_ref(15,_)  = 8.4*(1 - sex) + 14.3*sex;
    fm_ref(16,_)  = 8.8*(1 - sex) + 14.3*sex;
 NumericVector fm_ref_t(nind);
 int jmin;
 int jmax;
 double diff;
 for(int i=0;i<nind;i++){
  if(t(i)>=18.0){
     fm_ref_t(i)=fm_ref(16,i);
  }else{
   jmin=floor(t(i));
   jmin=std::max(jmin,2);
   jmin=jmin-2;
   jmax= std::min(jmin+1,17);
   diff= t(i)-floor(t(i));
   fm_ref_t(i)=fm_ref(jmin,i)+diff*(fm_ref(jmax,i)-fm_ref(jmin,i));
  } 
}
  return fm_ref_t;
}

NumericVector Child::IntakeReference(NumericVector t){
    NumericVector EB      = EB_impact(t);
    NumericVector FFMref  = FFMReference(t);
    NumericVector FMref   = FMReference(t);
    NumericVector delta   = Delta(t);
    NumericVector growth  = Growth_dynamic(t);
    NumericVector p       = cP(FFMref, FMref);
    NumericVector rhoFFM  = cRhoFFM(FFMref);
    return EB + K + (22.4 + delta)*FFMref + (4.5 + delta)*FMref +
                230.0/rhoFFM*(p*EB + growth) + 180.0/rhoFM*((1-p)*EB-growth);
}

NumericVector Child::Expenditure(NumericVector t, NumericVector FFM, NumericVector FM){
    NumericVector delta     = Delta(t);
    NumericVector Iref      = IntakeReference(t);
    NumericVector Intakeval = Intake(t);
    NumericVector DeltaI    = Intakeval - Iref;
    NumericVector p         = cP(FFM, FM);
    NumericVector rhoFFM    = cRhoFFM(FFM);
    NumericVector growth    = Growth_dynamic(t);
    NumericVector Expend    = K + (22.4 + delta)*FFM + (4.5 + delta)*FM +
                                0.24*DeltaI + (230.0/rhoFFM *p + 180.0/rhoFM*(1.0-p))*Intakeval +
                                growth*(230.0/rhoFFM -180.0/rhoFM);
    return Expend/(1.0+230.0/rhoFFM *p + 180.0/rhoFM*(1.0-p));
}

//Rungue Kutta 4 method for Adult
List Child::rk4 (double days){
    
    //Initial time
    NumericMatrix k1, k2, k3, k4;
    
    //Estimate number of elements to loop into
    int nsims = floor(days/dt);
    
    //Create array of states
    NumericMatrix ModelFFM(nind, nsims + 1); //in rcpp
    NumericMatrix ModelFM(nind, nsims + 1); //in rcpp
    NumericMatrix ModelBW(nind, nsims + 1); //in rcpp
    NumericMatrix AGE(nind, nsims + 1); //in rcpp
    NumericVector TIME(nsims + 1); //in rcpp
    
    //Create initial states
    ModelFFM(_,0) = FFM;
    ModelFM(_,0)  = FM;
    ModelBW(_,0)  = FFM + FM;
    TIME(0)  = 0.0;
    AGE(_,0)  = age;
    
    //Loop through all other states
    bool correctVals = true;
    for (int i = 1; i <= nsims; i++){
        /*if (check){
            for (int k = 0; k < nind; k++){
                //Need to correct in windows there is no isfinite.
                if(ModelFFM(k,i-1)<=0|| !isfinite(ModelFFM(k,i-1)) || ModelFM(k,i-1)<=0|| !isfinite(ModelFM(k,i-1))){
                    Rcout << "First error in person "<< k+1 <<std::endl;
                    correctVals = false;
                    break;
                }
            }
        }
        
        if (!correctVals) {
            break;
        }*/
        
        
        //Rungue kutta 4 (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
        k1 = dMass(AGE(_,i-1), ModelFFM(_,i-1), ModelFM(_,i-1));
        k2 = dMass(AGE(_,i-1) + 0.5 * dt/365.0, ModelFFM(_,i-1) + 0.5 * k1(0,_), ModelFM(_,i-1) + 0.5 * k1(1,_));
        k3 = dMass(AGE(_,i-1) + 0.5 * dt/365.0, ModelFFM(_,i-1) + 0.5 * k2(0,_), ModelFM(_,i-1) + 0.5 * k2(1,_));
        k4 = dMass(AGE(_,i-1) + dt/365.0, ModelFFM(_,i-1) + k3(0,_), ModelFM(_,i-1) +  k3(1,_));
        
        //Update of function values
        //Note: The dt is factored from the k1, k2, k3, k4 defined on the Wikipedia page and that is why
        //      it appears here.
        ModelFFM(_,i) = ModelFFM(_,i-1) + dt*(k1(0,_) + 2.0*k2(0,_) + 2.0*k3(0,_) + k4(0,_))/6.0;        //ffm
        ModelFM(_,i)  = ModelFM(_,i-1)  + dt*(k1(1,_) + 2.0*k2(1,_) + 2.0*k3(1,_) + k4(1,_))/6.0;        //fm
        
        //Update weight
        ModelBW(_,i) = ModelFFM(_,i) + ModelFM(_,i);
        
        //Update TIME(i-1)
        TIME(i) = TIME(i-1) + dt; // Currently time counts the time (days) passed since start of model
        
        //Update AGE variable
        AGE(_,i) = AGE(_,i-1) + dt/365.0; //Age is variable in years
    }
    
    return List::create(Named("Time") = TIME,
                        Named("Age") = AGE,
                        Named("Fat_Free_Mass") = ModelFFM,
                        Named("Fat_Mass") = ModelFM,
                        Named("Body_Weight") = ModelBW,
                        Named("Correct_Values")=correctVals,
                        Named("Model_Type")="Children");


}

NumericMatrix  Child::dMass (NumericVector t, NumericVector FFM, NumericVector FM){
    
    NumericMatrix Mass(2, nind); //in rcpp;
    NumericVector rhoFFM    = cRhoFFM(FFM);
    NumericVector p         = cP(FFM, FM);
    NumericVector growth    = Growth_dynamic(t);
    NumericVector expend    = Expenditure(t, FFM, FM);
    Mass(0,_)               = (1.0*p*(Intake(t) - expend) + growth)/rhoFFM;    // dFFM
    Mass(1,_)               = ((1.0 - p)*(Intake(t) - expend) - growth)/rhoFM; //dFM
    return Mass;
    
}

void Child::getParameters(void){
    
    //General constants
    rhoFM    = 9.4*1000.0;
    deltamin = 10.0;
    P        = 12.0;
    h        = 10.0;
    
    //Number of individuals
    nind     = age.size();
    
    //Sex specific constants
    ffm_beta0 = 2.9*(1 - sex)  + 3.8*sex;
    ffm_beta1 = 2.9*(1 - sex)  + 2.3*sex;
    fm_beta0  = 1.2*(1 - sex)  + 0.56*sex;
    fm_beta1  = 0.41*(1 - sex) + 0.74*sex;
    K         = 800*(1 - sex)  + 700*sex;
    deltamax  = 19*(1 - sex)   + 17*sex;
    A         = 3.2*(1 - sex)  + 2.3*sex;
    B         = 9.6*(1 - sex)  + 8.4*sex;
    D         = 10.1*(1 - sex) + 1.1*sex;
    tA        = 4.7*(1 - sex)  + 4.5*sex;       //years
    tB        = 12.5*(1 - sex) + 11.7*sex;      //years
    tD        = 15.0*(1-sex)   + 16.2*sex;      //years
    tauA      = 2.5*(1 - sex)  + 1.0*sex;       //years
    tauB      = 1.0*(1 - sex)  + 0.9*sex;       //years
    tauD      = 1.5*(1 - sex)  + 0.7*sex;       //years
    A_EB      = 7.2*(1 - sex)  + 16.5*sex;
    B_EB      = 30*(1 - sex)   + 47.0*sex;
    D_EB      = 21*(1 - sex)   + 41.0*sex;
    tA_EB     = 5.6*(1 - sex)  + 4.8*sex;
    tB_EB     = 9.8*(1 - sex)  + 9.1*sex;
    tD_EB     = 15.0*(1 - sex) + 13.5*sex;
    tauA_EB   = 15*(1 - sex)   + 7.0*sex;
    tauB_EB   = 1.5*(1 -sex)   + 1.0*sex;
    tauD_EB   = 2.0*(1 - sex)  + 1.5*sex;
    A1        = 3.2*(1 - sex)  + 2.3*sex;
    B1        = 9.6*(1 - sex)  + 8.4*sex;
    D1        = 10.0*(1 - sex) + 1.1*sex;
    tA1       = 4.7*(1 - sex)  + 4.5*sex;
    tB1       = 12.5*(1 - sex) + 11.7*sex;
    tD1       = 15.0*(1 - sex) + 16.0*sex;
    tauA1     = 1.0*(1 - sex)  + 1.0*sex;
    tauB1     = 0.94*(1 - sex) + 0.94*sex;
    tauD1     = 0.69*(1 - sex) + 0.69*sex;
}


//Intake in calories
NumericVector Child::Intake(NumericVector t){
    if (generalized_logistic) {
        return A_logistic + (K_logistic - A_logistic)/pow(C_logistic + Q_logistic*exp(-B_logistic*t), 1/nu_logistic); //t in years
    } else {
        int timeval = floor(365.0*(t(0) - age(0))/dt); //Example: Age: 6 and t: 7.1 => timeval = 401 which corresponds to the 401 entry of matrix
        return EIntake(timeval,_);
    }
    
}


