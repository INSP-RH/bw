//
//  Child.cpp
//  AdultHall
//
//  Created by Rodrigo Zepeda Tello on 26/05/17.
//  Copyright © 2017 Instituto Nacional de Salud Publica. All rights reserved.
//

#include "child_weight.h"
#include <stdio.h>     
#include <math.h>
#include <vector>
#include <iostream>   
#include <algorithm>  
using namespace std;

Child::Child(NumericVector input_age, NumericVector input_sex, NumericVector input_FFM, NumericVector input_FM, NumericMatrix input_EIntake, bool checkValues){
    age   = input_age;
    sex   = input_sex;
    FM    = input_FM;
    FFM   = input_FFM;
    EIntake = input_EIntake;
    check = checkValues;
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

NumericVector Child::FFMReference(NumericVector t){ //Linear model may be wrong
        return ffm_beta0 + ffm_beta1*t;
}

NumericVector Child::FMReference(NumericVector t){
        return fm_beta0 + fm_beta1*t;
}

/*NumericVector Child::IntakeReference(NumericVector t){
    NumericVector EB      = EB_impact(t);
    NumericVector FFMref  = FFMReference(t);
    NumericVector FMref   = FMReference(t);
    NumericVector delta   = Delta(t);
    NumericVector growth  = Growth_dynamic(t);
    NumericVector p       = cP(FFMref, FMref);
    NumericVector rhoFFM  = cRhoFFM(FFMref);
    return EB + K + (22.4 + delta)*FFMref + (4.5 + delta)*FMref +
                230.0/rhoFFM*(p*EB + growth) + 180.0/rhoFM*((1-p)*EB-growth);
}*/
NumericVector Child::IntakeReference(NumericVector t){
    double req_t[nind];
    int j;
    for(int i=0;i<=nind;i++){
        j=std::min(floor(t(i)),18.0);
        j=std::max(j,1.0);
        j=j-1;
        req_t[i]=req(j,i);
    }
   //return req(std::min(floor(t),18))
   //double req[]  = {7, 5, 16, 8};
   return req_t;

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
                                0.24*DeltaI + (230.0/rhoFFM *p + 180/rhoFM*(1-p))*Intakeval +
                                growth*(230.0/rhoFFM -180.0/rhoFM);
    
    return Expend/(1+230/rhoFFM *p + 180/rhoFM*(1-p));
}

//Rungue Kutta 4 method for Adult
List Child::rk4 (double days){
    
    //Set dt to 1
    double dt = 1.0;
    //double dt = 1.0/365.0;
    //Initial time
    NumericMatrix k1, k2, k3, k4;
    
    //Estimate number of elements to loop into
    int nsims = days; //Need to control here for > 18 yrs
    
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
        if (check){
            for (int k = 0; k < nind; k++){
                /* //Need to correct in windows there is no isfinite.
                if(ModelFFM(k,i-1)<=0|| !isfinite(ModelFFM(k,i-1)) || ModelFM(k,i-1)<=0|| !isfinite(ModelFM(k,i-1))){
                    Rcout << "First error in person "<< k+1 <<std::endl;
                    correctVals = false;
                    break;
                }*/
            }
        }
        
        if (!correctVals) {
            break;
        }
        
        
        //Rungue kutta 4 (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
        k1 = dMass(AGE(_,i-1), ModelFFM(_,i-1), ModelFM(_,i-1));
        k2 = dMass(AGE(_,i-1) + 0.5 * dt, ModelFFM(_,i-1) + 0.5 * dt * k1(0,_), ModelFM(_,i-1) + 0.5 * dt * k1(1,_));
        k3 = dMass(AGE(_,i-1) + 0.5 * dt, ModelFFM(_,i-1) + 0.5 * dt * k2(0,_), ModelFM(_,i-1) + 0.5 * dt * k2(1,_));
        k4 = dMass(AGE(_,i-1) + dt, ModelFFM(_,i-1) +dt * k3(0,_), ModelFM(_,i-1) + dt * k3(1,_));
        
        //Update FFM and FM
        ModelFFM(_,i) = ModelFFM(_,i-1) + dt * (k1(0,_) + 2.0*k2(0,_) + 2.0*k3(0,_) + k4(0,_))/6.0;        //ffm
        ModelFM(_,i)  = ModelFM(_,i-1) + dt * (k1(1,_) + 2.0*k2(1,_) + 2.0*k3(1,_) + k4(1,_))/6.0;        //fm
        
        //Update weight
        ModelBW(_,i) = ModelFFM(_,i) + ModelFM(_,i);

        //Update TIME(i-1)
        TIME(i) = TIME(i-1) + 1;
        
        //Update AGE variable
        AGE(_,i) = AGE(_,i-1) + dt/365;
        //AGE(_,i) = AGE(_,i-1) + dt;
    }
    
    return List::create(Named("Time") = TIME,
                        Named("Age") = AGE,
                        Named("Fat_Free_Mass") = ModelFFM,
                        Named("Fat_Mass") = ModelFM,
                        Named("Body_Weight") = ModelBW,
                        Named("Correct_Values")=correctVals);


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
    rhoFM    = 9.4*1000;
    deltamin = 10.0;
    P        = 12.0;
    h        = 10.0;
    
    //Number of individuals
    nind     = age.size();
    
    //Sex specific constants
    ffm_beta0 = 2.9*(1 - sex) + 3.8*sex;
    ffm_beta1 = 2.9*(1 - sex) + 2.3*sex;
    fm_beta0  = 1.2*(1 - sex) + 0.56*sex;
    fm_beta1  = 0.41*(1 - sex) + 0.74*sex;
    K         = 800*(1 - sex) + 700*sex;
    deltamax  = 19*(1 - sex) + 17*sex;
    A         = 3.2*(1 - sex) + 2.3*sex;
    B         = 9.6*(1 - sex) + 8.4*sex;
    D         = 10.1*(1 - sex) + 1.1*sex;
    tA        = 4.7*(1 - sex) + 4.5*sex;
    tB        = 12.5*(1 - sex) + 11.7*sex;
    tD        = 15.0*(1-sex) + 16.2*sex;
    tauA      = 2.5*(1 - sex) + 1.0*sex;
    tauB      = 1.0*(1 - sex) + 0.9*sex;
    tauD      = 1.5*(1 - sex) + 0.7*sex;
    A_EB      = 7.2*(1 - sex) + 16.5*sex;
    B_EB      = 30*(1 - sex) + 47.0*sex;
    D_EB      = 21*(1 - sex) + 41.0*sex;
    tA_EB     = 5.6*(1 - sex) + 4.8*sex;
    tB_EB     = 9.8*(1 - sex) + 9.1*sex;
    tD_EB     = 15.0*(1 - sex) + 13.5*sex;
    tauA_EB   = 15*(1 - sex) + 7.0*sex;
    tauB_EB   = 1.5*(1 -sex) + 1.0*sex;
    tauD_EB   = 2.0*(1 - sex) + 1.5*sex;
    A1        = 3.2*(1 - sex) + 2.3*sex;
    B1        = 9.6*(1 - sex) + 8.4*sex;
    D1        = 10.0*(1 - sex) + 1.1*sex;
    tA1       = 4.7*(1 - sex) + 4.5*sex;
    tB1       = 12.5*(1 - sex) + 11.7*sex;
    tD1       = 15.0*(1 - sex) + 16.0*sex;
    tauA1     = 1.0*(1 - sex) + 1.0*sex;
    tauB1     = 0.94*(1 - sex) + 0.94*sex;
    tauD1     = 0.69*(1 - sex) + 0.69*sex;
    
    req(0,_)   = 948*(1-sex)+865.0*sex;
    req(1,_)   =1129.0*(1 - sex) + 1047.0*sex;
    req(2,_)   =1252.0*(1 - sex) + 1156.0*sex;
    req(3,_)   =1360.0*(1 - sex) + 1241.0*sex;
    req(4,_)   =1467.0*(1 - sex) + 1330.0*sex;
    req(5,_)   =1573.0*(1 - sex) + 1428.0*sex;
    req(6,_)   =1692.0*(1 - sex) + 1554.0*sex;
    req(7,_)   =1830.0*(1 - sex) + 1698.0*sex;
    req(8,_)   =1978.0*(1 - sex) + 1854.0*sex;
    req(9,_)   =2150.0*(1 - sex) + 2006.0*sex;
    req(10,_)  =2341.0*(1 - sex) + 2149.0*sex;
    req(11,_)  =2548.0*(1 - sex) + 2276.0*sex;
    req(12,_)  =2770.0*(1 - sex) + 2379.0*sex;
    req(13,_)  =2990.0*(1 - sex) + 2249.0*sex;
    req(14,_)  =3178.0*(1 - sex) + 2491.0*sex;
    req(15,_)  =3322.0*(1 - sex) + 2503.0*sex;
    req(16,_)  =3410.0*(1 - sex) + 2503.0*sex;
}


//Intake in calories
NumericVector Child::Intake(NumericVector t){
    double timeval = t(0) - age(0);
    return EIntake(floor(timeval),_);
}

//This needs to be an R inputed function
/*NumericVector Child::Intake(NumericVector t){
    return IntakeReference(t);
}*/

