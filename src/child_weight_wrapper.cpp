//
//  child_weight_wrapper.cpp
//
//  This is a function that uses Rcpp to return
//  weight change for Children using the dynamic
//  weight model by Kevin D. Hall et al.
//
//  Input:
//  age  .-  Years since individual first arrived to Earth
//  sex  .-  Either 1 = "female" or 0 = "male"
//  FFM  .-  Fat Free Mass (kg) of the individual
//  FM   .-  Fat Mass (kg) of the individual
//  days .-  Days to model (integer)
//
//  Note:
//  Weight = FFM + FM. No extracellular fluid or glycogen is considered
//  Please see child_weight.hpp for additional information
//
//  Authors:
//  Dalia Camacho-García-Formentí
//  Rodrigo Zepeda-Tello
//
//  Copyright: Instituto Nacional de Salud Pública de México

#include <Rcpp.h>
#include "child_weight.h"

// [[Rcpp::export]]
List child_weight_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, NumericMatrix input_EIntake, double days, double dt, bool checkValues){
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM, input_EIntake, dt, checkValues);
    
    //Run model using RK4
    return Person.rk4(days);
    
}

// [[Rcpp::export]]
List child_weight_wrapper_richardson(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, double K, double Q, double A, double B, double nu, double C, double days, double dt, bool checkValues){
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM, K, Q, A, B, nu, C, dt, checkValues);
    
    //Run model using RK4
    return Person.rk4(days);
    
}

// [[Rcpp::export]]
NumericMatrix intake_reference_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, double days,  double dt){
    
    //Energy intake input empty matrix
    NumericMatrix EI(1,1);
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM, EI, dt, false);
    
    //Energy matrix
    NumericMatrix EnergyIntake(age.size(), floor(days/dt) + 1);
    
    //Get energy matrix
    for (double i = 0; i < floor(days/dt) + 1; i++){
        EnergyIntake(_,i) = Person.IntakeReference(age + dt*i/365.0);
    }
    
    return EnergyIntake;
    
}

// [[Rcpp::export]]
List mass_reference_wrapper(NumericVector age, NumericVector sex){
    
    //Input empty matrices
    NumericMatrix EI(1,1);
    NumericMatrix inputFM(1,1);
    NumericMatrix inputFFM(1,1);
    
    //Create new adult with characteristics
    Child Person (age,  sex, inputFFM, inputFM, EI, 1.0, false);
    
    //Energy matrix
    NumericVector FM  = Person.FMReference(age);
    NumericVector FFM = Person.FFMReference(age);
    
    return List::create(Named("FM")  = FM,
                        Named("FFM") = FFM);
    
}


