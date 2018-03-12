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
List child_weight_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, NumericMatrix input_EIntake, double days, bool checkValues){
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM, input_EIntake, checkValues);
    
    //Run model using RK4
    return Person.rk4(days);
    
}

// [[Rcpp::export]]
NumericMatrix intake_reference_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, double days){
    
    //Energy intake input empty matrix
    NumericMatrix EI(1,1);
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM, EI, false);
    
    //Energy matrix
    NumericMatrix EnergyIntake(age.size(), (days));
    
    //Get energy matrix
    for (double i = 0; i <= (days); i++){
        EnergyIntake(_,i) = Person.IntakeReference(age + i/365.0);
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
    Child Person (age,  sex, inputFFM, inputFM, EI, false);
    
    //Energy matrix
    NumericVector FM  = Person.FMReference(age);
    NumericVector FFM = Person.FFMReference(age);
    
    return List::create(Named("FM")  = FM,
                        Named("FFM") = FFM);
    
}


