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




#include <Rcpp.h>
#include "child_weight.h"

// [[Rcpp::export]]
List child_weight_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, NumericMatrix input_EIntake, double days, double dt, bool checkValues){
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM, input_EIntake, dt, checkValues);
    
    //Run model using RK4
    return Person.rk4(days - 1); //days - 1 to account for extra day (c++ indexing starts in 0; R in 1)
    
}

// [[Rcpp::export]]
List child_weight_wrapper_richardson(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, double K, double Q, double A, double B, double nu, double C, double days, double dt, bool checkValues){
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM, K, Q, A, B, nu, C, dt, checkValues);
    
    //Run model using RK4
    return Person.rk4(days - 1); //days - 1 to account for extra day (c++ indexing starts in 0; R in 1)
    
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


