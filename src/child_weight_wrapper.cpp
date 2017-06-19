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
#include "child_weight.hpp"

// [[Rcpp::export]]
List child_weight_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, double days){
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM);
    
    //Run model using RK4
    return Person.rk4(days);
    
}

