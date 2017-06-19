//
//  adult_weight_wrapper.cpp
//  
//
//  Created by Rodrigo Zepeda Tello on 06/06/17.
//
//
#include <Rcpp.h>
#include "adult_weight.hpp"

// [[Rcpp::export]]
List adult_weight_wrapper(NumericVector bw, NumericVector ht, NumericVector age,
                          NumericVector sex, NumericMatrix EIchange,
                          NumericMatrix NAchange, NumericVector PAL,
                          NumericVector pcarb_base, NumericVector pcarb,
                          double days, bool checkValues){
    
    //Create new adult with characteristics
    Adult Person (bw,  ht, age, sex, EIchange, NAchange, PAL, pcarb,  pcarb_base, true);
    
    //Run model using RK4
    return Person.rk4(days);
    
}

// [[Rcpp::export]]
List adult_weight_wrapper_EI(NumericVector bw, NumericVector ht, NumericVector age,
                          NumericVector sex, NumericMatrix EIchange,
                          NumericMatrix NAchange, NumericVector PAL,
                          NumericVector pcarb_base, NumericVector pcarb,
                             NumericVector input_EIntake, double days, bool checkValues){
    
    //Create new adult with characteristics
    Adult Person (bw,  ht, age, sex, EIchange, NAchange, PAL, pcarb,  pcarb_base, input_EIntake, true);
    
    //Run model using RK4
    return Person.rk4(days);
    
}
