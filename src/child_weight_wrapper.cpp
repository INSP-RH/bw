//
//  child_weight_wrapper.cpp
//  
//
//  Created by Rodrigo Zepeda Tello on 06/06/17.
//
//
#include <Rcpp.h>
#include "child_weight.hpp"

// [[Rcpp::export]]
List child_weight_wrapper(NumericVector age, NumericVector sex, NumericVector FFM, NumericVector FM, double days){
    
    //Create new adult with characteristics
    Child Person (age,  sex, FFM, FM);
    
    //Run model using RK4
    return Person.rk4(days);
    
}

