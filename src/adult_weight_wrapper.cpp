//
//  adult_weight_wrapper.cpp
//  
//  This is a function that uses Rcpp to return
//  weight change for adults using the dynamic
//  weight model by Kevin D. Hall et al.
//
//  Input:
//  bw              .-  Body weight (kg).
//  ht              .-  Height (m).
//  age             .-  Years since individual first arrived to Earth.
//  sex             .-  Either 1 = "female" or 0 = "male".
//  EIchange        .-  Change in energy intake (kcal).
//  NAchange        .-  Change in sodium consumption (mg).
//  PAL             .-  Physical activity level. (Between 1.4 and 2.4)
//  pcarb           .-  Proportion of carbohydrates from diet throughout the time the model runs.
//  pcarb_baseline  .-  Proportion of carbohydrates from diet at baseline.
//  dt              .-  Time step used to solve the ODE system numerically.
//  extradata       .-  Energy intake at baseline (kcal).
//  checkValues     .-  Verify values of lean and fat mass are possible values. (Not infinite nor NA)
//  isEnergy        .-  Boolean to determine if energy intake at baseline is given
//  input_EI        .-  Energy intake (kcal). 
//  input_fat       .-  Fat Mass (kg) of the individual.
//
//  Authors:
//  Dalia Camacho-García-Formentí
//  Rodrigo Zepeda-Tello
//
// References:
//
//  Chow, Carson C, and Kevin D Hall. 2008. “The Dynamics of Human Body Weight Change.” PLoS Comput Biol 4 (3):e1000045.
//
//  Hall, Kevin D. 2010. “Predicting Metabolic Adaptation, Body Weight Change, and Energy Intake in Humans.”
//      American Journal of Physiology-Endocrinology and Metabolism 298 (3). Am Physiological Soc: E449–E466.
//
//  Hall, Kevin D, and Peter N Jordan. 2008. “Modeling Weight-Loss Maintenance to Help Prevent Body Weight Regain.”
//      The American Journal of Clinical Nutrition 88 (6). Am Soc Nutrition: 1495–1503.
//
//  Hall, Kevin D, Gary Sacks, Dhruva Chandramohan, Carson C Chow, Y Claire Wang, Steven L Gortmaker, and Boyd A Swinburn. 2011.
//      “Quantification of the Effect of Energy Imbalance on Bodyweight.” The Lancet 378 (9793). Elsevier: 826–37.
//
//  Mifflin, Mark D, Sachiko T St Jeor, Lisa A Hill, Barbara J Scott, Sandra A Daugherty, and YO Koh. 1990.
//      “A New Predictive Equation for Resting Energy Expenditure in Healthy Individuals.” The American Journal of Clinical Nutrition 51 (2).
//      Am Soc Nutrition: 241–47.
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
#include "adult_weight.h"

// [[Rcpp::export]]
List adult_weight_wrapper(NumericVector bw, NumericVector ht, NumericVector age,
                          NumericVector sex, NumericMatrix EIchange,
                          NumericMatrix NAchange, NumericVector PAL,
                          NumericVector pcarb_base, NumericVector pcarb, double dt,
                          double days, bool checkValues){
    
    //Create new adult with characteristics
    Adult Person (bw,  ht, age, sex, EIchange, NAchange, PAL, pcarb,  pcarb_base, dt, checkValues);
    
    //Run model using RK4
    return Person.rk4(days);
    
}

// [[Rcpp::export]]
List adult_weight_wrapper_EI(NumericVector bw, NumericVector ht, NumericVector age,
                          NumericVector sex, NumericMatrix EIchange,
                          NumericMatrix NAchange, NumericVector PAL,
                          NumericVector pcarb_base, NumericVector pcarb, double dt,
                             NumericVector extradata, double days, bool checkValues, bool isEnergy){
    
    //Create new adult with characteristics
    Adult Person (bw,  ht, age, sex, EIchange, NAchange, PAL, pcarb,  pcarb_base, dt, extradata, checkValues, isEnergy);
    
    //Run model using RK4
    return Person.rk4(days);
    
}

// [[Rcpp::export]]
List adult_weight_wrapper_EI_fat(NumericVector bw, NumericVector ht, NumericVector age,
                             NumericVector sex, NumericMatrix EIchange,
                             NumericMatrix NAchange, NumericVector PAL,
                             NumericVector pcarb_base, NumericVector pcarb, double dt,
                             NumericVector input_EI, NumericVector input_fat,
                                 double days, bool checkValues){
    
    //Create new adult with characteristics
    Adult Person (bw,  ht, age, sex, EIchange, NAchange, PAL, pcarb,  pcarb_base, dt, input_EI, input_fat, checkValues);
    
    //Run model using RK4
    return Person.rk4(days);
    
}
