//
//  child_weight_wrapper.cpp
//
//  This is a function that uses Rcpp to return
//  weight change for children using the dynamic
//  weight model by Kevin D. Hall et al.
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


