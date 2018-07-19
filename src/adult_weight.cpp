//
//  adult_weight.cpp
//  
//
//  This is a function that calculates weight change for adults using the dynamic
//  weight model by Kevin D. Hall et al. and Runge Kutta method to solve the ODE system.
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

#include "adult_weight.h"

//Default Constructor for an Adult.
Adult::Adult(NumericVector weight, NumericVector height, NumericVector age_yrs,
             NumericVector sexstring, NumericMatrix input_EIchange,
             NumericMatrix input_NAchange, NumericVector physicalactivity,
             NumericVector percentc, NumericVector percentb, double input_dt, bool checkValues){
    
    //Build model from parameters
    build(weight, height, age_yrs, sexstring, input_EIchange, input_NAchange,
          physicalactivity, percentc, percentb, input_dt, checkValues);
    
}

//Constructor with energy intake vector or fat vector
Adult::Adult(NumericVector weight, NumericVector height, NumericVector age_yrs,
             NumericVector sexstring, NumericMatrix input_EIchange,
             NumericMatrix input_NAchange, NumericVector physicalactivity,
             NumericVector percentc, NumericVector percentb, double input_dt, NumericVector extradata,
             bool checkValues, bool isEnergy){
    
    
    //Build model from parameters
    build(weight, height, age_yrs, sexstring, input_EIchange, input_NAchange,
          physicalactivity, percentc, percentb, input_dt, extradata, checkValues, isEnergy);
    
}

//Constructor with energy intake vector and fat vector
Adult::Adult(NumericVector weight, NumericVector height, NumericVector age_yrs,
             NumericVector sexstring, NumericMatrix input_EIchange,
             NumericMatrix input_NAchange, NumericVector physicalactivity,
             NumericVector percentc, NumericVector percentb, double input_dt, NumericVector input_EI,
             NumericVector input_fat, bool checkValues){
    
    
    //Build model from parameters
    build(weight, height, age_yrs, sexstring, input_EIchange, input_NAchange,
          physicalactivity, percentc, percentb, input_dt ,input_EI, input_fat, checkValues);
    
}

//Function to build a new Adult
void Adult::build(NumericVector weight, NumericVector height, NumericVector age_yrs,
                  NumericVector sexstring, NumericMatrix input_EIchange,
                  NumericMatrix input_NAchange, NumericVector physicalactivity,
                  NumericVector percentc, NumericVector percentb, double input_dt, bool checkValues){
    
    //Assign parameters
    dt         = input_dt; //Time step set to 1 because of matrix use (each time is a row in EIchange)
    bw         = weight;
    ht         = height;
    age        = age_yrs;
    sex        = sexstring;
    EIchange   = input_EIchange;
    NAchange   = input_NAchange;
    PAL        = physicalactivity;
    pcarb      = percentc;
    pcarb_base = percentb;
    check      = checkValues;
    
    //Get energy
    getParameters();
    getRMR();
    getATinit();
    getECFinit();
    getBaselineMass();
    getCaloricSteadyState();
    getEnergy();
    getDelta();
    getK();
    getCarbConstants();
}

//Function to build a new Adult when input_EIintake and fat are included
void Adult::build(NumericVector weight, NumericVector height, NumericVector age_yrs,
                  NumericVector sexstring, NumericMatrix input_EIchange,
                  NumericMatrix input_NAchange, NumericVector physicalactivity,
                  NumericVector percentc, NumericVector percentb, double input_dt,
                  NumericVector extradata, bool checkValues, bool isEnergy){
    
    //Assign parameters
    dt         = input_dt; //For rk4
    bw         = weight;
    ht         = height;
    age        = age_yrs;
    sex        = sexstring;
    EIchange   = input_EIchange;
    NAchange   = input_NAchange;
    PAL        = physicalactivity;
    pcarb      = percentc;
    pcarb_base = percentb;
    check      = checkValues;
    
    //Get additional information
    getParameters();
    getRMR();
    getATinit();
    getECFinit();
    
    if (isEnergy){
        //Get energy
        EI     = extradata;
        
        //Get bw
        getBaselineMass();
    } else {
        //Get energy
        getCaloricSteadyState();
        getEnergy();
        
        //Get bw
        fat  = extradata;
        lean = bw - (ecfinit + fat + 3.7*G_base);
    }
    
    getDelta();
    getK();
    getCarbConstants();
}

//Function to build a new Adult when input_EIintake is included
void Adult::build(NumericVector weight, NumericVector height, NumericVector age_yrs,
                  NumericVector sexstring, NumericMatrix input_EIchange,
                  NumericMatrix input_NAchange, NumericVector physicalactivity,
                  NumericVector percentc, NumericVector percentb, double input_dt,
                  NumericVector input_EI, NumericVector input_fat, bool checkValues){
    
    //Assign parameters
    dt         = input_dt; //For rk4
    bw         = weight;
    ht         = height;
    age        = age_yrs;
    sex        = sexstring;
    EIchange   = input_EIchange;
    NAchange   = input_NAchange;
    PAL        = physicalactivity;
    pcarb      = percentc;
    pcarb_base = percentb;
    check      = checkValues;
    
    //Get additional information
    getParameters();
    getRMR();
    getATinit();
    getECFinit();
    
    //Inputted ei and fat
    EI     = input_EI;
    fat    = input_fat;
    lean    = bw - (ecfinit + fat + 3.7*G_base);

    getDelta();
    getK();
    getCarbConstants();
}

//Destroyer
Adult::~Adult(void){
    
}

void Adult::getParameters(void){
    
    //Get size of model
    nind    = bw.size();
    
    //Set to true
    
    //Initialize other variables that are not dependent on the individual
    roG     = 4206.501; // 1000*17.6*0.23900573614 #Changed from kjoules to kcals
    Na      = 3220;     // (1000*3.22)#Sodium
    zetaNa  = 3000;
    zetaCI  = 4000;
    roF     = 9440.727; // 1000*39.5*0.23900573614 #Changed from kjoules to kcals
    roL     = 1816.444; // 1000*7.6*0.23900573614  #Changed from kjoules to kcals
    gammaF  = 3.107075; // 13*0.23900573614        #Changed from kjoules to kcals
    gammaL  = 21.98853; // 92*0.23900573614        #Changed from kjoules to kcals
    etaF    = 179.2543; // 750*0.23900573614       #Changed from kjoules to kcals
    etaL    = 229.4455; // 960*0.23900573614       #Changed from kjoules to kcals
    betaTEF = 0.1;
    betaAT  = 0.14;
    tauAT   = 14.0;
    C       = 10.4*(roL/roF);
    alfa1  = -(1 + etaL/roL)*C;     //Auxiliary functions from Pablo
    alfa2  = -(1 + etaF/roF);       //Auxiliary functions from Pablo
    rmrbw  = 9.99;        //Linear regression coefficient for rmr estimation
    rmrage = 4.92;        //Linear regression coefficient for rmr estimation
    rmrht  = 625.0;       //Linear regression coefficient for rmr estimation
    rmr_m  = 5.0;         //Linear regression coefficient for rmr estimation (men)
    rmr_f  = 161.0;       //Linear regression coefficient for rmr estimation (women)
    G_base = NumericVector(nind, 0.5);
}

//Estimation of Resting Metabolic Rate (rmr) in kcal
void Adult::getRMR(void){
    //These equations come from Miffin & St.Jeor
    //Recall that sex = 0 => "male" and sex = 1 => "female"
    rmr = (rmrbw*bw + rmrht*ht - rmrage*age + rmr_m)*(1-sex) +
    (rmrbw*bw + rmrht*ht - rmrage*age - rmr_f)*sex;
}

//Estimation of calories at baseline
void Adult::getCaloricSteadyState(void){
    //These estimation assumes Energy Intake = Energy Expenditure.
    //Energy is returned in kcal
    steadystate = rmr*PAL;
}

void Adult::getATinit(void){
    //Personal communication with Hall: Yes, since the model starts in a state of
    //energy balance, AT(0) = 0.
    atinit = NumericVector(nind, 0.0);
}

//Add energy intake
void Adult::getEnergy(void){
    EI = steadystate;
}

//Calculate parameter delta
void Adult::getDelta(void){
    delta =  ((1.0 - betaTEF)*PAL - 1.0)*rmr/bw;
}

//Get extracellular water by Silva's equation
void Adult::getECFinit(void){
    ecfinit = (0.025*age + 9.57*ht + 0.191*bw - 12.4)*(1.0-sex) + (-4.0 + 5.98*ht + 0.167*bw)*sex;
}


//Estimation of initial fat and lean masses
void Adult::getBaselineMass(void){
    
    fat =  (bw * (0.14 * age + 37.31 * log(bw/( pow (ht,2.0))) - 103.94)/100.0)*(1-sex) +
    (bw * (0.14 * age + 39.96 * log(bw/( pow (ht,2.0))) - 102.01)/100.0)*sex;
    
    
    //Get lean mass:
    //“The initial lean body mass is simply the difference between the initial BW,
    //the initial F, the initial ECF, and the initial G and its associated water.”
    lean = bw - (ecfinit + fat + 3.7*G_base);
}

//Thermal effect of feeding
NumericVector Adult::TEF(double t){
    return betaTEF*deltaEI(t);
}

//Glycogen
NumericVector Adult::dG(double t, NumericVector G){
    return (CI(t) - kG*pow(G, 2.0))/roG;
}

//Adaptive Thermogenesis derivative
NumericVector Adult::dAT(double t, NumericVector AT){
    return (betaAT *deltaEI(t) - AT)*(1.0 /tauAT);
}

//Extracellular fluid derivative
NumericVector Adult::dECF(double t, NumericVector ECF){
    return ( deltaNA(t) - zetaNa*(ECF - ecfinit) - zetaCI*(1.0 - CI(t)/CIb) )/Na;
}

//Carbohydrate constants
void Adult::getCarbConstants(void){
    CIb = pcarb_base * EI;
    kG  = CIb/( pow (G_base, 2.0) );
}


//Carbohydrate intake
NumericVector Adult::CI(double t){
    return pcarb * TotalIntake(t);
}

//Total energy intake
NumericVector Adult::TotalIntake (double t){
    return EI + deltaEI(t);
}

//Get K constant
void Adult::getK(){
    /*
     Hall personnal communication:
     You are definitely on the right track here. The energy expenditure in the baseline energy balanced
     state is given my EE = PAL*RMR, where PAL is the specified baseline physical activity level and the
     RMR is from the Mifflin-St Jeor equation. The physical activity pa- rameter, delta, at the baseline
     steady state is determined by equation 8 and therefore you can solve for K.
     */
    K = (rmr * PAL) - gammaL * lean - gammaF * fat - delta * bw;
}

//Get fat mass as function of lean tissue
NumericVector Adult::fatMass(NumericVector L){
    return fat * exp(roL * (L - lean)/(roF * C));
}

//Lean tissue derivative
NumericVector Adult::dL(double t, NumericVector L, NumericVector G,
                        NumericVector AT, NumericVector ECF){
    return R(t, L, G, AT, ECF)*(C/roL);
}

//R helper for Lean derivative
NumericVector Adult::R(double t, NumericVector L, NumericVector G,
                       NumericVector AT, NumericVector ECF){
    NumericVector F      = fatMass(L);
    NumericVector weight = L + F + ECF + 3.7*(G);
    NumericVector R3     = K + delta*weight + TEF(t) + AT - TotalIntake(t) + dG(t, G);
    return (R3 + gammaL*L + gammaF*F)/(alfa1 + alfa2*F);
}

//Classifier for bMI
StringVector Adult::BMIClassifier(NumericVector BMI){
    StringVector classification(BMI.size());
    /*for(int i = 0; i < BMI.size(); i++){
        classification(i) = "Unknown";
        if (BMI(i) < 16){
            classification(i) = "Severe Thinness";
        } else if (BMI(i) >= 16 && BMI(i) < 17){
            classification(i) = "Moderate Thinness";
        } else if (BMI(i) >= 17 && BMI(i) < 18.5){
            classification(i) = "Mild Thinness";
        } else if (BMI(i) >= 18.5 && BMI(i) < 25){
            classification(i) = "Normal";
        } else if (BMI(i) >= 25 && BMI(i) < 30){
            classification(i) = "Pre-Obese";
        } else if (BMI(i) >= 30 && BMI(i) < 35){
            classification(i) = "Obese class I";
        } else if (BMI(i) >= 35 && BMI(i) < 40){
            classification(i) = "Obese class II";
        } else if (BMI(i) >= 40){
            classification(i) = "Obese class III";
        }
    }*/
    for(int i = 0; i < BMI.size(); i++){
        classification(i) = "Unknown";
        if (BMI(i) < 18.5){
            classification(i) = "Underweight";
        } else if (BMI(i) >= 18.5 && BMI(i) < 25){
            classification(i) = "Normal";
        } else if (BMI(i) >= 25 && BMI(i) < 30){
            classification(i) = "Pre-Obese";
        } else if (BMI(i) >= 30){
            classification(i) = "Obese";
        }
    }
    return classification;
}


//Rungue Kutta 4 method for Adult
List Adult::rk4(double days){
    
    //Initial TIME(i-1)
    NumericVector k1, k2, k3, k4;
    
    //Estimate number of elements to loop into
    const int nsims = std::min(ceil(days/dt), EIchange.nrow() - 1.0);
    
    NumericMatrix AT(nind, nsims + 1); //in rcpp
    NumericMatrix ECF(nind, nsims + 1); //in rcpp
    NumericMatrix GLY(nind, nsims + 1); //in rcpp
    NumericMatrix L(nind, nsims + 1); //in rcpp
    NumericMatrix F(nind, nsims + 1); //in rcpp
    NumericMatrix BW(nind, nsims + 1); //in rcpp
    NumericMatrix BMI(nind, nsims + 1); //in rcpp
    NumericMatrix TEI(nind, nsims + 1); //in rcpp
    NumericMatrix AGE(nind, nsims + 1); //in rcpp
    StringMatrix CAT(nind, nsims + 1); //in rcpp
    
    NumericVector TIME(nsims + 1); //in rcpp
    
    //Create initial states in rcpp
    AT(_,0)  = atinit;
    ECF(_,0) = ecfinit;
    GLY(_,0) = G_base;
    L(_,0)   = lean;
    F(_,0)   = fatMass(lean);
    BW(_,0)  = bw;
    BMI(_,0) = bw/pow(ht,2.0);
    CAT(_,0) = BMIClassifier(BMI(_,0));
    TEI(_,0) = EI;
    TIME(0)  = 0.0;
    AGE(_,0) = age;
    
    
    //Loop through all other states
    bool correctVals = true;
    for (int i = 1; i <= nsims; i++){
        
        
        //Adaptive thermogenesis
        k1 = dAT(TIME(i-1), AT(_, i-1)); // f(t_n , y_n)
        k2 = dAT(TIME(i-1) + 0.5 * dt, AT(_, i-1) + 0.5 * dt * k1); // f(t_n + h/2, y_n + h/2 k1)
        k3 = dAT(TIME(i-1) + 0.5 * dt, AT(_, i-1) + 0.5 * dt * k2); // f(t_n + h/2, y_n + h/2 k2)
        k4 = dAT(TIME(i-1) + dt, AT(_, i-1) + dt * k3); // f(t_n + h, y_n + h k3)
        
        //Update AT
        AT(_,i) = AT(_, i-1) + dt * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
        
        //Extracellular fluid
        k1 = dECF(TIME(i-1), ECF(_, i-1)); // f(t_n , y_n)
        k2 = dECF(TIME(i-1) + 0.5 * dt, ECF(_, i-1) + 0.5 * dt * k1); // f(t_n + h/2, y_n + h/2 k1)
        k3 = dECF(TIME(i-1) + 0.5 * dt, ECF(_, i-1) + 0.5 * dt * k2); // f(t_n + h/2, y_n + h/2 k2)
        k4 = dECF(TIME(i-1) + dt, ECF(_, i-1) + dt * k3); // f(t_n + h, y_n + h k3)
        
        //Update ECF
        ECF(_,i) = ECF(_, i-1) + dt * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
        
        //Glycogen
        k1 = dG(TIME(i-1), GLY(_, i-1));
        k2 = dG(TIME(i-1) + 0.5 * dt, GLY(_, i-1) + 0.5 * dt * k1);
        k3 = dG(TIME(i-1) + 0.5 * dt, GLY(_, i-1) + 0.5 * dt * k2);
        k4 = dG(TIME(i-1) + dt, GLY(_, i-1) + dt * k3);
        
        //Update Glycogen
        GLY(_,i) = GLY(_, i-1) + dt * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
        
        //Lean Mass
        k1 = dL(TIME(i-1) , L(_, i-1), GLY(_, i-1), AT(_, i-1), ECF(_, i-1));
        k2 = dL(TIME(i-1) + 0.5 * dt, L(_, i-1) + 0.5 * dt * k1, 0.5*(GLY(_, i) + GLY(_, i-1)),
                0.5*(AT(_,i) + AT(_, i-1)), 0.5*(ECF(_,i) +ECF(_,i-1)));
        k3 = dL(TIME(i-1) + 0.5 * dt, L(_, i-1) + 0.5 * dt * k2, 0.5*(GLY(_,i) + GLY(_, i-1)),
                0.5*(AT(_,i) + AT(_, i-1)), 0.5*(ECF(_,i) + ECF(_,i-1)));
        k4 = dL(TIME(i-1) + dt, L(_, i-1) + dt*k3, GLY(_,i), AT(_,i), ECF(_,i));
        
        //Update L
        L(_,i) = L(_, i-1) + dt * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
        
        //Update F
        F(_,i) = fatMass(L(_,i));
        
        //Update bw
        BW(_,i) = F(_,i) + L(_,i) + ECF(_,i) + 3.7*GLY(_,i);
        
        //Update BMI
        BMI(_,i) = BW(_,i)/pow(ht,2.0);
        
        //Classify BMI
        CAT(_,i) = BMIClassifier(BMI(_,i));
        
        //Update TIME(i-1)
        TIME(i) = TIME(i-1) + dt;
        
        //Update age
        AGE(_,i) = AGE(_,i-1) + dt/365.0;
        
        //Get energy intake
        TEI(_,i) = TotalIntake(TIME(i));
        
    }
    
    return List::create(Named("Time") = TIME,
                        Named("Age") = AGE,
                        Named("Adaptive_Thermogenesis") = AT,
                        Named("Extracellular_Fluid") = ECF,
                        Named("Glycogen") = GLY,
                        Named("Fat_Mass") = F,
                        Named("Lean_Mass")   = L,
                        Named("Body_Weight") = BW,
                        Named("Body_Mass_Index") = BMI,
                        Named("BMI_Category") = CAT,
                        Named("Energy_Intake") = TEI,
                        Named("Correct_Values")=correctVals,
                        Named("Model_Type")="Adult");
    
}



//Change in calories
NumericVector Adult::deltaEI(double t){
    return EIchange(floor(t/dt),_);
}

//Change in sodiumxs
NumericVector Adult::deltaNA(double t){
    return NAchange(floor(t/dt),_);
}

