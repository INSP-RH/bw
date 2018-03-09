//
//  Child.hpp
//  AdultHall
//
//  Created by Rodrigo Zepeda Tello on 26/05/17.
//  Copyright © 2017 Instituto Nacional de Salud Publica. All rights reserved.
//

#ifndef child_weight_h
#define child_weight_h

#include <math.h>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

//Create a Adult class to contain individual parameters
//--------------------------------------------------------------------------------
class Child {
public:
    
    //Constructor and destroyer
    Child(NumericVector input_age, NumericVector input_sex, NumericVector input_FFM, NumericVector input_FM, NumericMatrix input_EIntake, bool checkValues);
    ~Child(void);
    
    //Constants
    NumericVector age;  //Age (yrs)
    NumericVector sex;  //0 = "male"; 1 = "female"
    NumericVector FFM;  //Fat Free Mass (kg)
    NumericVector FM;   //Fat Mass (kg)
    NumericMatrix EIntake;
    bool          check; // Check values are correct
    
    //Functions
    //---------------------------------------------------------------------------
    List rk4(double days);
    
    //Reference functions for reference children
    NumericVector IntakeReference(NumericVector t);
    NumericVector FFMReference(NumericVector t);
    NumericVector FMReference(NumericVector t);
    
private:
    
    //Private unchanging constants
    double rhoFM; //kcals/g -> kcals/kg
    double deltamin;
    double P;
    double h;
    
    //Number of individuals
    int nind;
    
    //Constants additional
    NumericVector K;
    NumericVector deltamax;
    
    //Constants for g FROM DYNAMICS PAPER
    NumericVector A;
    NumericVector tA;
    NumericVector tauA;
    NumericVector B;
    NumericVector tB;
    NumericVector tauB;
    NumericVector D;
    NumericVector tD;
    NumericVector tauD;
    
    //Constants for g FROM IMPACT PAPER
    NumericVector A1;
    NumericVector tA1;
    NumericVector tauA1;
    NumericVector B1;
    NumericVector tB1;
    NumericVector tauB1;
    NumericVector D1;
    NumericVector tD1;
    NumericVector tauD1;
    
    //Constants for EB FROM IMPACT PAPER
    NumericVector A_EB;
    NumericVector tA_EB;
    NumericVector tauA_EB;
    NumericVector B_EB;
    NumericVector tB_EB;
    NumericVector tauB_EB;
    NumericVector D_EB;
    NumericVector tD_EB;
    NumericVector tauD_EB;
    
    //Constants for linear coefficients of ffm and fm regressions
    NumericVector ffm_beta0;
    NumericVector ffm_beta1;
    NumericVector fm_beta0;
    NumericVector fm_beta1;
    
    //WHO energy requirements
    NumericVector req;

    //
    void build(void);
    void getParameters();
    NumericVector Growth_dynamic(NumericVector t); //Growth function from Dynamics...
    NumericVector Growth_impact(NumericVector t);   //Growth function from Impact...
    NumericVector EB_impact(NumericVector t);   //Energy Balance function from Impact...
    NumericVector cRhoFFM(NumericVector input_FFM); //Crho function
    NumericVector general_ode(NumericVector t, NumericVector input_A, NumericVector input_B,
                              NumericVector input_D, NumericVector input_tA,
                              NumericVector input_tB, NumericVector input_tD,
                              NumericVector input_tauA, NumericVector input_tauB,
                              NumericVector input_tauD);
    NumericVector cP(NumericVector FFM, NumericVector FM);
    NumericVector Delta(NumericVector t);
    NumericVector Expenditure(NumericVector t, NumericVector FFM, NumericVector FM);
    NumericVector Intake(NumericVector t);
    NumericMatrix dMass (NumericVector time, NumericVector FFM, NumericVector FM);
};


#endif /* Child_h */
