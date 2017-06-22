//
//  adult_weight.hpp
//  
//
//  Created by Rodrigo Zepeda Tello on 06/06/17.
//
//

#ifndef adult_weight_h
#define adult_weight_h

#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

//Create a Adult class to contain individual parameters
//--------------------------------------------------------------------------------
class Adult {
public:
    
    //Constructor for when initial energy intake is estimated by the model
    Adult(NumericVector weight, NumericVector height, NumericVector age_yrs,
          NumericVector sexstring, NumericMatrix input_EIchange,
          NumericMatrix input_NAchange, NumericVector physicalactivity,
          NumericVector percentc, NumericVector percentb, bool checkValues);
    
    //Constructor for when initial energy or initial fat intake is added by user
    Adult(NumericVector weight, NumericVector height, NumericVector age_yrs,
          NumericVector sexstring, NumericMatrix input_EIchange,
          NumericMatrix input_NAchange, NumericVector physicalactivity,
          NumericVector percentc, NumericVector percentb, NumericVector extradata, bool checkValues, bool isEnergy);
    
    //Constructor for when initial energy intake and initial fat is added by user
    Adult(NumericVector weight, NumericVector height, NumericVector age_yrs,
          NumericVector sexstring, NumericMatrix input_EIchange,
          NumericMatrix input_NAchange, NumericVector physicalactivity,
          NumericVector percentc, NumericVector percentb, NumericVector input_EI,
          NumericVector input_fat, bool checkValues);
    
    //Destroyer
    ~ Adult();
    
    //Constants depending on the Adult
    //---------------------------------------------------------------------------
    NumericVector bw;              //Weight (kg)
    NumericVector ht;              //Height (m)
    NumericVector age;             //Age (yrs)
    NumericVector sex;             //0 = "male"; 1 = "female"
    NumericVector EI;              //Energy intake (kcal)
    NumericVector PAL;             //Physical Activity Level PAL
    NumericVector fat;             //Fat mass at baseline (kg)
    NumericVector lean;            //Lean mass at baseline (kg)
    NumericVector steadystate;     //Steady state of energy intake for no weight change according to Miffin % St Jeor (kcal)
    NumericVector G_base;          //Glycogen at baseline (kg)
    NumericVector ecfinit;         //Initial extracellular fluid (kg)
    NumericVector CIb;             //Carbohydrate intake at baseline (kcal)
    NumericVector pcarb;           //% carbohydrates after change
    NumericVector pcarb_base;      //% carbohydrates at baseline
    
    //Numeric vectors containing EI and NA changes
    NumericMatrix EIchange;
    NumericMatrix NAchange;
    

    
    //Functions
    //---------------------------------------------------------------------------
    List rk4(double days); //in Rcpp:
    
private:
    
    //Constants depending on the Adult
    //---------------------------------------------------------------------------
    NumericVector kG;              //Constant
    NumericVector K;               //Energy balance constant at baseline
    NumericVector rmr;             //Resting Metabolic Rate (kcal)
    NumericVector delta;           //Delta parameter of activity
    NumericVector atinit;          //Initial Adaptive Thermogenesis
    
    //Pre-defined parameters applicable to the whole population
    //---------------------------------------------------------------------------
    double roG;     //1000*17.6*0.23900573614 #Changed from kjoules to kcals
    double Na ;     // (1000*3.22)#Sodium
    double zetaNa;
    double zetaCI;
    double roF;     //1000*39.5*0.23900573614 #Changed from kjoules to kcals
    double roL;     //1000*7.6*0.23900573614  #Changed from kjoules to kcals
    double gammaF;  //13*0.23900573614        #Changed from kjoules to kcals
    double gammaL;  // 92*0.23900573614       #Changed from kjoules to kcals
    double etaF;    // 750*0.23900573614      #Changed from kjoules to kcals
    double etaL;    // 960*0.23900573614      #Changed from kjoules to kcals
    double betaTEF;
    double betaAT;
    double tauAT;
    double C;       //10.4*(roL/roF)
    double alfa1;   //Auxiliary functions from Pablo
    double alfa2;   //Auxiliary functions from Pablo
    double rmrbw;
    double rmrage;
    double rmrht;
    double rmr_m;
    double rmr_f;
    int    nind; //Number of individuals in model
    double dt;   //Delta t for Rungue Kutta 4
    bool check;
    
    //Auxiliary functions
    void getRMR(void);
    void getParameters(void);
    void getBaselineMass(void);
    void getCaloricSteadyState(void);
    void getEnergy(void);
    void getDelta(void);
    void getK(void);
    void getCarbConstants(void);
    void getATinit(void);
    void getECFinit(void);
    void build(NumericVector weight, NumericVector height, NumericVector age_yrs,
               NumericVector sexstring, NumericMatrix input_EIchange,
               NumericMatrix input_NAchange, NumericVector physicalactivity,
               NumericVector percentc, NumericVector percentb, bool checkValues);
    void build(NumericVector weight, NumericVector height, NumericVector age_yrs,
               NumericVector sexstring, NumericMatrix input_EIchange,
               NumericMatrix input_NAchange, NumericVector physicalactivity,
               NumericVector percentc, NumericVector percentb, NumericVector extradata,
               bool checkValues, bool isEnergy);
    void build(NumericVector weight, NumericVector height, NumericVector age_yrs,
               NumericVector sexstring, NumericMatrix input_EIchange,
               NumericMatrix input_NAchange, NumericVector physicalactivity,
               NumericVector percentc, NumericVector percentb, NumericVector input_EI,
               NumericVector input_fat,bool checkValues);
    NumericVector TotalIntake (double t);
    StringVector  BMIClassifier(NumericVector BMI);
    NumericVector CI(double t);
    NumericVector R(double t, NumericVector L, NumericVector G,
                    NumericVector AT, NumericVector ECF);
    NumericVector fatMass(NumericVector L);
    NumericVector deltaEI(double t);
    NumericVector deltaNA(double t);
    NumericVector TEF(double t);
    NumericVector dAT(double t, NumericVector AT);
    NumericVector dECF(double t, NumericVector ECF);
    NumericVector dG(double t, NumericVector G);
    NumericVector dL(double t, NumericVector L, NumericVector G,
                     NumericVector AT, NumericVector ECF);
    
    
};

#endif /* adult_weight_h */
