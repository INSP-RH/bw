#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix EnergyBuilder(NumericMatrix Energy, NumericVector Time, 
                            std::string interpol){
  
  //INPUT:
  //Energy - includes measurements of energy comsumption where each row is an individual
  //and each column is a time in Time in which it was measured
  //Time   - Vector of measurement times of the energy. First element of Time should be a 0
  //otherwise the model does not make any sense.
  //interpol .- Interpolation mode: linear, exponential, stepwise_r, stepwise_l, 
  //brownian and logarihmmic
  
  //Number of times to calculate
  int days = floor(Time(Time.size()-1));
  int j    = 0; //indicator of time value we are taking
  
  //Numeric matrix to return
  NumericMatrix Evalues(Energy.nrow(), days + 1);
  
  double K = 5000; //To avoid logarithm starting at 0 we displace the exponential to let for a maximum y2 - y1 of 1000.
  
  //Brownian bridge
  if (interpol.compare("Brownian") == 0){
   
   for (int j = 0; j < (Time.size()-1); j++){
     
     //Get times
     double T = Time(j+1);
     double t = Time(j);
     
     //Simulate W brownian path
     NumericMatrix W(Energy.nrow(), (T - t) + 1); //By default W(_, 0) = 0;
     for (int i = 1; i < (T - t + 1); i++){
       W(_, i) = W(_,i-1) + rnorm(Energy.nrow());
     }
     
     //Get brownian bridge
     for (int i = 0 ; i < (T - t + 1); i++){
       Evalues(_,i + t) = Energy(_,j)*( (T - t) - i )/(T - t) + Energy(_,j+1)*i/(T-t) + 
         W(_, i) -  (i/(T-t))*W(_, (T-t));  
     }
   }
   
  } else {
    
    //Case linear
    for (int i = 0; i < days; i++){
      
      //Linear case
      if (interpol.compare("Linear") == 0){
        Evalues(_,i) =  (Energy(_, j + 1) - Energy(_, j))/(Time(j+1) - Time(j))*(i - Time(j)) + 
          Energy(_, j);
      } else if (interpol.compare("Stepwise_L") == 0){
        Evalues(_,i) = Energy(_,j);
      } else if (interpol.compare("Stepwise_R") == 0){
        Evalues(_,i) = Energy(_,j+1);
      } else if (interpol.compare("Exponential") == 0){
        Evalues(_,i) =  exp((log(Energy(_, j + 1) - Energy(_, j) + K) - log(K))/(Time(j+1) - Time(j))*(i - Time(j)) + log(K)) - K + Energy(_, j);
      } else if (interpol.compare("Logarithmic") == 0){
        Evalues(_,i) =  1000*log( (exp( (Energy(_, j + 1) - Energy(_, j))/1000) -1)/(Time(j+1) - Time(j))*(i - Time(j)) + 1) + Energy(_,j);
        //Evalues(_,i) =  log( (exp(Energy(_,j+1) - Energy(_,j)) - 1)/(Time(j+1) - Time(j))*(i - Time(j)) + 1) + Energy(_,j);
      } 
      
      //Update to next time
      if (i + 1 >= Time(j + 1)){
          j = j + 1;  
      }
      
    }
    
    //Last day
    Evalues(_, Evalues.ncol() - 1) = Energy(_,Energy.ncol() - 1);
    
  }
  
  
  return Evalues;
}
