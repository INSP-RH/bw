#' @title Dynamic Children Weight Change Model
#'
#' @description Estimates weight given age, sex, fat mass, and fat free mass, 
#'
#' @param age      (vector) Age of individual (yrs)
#' @param sex      (vector) Sex either \code{"female"} or \code{"male"}
#' @param FM       (vector) Fat Mass at Baseline
#' @param FFM      (vector) Fat Free Mass at Baseline
#' @param EI       (matrix) Numeric Matrix with energy intake
#' @param richardsonparams (list) List of parameters for Richardson's curve for energy. See details.
#' 
#' \strong{ Optional }
#' @param days     (numeric) Days to run the model.
#' @param checkValues (boolean) Checks whether values of fat mass and free fat mass are possible
#' @param dt       (double) Time step for Rungue-Kutta method
#' 
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @details \code{richardsonparams} is a named list of parameters:
#' \code{K}, \code{A}, \code{Q}, \code{C}, \code{B}, \code{nu}
#' which result in Richardon's curve:
#' \deqn{A + \frac{K-A}{(C + Q exp(-B*t))^{1/nu}}}
#' The Richardson's curve is another option for modelling the energy
#' intake for a child: by specifying the parameters no energy input
#' is needed; instead Energy is assumed to follow the equation:
#' \deqn{EI(t) = A + \frac{K-A}{(C + Q exp(-B*t))^{1/nu}}}
#' 
#' @useDynLib bw
#' @import compiler
#' @importFrom Rcpp evalCpp 
#' 
#' @references Hall, K. D., Butte, N. F., Swinburn, B. A., & Chow, C. C. (2013). 
#' \emph{Dynamics of childhood growth and obesity: development and validation of a 
#' quantitative mathematical model}. The Lancet Diabetes & Endocrinology, 1(2), 97-105.
#' 
#' Haschke, F. (1989). \emph{Body Composition During Adolescence.} 
#' Body Composition Measurements in Infants and Children.
#' Ross Laboratories Columbus, OH, 76–83.
#' 
#' Fomon, Samuel J, Ferdinand Haschke, Ekhard E Ziegler, and Steven E Nelson. 1982. 
#' \emph{Body Composition of Reference Children from Birth to Age 10 Years.}
#' The American Journal of Clinical Nutrition 35 (5). Am Soc Nutrition: 1169–75.
#'  
#' Ellis, Kenneth J, Roman J Shypailo, Steven A Abrams, and William W Wong. 2000. 
#' \emph{The Reference Child and Adolescent Models of Body Composition: A Contemporary Comparison.} 
#' Annals of the New York Academy of Sciences 904 (1). Wiley Online Library: 374–82.
#'  
#' Deurenberg, Paul, Jan A Weststrate, and Jaap C Seidell. 1991. 
#' \emph{Body Mass Index as a Measure of Body Fatness: Age-and Sex-Specific Prediction Formulas.} 
#' British Journal of Nutrition 65 (2). Cambridge University Press: 105–14.
#' 
#' Katan, Martijn B, Janne C De Ruyter, Lothar DJ Kuijper, Carson C Chow, Kevin D Hall, and Margreet R Olthof. 2016.
#' \emph{Impact of Masked Replacement of Sugar-Sweetened with Sugar-Free Beverages on Body Weight Increases with Initial Bmi:
#' Secondary Analysis of Data from an 18 Month Double–Blind Trial in Children.} 
#' PloS One 11 (7). Public Library of Science: e0159771.
#' 
#' @seealso @\code{\link{adult_weight}} for the weight change model for adults;
#' \code{\link{model_plot}} for plotting the results and 
#' \code{\link{model_mean}} for aggregate data estimation. 
#' 
#' @examples 
#' #EXAMPLE 1: INDIVIDUAL MODELLING
#' #--------------------------------------------------------
#' #For one child with default energy intake
#' child_weight(6,"male")
#' 
#' #For a child with specific energy intake
#' child_weight(6,"male",2.5, 16, as.matrix(rep(2000, 365)), days = 365)
#' 
#' #Using Richardson's energy
#' girl <- child_weight(6,"female", days=365, dt = 5, 
#'                      richardsonparams = list(K = 2700, Q = 10, 
#'                      B = 12, A = 3, nu = 4, C = 1))
#' plot(girl$Body_Weight[1,])
#' 
#' #EXAMPLE 2: DATASET MODELLING
#' #--------------------------------------------------------
#' #Antropometric data
#' FatFree <- c(32, 17.2, 18.8, 20, 24.1)
#' Fat     <- c(4.30, 2.02, 3.07, 1.12, 2.93)
#' ages    <- c(10, 6.2, 5.4, 4, 4.1)
#' sexes   <- c("male", "female", "female", "male", "male") 
#' 
#' #With specific energy intake
#' eintake <- matrix(rep(2000, 365*5), ncol = 5)
#' 
#' #Returns a weight change matrix and other matrices
#' model_weight <- child_weight(ages, sexes, Fat, FatFree, eintake)
#' 
#' model_weight_2 <- child_weight(ages, sexes, Fat, FatFree, 
#'                     richardsonparams = list(K = 2700, Q = 10, 
#'                     B = 12, A = 3, nu = 4, C = 1))
#'          
#' @export
#'

child_weight <- function(age, sex, FM = child_reference_FFMandFM(age, sex)$FM, 
                         FFM = child_reference_FFMandFM(age, sex)$FFM, 
                         EI = NA, 
                         richardsonparams = list(K = NA, Q = NA, B = NA, A = NA, nu = NA, C = NA),
                         days = 365, dt = 1, checkValues = TRUE){
  
  #Check all variables are positive
  if (any(age < 0) || any(FM < 0) || any(FFM < 0)){
    stop("Cannot handle negative values for age, FM and FFM.")
  }
  
  #Check days > 0
  if (days <= 0){
    stop("Don't know how to handle negative time scales.Please make sure days > 0.")
  }
  
  #Check dimensions of inputs
  if (length(age) != length(sex) || length(age) != length(FM) 
      || length(age) != length(FFM)){
    stop("Dimension mismatch: age, sex, FM and FFM must have same length.")
  }
  
  #Check sex is "male" and "female"
  if (length(which(!(sex %in% c("male","female")))) > 0){
    stop(paste0("Invalid sex. Please specify either 'male' of 'female'"))
  }
  
  #Check that we don't go over 18 yrs where we have no data
  if (max(age) + days/365 > 18){
    warning(paste0("Some individuals reach age > 18 before model stops. Results",
                   " might not be accurate for adults. Please use adult_weight",
                   " instead."))
  }
  
  #Check that dt is > 0
  if (dt < 0 || dt > days){
    stop(paste0("Invalid time step dt; please choose 0 < dt < days"))
  }
  
  #Check if is na logistic and params
  if (is.na(EI[1]) & (is.na(richardsonparams$K) || is.na(richardsonparams$Q) || 
                   is.na(richardsonparams$A) || is.na(richardsonparams$B) || 
                   is.na(richardsonparams$nu) || is.na(richardsonparams$C))){
    message("Creating default energy intake for healthy child.")
    EI <- child_reference_EI(age, sex, FM, FFM, days, dt) 
  }
  
  #Change sex to numeric for c++
  newsex                         <- rep(0, length(sex))
  newsex[which(sex == "female")] <- 1
  
  #Choose between richardson curve or given energy intake
  if (!is.na(EI[1])){
    message("Using user's energy intake")
    wt <- child_weight_wrapper(age, newsex, FFM, FM, as.matrix(EI), days, dt, checkValues)  
  } else {
    message("Using Richardson's function")
    wt <- child_weight_wrapper_richardson(age, newsex, FFM, FM, richardsonparams$K, 
                               richardsonparams$Q, richardsonparams$A, 
                               richardsonparams$B, richardsonparams$nu, 
                               richardsonparams$C, days, dt, checkValues)
  }
  
  
  return(wt)
  
  
}
