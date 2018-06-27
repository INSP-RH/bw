#' @title Dynamic Adult Weight Change Model
#'
#' @description Estimates weight change given energy and sodium intake changes at 
#' individual level.
#'
#' @param bw       (vector) Body weight for model (kg)
#' @param ht       (vector) Height for model (m)
#' @param age      (vector) Age of individual (yrs)
#' @param sex      (vector) Sex either \code{"female"} or \code{"male"}
#' @param EIchange (matrix) Matrix of caloric intake change (kcals)
#' @param NAchange (matrix) Vector of sodium intake change (mg)
#'
#' \strong{ Optional }
#' @param EI          (vector) Energy Intake at Baseline.
#' @param fat         (vector) Vector containing fat mass. Recall that 
#' @param PAL         (vector) Physical activity level.
#' @param pcarb       (vector) Percent carbohydrates after intake change.
#' @param pcarb_base  (vector) Percent carbohydrates at baseline.
#' @param days        (double) Days to run the model.
#' @param dt          (double) Time step for model; default 1 day (\code{dt = 1})
#' @param checkValues (boolean) Check whether the values from the model are biologically feasible.
#' 
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' 
#' @details \code{EIchange} and \code{NAchange} must be consumption
#' change matrices. Each row should represent consumption at each
#' day. That is, each row of \code{EIchange} and \code{NAchange} 
#' represents a day in consumption change since baseline. Consumption
#' change is non-cummulative and it's all from baseline. 
#' As an example, \code{EIchange <- rep(-100, 50)} represents that 
#' each day \code{-100} kcals are reduced from consumption. 
#' 
#' 
#' @useDynLib bw
#' @import compiler
#' @importFrom Rcpp evalCpp 
#' 
#' @references Chow, Carson C, and Kevin D Hall. 2008. \emph{The Dynamics of Human Body Weight Change.} PLoS Comput Biol 4 (3):e1000045.
#'
#' Hall, Kevin D. 2010. \emph{Predicting Metabolic Adaptation, Body Weight Change, and Energy Intake in Humans.}
#'    American Journal of Physiology-Endocrinology and Metabolism 298 (3). Am Physiological Soc: E449–E466.
#'
#' Hall, Kevin D, and Peter N Jordan. 2008. \emph{Modeling Weight-Loss Maintenance to Help Prevent Body Weight Regain.}
#'     The American Journal of Clinical Nutrition 88 (6). Am Soc Nutrition: 1495–1503.
#'     
#' Hall, Kevin D, Gary Sacks, Dhruva Chandramohan, Carson C Chow, Y Claire Wang, Steven L Gortmaker, and Boyd A Swinburn. 2011.
#' \emph{Quantification of the Effect of Energy Imbalance on Bodyweight.} The Lancet 378 (9793). Elsevier: 826–37.
#´
#´ Mifflin, Mark D, Sachiko T St Jeor, Lisa A Hill, Barbara J Scott, Sandra A Daugherty, and YO Koh. 1990.
#' \emph{A New Predictive Equation for Resting Energy Expenditure in Healthy Individuals.} 
#' The American Journal of Clinical Nutrition 51 (2). Am Soc Nutrition: 241–47.
#' 
#' @seealso \code{\link{model_plot}} for plotting the results and 
#' \code{\link{model_mean}} for aggregate data estimation. \code{\link{child_weight}} 
#' implements a similar model for children.
#' 
#' @examples 
#' #EXAMPLE 1: INDIVIDUAL MODELLING
#' #--------------------------------------------------------
#' #For one female in a diet of 100 kcal reduction. 
#' adult_weight(80, 1.8, 40, "female", rep(-100, 365))
#' 
#' #Same female also reducing sodium in -25mg
#' adult_weight(80, 1.8, 40, "female", rep(-100, 365), rep(-25, 365))
#' 
#' #Same female modelled for 400 days
#' adult_weight(80, 1.8, 40, "female", rep(-100, 400), rep(-25, 400), days = 400)
#' 
#' #Same female reducing -50 kcals per 100 days and not reducing sodium
#' kcalvec <-c(rep(-50, 100), rep(-100, 100), rep(-150, 100), rep(-200, 100))
#' adult_weight(80, 1.8, 40, "female", kcalvec, days = 400)
#' 
#' #Same female with known energy intake
#' adult_weight(80, 1.8, 40, "female", rep(-100, 365), rep(-25, 365), EI = 2000)
#' 
#' #Same female with known fat mass
#' adult_weight(80, 1.8, 40, "female", rep(-100, 365), rep(-25, 365), fat = 32)
#' 
#' #Same female with known fat mass and known energy consumption
#' adult_weight(80, 1.8, 40, "female", rep(-100, 365), rep(-25, 365), EI = 2000, fat = 32)
#' 
#' #EXAMPLE 2: DATASET MODELLING
#' #--------------------------------------------------------
#' 
#' #Antropometric data
#' weights <- c(45, 67, 58, 92, 81)
#' heights <- c(1.30, 1.73, 1.77, 1.92, 1.73)
#' ages    <- c(45, 23, 66, 44, 23)
#' sexes   <- c("male", "female", "female", "male", "male") 
#' 
#' #Matrix of energy consumption reduction: 
#' EIchange <- rbind(rep(-100, 365), rep(-200, 365), rep(-200, 365), 
#'                   rep(-123, 365), rep(-50, 365))
#' 
#' #Returns a weight change matrix and other matrices
#' model_weight <- adult_weight(weights, heights, ages, sexes, 
#'                              EIchange)["Body_Weight"][[1]]
#' 
#' @export


adult_weight <- function(bw, ht, age, sex, 
                         EIchange = matrix(0, ncol = abs(ceiling(days/dt)), nrow = length(bw)), 
                         NAchange = matrix(0, ncol = abs(ceiling(days/dt)), nrow = length(bw)), 
                         EI = NA, fat = rep(NA, length(bw)),
                         PAL = rep(1.5, length(bw)), 
                         pcarb_base = rep(0.5, length(bw)), 
                         pcarb = pcarb_base,  days = 365, dt = 1,
                         checkValues = TRUE){
  
  #Check that EIchange and Nachange are matrices
  if (is.vector(EIchange)){
    EIchange <- matrix(EIchange, nrow = 1)
  }  
  if (is.vector(NAchange)){
    NAchange <- matrix(NAchange, nrow = 1)
  }
  
  if (any(dim(EIchange) != dim(NAchange))){
    stop("Dimension mismatch. NAchange and EIchange don't have the same dimensions.")
  }
  
  #Check that all parameters have same length
  if (length(bw) != length(ht)  || length(bw) != length(age) || 
      length(bw) != length(sex) || length(bw) != length(PAL) || 
      length(bw) != length(pcarb_base) || length(bw) != length(pcarb) ||
      length(bw) != length(fat)){
    stop(paste0("Dimension mismatch. bw, ht, age, sex, PAL, fat, pcarb_base", 
                "and pcarb don't have the same length"))
  }
  
  
  #Check that dt is > 0
  if (dt < 0 || dt > days){
    stop(paste0("Invalid time step dt; please choose 0 < dt < days"))
  }
  
  #Check that EIchange has the same number of rows as the length of bw
  if (nrow(EIchange) != length(bw)){
    stop(paste("Dimension mismatch. EIchange must have the", 
               "same amount of rows as individuals."))
  }
  
  #Check that they have as many columns as days
  if (ncol(EIchange) != ceiling(days/dt)){
    warning(paste("Dimension mismatch. EIchange and NAchange must have", 
                  ceiling(days/dt), "columns"))
  }

  
  #Check that age, bw and height are positive
  if (any(bw <= 0) || any(ht <= 0) || any(age < 0)){
    stop(paste0("Don't know how to handle negative or zero values ",
                "in bw and ht. Nor  negative values in age."))
  }
  
  #Check that age is greater than 18
  if(any(age<18)){
    warning(paste0("Some individuals' age is less than 18. Results",
                   " are not be accurate for children and adolescents.",
                   "Please use child_weight instead."))
  }
  
  # Check pcarb and pcarb_base are between 0 and 1
  if(any(pcarb_base > 1) || any(pcarb_base<0) || any(pcarb > 1) || any(pcarb<0)){
    stop(paste0("The variables pcarb and pcarb_base are ",
                "the proportion of carbohydrates consumed.",
                "Therefore they must take values between 0 and 1."))
  }
  # Check PAL values
  if(any(PAL <=0)){
    stop("PAL must have a positive value")
  }else if(any(PAL<1.4)){
    warning(paste("Some individuals have a PAL less than 1.4, which is only",
                  "viable in extreme cases such as: elderly mental patients,",
                  "adolescents with cerebral palsy or myelodysplasia",
                  "and resting adults confined to a whole body calorimeter." ,
                  "(WHO, ENERGY REQUIREMENTS OF ADULTS)"))
  }else if(any(PAL > 2.4)){
    warning(paste("Some individuals have a PAL greater than 2.4, which rarely occurs",
                 "and is not sustainable in the long term.",
                 "(WHO, ENERGY REQUIREMENTS OF ADULTS)"))
  }
  
  #Check days > 0
  if (days <= 0){
    stop("Don't know how to handle negative time scales.Please make sure days > 0.")
  }
  
  #Check sex is "male" and "female"
  if (length(which(!(sex %in% c("male","female")))) > 0){
    stop(paste0("Invalid sex. Please specify either 'male' of 'female'"))
  }
  
  
  #Change sex to numeric for c++
  newsex                         <- rep(0, length(sex))
  newsex[which(sex == "female")] <- 1
  
  #Check fat/energy are inputted
  isfat <- any(is.na(fat))
  isEI  <- any(is.na(EI))
  
  #Change because c++ takes them as transpose
  EIchange <- t(EIchange)
  NAchange <- t(NAchange)
  
  #Run C++ program to estimate weight there are 3 constructors depending
  #on if you have energy intake or fat intake or not.
  if (isfat && isEI){
    wl <- adult_weight_wrapper(bw, ht, age, newsex, EIchange, NAchange,
                               PAL, pcarb_base, pcarb, dt, ceiling(days), checkValues)  
  } else if (!isEI && isfat) {
    wl <- adult_weight_wrapper_EI(bw, ht, age, newsex, EIchange, NAchange,
                                  PAL, pcarb_base, pcarb, dt, EI, ceiling(days), checkValues, TRUE)  
  } else if (isEI && !isfat) {
    wl <- adult_weight_wrapper_EI(bw, ht, age, newsex, EIchange, NAchange,
                                  PAL, pcarb_base, pcarb, dt, fat, ceiling(days), checkValues, FALSE)  
  } else if (!isEI && !isfat){
    wl <- adult_weight_wrapper_EI_fat(bw, ht, age, newsex, EIchange, NAchange,
                                      PAL, pcarb_base, pcarb, dt, EI, fat, ceiling(days), checkValues)  
  }
  if(wl$Correct_Values[1]==FALSE){
    stop("One of the variables takes either negative values, or NaN, NA or infinity")
  }
  return(wl)
  
  
}