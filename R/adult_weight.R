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
#' @param PAL         (vector) Physical activity level.
#' @param pcarb       (vector) Percent carbohydrates after intake change.
#' @param pcarb_base  (vector) Percent carbohydrates at baseline.
#' @param days        (double) Days to run the model.
#' @param checkValues (boolean) Check whether the values from the model are biologically feasible.
#'
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @details \code{EIchange} and \code{NAchange} must be consumption
#' change matrices. Each row should represent consumption at each
#' day. That is, each row of \code{EIchange} and \code{NAchange} 
#' represents a day in consumption change since baseline. Consumption
#' change is non-cummulative and it's all from baseline. 
#' As an example, \code{EIchange <- rep(-100, 50)} represents that 
#' each day \code{-100} kcals are reduced from consumption. 
#' 
#' @note I still haven't solved how to be able to change de
#' Rungue-Kutta \code{dt} parameter. For now it is a day; but for faster
#' models it should be changeable. 
#' 
#' @note For faster modelling code should be computed in paralell.
#' \url{https://rcppcore.github.io/RcppParallel/}
#' 
#' @useDynLib bw
#' @import compiler
#' @importFrom Rcpp evalCpp 
#' 
#' @references HALL, Kevin D., et al. \emph{Quantification of the effect of energy 
#' imbalance on bodyweight}. The Lancet, 2011, vol. 378, no 9793, p. 826-83
#' 
#' @seealso \code{\link{model_plot}} for plotting the results and 
#' \code{\link{model_mean}} for aggregate data estimation. 
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
#' EIchange <- cbind(rep(-100, 365), rep(-200, 365), rep(-200, 365), 
#'                   rep(-123, 365), rep(-50, 365))
#' 
#' #Returns a weight change matrix and other matrices
#' model_weight <- adult_weight(weights, heights, ages, sexes, 
#'                              EIchange)["Body_Weight"][[1]]
#' 
#' @export


adult_weight <- function(bw, ht, age, sex, 
                         EIchange = matrix(0, ncol = length(bw), nrow = ceiling(days)), 
                         NAchange = matrix(0, ncol = length(bw), nrow = ceiling(days)), 
                         EI = NA,
                         PAL = rep(1.5, length(bw)), 
                         pcarb_base = rep(0.5, length(bw)), 
                         pcarb = pcarb_base,  days = 365,
                         checkValues = TRUE){
  
  #Check that dimension of EIchange and Nachange match
  EIchange <- as.matrix(EIchange)
  NAchange <- as.matrix(NAchange)
  
  if (any(dim(EIchange) != dim(NAchange))){
    stop("Dimension mismatch. NAchange and EIchange don't have the same dimensions.")
  }
  
  #Check that all parameters have same length
  if (length(bw) != length(ht)  || length(bw) != length(age) || 
      length(bw) != length(sex) || length(bw) != length(PAL) || 
      length(bw) != length(pcarb_base) || length(bw) != length(pcarb)){
    stop(paste0("Dimension mismatch. bw, ht, age, sex, PAL, pcarb_base", 
                "and pcarb don't have the same length"))
  }
  
  #Check that they have as many rows as days
  if (nrow(EIchange) != ceiling(days)){
    warning(paste("Dimension mismatch. EIchange and NAchange must have", 
                  ceiling(days), "rows"))
  }
  
  #Check that age, bw and height are positive
  if (any(bw < 0) || any(ht < 0) || any(age < 0)){
    stop("Don't know how to handle negative values in bw, ht or age")
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
  
  #Run C program to estimate the loop
  if (any(is.na(EI))){
    wl <- adult_weight_wrapper(bw, ht, age, newsex, EIchange, NAchange,
                               PAL, pcarb_base, pcarb, ceiling(days), checkValues)  
  } else {
    wl <- adult_weight_wrapper_EI(bw, ht, age, newsex, EIchange, NAchange,
                                PAL, pcarb_base, pcarb, EI, ceiling(days), checkValues)  
  }
  if(wl$Correct_Values[1]==FALSE){
    stop("One of the variables takes either negative values, or NaN, NA or infinity")
  }
  return(wl)
  
  
}