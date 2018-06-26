#' @title Energy Matrix Interpolating Function
#'
#' @description Creates a matrix interpolating energy consumption
#' from measurements at specific moments in time. 
#'
#' @param energy   (matrix) Matrix with each row representing an individual and each column
#' a moment in time in which energy was measured. Energy is assumed to be measured at time 0 
#' initially.
#' 
#' @param time     (vector) Vector of times at which the measurements (columns of energy) 
#' were made. \strong{Note} that first element of time most always be \code{0}. 
#' 
#' \strong{ Optional }
#' @param interpolation (string) Way to interpolate the values between measurements. Currently
#' supporting \code{"Linear"}, \code{"Exponential"}, \code{"Stepwise_R"},  \code{"Stepwise_L"},
#' \code{"Logarithmic"} and \code{"Brownian"}.
#' 
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' 
#' @useDynLib bw
#' @import compiler
#' @importFrom Rcpp evalCpp 
#' 
#' @seealso \code{\link{adult_weight}} for weight change in adults and
#' \code{\link{child_weight}} for children weight change. 
#' 
#' @examples 
#' #EXAMPLE 1: INDIVIDUAL MODELLING
#' #--------------------------------------------------------
#' 
#' #Get energy consumption
#' myconsumption <- energy_build(c(0, 200, -500), c(0, 365*2, 365*4), "Linear")
#' plot(1:(365*4), myconsumption, type = "l")
#' 
#' #Change interpolation to exponential
#' myexponential <- energy_build(c(0, 200, -500), c(0, 365*2, 365*4), "Exponential")
#' lines(1:(365*4), myexponential, type = "l", col = "red")
#' 
#' mystepwise    <- energy_build(c(0, 200, -500), c(0, 365*2, 365*4), "Stepwise_R")
#' lines(1:(365*4), mystepwise, type = "l", col = "blue")
#' 
#' mystepwise2    <- energy_build(c(0, 200, -500), c(0, 365*2, 365*4), "Stepwise_L")
#' lines(1:(365*4), mystepwise2, type = "l", col = "green")
#' 
#' mylogarithmic <- energy_build(c(0, 200, -500), c(0, 365*2, 365*4), "Logarithmic")
#' lines(1:(365*4), mylogarithmic, type = "l", col = "purple")
#' 
#' mybrownian    <- energy_build(c(0, 200, -500), c(0, 365*2, 365*4), "Brownian")
#' lines(1:(365*4), mybrownian, type = "l", col = "forestgreen")
#' 
#' #EXAMPLE 2: GROUP MODELLING
#' #--------------------------------------------------------
#' 
#' #Get energy consumption
#' multiple <- energy_build(cbind(runif(10,1000,2000), 
#'                                  runif(10,1000,2000), 
#'                                  runif(10,1000,2000)), c(0, 142, 365),
#'                                  "Brownian")
#' matplot(1:365, t(multiple), type = "l")
#' @export
#'

energy_build <- function(energy, time, interpolation = "Brownian"){
  
  #Set energy as matrix
  if (is.vector(energy)){
    energy <- matrix(energy, nrow = 1)
  }
  
  #Check that time is a vector
  if(is.vector(time)==FALSE){
    stop("Variable time should be a vector. Time values are the same for all individuals")
  }
  
  #Check that energy has same columns as time length
  if (ncol(energy) != length(time)){
    stop(paste0("energy matrix has different number of columns than length(time).",
                "Recall that each column of energy matrix represents a measurement",
                "in time."))
  }
  
  # Check that time values are non negative
  if(any(time<0)){
    stop("Values in time should be positive.")
  }
  
  # Check that time values are integers
  if(any((round(time)!=time))){
    stop("Values in time should be integer.")
  }
  
  #Check that first time element is 0
  if (time[1] != 0){
    stop("First element of time, time[1], must always equal 0")
  }
  
  # Check that time is increasing
  
  if(any(time!=time[order(time)]) || length(time)!= length(unique(time))){
    stop("Time values should be increasing, and there should be no repeated values.")
  }
  
  #Check that interpolation in list
  if (!(interpolation %in% c("Linear","Exponential","Logarithmic",
                             "Stepwise_L","Stepwise_R","Brownian"))){
    stop(paste0("Invalid interpolation. Please choose one of the following:",
                "\n - 'Linear' \n - 'Exponential' \n - 'Logarithmic' \n - 'Stepwise_L'",
                "\n - 'Stepwise_R' \n - 'Brownian'"))
  }
  
  #Run energy builder
  return( EnergyBuilder(energy, time, interpolation)[,-1] )
  
}