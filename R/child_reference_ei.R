#' @title Energy Intake Matrix
#'
#' @description Estimates weight given age, sex, fat mass, and fat free mass, 
#'
#' @param age      (vector) Age of individual (yrs)
#' @param sex      (vector) Sex either \code{"female"} or \code{"male"}
#' @param FM       (vector) Fat Mass at Baseline
#' @param FFM      (vector) Fat Free Mass at Baseline
#' 
#' \strong{ Optional }
#' @param days     (numeric) Days to run the model.
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @useDynLib bw
#' @import compiler
#' @importFrom Rcpp evalCpp 
#' 
#' @references Hall, K. D., Butte, N. F., Swinburn, B. A., & Chow, C. C. (2013). 
#' \emph{Dynamics of childhood growth and obesity: development and validation of a 
#' quantitative mathematical model}. The Lancet Diabetes & Endocrinology, 1(2), 97-105.
#' 
#' @seealso \code{\link{model_plot}} for plotting the results and 
#' \code{\link{model_mean}} for aggregate data estimation. 
#' 
#' @examples
#' #One child
#' child_reference_EI(6, "male", 2, 4, 10)
#' 
#' #Several children
#' child_reference_EI(sample(6:12, 10, replace = TRUE), 
#'                    sample(c("male","female"), 10, replace = TRUE), 
#'                    sample(2:10, 10, replace = TRUE), 
#'                    sample(2:10, 10, replace = TRUE),
#'                    365)
#' 
#' @keywords internal
#' @export

child_reference_EI <- function(age, sex, FM, FFM, days){
  
  #Change sex to numeric for c++
  newsex                         <- rep(0, length(sex))
  newsex[which(sex == "female")] <- 1
  
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
  
  #Get c++ function
  t(intake_reference_wrapper(age, newsex, FM, FFM, days))
}