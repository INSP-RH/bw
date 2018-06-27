#' @title FFM and FM reference
#'
#' @description Estimates FFM and FM reference given age, sex.
#'
#' @param age      (vector) Age of individual (yrs)
#' @param sex      (vector) Sex either \code{"female"} or \code{"male"}
#' 
#' \strong{ Optional }
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
#' @seealso \code{\link{model_plot}} for plotting the results and 
#' \code{\link{model_mean}} for aggregate data estimation. 
#' 
#' @examples
#' #One child
#' child_reference_FFMandFM(6, "male")
#' 
#' #Several children
#' child_reference_FFMandFM(sample(6:12, 10, replace = TRUE), 
#'                    sample(c("male","female"), 10, replace = TRUE))
#' 
#' @keywords internal
#' @export

child_reference_FFMandFM <- function(age, sex){
  
  #Change sex to numeric for c++
  newsex                         <- rep(0, length(sex))
  newsex[which(sex == "female")] <- 1
  
  #Check all variables are positive
  if (any(age < 0) ){
    stop("Cannot handle negative values for age")
  }
  
  #Check dimensions of inputs
  if (length(age) != length(sex)){
    stop("Dimension mismatch: age and sex must have same length.")
  }
  
  #Check sex is "male" and "female"
  if (length(which(!(sex %in% c("male","female")))) > 0){
    stop(paste0("Invalid sex. Please specify either 'male' of 'female'"))
  }
  
  #Check that we don't go over 18 yrs where we have no data
  if (max(age) > 18){
    warning(paste0("Some individuals reach age > 18 before model stops. Results",
                   " might not be accurate for adults. Please use adult_weight",
                   " instead."))
  }
  
  #Get c++ function
  mass_reference_wrapper(age, newsex)
}