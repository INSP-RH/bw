#' @title Get Mean results from Adult Weight Change Model
#'
#' @description Gets survey means \code{\link[survey]{svymean}}, standard error and
#' confidence interval estimates of \code{\link{adult_weight}} or \code{\link{child_weight}}.
#'
#' @param weight     (list) List from \code{\link{adult_weight}} or \code{\link{adult_weight}}.
#'
#' \strong{ Optional }
#' @param design A \code{survey.design} object. See \code{\link[survey]{svydesign}} 
#' for additional information on design objects. 
#' @param meanvars (vector) Strings indicating which variables are required to
#' estimate the mean. 
#' @param days   (vector) Vector of days in which to compute the estimates
#' @param confidence (numeric) Confidence level (\code{default = 0.95})
#' @param group (vector) Variable in which to group the results.
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' 
#' @details The default \code{design} is that of simple random sampling.
#' 
#' @importFrom survey svyby
#' @importFrom survey svymean
#' @importFrom survey svyvar
#' @importFrom stats update
#' @importFrom stats coef
#' @importFrom stats confint
#' @importFrom survey SE
#' 
#' @examples 
#' #EXAMPLE 1A: RANDOM SAMPLE MODELLING FOR ADULTS
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
#' #Create weight change model
#' model_weight <- adult_weight(weights, heights, ages, sexes, 
#'                              EIchange)
#'                              
#' #Calculate survey mean and variance for 25 days
#' aggregate_data <- model_mean(model_weight)
#' 
#' #You can plot the mean with ci
#' if(require(ggplot2)){
#' ggplot(subset(aggregate_data, variable == "Body_Weight")) + 
#'     geom_line(aes(x = time, y = mean)) +
#'     geom_line(aes(x = time, y = Lower_CI_mean), linetype = "dashed") +
#'     geom_line(aes(x = time, y = Upper_CI_mean), linetype = "dashed") +
#'     theme_classic() + xlab("Days") + ylab("Mean Body Weight (kg)")
#' }
#' 
#' #EXAMPLE 1C: RANDOM SAMPLE MODELLING FOR CHILDREN
#' #--------------------------------------------------------
#' #Antropometric data
#' FatFree <- c(32, 17.2, 18.8, 20, 24.1)
#' Fat     <- c(4.30, 2.02, 3.07, 1.12, 2.93)
#' ages    <- c(10, 6.2, 5.4, 4, 4.1)
#' sexes   <- c("male", "female", "female", "male", "male") 
#' 
#' #Returns a weight change matrix and other matrices
#' model_weight <- child_weight(ages, sexes, Fat, FatFree)
#' 
#' #Calculate survey mean and variance for 25 days
#' aggregate_data <- model_mean(model_weight)
#' 
#' #You can plot the mean with ci
#' if(require(ggplot2)){
#' ggplot(subset(aggregate_data, variable == "Body_Weight")) + 
#'     geom_line(aes(x = time, y = mean)) +
#'     geom_line(aes(x = time, y = Lower_CI_mean), linetype = "dashed") +
#'     geom_line(aes(x = time, y = Upper_CI_mean), linetype = "dashed") +
#'     theme_classic() + xlab("Days") + ylab("Mean Body Weight (kg)")
#' }
#'   
#' #EXAMPLE 2A: SURVEY DATA FOR ADULTS
#' #-------------------------------------------------------
#' #Data frame for use in survey
#' probs   <- runif(10, 20, 60)
#' datasvy <- data.frame(
#'   id    = 1:10,
#'   bw    = runif(10,60,90),
#'   ht    = runif(10, 1.5, 2),
#'   age   = runif(10, 18, 80),
#'   sex   = sample(c("male","female"),10, replace = TRUE),
#'   kcal  = runif(10, 2000, 3000),
#'   group = sample(c(0,1), 10, replace = TRUE),
#'   svyw  = probs/sum(probs))
#' 
#' #Days
#' days <- 365
#' 
#' #Energy intake matrix
#' EIchange <- matrix(NA, nrow = days, ncol = 0)
#' for(i in 1:nrow(datasvy)){
#'     EIchange <- cbind(EIchange, rep(datasvy$kcal[i], days))
#' }
#' 
#' #Calculate weight change                   
#' svyweight <- adult_weight(datasvy$bw, datasvy$ht, datasvy$age, 
#'                           datasvy$sex, EIchange)
#'                           
#' #Create survey design using survey package                           
#' design <- survey::svydesign(id = ~id, weights = datasvy$svyw, 
#' data = datasvy)
#' 
#' #Group to calculate means
#' group  <- datasvy$group     
#'     
#' #Calculate survey mean and variance for 25 days
#' aggregate_data <- model_mean(svyweight, design = design, group = group)
#' 
#' #You can plot the mean with ci
#' if(require(ggplot2)){
#' ggplot(subset(aggregate_data, variable == "Body_Weight")) + 
#'     geom_ribbon(aes(x = time, ymin = Lower_CI_mean, ymax = Upper_CI_mean,
#'     fill = factor(group)), alpha = 0.25) +
#'     geom_line(aes(x = time, y = mean, color = factor(group)), size = 2) +
#'     theme_classic() + xlab("Days") + ylab("Mean Body Weight (kg)") 
#' }                     
#' 
#' #EXAMPLE 2A: SURVEY DATA FOR CHILDREN
#' #-------------------------------------------------------
#' #Data frame for use in survey
#' probs   <- runif(10, 20, 60)
#' datasvy <- data.frame(
#'   id      = 1:10,
#'   age     = runif(10, 2, 12),
#'   sex     = sample(c("male","female"),10, replace = TRUE),
#'   fat     = runif(10, 2, 10),
#'   fatfree = runif(10, 8, 15),
#'   group   = sample(c(0,1), 10, replace = TRUE),
#'   svyw    = probs/sum(probs))
#' 
#' #Days
#' days <- 365
#' 
#' #Calculate weight change                   
#' svyweight <- child_weight(datasvy$age, datasvy$sex, datasvy$fat, datasvy$fatfree)
#'                           
#' #Create survey design using survey package                           
#' design <- survey::svydesign(id = ~id, weights = datasvy$svyw, 
#' data = datasvy)
#' 
#' #Group to calculate means
#' group  <- datasvy$group     
#'     
#' #Calculate survey mean and variance for 25 days
#' aggregate_data <- model_mean(svyweight, design = design, group = group)
#' 
#' #You can plot the mean with ci
#' if(require(ggplot2)){
#' ggplot(subset(aggregate_data, variable == "Body_Weight")) + 
#'     geom_ribbon(aes(x = time, ymin = Lower_CI_mean, ymax = Upper_CI_mean,
#'     fill = factor(group)), alpha = 0.25) +
#'     geom_line(aes(x = time, y = mean, color = factor(group)), size = 2) +
#'     theme_classic() + xlab("Days") + ylab("Mean Body Weight (kg)") 
#' }                     
#'                                                             
#' @export

model_mean <- function(weight, 
                       meanvars = names(weight)[-which(names(weight) %in% c("Time", "BMI_Category"))], 
                       days     = seq(0, length(weight[["Time"]]) - 1, length.out = 25),
                       group    = rep(1,nrow(weight[[meanvars[1]]])),
                       design   = svydesign(ids=~1, weights = rep(1,nrow(weight[[meanvars[1]]])),
                                            data = as.data.frame(weight[meanvars[1]])),
                       confidence = 0.95){
  
  #Throw warning that it will take time
  if (length(days) > 50){
    warning("This process will take some time")
  }
  
  #Check confidence
  if(confidence > 1 || confidence <= 0){
    stop("Invalid confidence level. Confidence must be between 0 and 1")
  }
  
  #Check that BMI_Category not in meanvars
  if ("BMI_Category" %in% meanvars){
    stop("Cannot estimate BMI_Category mean. Please use 'adult_bmi' function instead.")
  }
  
  #Check that meanvars are in names(weight)
  if (!all(meanvars %in% names(weight))){
    stop(paste0("Not all variables specified in meanvars are available ",
                "in weight. You must use one of the following: '", 
                paste0(names(weight)[-which(names(weight) %in% c("Time", "BMI_Category", "Age"))], collapse = "', '"),"'."))
  }
  
  #Check that time is part of weight
  if (!("Time" %in% names(weight))){
    stop("Invalid weight parameter. Weight must include vector 'Time'.")
  }
  
  #Set time to integers
  days <- which(weight[["Time"]] %in% floor(days))
  
  #Get number of variables to plot
  nvars <- length(meanvars) 
  
  #Create empty data frame
  weightdata <- data.frame(matrix(NA, nrow = 0, ncol = 11))
  
  #Update design to add group
  design <- update(design, group = group)
  
  #Loop through every day
  for(t in 1:length(days)){
    
    #Loop through each of the variables under consideration
    for (var in 1:nvars){
      
      #Weight update to add variable of interest
      design <- update(design, myvar = weight[[meanvars[var]]][,days[t]])
      
      #Get mean and var
      mymean <- svyby(~myvar, ~group, design, svymean)
      myvar  <- svyby(~myvar, ~group, design, svyvar)
      
      #Get confidence intervals for each
      confmean <- confint(mymean, level = confidence)
      confvar  <- confint(myvar, level = confidence)
      
      #Get into data frame
      thisdata <- data.frame(days[t], meanvars[var], mymean$group, 
                        coef(mymean), SE(mymean), confmean, 
                        coef(myvar), SE(myvar), confvar)
      
      #Bind together
      weightdata <- rbind(weightdata, thisdata)
      
    }
  }
  
  #Add column names
  colnames(weightdata) <- c("time", "variable", "group", "mean", "SE_mean",
                            "Lower_CI_mean", "Upper_CI_mean", "variance",
                            "SE_variance", "Lower_CI_variance", 
                            "Upper_CI_variance")
  
  #Remove rownames
  rownames(weightdata) <- c()
  
  #Return data frame
  return(weightdata)
  
}