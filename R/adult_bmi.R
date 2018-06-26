#' @title Get BMI prevalence results from Adult Weight Change Model
#'
#' @description Gets survey proportions \code{\link[survey]{svytable}}, standard error and
#' confidence interval estimates of BMI from \code{\link{adult_weight}}. 
#'
#' @param weight     (list) List from \code{\link{adult_weight}}
#'
#' \strong{ Optional }
#' @param design A \code{survey.design} object. See \code{\link[survey]{svydesign}} 
#' for additional information on design objects. 
#' @param days   (vector) Vector of days in which to compute the estimates
#' @param confidence (numeric) Confidence level (\code{default = 0.95})
#' @param group (vector) Variable in which to group the results.
#' 
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' 
#' @details The default \code{design} is that of simple random sampling.
#' 
#' @importFrom survey svyby
#' @importFrom stats update
#' @importFrom survey svymean
#' @importFrom survey svydesign
#' @importFrom stats coef
#' @importFrom stats confint
#' 
#' @examples 
#' #EXAMPLE 1: RANDOM SAMPLE MODELLING
#' #--------------------------------------------------------
#' 
#' #Antropometric data
#' weights <- c(45, 67, 58, 67, 81)
#' heights <- c(1.30, 1.73, 1.77, 1.92, 1.73)
#' ages    <- c(45, 23, 66, 44, 23)
#' sexes   <- c("male", "female", "female", "male", "male") 
#' 
#' #Matrix of energy consumption reduction: 
#' EIchange <- rbind(rep(-100, 365), rep(-200, 365), rep(-200, 365), 
#'                   rep(-123, 365), rep(-50, 365))
#' 
#' #Create weight change model
#' model_weight <- adult_weight(weights, heights, ages, sexes, 
#'                              EIchange)
#'                              
#' #Calculate proportions
#' adult_bmi(model_weight)
#' 
#' #EXAMPLE 2: Survey data
#' #-------------------------------------------------------
#' set.seed(7423)
#' 
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
#' #Days to model
#' days <- 365
#' 
#' #Energy intake matrix
#' EIchange <- matrix(NA, ncol = days, nrow = 0)
#' for(i in 1:nrow(datasvy)){
#'     EIchange <- rbind(EIchange, rep(datasvy$kcal[i], days))
#' }
#' 
#' #Calculate weight change                   
#' weight <- adult_weight(datasvy$bw, datasvy$ht, datasvy$age, 
#'                           datasvy$sex, EIchange)
#' 
#' 
#' #Create survey design using survey package                           
#' design <- survey::svydesign(id = ~id, weights = datasvy$svyw, 
#' data = datasvy)
#'    
#' #' #Group to calculate means
#' group  <- datasvy$group     
#' 
#' #Calculate survey mean and variance for 25 days
#' adult_bmi(weight, design = design, group = group)
#'  
#' @export

adult_bmi  <- function(weight, 
                       days   = seq(0, length(weight[["Time"]])-1, length.out = 25),
                       group  = rep(1,nrow(weight[["BMI_Category"]])),
                       design = svydesign(ids=~1, weights = rep(1,nrow(weight[["BMI_Category"]])),
                                          data = as.data.frame(weight[["BMI_Category"]])),
                       confidence = 0.95){
  
  #Throw message that it will take time
  if (length(days) > 50){
    message("This process will take some time...")
  }
  
  # Check time values
  if(any(unlist(lapply(days, function(day){
    if(day %in% weight$Time){TRUE}
    else{FALSE}})))==FALSE){
    stop("Some time values are not available in object weight")
  }
  
  # Check groups do not exceed individuals
  if(length(group)>1 & length(group)!=nrow(weight$BMI_Category)){
    stop(paste("Dimension mismatch.",
               "Group must be defined for every individual or a unique for all individuals."))
  }
  
  #Check confidence
  if(confidence > 1 || confidence <= 0){
    warning("Invalid confidence level. Confidence must be between 0 and 1")
  }
  
  #Create empty data frame
  mydata <- as.data.frame(matrix(NA, ncol = 5, nrow = 0))
  
  #Set time to integers
  days <- which(weight[["Time"]] %in% floor(days))
  
  #Update design to add group
  design <- update(design, group = group)
  
  #Loop through every day
  for(t in 1:length(days)){
    
    #Weight update to add variable of interest
    myvar  <- as.factor(weight[["BMI_Category"]][,days[t]])
    design <- update(design, bmi_ = myvar)
    
    #Get mean and ci
    if (length(levels(myvar)) < 2){
      mu        <- 1
      confmean  <- matrix(NA, ncol = 2)
      names(mu) <- levels(myvar)
    } else {
      mymean    <- svyby(~bmi_, ~group, design, svymean)
      confmean  <- confint(mymean, level = confidence)
      mu        <- coef(mymean)
    }
    
    
    
    #Add to same data frame
    varnames <- unlist(lapply(names(mu), function(x){gsub(".*bmi_","",x)}))
    
    #Empty names to allow for data frame to work and bind
    names(mu)          <- c()
    colnames(confmean) <- c()
    
    #Create data frame
    if(exists("mymean")){
      today    <- data.frame(Day = weight[["Time"]][days[t]], 
                             Group = mymean$group, 
                             BMI_Category = varnames, 
                             Mean = mu,
                             confmean)
    }else{
      today    <- data.frame(Day = weight[["Time"]][days[t]], 
                             Group = NA, 
                             BMI_Category = varnames, 
                             Mean = mu,
                             confmean)
    }

    
    #Bind to previous data
    mydata <- rbind(mydata, today)
    
    
  }
  
  #Clear rownames
  rownames(mydata) <- c()
  
  #Add colnames
  colnames(mydata) <- c(colnames(mydata)[1:(ncol(mydata) - 2)], colnames(confmean))
  
  
  #Return data frame
  return(mydata)
  
}