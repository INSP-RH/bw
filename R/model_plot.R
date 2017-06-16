#' @title Plot Results from Weight Change Model
#'
#' @description Generates a plot for list from \code{\link{adult_weight}} or
#' \code{\link{child_weight}}. 
#'
#' @param weight     (list) List from \code{\link{adult_weight}} or \code{\link{child_weight}}
#'
#' \strong{ Optional }
#' @param plotvars   (vector) String vector of the plots to generate (default generates all)
#' @param title      (string) Title of plot collection
#' @param ncol       (string) Number of columns to include in plot
#' @param timevar    (string) String indicating which of the variables in weight list indicates time.
#'
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @details It returns a grid object
#' 
#' @import ggplot2
#' @import gridExtra
#' @importFrom reshape2 melt
#' 
#' @examples 
#' #EXAMPLE 1A: INDIVIDUAL MODELLING FOR ADULTS
#' #--------------------------------------------------------
#' myweight <- adult_weight(80, 1.8, 40, "female", rep(-100, 365))
#' 
#' #You can plot all the variables
#' model_plot(myweight)
#' 
#' #Or only one of them
#' model_plot(myweight, "Body_Weight", ncol = 1)
#' 
#' #EXAMPLE 1C: INDIVIDUAL MODELLING FOR CHILDREN
#' #--------------------------------------------------------
#' myweight <- child_weight(5, "female", 12, 4)
#' 
#' #You can plot all the variables
#' model_plot(myweight)
#' 
#' #Or only one of them
#' model_plot(myweight, "Body_Weight", ncol = 1)
#' 
#' #EXAMPLE 2A: DATASET MODELLING FOR ADULTS
#' #--------------------------------------------------------
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
#'                              EIchange)
#'                              
#' #Create all plots                            
#' model_plot(model_weight)
#' 
#' #Plot Body Mass Index
#' model_plot(model_weight, "Body_Mass_Index")
#' 
#' #EXAMPLE 2C: DATASET MODELLING FOR CHILDREN
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
#' #Create all plots                            
#' model_plot(model_weight)
#' 
#' #Plot Body Mass Index
#' model_plot(model_weight, "Fat_Mass")
#'                             
#' @export

model_plot <- function(weight, 
                       plotvars = names(weight)[-which(names(weight) %in% c("Time", "BMI_Category", "Age"))], 
                       timevar  = "Time", title = "Hall's model results", ncol = 2){
  
  #Get time variable
  time <- weight[[timevar]]
  
  #Get all other variables and assign to each a plot
  nplots <- length(plotvars)
  
  #Set ncol to 1 if nplots = 1
  if (nplots == 1){ncol = 1}
  
  #Create plotlist
  plotlist <- list()
  
  #Create a plot for each variable
  for (i in 1:nplots){
    
    #Get data
    plot_data      <- as.data.frame(t(weight[[plotvars[i]]]))
    plot_data$id   <- 1:nrow(plot_data)
    assign(paste0("plot_data",i), melt(plot_data, id.var="id"))
    
    plotlist[[i]] <- 
           ggplot(get(paste0("plot_data",i))) + 
             geom_line(aes_string(x = "id", 
                                  y = "value", 
                                  group = "variable", 
                                  colour = "variable")) +
             xlab("Time") + ylab(gsub("_"," ", plotvars[i])) + theme_classic() + 
             theme(legend.position = "none") 
    
  }
  
  do.call("grid.arrange", c(plotlist, ncol= ncol, top = title))
  
}