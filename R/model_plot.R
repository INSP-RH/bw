#' @title Plot Results from Weight Change Model
#'
#' @description Generates a plot for list from \code{\link{adult_weight}} or
#' \code{\link{child_weight}}. 
#'
#' @param model     (list) List from \code{\link{adult_weight}} or \code{\link{child_weight}}
#'
#' \strong{ Optional }
#' @param plotvars   (vector) String vector of the plots to generate (default generates all)
#' @param title      (string) Title of plot collection
#' @param ncol       (string) Number of columns to include in plot
#' @param timevar    (string) String indicating which of the variables in model list indicates time.
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
#' mymodel <- adult_weight(80, 1.8, 40, "female", rep(-100, 365))
#' 
#' #You can plot all the variables
#' model_plot(mymodel)
#' 
#' #Or only one of them
#' model_plot(mymodel, "Body_Weight", ncol = 1)
#' 
#' #EXAMPLE 1C: INDIVIDUAL MODELLING FOR CHILDREN
#' #--------------------------------------------------------
#' mymodel <- child_weight(5, "female", 12, 4)
#' 
#' #You can plot all the variables
#' model_plot(mymodel)
#' 
#' #Or only one of them and specify by age
#' model_plot(mymodel, "Body_Weight", ncol = 1)
#' 
#' #EXAMPLE 2A: DATASET MODELLING FOR ADULTS
#' #--------------------------------------------------------
#' \donttest{
#' #Antropometric data
#' models <- c(45, 67, 58, 92, 81)
#' heights <- c(1.30, 1.73, 1.77, 1.92, 1.73)
#' ages    <- c(45, 23, 66, 44, 23)
#' sexes   <- c("male", "female", "female", "male", "male") 
#' 
#' #Matrix of energy consumption reduction: 
#' EIchange <- rbind(rep(-100, 365), rep(-200, 365), rep(-200, 365), 
#'                   rep(-123, 365), rep(-50, 365))
#' 
#' #Returns a model change matrix and other matrices
#' model_model <- adult_weight(models, heights, ages, sexes, 
#'                              EIchange)
#'                              
#' #Create all plots                            
#' model_plot(model_model)
#' 
#' #Plot Body Mass Index
#' model_plot(model_model, "Body_Mass_Index")
#' }
#' 
#' \donttest{
#' #EXAMPLE 2C: DATASET MODELLING FOR CHILDREN
#' #--------------------------------------------------------
#' #Antropometric data
#' FatFree <- c(32, 17.2, 18.8, 20, 24.1)
#' Fat     <- c(4.30, 2.02, 3.07, 1.12, 2.93)
#' ages    <- c(10, 6.2, 5.4, 4, 4.1)
#' sexes   <- c("male", "female", "female", "male", "male") 
#' 
#' #Returns a model change matrix and other matrices
#' model_model <- child_weight(ages, sexes, Fat, FatFree)
#'                              
#' #Create all plots                            
#' model_plot(model_model)
#' 
#' #Plot Body Mass Index
#' model_plot(model_model, "Fat_Mass")
#' }                           
#' @export

model_plot <- function(model, 
                       plotvars = names(model)[-which(names(model) %in% c("Time", "BMI_Category", "Age", "Correct_Values", "Model_Type"))], 
                       timevar  = "Time", title = "Hall's model results", ncol = 2){
  
  #Check object is list
  if (!is.list(model)){
    stop("Please input a model object comming from child_weight or adult_weight.")
  }
  
  #Check timevar makes sense
  if (!(timevar %in% c("Age","Time"))){
    stop(paste("Invalid timevar = ", timevar, "please select 'Time' or 'Age'."))
  }
  
  #Check that we have only one individual
  if (nrow(as.matrix(model[[timevar]])) > 1 && timevar == "Age"){
    warning("We cannot plot against age for several individuals. Defaulting timevar = 'Time'")
    timevar  <- "Time"
  }
  
  
  #Check title is string
  if(!is.character(title)){
    stop("Invalid plot title. Please input a string.")
  } 
  
  #Check number of columns is numeric
  if(!all.equal(ncol, as.integer(ncol))){
    stop(paste("Invalid ncol = ", ncol,
          "Input an integer value for the number of columns ncol."))
  }
  
  #Check that timevar is in list
  if(!(timevar %in% names(model))){
    stop(paste(timevar, " is not part of names(model):", paste0(names(model), collapse = ", ")))
  }
  
  #Get time variable
  time <- model[[timevar]]
  
  #Check for child model
  if ("Model_Type" %in% names(model)){
    if (model[["Model_Type"]] == "Adult"){
      if(!all(plotvars %in% c("Adaptive_Thermogenesis", "Extracellular_Fluid",
                           "Glycogen", "Fat_Mass", "Lean_Mass", "Body_Weight", 
                           "Body_Mass_Index", "Energy_Intake"))){
        warning(paste("Not all specified plotvars are in current model. For Adult \n",
                      "Adaptive_Thermogenesis, Extracellular_Fluid, Glycogen, Energy_Intake,",
                      "Fat_Mass, Lean_Mass, Body_Weight & Body_Mass_Index",
                      "are the valid variables."))}
    } else if (model[["Model_Type"]] == "Children"){
      if(!all(plotvars %in% c("Fat_Mass", "Fat_Free_Mass", "Body_Weight"))){
        warning(paste("Not all specified plotvars are in current model. For Children \n",
                      "Fat_Mass, Fat_Free_Mass, Body_Weight",
                      "are the valid variables."))}
    } else {
      warning(paste("Unknown Model_Type:",model[["Model_Type"]] ))
    }
  } else {
    warning("No 'Model_Type' found in names(model)")
  }
  
  #Loop to check if the plotvars are in model and numeric
  plotvars2 <- c()
  for (i in 1:length(plotvars)){
    if (!(plotvars[i] %in% names(model))){
      warning(paste0("Ignoring plotvar '", plotvars[i],"' as it is not in current model."))
    } else {
      plotvars2 <- c(plotvars2, plotvars[i])
    }
  }
  
  #Update
  plotvars <- plotvars2
  
  #Get all other variables and assign to each a plot
  nplots <- length(plotvars)
  
  if (nplots < 1){stop(paste("No valid variables for plotting."))}
  

  #Set ncol to 1 if nplots = 1
  if (nplots == 1){ncol = 1}
  
  #Create plotlist
  plotlist <- list()
  
  #Create a plot for each variable
  for (i in 1:nplots){
    
    #Get data
    plot_data      <- as.data.frame(t(model[[plotvars[i]]]))
    plot_data$id   <- time[1:length(time)]
    assign(paste0("plot_data",i), melt(plot_data, id.var="id"))
    
    plotlist[[i]] <- 
           ggplot(get(paste0("plot_data",i))) + 
             geom_line(aes_string(x = "id", 
                                  y = "value", 
                                  group = "variable", 
                                  colour = "variable")) +
             xlab(timevar) + ylab(gsub("_"," ", plotvars[i])) + theme_classic() + 
             theme(legend.position = "none") 
    
  }
  
  if (nplots == 1){ plotlist[[1]] + ggtitle(title) } else {
    do.call("grid.arrange", c(plotlist, ncol= ncol, top = title))  
  }
  
  
}