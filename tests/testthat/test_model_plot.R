context("Plotting of model")

test_that("Checking plot errors",{
  
  # Check that is list the object
  expect_error({
    model_plot(1)
  })
  
  # Check that the list contains the objects to plot
  expect_error({
    model_plot(list("Hello" = 12,  "Bye" = 43))
  })
  
  # Check that the list contains the objects to plot
  expect_error({
    model_plot(adult_weight(80,180,33,"male", days = 3), timevar = "Unicorn")
  })
  
  # Check that title is string
  expect_error({
    model_plot(adult_weight(80,180,33,"male", days = 3), title = function(x){x^2})
  })
  
  # Warning for plotvars of adult used in children
  expect_warning(expect_error({
    model_plot(child_weight(10,"male", days = 3), plotvars = c("Extracellular_Fluid"))
  }))
  
  # Check that ncol is integer
  expect_error({
    model_plot(adult_weight(80,180,33,"male", days = 3), ncol = pi)
  })
  
  # Check that ncol is integer
  expect_error({
    mymodel <- adult_weight(80,180,33,"male", days = 3)
    model_plot(mymodel[-1])
  })
  
  # Check that thwere is at least one plotvar
  expect_warning(expect_error({
    model_plot(adult_weight(80,180,33,"male", days = 3), plotvars = c("Hello"))
  }))
  
})

test_that("Checking plot warnings",{
  
  # Warning for plotvars not present
  expect_warning({
    model_plot(adult_weight(80,180,33,"male", days = 3), plotvars = c("Hello", "Lean_Mass"))
  })
  
  # Warning for invalid time for several individuals
  expect_warning({
    model_plot(adult_weight(c(70,80),rep(180,2),rep(33,2),rep("male",2), days = 3), 
               timevar = "Age")
  })
  
  # Warning for invalid time for several individuals
  expect_warning({
    mymodel <- adult_weight(80,180,33,"male", days = 3)
    mymodel$Model_Type <- "Martian"
    model_plot(mymodel)
  })
  
  # Warning for plotvars of adult used in children
  expect_warning({
    model_plot(child_weight(10,"male", days = 3), plotvars = c("Extracellular_Fluid","Body_Weight"))
  })
  
})

test_that("Test plot object",{
  
  #If only one colname is specified it should be a ggplot object
  expect_true(is.ggplot(model_plot(child_weight(10,"male", days = 3), "Body_Weight")))
  
  #Check plot labels for adult model
  expect_true({
    
    #Check that xlab and ylab match
    boolvec <- c()
    for (timevar in c("Time","Age")){
      for (var in c("Adaptive_Thermogenesis", "Extracellular_Fluid","Glycogen", 
                    "Fat_Mass", "Lean_Mass", "Body_Weight", "Body_Mass_Index")){
        myplot <- model_plot(adult_weight(80,180,22,"male", days = 3), var, timevar = timevar)    
        boolvec <- c(boolvec, all(c(myplot$labels$x == timevar, myplot$labels$y == gsub("_", " ", var))))
      }
    }
    
    all(boolvec)
    
  })
  
  #Check plot labels for child model
  expect_true({
    
    #Check that xlab and ylab match
    boolvec <- c()
    for (timevar in c("Time","Age")){
      for (var in c("Fat_Mass", "Fat_Free_Mass", "Body_Weight")){
        myplot <- model_plot(child_weight(7,"male", days = 3), var, timevar = timevar, title = "Title")    
        boolvec <- c(boolvec, all(c(myplot$labels$x == timevar, 
                                    myplot$labels$title == "Title", 
                                    myplot$labels$y == gsub("_", " ", var))))
      }
    }
    
    all(boolvec)
    
  })
  
  #Check ncol for child model
  expect_true({
    
    boolvec <- c()
    mymodel <- child_weight(7,"male", days = 3)
    for(i in 1:3){
      mygrid  <- model_plot(mymodel, ncol = i)
      boolvec <- c(boolvec, ncol(mygrid) == i)
    }
    
    all(boolvec)
    
  })
  
  #Check ncol for adult model
  expect_true({
    
    boolvec <- c()
    mymodel <- adult_weight(80,180,33,"male", days = 3)
    for(i in 1:8){
      mygrid  <- model_plot(mymodel, ncol = i)
      boolvec <- c(boolvec, ncol(mygrid) == i)
    }
    all(boolvec)
    
  })
  
})
