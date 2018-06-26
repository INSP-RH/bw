context("Mean of model")

test_that("Checking mean errors",{
  
  #Confidence > 1
  expect_error({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
      id      = 1:10,
      age     = runif(10,20,60),
      sex     = sample(c("male","female"),10, replace = TRUE),
      weight  = runif(10,60,80),
      height  = runif(10,160,180),
      group   = sample(c(0,1), 10, replace = TRUE),
      svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                        data = datasvy)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, 
                                 datasvy$age, datasvy$sex, days = 3)
    
    #Calculate survey mean and variance for 25 days
    model_mean(model_weight, design = design, days = 1, confidence = 2)
    
  })
  
  #Negative confidence 
  expect_error({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
      id      = 1:10,
      age     = runif(10,20,60),
      sex     = sample(c("male","female"),10, replace = TRUE),
      weight  = runif(10,60,80),
      height  = runif(10,160,180),
      group   = sample(c(0,1), 10, replace = TRUE),
      svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                        data = datasvy)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, 
                                 datasvy$age, datasvy$sex, days = 3)
    
    #Calculate survey mean and variance for 25 days
    model_mean(model_weight, design = design, days = 1, confidence = -2)
    
  })
  
  #Adult bmi to calculate bmi category
  expect_error({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
      id      = 1:10,
      age     = runif(10,20,60),
      sex     = sample(c("male","female"),10, replace = TRUE),
      weight  = runif(10,60,80),
      height  = runif(10,160,180),
      group   = sample(c(0,1), 10, replace = TRUE),
      svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                        data = datasvy)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, 
                                 datasvy$age, datasvy$sex, days = 3)
    
    #Calculate survey mean and variance for 25 days
    model_mean(model_weight, design = design, days = 1, meanvars = "BMI_Category")
    
  })
  
  #Different value  meanvar
  expect_error({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
      id      = 1:10,
      age     = runif(10,20,60),
      sex     = sample(c("male","female"),10, replace = TRUE),
      weight  = runif(10,60,80),
      height  = runif(10,160,180),
      group   = sample(c(0,1), 10, replace = TRUE),
      svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                        data = datasvy)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, 
                                 datasvy$age, datasvy$sex, days = 3)
    
    #Calculate survey mean and variance for 25 days
    model_mean(model_weight, design = design, days = 1, meanvars = "I_am_bored")
    
  })
  
  #Check that time is on list
  expect_error({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
      id      = 1:10,
      age     = runif(10,20,60),
      sex     = sample(c("male","female"),10, replace = TRUE),
      weight  = runif(10,60,80),
      height  = runif(10,160,180),
      group   = sample(c(0,1), 10, replace = TRUE),
      svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                        data = datasvy)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, 
                                 datasvy$age, datasvy$sex, days = 3)
    
    #Calculate survey mean and variance for 25 days
    timelist <- which(names(model_weight) == "Time")
    model_mean(model_weight[-timelist], design = design, days = 1)
    
  })
  
})

test_that("Checking mean warnings",{
  
  # Warning for taking time to calculate
  expect_warning({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
       id      = 1:10,
       age     = runif(10,20,60),
       sex     = sample(c("male","female"),10, replace = TRUE),
       weight  = runif(10,60,80),
       height  = runif(10,160,180),
       group   = sample(c(0,1), 10, replace = TRUE),
       svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                                data = datasvy)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, datasvy$age, datasvy$sex, days = 51)
    
    #Calculate survey mean and variance for 25 days
    model_mean(model_weight, design = design, days = 1:51)
    
  })
  
  expect_warning({
    
    #Antropometric data
    datasvy <- data.frame(
      id      = 1,
      age     = 30,
      sex     = "female",
      weight  = 60,
      height  = 1.89)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, datasvy$age, datasvy$sex, days = 1)
    
    #Calculate survey mean and variance for 25 days
    bw   <- model_mean(model_weight, days = 0)
    
  })
  
})

test_that("Test mean object",{
  
  #Check results are groupped
  expect_true({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
      id      = 1:10,
      age     = rep(30,10),
      sex     = sample(c("male","female"),10, replace = TRUE),
      weight  = runif(10,60,80),
      height  = runif(10,160,180),
      group   = rep(1, 10),
      svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                        data = datasvy)
    
    #Group to calculate means
    group  <- datasvy$group  
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, datasvy$age, datasvy$sex, days = 2)
    
    #Calculate survey mean and variance for 25 days
    all(model_mean(model_weight, design = design, group = group, days = 1:2)$group == 1)
    
  })
  
  #Check results are groupped
  expect_true({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
      id      = 1:10,
      age     = rep(30,10),
      sex     = sample(c("male","female"),10, replace = TRUE),
      weight  = runif(10,60,80),
      height  = runif(10,160,180),
      group   = c(rep(1, 5), rep(2, 5)),
      svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                        data = datasvy)
    
    #Group to calculate means
    group  <- datasvy$group  
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, datasvy$age, datasvy$sex, days = 2)
    
    #Calculate survey mean and variance for 25 days
    length(unique(model_mean(model_weight, design = design, group = group, days = 1:2)$group)) == 2
    
  })
  
  #Check results are correct
  expect_true({
    
    #Antropometric data
    probs   <- runif(10, 20, 60)
    datasvy <- data.frame(
      id      = 1:10,
      age     = rep(30,10),
      sex     = sample(c("male","female"),10, replace = TRUE),
      weight  = runif(10,60,80),
      height  = runif(10,160,180),
      svyw    = probs/sum(probs))
    
    #Create survey design using survey package                           
    design <- svydesign(id = ~id, weights = datasvy$svyw, 
                        data = datasvy)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, datasvy$age, datasvy$sex, days = 2)
    
    #Calculate survey mean and variance for 25 days
    bw   <- model_mean(model_weight, design = design, days = 0)
    hall <- bw[which(bw[,which(colnames(bw) == "variable")] == "Body_Weight"),which(colnames(bw) == "mean")]
    
    #Check that the mean is correct
    abs(weighted.mean(datasvy$weight, datasvy$svyw) - hall) < 0.0001
    
  })
  
  #Check results are correct for only 1 individual
  expect_warning(expect_true({
    
    #Antropometric data
    datasvy <- data.frame(
      id      = 1,
      age     = 30,
      sex     = "female",
      weight  = 60,
      height  = 1.89)
    
    #Returns a weight change matrix and other matrices
    model_weight <- adult_weight(datasvy$weight, datasvy$height, datasvy$age, datasvy$sex, days = 1)
    
    #Calculate survey mean and variance for 25 days
    bw   <- model_mean(model_weight, days = 0)
    hall <- bw[which(bw[,which(colnames(bw) == "variable")] == "Body_Weight"),which(colnames(bw) == "mean")]
    
    #Check that the mean is correct
    abs(hall - datasvy$weight) < 0.0001
    
  }))
  
})
