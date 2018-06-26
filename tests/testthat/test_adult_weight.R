context("Adult weight change function")

test_that("Checking adult_weight  errors",{
  
  # Check that age >= 0
  expect_error({
    adult_weight(bw = 76, ht = 1.73, age = -36, sex = "male")  
  })
  
  # Check that weight > 0 
  expect_error({
    adult_weight(bw = 0, ht = 1.73, age = 36, sex = "female")  
  })
  
  # Check that FFM >= 0 
  expect_error({
    adult_weight(bw = 76, ht = 0, age = 36, sex = "male")  
  })
  
  
  # Check that days > 0 
  expect_error({
    adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", days = 0)  
  })
  
  
  # Check that bw, ht, age, sex, PAL, pcarb_base, and p_carb have the same length
  # bw
  expect_error({
    adult_weight(bw = c(76), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"))  
  })
  
  # ht
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73,1.6,1.8), age = c(36,43),
                 sex = c("male", "female"))    })
  # age
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43,76,19),
                 sex = c("male", "female"))    })
  # sex
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male"))    
    })
  
  # PAL
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"), PAL = 1.4)    
  })
  
  # pcarb_base
  # sex
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"), pcarb_base = c(.5,.6,.8))    
  })
  # pcarb
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"), pcarb = c(.5,.6,.8))    
  })
  
  # Check that sex is "male" or "female"
  expect_error({
    adult_weight(bw = 76, ht = 1.73, age = 36, sex = "Mail")  
  })
  
  # Check that time step (dt) is less than time to run the model
  
  expect_error({
    adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", days = 30, dt=365)  
  })
  
  # Check that dt is positive
  expect_error({
    adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", days = 30, dt=-1)  
  })
  
  # Check change in energy intake and change in sodium have the same dimension
  expect_error({
    adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", 
                 EIchange = rep(-40, 366), NAchange = rep(-0.3,365))  
  })
  # Consumption matrix is defined by rows 
  # i.e Consumption for the i-th individual is in the i-th row
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"), 
                 EIchange = t(matrix(c(rep(-20,365), rep(-5, 365)),byrow = T,
                                   nrow = 2)),
                 NAchange = t(matrix(c(rep(-2,365), rep(-1, 365)),byrow = T,
                                    nrow = 2)), days = 365)      
  })
  
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"), 
                 EIchange = (matrix(c(rep(-20,365), rep(-5, 365)),byrow = T,
                                     nrow = 2)),
                 NAchange = (matrix(c(rep(-2,365), rep(-1, 365)),byrow = T,
                                     nrow = 2)), days = 365, 
                 pcarb = c(1.2, 0.6))      
  })
  
  expect_error({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"), 
                 EIchange = (matrix(c(rep(-20,365), rep(-5, 365)),byrow = T,
                                    nrow = 2)),
                 NAchange = (matrix(c(rep(-2,365), rep(-1, 365)),byrow = T,
                                    nrow = 2)), days = 365, 
                 pcarb = c(0.4, 0.6), PAL=c(1.6,0))      
  })
})

test_that("Checking adult_weight warning",{
  # Age less than 18
  expect_warning({
    adult_weight(bw = 76, ht = 1.73, age = 17, sex = "male")  
  })
  # EIchange and NAchange are reccomended to have the same dimensions as the 
  # number of days to run the model
  expect_warning({
    adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", 
                 EIchange = rep(-40, 363), NAchange = rep(-0.3,363))  
  })
  # Extreme PAL levels
  # Low PAL
  expect_warning({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"), 
                 EIchange = (matrix(c(rep(-20,365), rep(-5, 365)),byrow = T,
                                    nrow = 2)),
                 NAchange = (matrix(c(rep(-2,365), rep(-1, 365)),byrow = T,
                                    nrow = 2)), days = 365, 
                 pcarb = c(0.4, 0.6), PAL=c(1.6,1.2))      
  })
  
  # High PAL
  expect_warning({
    adult_weight(bw = c(76,54), ht = c(1.73, 1.6), age = c(36,43),
                 sex = c("male", "female"), 
                 EIchange = (matrix(c(rep(-20,365), rep(-5, 365)),byrow = T,
                                    nrow = 2)),
                 NAchange = (matrix(c(rep(-2,365), rep(-1, 365)),byrow = T,
                                    nrow = 2)), days = 365, 
                 pcarb = c(0.4, 0.6), PAL=c(1.6,2.7))      
  })
  
})

test_that("Checking adult_weight results",{
  # Check input gives same result to default
  # PAL
  expect_equal(adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male"),
               adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", PAL = 1.5))
  
  #Proportion of carbohydrates
  expect_equal(adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male"),
               adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male",
                            pcarb_base = 0.5, pcarb = 0.5))
  
  
  # Check similar results to Hall programmed model. Relative error less than 5%.

  #Example 1: Male, 76 kg, 1.73 m,  36 years old, PAL 1.5
  
  # Energy to maintain 70 kg (2360 kcal)
  expect_lt({
    result <- adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", PAL = 1.5,
                 EI=rep(2360,365))$Body_Weight[365]
    abs(result-70)/70
  }, 0.05)
  
  # Energy to reach 70 kg (2287 kcal)
  expect_lt({
    result <- adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", PAL = 1.5,
                           EI=rep(2287,365))$Body_Weight[365]
    abs(result-70)/70
  }, 0.05)
  
  # Energy to maintain 76 kg (2503 kcal)
  expect_lt({
    result <- adult_weight(bw = 76, ht = 1.73, age = 36, sex = "male", PAL = 1.5,
                           EI=rep(2503,365))$Body_Weight[365]
    abs(result-76)/76
  }, 0.05)
  
  # Example 2: Female 58 kg, 1.64m, 21 years old, PAL=1.7 
  
  # Energy to maintain 55 kg (2187)
  expect_lt({
    result <- adult_weight(bw = 58, ht = 1.64, age = 21, sex = "female", PAL = 1.7,
                           EI=rep(2187,365))$Body_Weight[365]
    abs(result-55)/55
  }, 0.05)
  
  # Energy to reach 55 kg (2187)
  expect_lt({
    result <- adult_weight(bw = 58, ht = 1.64, age = 21, sex = "female", PAL = 1.7,
                           EI=rep(2159,365))$Body_Weight[365]
    abs(result-55)/55
  }, 0.05)
  
  # Energy to maintain 58 kg (2278)
  expect_lt({
    result <- adult_weight(bw = 58, ht = 1.64, age = 21, sex = "female", PAL = 1.7,
                           EI=rep(2278,365))$Body_Weight[365]
    abs(result-58)/58
  }, 0.05)
 
})
