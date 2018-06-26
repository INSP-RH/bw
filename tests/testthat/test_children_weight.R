context("Child weight change function")

test_that("Checking child_weight  errors",{
  
  # Check that age >= 0
  expect_error({
    child_weight(age=-1, sex="female", FM=2.7, FFM = 16)
  })
  
  # Check that FM >= 0 
  expect_error({
    child_weight(age=5, sex="female", FM=-4, FFM = 16)
  })
  
  # Check that FFM >= 0 
  expect_error({
    child_weight(age=5, sex="female", FM=2.7, FFM = -6)
  })
  
  
  # Check that days > 0 
  expect_error({
    child_weight(age=5, sex="female", FM=2.7, FFM = 16, days = 0)
  })
  
  
  # Check that sex, age, FM and FFM have the same length
  expect_error({
    child_weight(age=c(5,2), sex="female", FM=2.7, FFM = 16)
  })
  
  expect_error({
    child_weight(age=4, sex=c("female", "male"), FM=2.7, FFM = 16)
  })
  
  expect_error({
    child_weight(age=4, sex="female", FM=c(2.7,3), FFM = 16)
  })
  
  expect_error({
    child_weight(age=4, sex="female", FM=2.7, FFM = c(16,2))
  })
  
  # Check that sex is "male" or "female"
  expect_error({
    child_weight(age=5, sex="Female", FM=2.7, FFM = 16)
  })
  
  # Check that time step is less than time to run the model
  
  expect_error({
    child_weight(age=5, sex="female", FM=2.7, FFM = 16, days=10, dt=12)
  })
  
  # Check consumption has correct dimensions
  expect_error({
    child_weight(age=5, sex="female", FM=2.7, FFM = 16, EI=rep(1635.656,35), days=365)
  })
  # Consumption matrix is defined by columns. 
  # i.e Consumption for the i-th individual is in the i-th column
  expect_error({
    child_weight(age=c(5,4), sex=c("female", "male"), FM=c(2.7,3), 
                 FFM = c(16,14), EI=matrix(rep(1635.656,365*2),nrow = 2))
  })
})

test_that("Checking child_weight warning",{
  expect_warning({
    child_weight(age=5, sex="female", FM=2.7, FFM = 16, days=365*20)
  })
})

test_that("Checking child_weight message",{
  
  expect_message({
    child_weight(age=5, sex="female", FM=2.7, FFM = 16)
  })
  
  expect_message({
    child_weight(age=5, sex="female", FM=2.7, FFM = 16, EI=rep(1635.656,365))
  })
})

test_that("Checking child_weight results",{
  # Check that if there are no specifications children are estimated with reference formulas
  expect_equal({child_weight(age = 5, sex = "female")},
               {child_weight(age = 5, sex = "female", 
                             FM=child_reference_FFMandFM(5, "female")$FM,
                             FFM=child_reference_FFMandFM(5, "female")$FFM)})
  
  expect_equal({child_weight(age = 5, sex = "female")},
               {child_weight(age = 5, sex = "female", 
                             EI=child_reference_EI(5, "female",
                                                   FM=child_reference_FFMandFM(5, "female")$FM,
                                                   FFM=child_reference_FFMandFM(5, "female")$FFM,
                                                   days = 365))})
  
  # Check reference children have similar values to the reference after 10 years
  # of running Hall model 
  # Age 5
 #Female fat mass
   expect_lt({
    FM      <- child_reference_FFMandFM(15,"female")$FM
    FMmodel <- child_weight(age=5, sex="female", days=365*10)$Fat_Mass[365*10]
    abs(FM-FMmodel)/FM
  }, 0.1)
   
  
   #Female fat free mass
   expect_lt({
     FFM      <- child_reference_FFMandFM(15,"female")$FFM
     FFMmodel <- child_weight(age=5, sex="female", days=365*10)$Fat_Free_Mass[365*10]
     abs(FFM-FFMmodel)/FFM
   }, 0.1)
   
  # Male fat mass
  expect_lt({
    FM      <- child_reference_FFMandFM(15,"male")$FM
    FMmodel <- child_weight(age=5, sex="male", days=365*10)$Fat_Mass[365*10]
    abs(FM-FMmodel)/FM
  }, 0.1)
  
  # Male fat free mass
  expect_lt({
    FFM      <- child_reference_FFMandFM(15,"male")$FFM
    FFMmodel <- child_weight(age=5, sex="male", days=365*10)$Fat_Free_Mass[365*10]
    abs(FFM-FFMmodel)/FFM
  }, 0.1)
  
  #Age 8
  #Female fat mass
  expect_lt({
    FM      <- child_reference_FFMandFM(18,"female")$FM
    FMmodel <- child_weight(age=8, sex="female", days=365*10)$Fat_Mass[365*10]
    abs(FM-FMmodel)/FM
  }, 0.1)
  
  #Female fat free mass
  expect_lt({
    FFM      <- child_reference_FFMandFM(18,"female")$FFM
    FFMmodel <- child_weight(age=8, sex="female", days=365*10)$Fat_Free_Mass[365*10]
    abs(FFM-FFMmodel)/FFM
  }, 0.1)
  
  # Male fat mass
  expect_lt({
    FM      <- child_reference_FFMandFM(18,"male")$FM
    FMmodel <- child_weight(age=8, sex="male", days=365*10)$Fat_Mass[365*10]
    abs(FM-FMmodel)/FM
  }, 0.1)
  
  # Male fat free mass
  expect_lt({
    FFM      <- child_reference_FFMandFM(18,"male")$FFM
    FFMmodel <- child_weight(age=8, sex="male", days=365*10)$Fat_Free_Mass[365*10]
    abs(FFM-FFMmodel)/FFM
  }, 0.1)
  
})

  