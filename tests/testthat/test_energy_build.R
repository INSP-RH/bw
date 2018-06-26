context("Energy build function")

test_that("Checking energy_build  errors",{
  
  # Check that the columns in energy are the same as the length of time vector
  expect_error({
    energy_build(energy = c(1220, 2600), time = c(0, 5*365, 10*365))  
  })
  
  #  Check that time is a vector
  expect_error({
    energy_build(energy = matrix(c(1220, 2600, 1500, 2500), byrow = TRUE, nrow = 2),
                 time =   matrix(c(0, 365, 0, 365), byrow = TRUE, nrow = 2))  
  })
  
  # Check that time values are nonnegative
  expect_error({
    energy_build(energy = matrix(c(1220, 2600, 1500, 2500), byrow = TRUE, nrow = 2),
                 time =   c(0, -365))
  })
  
  # Check that time values are integer
  expect_error({
    energy_build(energy = matrix(c(1220, 2600, 1500, 2500), byrow = TRUE, nrow = 2),
                 time =   c(0, 365.4))
  })
  
  # Check that the first time value is equal to zero
  expect_error({
    energy_build(energy = matrix(c(1220, 2600, 1500, 2500), byrow = TRUE, nrow = 2),
                 time =   c(1, 365))
  })
  
  # Check that time values are increasing
  expect_error({
    energy_build(energy = matrix(c(1220, 2600, 2400, 1500, 2500, 2000), byrow = TRUE, nrow = 2),
                 time =   c(0, 365,12))
  })
  
  # Check there are no repeated values
  expect_error({
    energy_build(energy = matrix(c(1220, 2600, 2400, 1500, 2500, 2000), byrow = TRUE, nrow = 2),
                 time =   c(0, 365,365))
  })
  
  # Check interpolation name
  expect_error({
    energy_build(energy = matrix(c(1220, 2600, 2400, 1500, 2500, 2000), byrow = TRUE, nrow = 2),
                 time =   c(0, 365, 365*2), interpolation = "Unknown")
  })
})

test_that("Checking energy_build different interpolation inputs.",{
  # Check thar energy interpolations work
  # Linear
  expect_silent({
    energy_build(energy = c(1220, 2600), time = c(0, 5*365), interpolation = "Linear") 
  })
  
  # Exponential
  expect_silent({
    energy_build(energy = c(1220, 2600), time = c(0, 5*365), interpolation = "Exponential") 
  })
  
  # Logarithmic
  expect_silent({
    energy_build(energy = c(1220, 2600), time = c(0, 5*365), interpolation = "Logarithmic") 
  })
  
  # Stepwise_L
  expect_silent({
    energy_build(energy = c(1220, 2600), time = c(0, 5*365), interpolation = "Stepwise_L") 
  })
  
  # Stepwise_R
  expect_silent({
    energy_build(energy = c(1220, 2600), time = c(0, 5*365), interpolation = "Stepwise_R") 
  })
  
  # Brownian
  expect_silent({
    energy_build(energy = c(1220, 2600), time = c(0, 5*365), interpolation = "Brownian") 
  })
  
  
})
