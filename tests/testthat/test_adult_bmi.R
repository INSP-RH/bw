context("Adult bmi function")

test_that("Checking adult_bmi errors",{
  W <-adult_weight(bw = c(76, 58, 65), ht = c(1.73, 1.64, 1.65), 
                   age = c(36,21, 56), sex = c("male", "female", "female"), PAL = c(1.5,1.7, 1.5),
                   EI=matrix(c(rep(2287,365), rep(2159, 365), rep(1762,365)), 
                             nrow =3, byrow = T))
  
  # Check that no results are given if time values are not present in weight object
  expect_error({
    adult_bmi(W, days = 370)
  })
  
  expect_error({
    adult_bmi(W, days = 3.9)
  })
  
  # Check dimension of groups is either one or the same length as number of individuals
  # More groups
  expect_error({
    adult_bmi(W, days = c(0, 10,40, 50), group = c("Sedentary", "Athletic", "Sedentary", "Sedentary"))
  })
  
  # Less groups (!=1)
  expect_error({
    adult_bmi(W, days = c(0, 10,40, 50), group = c("Sedentary", "Athletic"))
  })
})

# Check warning
test_that("Checking adult_bmi warnings",{
  W <-adult_weight(bw = c(76, 58, 65), ht = c(1.73, 1.64, 1.65), 
                   age = c(36,21, 56), sex = c("male", "female", "female"), PAL = c(1.5,1.7, 1.5),
                   EI=matrix(c(rep(2287,365), rep(2159, 365), rep(1762,365)), 
                             nrow =3, byrow = T))
  
  #Check confidence level is plausible
  expect_warning({
    adult_bmi(W, days = c(0, 10,40, 50), confidence = 1.1)
  })
})

# Check message
test_that("Checking adult_bmi messages",{
  W <-adult_weight(bw = c(76, 58, 65), ht = c(1.73, 1.64, 1.65), 
                   age = c(36,21, 56), sex = c("male", "female", "female"), PAL = c(1.5,1.7, 1.5),
                   EI=matrix(c(rep(2287,365), rep(2159, 365), rep(1762,365)), 
                             nrow =3, byrow = T))
  
  #Check confidence level is plausible
  expect_message({
    adult_bmi(W, days = seq(0,364,by=5))
  })
})

# Check results at baseline

test_that("Check bmi results at baseline",{
  bw  <- c(76, 58, 65, 88, 37, 82)
  ht  <- c(1.73, 1.64, 1.65, 1.70, 1.5, 1.8)
  
  W   <- adult_weight(bw = bw, ht = ht, age = c(36, 21, 56, 44, 28, 63), 
                   sex = c("male", "female", "female", "male", "female", "male"),
                   PAL = c(1.5,1.7, 1.5, 1.5, 1.6, 1.7),
                   EI=matrix(c(rep(2287,365), rep(2159, 365), 
                               rep(1762,365), rep(2600,365),
                               rep(1610,365),rep(2900,365)), 
                             nrow =5, byrow = T))
  result <- adult_bmi(W, days = 0)
  
  bmi <- bw/(ht^2) 
  under  <- length(which(bmi<18.5))/length(bw)
  norm   <- length(which(bmi>=18.5 & bmi<25))/length(bw)
  over   <- length(which(bmi>=25 & bmi<30))/length(bw)
  obese  <- length(which(bmi>=30))/length(bw)
  
  # Check underweight category
  expect_equal({
    result$Mean[which(result$BMI_Category=="Underweight")]
  }, under)
  
  # Check normal category
  expect_equal({
    result$Mean[which(result$BMI_Category=="Normal")]
  }, norm)
  
  # Check pre-obese category
  expect_equal({
    result$Mean[which(result$BMI_Category=="Pre-Obese")]
  }, over)
  
  # Check obese category
  expect_equal({
    result$Mean[which(result$BMI_Category=="Obese")]
  }, obese)
})
