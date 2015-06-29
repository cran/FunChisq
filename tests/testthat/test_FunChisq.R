library(testthat)
library(FunChisq)
context("Testing the FunChisq package")

test_that("Testing the exact functional test", {
  x1 <- matrix(c(12, 26, 18, 0, 8, 12), nrow=2, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x1), 8), 0.042556227)
  expect_equal(signif(exact.functional.test(t(x1)), 8), 0.027271581)
  
  x2 <- matrix(c(0,0,0,0,0,0,0,0,0), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x2), 1)
  
  x3 <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  expect_equal(signif(exact.functional.test(x3), 8), 0.002997003)
  expect_equal(signif(exact.functional.test(t(x3)), 8), 0.0065490065)
  
  if(0) { # This test case causes hang on windows. To be fixed.
    x4 <- matrix(rep(10,25), nrow=5)
    expect_equal(exact.functional.test(x4), 1)
  }  
  
  x5 <- matrix(c(4,0,0,0,4,0,0,0,4), nrow=3, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x5), 8), 0.00017316017)
  
  x6 <- matrix(c(2,0,0,2), nrow=2, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x6), 8), 
               signif(as.numeric(stats::fisher.test(x6)$p.value), 8))
  
  x7 <- matrix(c(2,2,2,2), nrow=2, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x7), 8), 
               signif(as.numeric(stats::fisher.test(x7)$p.value), 8))
    
  x8 <- matrix(c(0,10,15,20,5,0,25,0,0), nrow=3, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x8), 8), 
               signif(as.numeric(fun.chisq.test(x8)$p.value), 8))
  
  x9 <- matrix(c(1,1,1,1,1,1,1,1,1), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x9), 1)
})

test_that("Testing the functional chi-square test", {
  x <- matrix(c(20,0,20,0,20,0,5,0,5), 3)
  expect_equal(signif(as.numeric(fun.chisq.test(x)$p.value), 8), 8.5822345e-15)
  expect_equal(signif(as.numeric(fun.chisq.test(t(x))$p.value), 8), 3.6385174e-13)

  expect_equal(signif(as.numeric(fun.chisq.test(x, method="normalized")$p.value), 8), 
               5.1061401e-128)
  expect_equal(signif(as.numeric(fun.chisq.test(t(x), method="normalized")$p.value), 8), 
               4.1897164e-101)
  
})

test_that("Testing the comparative functional chi-square test", {
  x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  y <- t(x)
  z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
  data <- list(x,y,z)
  # expect_equal(as.numeric(cp.fun.chisq.test(data)$p.value), 0.0001876212-8.1134e-12)
  expect_equal(signif(as.numeric(cp.fun.chisq.test(data)$p.value), 8), 0.00018762119)
  expect_equal(signif(as.numeric(cp.fun.chisq.test(data, method="normalized")$p.value), 8), 
               1.0052639e-07)
})
