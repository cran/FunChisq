require(testthat)
require(FunChisq)
context("Test the tests in FunChisq package")

test_that("Test the exact functional test", {
  x1 <- matrix(c(12, 26, 18, 0, 8, 12), nrow=2, byrow = TRUE)
  expect_equal(exact.functional.test(x1), 0.042556228)
  expect_equal(exact.functional.test(t(x1)), 0.027271578)
  
  x2 <- matrix(c(0,0,0,0,0,0,0,0,0), nrow=3, byrow = TRUE)
  expect_identical(exact.functional.test(x2), 1)
  
  x3 <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  expect_equal(exact.functional.test(x3), 0.002997003)
  expect_equal(exact.functional.test(t(x3)), 0.006549005)
  
  x4 <- matrix(rep(10,25), nrow=5)
  expect_identical(exact.functional.test(x4), 1)
  
  x5 <- matrix(c(4,0,0,0,4,0,0,0,4), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x5), 0.00017316)
  
  x6 <- matrix(c(2,0,0,2), nrow=2, byrow = TRUE)
  expect_true(abs(exact.functional.test(x6) - fisher.test(x6)$p.value) <= 1e-6)
  
  x7 <- matrix(c(2,2,2,2), nrow=2, byrow = TRUE)
  expect_equal(exact.functional.test(x7), fisher.test(x7)$p.value)
  
  x8 <- matrix(c(0,10,15,20,5,0,25,0,0), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x8), as.numeric(fun.chisq.test(x8)$p.value))
})

test_that("Test the functional chi-square test", {
  x <- matrix(c(20,0,20,0,20,0,5,0,5), 3)
  expect_equal(as.numeric(fun.chisq.test(x)$p.value), 8.582e-15)
  expect_equal(as.numeric(fun.chisq.test(t(x))$p.value), 3.638517e-13)
  
})

test_that("Test the comparative functional chi-square test", {
  x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  y <- t(x)
  z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
  data <- list(x,y,z)
  expect_equal(as.numeric(cp.fun.chisq.test(data)$p.value), 0.0001876212-8.1134e-12)
})
