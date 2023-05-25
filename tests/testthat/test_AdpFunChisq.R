#Testing Adapted Functional Chi-squared test
#Created by: Sajal Kumar
#Date : 5/17/2021

library(testthat)
library(FunChisq)
library(dqrng)

context("Testing AdpFunChisq")

test_that("Testing the adapted functional chi-squared test", {

  # set seed
  dqset.seed(123)

  # Test 1 : Portrait functional table
  x = matrix(c(12, 1, 0,
               0, 15, 2,
               20, 1, 1,
               2, 18, 1,
               0, 0, 15), nrow=5, ncol=3, byrow=TRUE)

  xtoy = fun.chisq.test(x, method="adapted")
  ytox = fun.chisq.test(t(x), method="adapted")

  # better stats
  expect_equivalent(signif(xtoy$p.value, digits = 2), 1.1e-23)
  expect_equivalent(round(xtoy$statistic, digits = 2), 127.21)
  expect_equivalent(round(xtoy$estimate, digits = 2), 0.86)

  # worse stats
  expect_equivalent(signif(ytox$p.value, digits = 2), 5e-08)
  expect_equivalent(round(ytox$statistic, digits = 2), 49.53)
  expect_equivalent(round(ytox$estimate, digits = 2), 0.63)


  # Test 2 : Landscape functional table
  x = matrix(c(20, 1, 0, 1,
               0, 15, 1, 2,
               18, 1, 3, 4), nrow=3, ncol=4, byrow=TRUE)

  xtoy = fun.chisq.test(x, method="adapted")
  ytox = fun.chisq.test(t(x), method="adapted")

  # better stats
  expect_equivalent(signif(xtoy$p.value, digits = 2), 4.7e-10)
  expect_equivalent(round(xtoy$statistic, digits = 2), 53.48)
  expect_equivalent(round(xtoy$estimate, digits = 2), 0.73)

  # worse stats
  expect_equivalent(signif(ytox$p.value, digits = 2), 2.6e-08)
  expect_equivalent(round(ytox$statistic, digits = 2), 46.26)
  expect_equivalent(round(ytox$estimate, digits = 2), 0.6)


  # Test 3 : Independent table
  x = matrix(c(15, 16, 12, 13,
               11, 15, 12, 17,
               18, 14, 13, 14), nrow=3, ncol=4, byrow=TRUE)

  xtoy = fun.chisq.test(x, method="adapted")
  ytox = fun.chisq.test(t(x), method="adapted")

  # better stats
  expect_equivalent(signif(xtoy$p.value, digits = 2), 0.81)
  expect_equivalent(round(xtoy$statistic, digits = 2), 3)
  expect_equivalent(round(xtoy$estimate, digits = 2), 0.09)

  # worse stats
  expect_equivalent(signif(ytox$p.value, digits = 2), 0.89)
  expect_equivalent(round(ytox$statistic, digits = 2), 2.31)
  expect_equivalent(round(ytox$estimate, digits = 2), 0.08)

})
