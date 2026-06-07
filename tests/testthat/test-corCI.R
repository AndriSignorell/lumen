library(testthat)
library(lumen)

tol <- 1e-4

test_that("corCI: output has names cor/lci/uci", {
  res <- corCI(0.5, n = 30)
  expect_named(res, c("cor", "lci", "uci"))
})

test_that("corCI: cor element equals input rho", {
  expect_equal(unname(corCI(0.5, n = 30)["cor"]), 0.5)
})

test_that("corCI: CI contains rho", {
  res <- corCI(0.5, n = 30)
  expect_true(res["lci"] <= 0.5 && res["uci"] >= 0.5)
})

test_that("corCI: CI width decreases with larger n", {
  w30  <- diff(unname(corCI(0.5, n = 30)[c("lci","uci")]))
  w300 <- diff(unname(corCI(0.5, n = 300)[c("lci","uci")]))
  expect_true(w300 < w30)
})

test_that("corCI: CI widens with higher conf.level", {
  w95 <- diff(unname(corCI(0.5, n = 50, conf.level = 0.95)[c("lci","uci")]))
  w99 <- diff(unname(corCI(0.5, n = 50, conf.level = 0.99)[c("lci","uci")]))
  expect_true(w99 > w95)
})

test_that("corCI: rho=0 gives symmetric CI", {
  res <- corCI(0, n = 50)
  expect_equal(unname(res["lci"]), -unname(res["uci"]), tolerance = tol)
})

test_that("corCI: alternative='greater' sets lci < rho", {
  res <- corCI(0.5, n = 50, alternative = "greater")
  expect_true(res["lci"] < 0.5)
  expect_equal(unname(res["uci"]), 1, tolerance = tol)
})

test_that("corCI: alternative='less' sets uci > rho", {
  res <- corCI(0.5, n = 50, alternative = "less")
  expect_true(res["uci"] > 0.5)
  expect_equal(unname(res["lci"]), -1, tolerance = tol)
})

test_that("corCI: n < 3 throws error", {
  expect_error(corCI(0.5, n = 2))
})

test_that("corCI: |rho| > 1 throws error", {
  expect_error(corCI(1.1, n = 30))
})

test_that("corCI: known CI matches manual Fisher z computation", {
  # manual: z = atanh(0.5), se = 1/sqrt(30-3), z_alpha = qnorm(0.975)
  rho <- 0.5; n <- 30
  z   <- atanh(rho)
  se  <- 1 / sqrt(n - 3)
  lci <- tanh(z - qnorm(0.975) * se)
  uci <- tanh(z + qnorm(0.975) * se)
  res <- corCI(rho, n = n)
  expect_equal(unname(res["lci"]), lci, tolerance = tol)
  expect_equal(unname(res["uci"]), uci, tolerance = tol)
})
