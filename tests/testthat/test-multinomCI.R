library(testthat)
library(lumen)

x <- c(A = 20, B = 15, C = 25)

test_that("multinomCI: returns matrix with 3 columns", {
  res <- multinomCI(x)
  expect_true(is.matrix(res))
  expect_equal(ncol(res), 3L)
})

test_that("multinomCI: column names est/lci/uci", {
  expect_equal(colnames(multinomCI(x)), c("est", "lci", "uci"))
})

test_that("multinomCI: est sums to 1", {
  res <- multinomCI(x)
  expect_equal(sum(res[,"est"]), 1, tolerance = 1e-10)
})

test_that("multinomCI: lci <= est <= uci for all categories", {
  res <- multinomCI(x)
  expect_true(all(res[,"lci"] <= res[,"est"]))
  expect_true(all(res[,"est"] <= res[,"uci"]))
})

test_that("multinomCI: lci >= 0 and uci <= 1", {
  res <- multinomCI(x)
  expect_true(all(res[,"lci"] >= 0))
  expect_true(all(res[,"uci"] <= 1))
})

test_that("multinomCI: est = proportions", {
  res <- multinomCI(x)
  expect_equal(unname(res[,"est"]), unname(x / sum(x)), tolerance = 1e-10)
})

test_that("multinomCI: wider CI with higher conf.level", {
  w95 <- mean(multinomCI(x, conf.level = 0.95)[,"uci"] -
                multinomCI(x, conf.level = 0.95)[,"lci"])
  w99 <- mean(multinomCI(x, conf.level = 0.99)[,"uci"] -
                multinomCI(x, conf.level = 0.99)[,"lci"])
  expect_gt(w99, w95)
})

test_that("multinomCI: all methods return valid result", {
  methods <- c("sisonglaz","goodman","wald","waldcc","wilson")
  for (m in methods) {
    res <- multinomCI(x, method = m)
    expect_true(all(res[,"lci"] >= 0), label = paste(m, "lci>=0"))
    expect_true(all(res[,"uci"] <= 1), label = paste(m, "uci<=1"))
  }
})

test_that("multinomCI: sides='left' gives uci=1", {
  res <- multinomCI(x, sides = "left")
  expect_true(all(res[,"uci"] == 1))
})

test_that("multinomCI: sides='right' gives lci=0", {
  res <- multinomCI(x, sides = "right")
  expect_true(all(res[,"lci"] == 0))
})
