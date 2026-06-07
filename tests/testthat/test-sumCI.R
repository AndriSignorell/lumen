library(testthat)
library(lumen)

test_that("sumCI: returns named vector sum/lci/uci", {
  m <- matrix(c(10, 8, 12, 20, 17, 23), nrow = 2,
              dimnames = list(NULL, c("est","lci","uci")))
  res <- sumCI(m)
  expect_named(res, c("sum","lci","uci"))
})

test_that("sumCI: sum = sum of estimates", {
  m <- matrix(c(10, 20, 8, 17, 12, 23), nrow = 2)
  res <- sumCI(m)
  expect_equal(unname(res["sum"]), 30)
})

test_that("sumCI: lci <= sum <= uci", {
  m <- matrix(c(10, 20, 8, 17, 12, 23), nrow = 2)
  res <- sumCI(m)
  expect_true(res["lci"] <= res["sum"] && res["sum"] <= res["uci"])
})

test_that("sumCI: CI width = sqrt(sum of squared half-widths)", {
  m <- matrix(c(10, 20, 6, 16, 14, 24), nrow = 2)
  hw1 <- (14 - 6) / 2; hw2 <- (24 - 16) / 2
  expected_hw <- sqrt(hw1^2 + hw2^2)
  res <- sumCI(m)
  actual_hw <- (res["uci"] - res["lci"]) / 2
  expect_equal(unname(actual_hw), expected_hw, tolerance = 1e-10)
})

test_that("sumCI: single row gives trivial result", {
  m <- matrix(c(5, 3, 7), nrow = 1)
  res <- sumCI(m)
  expect_equal(unname(res["sum"]), 5)
  expect_equal(unname(res["lci"]), 3)
  expect_equal(unname(res["uci"]), 7)
})

test_that("sumCI: non-matrix input throws error", {
  expect_error(sumCI(c(1, 2, 3)))
})

test_that("sumCI: wrong number of columns throws error", {
  expect_error(sumCI(matrix(1:4, nrow = 2)))
})
