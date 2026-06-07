library(testthat)
library(lumen)

# NOTE: poissonCI() has a bug where .poissonCI_engine() has a parameter
# stdEst without a default value, causing NA results when called via
# .recycleApply(). The fix is: stdEst = NULL in .poissonCI_engine().
# These tests assume the fix is deployed (poissonCI.R updated).

tol <- 1e-4

test_that("poissonCI: default returns named vector with est/lci/uci", {
  res <- poissonCI(10, 1)
  expect_type(res, "double")
  expect_true(all(c("est", "lci", "uci") %in% names(res)))
})

test_that("poissonCI: est = x/n", {
  expect_equal(unname(poissonCI(10, 2)["est"]), 5)
  expect_equal(unname(poissonCI(0,  1)["est"]), 0)
})

test_that("poissonCI: lci <= est <= uci", {
  for (m in c("exact", "score", "wald", "byar")) {
    res <- poissonCI(10, 1, method = m)
    expect_false(is.na(res["lci"]), label = paste(m, "not NA"))
    expect_lte(unname(res["lci"]), unname(res["est"]), label = paste(m, "lci<=est"))
    expect_lte(unname(res["est"]), unname(res["uci"]), label = paste(m, "est<=uci"))
  }
})

test_that("poissonCI: lci >= 0 for x=0", {
  for (m in c("exact", "score", "byar")) {
    expect_gte(unname(poissonCI(0, 1, method = m)["lci"]), 0, label = m)
  }
})

test_that("poissonCI: exact matches poisson.test", {
  ref <- poisson.test(10, 1)$conf.int
  res <- poissonCI(10, 1, method = "exact")
  expect_equal(unname(res["lci"]), ref[1], tolerance = tol)
  expect_equal(unname(res["uci"]), ref[2], tolerance = tol)
})

test_that("poissonCI: wider CI with higher conf.level", {
  w95 <- diff(unname(poissonCI(10, 1, conf.level = 0.95)[c("lci","uci")]))
  w80 <- diff(unname(poissonCI(10, 1, conf.level = 0.80)[c("lci","uci")]))
  expect_gt(w95, w80)
})

test_that("poissonCI: sides=left gives uci=Inf", {
  expect_equal(unname(poissonCI(10, 1, sides = "left")["uci"]), Inf)
})

test_that("poissonCI: sides=right gives lci=0", {
  expect_equal(unname(poissonCI(10, 1, sides = "right")["lci"]), 0)
})

test_that("poissonCI: n>1 scales est correctly", {
  expect_equal(unname(poissonCI(10, 2)["est"]),
               unname(poissonCI(10, 1)["est"]) / 2, tolerance = tol)
})

test_that("poissonCI: multiple methods returns data.frame with est/lci/uci", {
  res <- poissonCI(10, 1, method = c("exact", "score", "wald", "byar"))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("est", "lci", "uci") %in% names(res)))
  expect_equal(nrow(res), 4L)
  expect_true(all(res[["lci"]] <= res[["est"]]))
  expect_true(all(res[["est"]] <= res[["uci"]]))
})
