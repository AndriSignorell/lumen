library(testthat)
library(lumen)

set.seed(1)
x <- rnorm(50)

test_that("madCI: returns named vector est/lci/uci", {
  res <- madCI(x)
  expect_true(all(c("est","lci","uci") %in% names(res)))
})

test_that("madCI: est equals mad(x)", {
  res <- madCI(x)
  expect_equal(unname(res["est"]), mad(x), tolerance = 1e-10)
})

test_that("madCI: lci <= est <= uci", {
  res <- madCI(x)
  expect_true(res["lci"] <= res["est"] && res["est"] <= res["uci"])
})

test_that("madCI: lci >= 0", {
  expect_gte(unname(madCI(abs(x))["lci"]), 0)
})

test_that("madCI: wider CI with higher conf.level", {
  w95 <- diff(unname(madCI(x, conf.level = 0.95)[c("lci","uci")]))
  w99 <- diff(unname(madCI(x, conf.level = 0.99)[c("lci","uci")]))
  expect_gt(w99, w95)
})

test_that("madCI: sides='left' gives uci=Inf", {
  res <- madCI(x, sides = "left")
  expect_equal(unname(res["uci"]), Inf)
})

test_that("madCI: sides='right' gives lci=0 or -Inf", {
  res <- madCI(x, sides = "right")
  expect_true(res["lci"] <= 0)
})

test_that("madCI: na.rm=TRUE handles NAs", {
  xna <- c(x, NA, NA)
  expect_equal(unname(madCI(xna, na.rm = TRUE)["est"]),
               mad(x), tolerance = 1e-10)
})

test_that("madCI: empty vector throws error", {
  expect_error(madCI(numeric(0)))
})
