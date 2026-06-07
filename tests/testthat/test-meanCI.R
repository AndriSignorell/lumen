library(testthat)
library(lumen)

set.seed(1)
x <- rnorm(50, mean = 5, sd = 2)

test_that("meanCI: returns named vector est/lci/uci", {
  res <- meanCI(x)
  expect_true(all(c("est","lci","uci") %in% names(res)))
})

test_that("meanCI: est equals mean(x)", {
  expect_equal(unname(meanCI(x)["est"]), mean(x), tolerance = 1e-10)
})

test_that("meanCI: lci <= est <= uci", {
  res <- meanCI(x)
  expect_true(res["lci"] <= res["est"] && res["est"] <= res["uci"])
})

test_that("meanCI: matches t.test CI", {
  ref <- t.test(x)$conf.int
  res <- meanCI(x)
  expect_equal(unname(res["lci"]), ref[1], tolerance = 1e-6)
  expect_equal(unname(res["uci"]), ref[2], tolerance = 1e-6)
})

test_that("meanCI: wider CI with higher conf.level", {
  w95 <- diff(unname(meanCI(x, conf.level = 0.95)[c("lci","uci")]))
  w99 <- diff(unname(meanCI(x, conf.level = 0.99)[c("lci","uci")]))
  expect_gt(w99, w95)
})

test_that("meanCI: narrower CI with larger n", {
  set.seed(1)
  x50  <- rnorm(50)
  x500 <- rnorm(500)
  w50  <- diff(unname(meanCI(x50)[c("lci","uci")]))
  w500 <- diff(unname(meanCI(x500)[c("lci","uci")]))
  expect_gt(w50, w500)
})

test_that("meanCI: sides='left' gives uci=Inf", {
  expect_equal(unname(meanCI(x, sides = "left")["uci"]), Inf)
})

test_that("meanCI: sides='right' gives lci=-Inf", {
  expect_equal(unname(meanCI(x, sides = "right")["lci"]), -Inf)
})

test_that("meanCI: sd argument uses known variance (z-interval)", {
  # with known sd, CI based on normal (narrower for large n)
  res_known <- meanCI(x, sd = 2)
  res_est   <- meanCI(x)
  # both should contain the mean
  expect_true(res_known["lci"] <= mean(x) && mean(x) <= res_known["uci"])
})

test_that("meanCI: na.rm=TRUE handles NAs", {
  xna <- c(x, NA)
  expect_equal(unname(meanCI(xna, na.rm = TRUE)["est"]),
               mean(x), tolerance = 1e-10)
})
