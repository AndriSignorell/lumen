library(testthat)
library(lumen)

set.seed(1); x <- rnorm(100)

test_that("quantileCI: returns matrix with est/lci/uci", {
  res <- quantileCI(x)
  expect_true(is.matrix(res))
  expect_true(all(c("est","lci","uci") %in% colnames(res)))
})

test_that("quantileCI: number of rows = length(probs)", {
  res <- quantileCI(x, probs = c(0.25, 0.5, 0.75))
  expect_equal(nrow(res), 3L)
})

test_that("quantileCI: est equals quantile(x, probs)", {
  probs <- c(0.25, 0.5, 0.75)
  res   <- quantileCI(x, probs = probs)
  expect_equal(unname(res[,"est"]), unname(quantile(x, probs)), tolerance = 1e-10)
})

test_that("quantileCI: lci <= est <= uci", {
  res <- quantileCI(x, probs = c(0.1, 0.5, 0.9))
  expect_true(all(res[,"lci"] <= res[,"est"]))
  expect_true(all(res[,"est"] <= res[,"uci"]))
})

test_that("quantileCI: wider CI with higher conf.level", {
  r95 <- quantileCI(x, probs = 0.5, conf.level = 0.95)
  r80 <- quantileCI(x, probs = 0.5, conf.level = 0.80)
  expect_gte(r95[1,"uci"] - r95[1,"lci"], r80[1,"uci"] - r80[1,"lci"])
})

test_that("quantileCI: na.rm=TRUE handles NAs", {
  xna <- c(x, NA)
  expect_true(is.matrix(quantileCI(xna, na.rm = TRUE)))
})

test_that("quantileCI: sides='left' gives uci=Inf", {
  res <- quantileCI(x, probs = 0.5, sides = "left")
  expect_equal(unname(res[1,"uci"]), Inf)
})

test_that("quantileCI: sides='right' gives lci=-Inf", {
  res <- quantileCI(x, probs = 0.5, sides = "right")
  expect_equal(unname(res[1,"lci"]), -Inf)
})
