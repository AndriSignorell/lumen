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


# quantileCI() ---------------------------------------------------------------

# coverage of [X_(l), X_(u)] for the prob-quantile, l/u in 0..n+1
covOS <- function(l, u, n, p) pbinom(u - 1, n, p) - pbinom(l - 1, n, p)

test_that("exact two-sided: order statistics and achieved coverage", {
  x <- c(8.1, 3.4, 5.9, 12.0, 7.7, 1.5, 9.8, 4.4, 6.6, 10.3,
         2.8, 11.1, 5.2, 7.0, 9.1, 3.9, 6.1, 8.8, 4.9, 10.9)
  r <- quantileCI(x, probs = 0.5)
  s <- sort(x)
  l <- match(r[, "lci"], s)
  u <- match(r[, "uci"], s)
  
  expect_identical(colnames(r), c("est", "lci", "uci"))
  expect_equal(r[, "est"], median(x), ignore_attr = TRUE)
  expect_equal(attr(r, "conf.level"), covOS(l, u, 20, 0.5))
  expect_gte(attr(r, "conf.level"), 0.95)
})

test_that("exact: smallest coverage not below the level (n = 100, median)", {
  # candidates [40, 60] (0.954) and the symmetric [40, 61] (0.965):
  # the rule takes the one closer to 0.95
  r <- quantileCI(1:100, probs = 0.5)
  expect_equal(unname(r[1, c("lci", "uci")]), c(40, 60))
  expect_equal(attr(r, "conf.level"), covOS(40, 60, 100, 0.5))
})

test_that("exact: several probs, one coverage each", {
  r <- quantileCI(1:100, probs = c(0.25, 0.75, 0.9))
  expect_identical(dim(r), c(3L, 3L))
  expect_length(attr(r, "conf.level"), 3)
  expect_true(all(attr(r, "conf.level") >= 0.95))
  expect_true(all(r[, "lci"] <= r[, "est"] & r[, "est"] <= r[, "uci"]))
})

test_that("exact one-sided bounds and their coverage", {
  n <- 20; x <- seq_len(n)
  
  lft <- quantileCI(x, probs = 0.5, sides = "left")
  l <- lft[, "lci"]
  expect_identical(unname(lft[, "uci"]), Inf)
  expect_equal(attr(lft, "conf.level"), 1 - pbinom(l - 1, n, 0.5))
  expect_gte(attr(lft, "conf.level"), 0.95)
  
  rgt <- quantileCI(x, probs = 0.5, sides = "right")
  u <- rgt[, "uci"]
  expect_identical(unname(rgt[, "lci"]), -Inf)
  # X_(u) >= q  <=>  Bin(n, p) <= u - 1
  expect_equal(attr(rgt, "conf.level"), pbinom(u - 1, n, 0.5))
  expect_gte(attr(rgt, "conf.level"), 0.95)
  
  # symmetric problem, symmetric answer (regression: 'right' was one
  # order statistic short and fell below the level)
  expect_equal(unname(u), n + 1 - unname(l))
  expect_equal(attr(rgt, "conf.level"), attr(lft, "conf.level"))
})

test_that("exact one-sided bounds cover at the reported level", {
  set.seed(1)
  n <- 20
  hit <- replicate(4000, {
    z <- rnorm(n)
    c(quantileCI(z, probs = 0.5, sides = "right")[, "uci"] >= 0,
      quantileCI(z, probs = 0.5, sides = "left")[, "lci"] <= 0)
  })
  expect_gt(mean(hit[1, ]), 0.96)
  expect_gt(mean(hit[2, ]), 0.96)
})

test_that("exact: too few observations give open bounds, not NA/errors", {
  r <- quantileCI(1:3, probs = 0.5)
  expect_equal(unname(r[1, c("lci", "uci")]), c(-Inf, Inf))
  expect_equal(attr(r, "conf.level"), 1)
  
  r <- quantileCI(1:3, probs = 0.5, sides = "left")
  expect_equal(unname(r[1, c("lci", "uci")]), c(-Inf, Inf))
  r <- quantileCI(1:3, probs = 0.5, sides = "right")
  expect_equal(unname(r[1, c("lci", "uci")]), c(-Inf, Inf))
})

test_that("exact: probs 0 and 1 give no NA", {
  r <- quantileCI(1:50, probs = c(0, 1))
  expect_false(anyNA(r))
})

test_that("boot method: bounds from boot.ci(), sides", {
  x <- mtcars$mpg
  set.seed(2)
  r <- quantileCI(x, probs = 0.5, method = "boot", type = "perc", R = 199)
  set.seed(2)
  b <- boot::boot(x, function(x, d) quantile(x[d], 0.5), R = 199)
  expect_equal(unname(r[1, c("lci", "uci")]),
               unname(boot::boot.ci(b, type = "perc")$percent[4:5]))
  expect_null(attr(r, "conf.level"))
  
  set.seed(2)
  n <- quantileCI(x, probs = 0.5, method = "boot", type = "norm", R = 199)
  set.seed(2)
  b <- boot::boot(x, function(x, d) quantile(x[d], 0.5), R = 199)
  expect_equal(unname(n[1, c("lci", "uci")]),
               unname(boot::boot.ci(b, type = "norm")$normal[2:3]))
  
  set.seed(2)
  l <- quantileCI(x, probs = c(0.25, 0.5), method = "boot", sides = "left",
                  type = "perc", R = 199)
  expect_true(all(l[, "uci"] == Inf))
  set.seed(2)
  rr <- quantileCI(x, probs = c(0.25, 0.5), method = "boot", sides = "right",
                   type = "perc", R = 199)
  expect_true(all(rr[, "lci"] == -Inf))
})

test_that("boot method passes invalid bootstrap arguments to the validator", {
  expect_error(quantileCI(1:20, method = "boot", type = "all"),
               "'type' must be one of")
})

test_that("input checks", {
  expect_error(quantileCI(c(1, NA, 3)), "missing values")
  expect_equal(quantileCI(c(1:20, NA), probs = 0.5, na.rm = TRUE),
               quantileCI(1:20, probs = 0.5))
  expect_error(quantileCI(letters), "must be numeric")
  expect_error(quantileCI(1), "at least two")
  expect_error(quantileCI(1:10, probs = 1.2), "'probs'")
  expect_error(quantileCI(1:10, probs = NA), "'probs'")
  expect_error(quantileCI(1:10, sides = "foo"))
  expect_error(quantileCI(1:10, method = "foo"))
})
