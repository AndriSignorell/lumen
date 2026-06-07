library(testthat)
library(lumen)

# -----------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------

set.seed(1)
n    <- 30
x1   <- rnorm(n)
x2   <- rnorm(n)
y    <- 1 + 2 * x1 - x2 + rnorm(n)
df   <- data.frame(y = y, x1 = x1, x2 = x2)
fit  <- lm(y ~ x1 + x2, data = df)

set.seed(42)
x_alt  <- rep(c(-1, 1), 50)
err1   <- rnorm(100)
err2   <- as.numeric(stats::filter(err1, 0.9, method = "recursive"))
df_iid <- data.frame(y = 1 + x_alt + err1, x = x_alt)
df_ar1 <- data.frame(y = 1 + x_alt + err2, x = x_alt)

# -----------------------------------------------------------------------
# pan() called as in durbinWatsonTest
# -----------------------------------------------------------------------

test_that("pan: called as in durbinWatsonTest gives value in [0,1]", {
  res <- lm.fit(model.matrix(fit), model.response(model.frame(fit)))$residuals
  dw  <- sum(diff(res)^2) / sum(res^2)
  n_  <- length(res)
  X_  <- model.matrix(fit)
  A_  <- diag(c(1, rep(2, n_ - 2), 1))
  A_[abs(row(A_) - col(A_)) == 1] <- -1
  MA  <- (diag(n_) - X_ %*% chol2inv(qr.R(qr(X_))) %*% t(X_)) %*% A_
  ev  <- Re(eigen(MA, only.values = TRUE)$values)
  ev  <- ev[ev > 1e-10]
  p   <- pan(c(dw, ev), length(ev), 0, 15)
  expect_gte(p, 0)
  expect_lte(p, 1)
})

# -----------------------------------------------------------------------
# S3 dispatch
# -----------------------------------------------------------------------

test_that("durbinWatsonTest.lm: returns htest with named DW in (0, 4)", {
  res <- durbinWatsonTest(fit, exact = TRUE)
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "DW")
  expect_gt(unname(res$statistic), 0)
  expect_lt(unname(res$statistic), 4)
})

test_that("durbinWatsonTest.formula: same statistic and p.value as lm method", {
  res_lm  <- durbinWatsonTest(fit,                    exact = TRUE)
  res_fml <- durbinWatsonTest(y ~ x1 + x2, data = df, exact = TRUE)
  expect_equal(unname(res_lm$statistic), unname(res_fml$statistic), tolerance = 1e-10)
  expect_equal(res_lm$p.value,           res_fml$p.value,           tolerance = 1e-8)
})

test_that("durbinWatsonTest.numeric: returns htest with p in [0,1]", {
  e_t <- c(-32.33, -26.603, 2.215, -16.967, -1.148, -2.512, -1.967, 11.669,
           -0.513,  27.032, -4.422, 40.032, 23.577,  33.94, -2.787, -8.606,
            0.575,   6.848,-18.971, -29.063)
  res <- durbinWatsonTest(e_t, exact = TRUE)
  expect_s3_class(res, "htest")
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("durbinWatsonTest.default: unsupported class gives error", {
  expect_error(durbinWatsonTest(list(a = 1)), "no applicable method")
})

# -----------------------------------------------------------------------
# p-value algebra
# -----------------------------------------------------------------------

test_that("durbinWatsonTest: greater + less = 1 (exact)", {
  r_gt <- durbinWatsonTest(fit, alternative = "greater", exact = TRUE)
  r_lt <- durbinWatsonTest(fit, alternative = "less",    exact = TRUE)
  expect_equal(r_gt$p.value + r_lt$p.value, 1, tolerance = 1e-6)
})

test_that("durbinWatsonTest: two.sided <= 2 * min(greater, less)", {
  r_two <- durbinWatsonTest(fit, alternative = "two.sided", exact = TRUE)
  r_gt  <- durbinWatsonTest(fit, alternative = "greater",   exact = TRUE)
  r_lt  <- durbinWatsonTest(fit, alternative = "less",      exact = TRUE)
  expect_lte(r_two$p.value, 2 * min(r_gt$p.value, r_lt$p.value) + 1e-10)
})

test_that("durbinWatsonTest: exact and approx agree for n=30", {
  r_exact  <- durbinWatsonTest(fit, alternative = "greater", exact = TRUE)
  r_approx <- durbinWatsonTest(fit, alternative = "greater", exact = FALSE)
  expect_equal(r_exact$p.value, r_approx$p.value, tolerance = 0.05)
})

# -----------------------------------------------------------------------
# Power / DW interpretation
# -----------------------------------------------------------------------

test_that("durbinWatsonTest: iid errors -> DW near 2, large p (greater)", {
  res <- durbinWatsonTest(y ~ x, data = df_iid, alternative = "greater")
  expect_equal(unname(res$statistic), 2, tolerance = 0.5)
  expect_gt(res$p.value, 0.05)
})

test_that("durbinWatsonTest: AR(1) errors -> DW < 2, small p (greater)", {
  res <- durbinWatsonTest(y ~ x, data = df_ar1, alternative = "greater")
  expect_lt(unname(res$statistic), 2)
  expect_lt(res$p.value, 0.05)
})

# -----------------------------------------------------------------------
# orderBy
# -----------------------------------------------------------------------

test_that("durbinWatsonTest: orderBy changes DW statistic", {
  set.seed(3)
  df3  <- cbind(df, z = rnorm(n))
  res1 <- durbinWatsonTest(fit, exact = TRUE)
  res2 <- durbinWatsonTest(y ~ x1 + x2, data = df3, orderBy = ~z, exact = TRUE)
  expect_false(isTRUE(all.equal(unname(res1$statistic), unname(res2$statistic))))
})
