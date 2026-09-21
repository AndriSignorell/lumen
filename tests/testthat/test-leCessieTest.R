library(testthat)
library(lumen)

set.seed(111)
x1  <- factor(sample(1:3, 50, replace = TRUE))
x2  <- rnorm(50)
obs <- sample(c(0, 1), 50, replace = TRUE)
fit <- glm(obs ~ x1 + x2, family = binomial)
f   <- fitted(fit)
X   <- model.matrix(fit)   # full design matrix, incl. intercept

test_that("leCessieTest: returns htest / LeCessieTest", {
  res <- leCessieTest(x = f, obs = obs, X = X)
  expect_s3_class(res, "htest")
  expect_s3_class(res, "LeCessieTest")
})

test_that("leCessieTest: statistic named Z", {
  res <- leCessieTest(x = f, obs = obs, X = X)
  expect_named(res$statistic, "Z")
})

test_that("leCessieTest: p.value in [0, 1]", {
  res <- leCessieTest(x = f, obs = obs, X = X)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("leCessieTest: sse, expected, sd are positive scalars", {
  res <- leCessieTest(x = f, obs = obs, X = X)
  expect_gt(res$sse,      0)
  expect_gt(res$expected, 0)
  expect_gt(res$sd,       0)
})

test_that("leCessieTest: well-specified model gives large p", {
  set.seed(1)
  n   <- 500
  x   <- rnorm(n)
  eta <- -1 + 2 * x
  y   <- rbinom(n, 1, plogis(eta))
  g   <- glm(y ~ x, family = binomial)
  res <- leCessieTest(
    x   = fitted(g), obs = y,
    X   = model.matrix(g)   # full design matrix, incl. intercept
  )
  expect_gt(res$p.value, 0.05)
})

test_that("leCessieTest: misspecified model gives small p", {
  set.seed(42)
  n   <- 1000
  x1  <- rnorm(n)
  x2  <- rnorm(n)
  # true model has interaction, fitted model misses it
  eta <- -1 + 2 * x1 + 2 * x2 + 5 * x1 * x2
  y   <- rbinom(n, 1, plogis(eta))
  g   <- glm(y ~ x1 + x2, family = binomial)   # missing interaction
  res <- leCessieTest(
    x   = fitted(g), obs = y,
    X   = model.matrix(g)   # full design matrix, incl. intercept
  )
  expect_lt(res$p.value, 0.05)
})

test_that("leCessieTest: input validation - length mismatch", {
  expect_error(leCessieTest(f[-1], obs, X), "same length")
})

test_that("leCessieTest: input validation - fit out of [0,1]", {
  bad <- f; bad[1] <- 1.5
  expect_error(leCessieTest(bad, obs, X), "probabilities")
})

test_that("leCessieTest: input validation - non-binary obs", {
  bad <- obs; bad[1] <- 2
  expect_error(leCessieTest(f, bad, X), "binary")
})

test_that("leCessieTest: input validation - X row mismatch", {
  expect_error(leCessieTest(f, obs, X[-1, ]), "same number of rows")
})

test_that("leCessieTest: print method runs without error", {
  res <- leCessieTest(x = f, obs = obs, X = X)
  expect_output(print(res))
})


test_that("leCessieTest.glm: matches the default method with the full design matrix", {
  res_default <- leCessieTest(x = f, obs = obs, X = X)
  res_glm     <- leCessieTest(fit)

  expect_equal(unname(res_glm$statistic), unname(res_default$statistic))
  expect_equal(res_glm$p.value, res_default$p.value)
})

test_that("leCessieTest.glm: data.name reflects the model formula", {
  res <- leCessieTest(fit)
  expect_equal(res$data.name, "obs ~ x1 + x2")
})

test_that("leCessieTest: omitting the intercept column warns and changes the result", {
  # regression test: earlier documentation instructed passing
  # model.matrix(fit)[, -1] (no intercept), which silently gives a
  # materially different, incorrect result
  Xno <- model.matrix(fit)[, -1, drop = FALSE]

  expect_warning(
    res_no <- leCessieTest(f, obs, Xno),
    "intercept"
  )
  res_full <- leCessieTest(f, obs, X)

  expect_false(isTRUE(all.equal(unname(res_no$statistic),
                                unname(res_full$statistic))))
})

test_that("leCessieTest.glm: rejects a non-binomial glm", {
  g <- glm(x2 ~ x1, family = gaussian)
  expect_error(leCessieTest(g), "binomial")
})


# -- added --------------------------------------------------------------------

set.seed(8)
lc <- data.frame(x = rnorm(300), z = runif(300))
lc$y <- rbinom(300, 1, plogis(-0.2 + 0.9 * lc$x))
lfit <- glm(y ~ x, family = binomial, data = lc)

# direct matrix form: Var(SSE) = d'(W - W X (X'WX)^-1 X'W) d, d = 1 - 2p
lcRef <- function(p, y, X) {
  W <- diag(p * (1 - p))
  d <- 1 - 2 * p
  V <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  (sum((y - p)^2) - sum(p * (1 - p))) / sqrt(drop(t(d) %*% V %*% d))
}

test_that("statistic equals the closed matrix form", {
  r <- leCessieTest(lfit)
  ref <- lcRef(fitted(lfit), lc$y, model.matrix(lfit))
  expect_equal(unname(r$statistic), ref)
  expect_equal(r$p.value, 2 * pnorm(-abs(ref)))
  expect_equal(r$sse, sum((lc$y - fitted(lfit))^2))
  expect_equal(r$expected, sum(fitted(lfit) * (1 - fitted(lfit))))
})

test_that("closed form also with a factor and an interaction", {
  d <- lc
  d$g <- gl(3, 100)
  f <- glm(y ~ x * g, family = binomial, data = d)
  expect_equal(unname(leCessieTest(f)$statistic),
               lcRef(fitted(f), d$y, model.matrix(f)))
})

test_that("size under a correctly specified model", {
  set.seed(4)
  p <- replicate(800, {
    x <- rnorm(200)
    y <- rbinom(200, 1, plogis(0.3 + x))
    leCessieTest(glm(y ~ x, family = binomial))$p.value
  })
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.025)
})

test_that("glm with na.exclude", {
  d <- lc
  d$x[c(3, 50, 99)] <- NA
  f1 <- glm(y ~ x, family = binomial, data = d, na.action = na.exclude)
  f2 <- glm(y ~ x, family = binomial, data = d)
  # regression: fitted() padded the excluded rows with NA
  expect_equal(leCessieTest(f1)$statistic, leCessieTest(f2)$statistic)
})

test_that("factor, logical and quasibinomial responses", {
  d <- lc
  d$yf <- factor(ifelse(d$y == 1, "yes", "no"))
  d$yl <- d$y == 1
  ref <- leCessieTest(lfit)$statistic
  expect_equal(leCessieTest(glm(yf ~ x, binomial, d))$statistic, ref)
  expect_equal(leCessieTest(glm(yl ~ x, binomial, d))$statistic, ref)
  expect_equal(leCessieTest(glm(y ~ x, quasibinomial, d))$statistic, ref)
})

test_that("glm method rejects unsupported models", {
  d <- lc
  d$n1 <- d$y; d$n0 <- 1 - d$y
  expect_error(leCessieTest(glm(cbind(n1, n0) ~ x, binomial, d)),
               "matrix responses")
  expect_error(leCessieTest(glm(y ~ x, binomial, d, weights = rep(1:3, length.out = 300))),
               "weighted")
  expect_error(leCessieTest(glm(y ~ x, binomial, d, weights = rep(2, 300))),
               "weighted")
})

test_that("default method: further checks", {
  p <- fitted(lfit); X <- model.matrix(lfit)
  expect_error(leCessieTest(as.character(p), lc$y, X), "numeric")
  expect_error(leCessieTest(replace(p, 1, NA), lc$y, X), "'x' must not contain")
  expect_error(leCessieTest(p, replace(lc$y, 1, NA), X), "'obs' must not contain")
  expect_error(leCessieTest(p, lc$y, as.data.frame(X)), "numeric matrix")
  expect_error(leCessieTest(p, lc$y, replace(X, 1, NA)), "'X' must not contain")
  # p = 0.5 everywhere: d = 0, the standard deviation vanishes
  expect_error(leCessieTest(rep(0.5, 300), lc$y, X), "standard deviation is zero")
})

test_that("print shows SSE and its expectation", {
  expect_output(print(leCessieTest(lfit)), "Sum of squared errors")
  expect_invisible(print(leCessieTest(lfit)))
})
