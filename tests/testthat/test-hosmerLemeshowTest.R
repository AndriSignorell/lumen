library(testthat)
library(lumen)

# n = 300 rather than 50: with 10 groups of 5, nearly every call warned that
# expected counts are below 5, which is the test's own sparsity, not a finding
set.seed(111)
x1  <- factor(sample(1:3, 300, replace = TRUE))
x2  <- rnorm(300)
obs <- sample(c(0, 1), 300, replace = TRUE)
fit <- glm(obs ~ x1 + x2, family = binomial)
f   <- fitted(fit)

test_that("hosmerLemeshowTest: returns htest / HosmerLemeshowTest (type C)", {
  res <- hosmerLemeshowTest(x = f, obs = obs, type = "C")
  expect_s3_class(res, "htest")
  expect_s3_class(res, "HosmerLemeshowTest")
})

test_that("hosmerLemeshowTest: returns htest / HosmerLemeshowTest (type H)", {
  # fixed [0, 1] bins: fitted values near 0.5 leave most bins empty or
  # sparse, so warnings are expected here (tested in their own block below)
  res <- suppressWarnings(hosmerLemeshowTest(x = f, obs = obs, type = "H"))
  expect_s3_class(res, "htest")
  expect_s3_class(res, "HosmerLemeshowTest")
})

test_that("hosmerLemeshowTest: statistic named X-squared", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_named(res$statistic, "X-squared")
})

test_that("hosmerLemeshowTest: parameter named df equals nGroups - 2", {
  res <- hosmerLemeshowTest(x = f, obs = obs, nGroups = 10)
  expect_named(res$parameter, "df")
  expect_equal(unname(res$parameter), res$nGroups - 2L)
})

test_that("hosmerLemeshowTest: p.value in [0, 1]", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("hosmerLemeshowTest: observed and expected are matrices", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_true(is.matrix(res$observed))
  expect_true(is.matrix(res$expected))
  expect_equal(colnames(res$observed), c("0s", "1s"))
  expect_equal(colnames(res$expected), c("0s", "1s"))
})

test_that("hosmerLemeshowTest: well-specified model gives large p", {
  set.seed(1)
  n   <- 500
  x   <- rnorm(n)
  eta <- -1 + 2 * x
  y   <- rbinom(n, 1, plogis(eta))
  g   <- glm(y ~ x, family = binomial)
  # the tail deciles of a well-fitting model have few expected events
  res <- suppressWarnings(hosmerLemeshowTest(x = fitted(g), obs = y))
  expect_gt(res$p.value, 0.05)
})

test_that("hosmerLemeshowTest: nGroups respected", {
  res <- hosmerLemeshowTest(x = f, obs = obs, nGroups = 5)
  expect_lte(res$nGroups, 5L)
  expect_gte(res$nGroups, 3L)
})

test_that("hosmerLemeshowTest: input validation - length mismatch", {
  expect_error(hosmerLemeshowTest(f[-1], obs), "same length")
})

test_that("hosmerLemeshowTest: input validation - fit out of [0,1]", {
  bad <- f; bad[1] <- -0.1
  expect_error(hosmerLemeshowTest(bad, obs), "probabilities")
})

test_that("hosmerLemeshowTest: input validation - non-binary obs", {
  bad <- obs; bad[1] <- 3
  expect_error(hosmerLemeshowTest(f, bad), "binary")
})

test_that("hosmerLemeshowTest: input validation - nGroups < 3", {
  expect_error(hosmerLemeshowTest(f, obs, nGroups = 2), "nGroups")
})

test_that("hosmerLemeshowTest: print method runs without error", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_output(print(res))
})

test_that("hosmerLemeshowTest: print with details runs without error", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_output(print(res, details = TRUE))
})


test_that("hosmerLemeshowTest.glm: matches the default method on the same data", {
  res_default <- hosmerLemeshowTest(x = f, obs = obs, type = "C")
  res_glm     <- hosmerLemeshowTest(fit)

  expect_equal(unname(res_glm$statistic), unname(res_default$statistic))
  expect_equal(res_glm$p.value, res_default$p.value)
})

test_that("hosmerLemeshowTest.glm: data.name reflects the model formula", {
  res <- hosmerLemeshowTest(fit)
  expect_equal(res$data.name, "obs ~ x1 + x2")
})

test_that("hosmerLemeshowTest.glm: rejects a non-binomial glm", {
  g <- glm(x2 ~ x1, family = gaussian)
  expect_error(hosmerLemeshowTest(g), "binomial")
})

test_that("hosmerLemeshowTest: type = 'H' uses fixed [0,1] bins and can drop empty groups", {
  set.seed(4)
  fit_narrow <- runif(80, 0.35, 0.65)   # occupies only a few of 10 [0,1] bins
  obs_narrow <- rbinom(80, 1, fit_narrow)

  # two warnings: dropped groups, and sparse expected counts in the rest
  expect_warning(
    expect_warning(
      res <- hosmerLemeshowTest(fit_narrow, obs_narrow, type = "H"),
      "empty group"
    ),
    "expected counts"
  )
  expect_lt(res$nGroups, 10L)
  expect_equal(unname(res$parameter), res$nGroups - 2L)
})


# -- added --------------------------------------------------------------------

set.seed(7)
hl <- data.frame(x = rnorm(400), z = runif(400))
hl$y <- rbinom(400, 1, plogis(-0.3 + 1.2 * hl$x))
hfit <- glm(y ~ x, family = binomial, data = hl)

# test_that("type C is identical to ResourceSelection::hoslem.test", {
#   skip_if_not_installed("ResourceSelection")
#   for (g in c(5, 10, 12)) {
#     a <- suppressWarnings(hosmerLemeshowTest(hfit, nGroups = g))
#     b <- ResourceSelection::hoslem.test(hfit$y, fitted(hfit), g = g)
#     expect_equal(unname(a$statistic), unname(b$statistic), info = g)
#     expect_equal(unname(a$parameter), unname(b$parameter), info = g)
#     expect_equal(a$p.value, b$p.value, info = g)
#   }
# })

test_that("statistic by hand from observed and expected counts", {
  a <- suppressWarnings(hosmerLemeshowTest(hfit))
  expect_equal(unname(a$statistic),
               sum((a$observed - a$expected)^2 / a$expected))
  expect_equal(sum(a$observed), 400)
  expect_equal(sum(a$expected), 400)
  # glm with intercept: expected 1s sum to observed 1s
  expect_equal(sum(a$expected[, "1s"]), sum(hl$y))
})

test_that("type H: groups are fixed deciles of [0, 1]", {
  a <- suppressWarnings(hosmerLemeshowTest(hfit, type = "H"))
  p <- fitted(hfit)
  grp <- cut(p, seq(0, 1, 0.1), include.lowest = TRUE)
  expect_equal(unname(a$observed[, "1s"]),
               as.vector(tapply(hl$y, grp, sum))[table(grp) > 0])
  expect_identical(a$method, "Hosmer-Lemeshow H statistic")
})

test_that("glm with na.exclude", {
  d <- hl
  d$x[c(3, 50, 99)] <- NA
  f1 <- glm(y ~ x, family = binomial, data = d, na.action = na.exclude)
  f2 <- glm(y ~ x, family = binomial, data = d)
  # regression: fitted() padded the excluded rows with NA
  expect_equal(suppressWarnings(hosmerLemeshowTest(f1))$statistic,
               suppressWarnings(hosmerLemeshowTest(f2))$statistic)
})

test_that("factor and logical responses", {
  d <- hl
  d$yf <- factor(ifelse(d$y == 1, "yes", "no"))
  d$yl <- d$y == 1
  ref <- suppressWarnings(hosmerLemeshowTest(hfit))$statistic
  expect_equal(suppressWarnings(hosmerLemeshowTest(
    glm(yf ~ x, binomial, d)))$statistic, ref)
  expect_equal(suppressWarnings(hosmerLemeshowTest(
    glm(yl ~ x, binomial, d)))$statistic, ref)
})

test_that("glm method rejects unsupported models", {
  d <- hl
  d$n1 <- d$y; d$n0 <- 1 - d$y
  expect_error(hosmerLemeshowTest(glm(cbind(n1, n0) ~ x, binomial, d)),
               "matrix responses")
  expect_error(hosmerLemeshowTest(glm(y ~ x, binomial, d, weights = rep(1:3, length.out = 400))),
               "weighted")
  expect_equal(suppressWarnings(hosmerLemeshowTest(
    glm(y ~ x, binomial, d, weights = rep(1, 400))))$statistic,
    suppressWarnings(hosmerLemeshowTest(hfit))$statistic)
})

test_that("few distinct fitted values: fewer groups, warning or error", {
  p <- rep(c(0.2, 0.4, 0.6, 0.8), each = 50)
  set.seed(1)
  o <- rbinom(200, 1, p)
  w <- character()
  r <- withCallingHandlers(hosmerLemeshowTest(p, o), warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("distinct groups", w)))
  # quantile() interpolates a break at 0.5 between the tied values 0.4 and
  # 0.6; the resulting empty group is dropped, 3 groups remain
  expect_true(any(grepl("empty group", w)))
  expect_identical(r$nGroups, 3L)
  expect_equal(unname(r$parameter), 1L)
  expect_error(hosmerLemeshowTest(rep(c(0.3, 0.6), 50), rbinom(100, 1, 0.5)),
               "at least 3 groups")
})

test_that("type H: fewer than 3 non-empty groups is an error", {
  p <- rep(c(0.31, 0.35), 50)
  expect_error(suppressWarnings(hosmerLemeshowTest(p, rbinom(100, 1, 0.3),
                                                   type = "H")),
               "fewer than 3")
})

test_that("argument checks", {
  p <- runif(50); o <- rbinom(50, 1, 0.5)
  for (bad in list(Inf, NA, 3.5, "10", c(5, 6)))
    expect_error(hosmerLemeshowTest(p, o, nGroups = bad), "'nGroups'",
                 info = format(bad))
  expect_error(hosmerLemeshowTest(c(p[-1], NA), o), "missing values")
  expect_error(hosmerLemeshowTest(p, o, type = "X"))
  expect_error(hosmerLemeshowTest(as.character(p), o), "numeric")
})

test_that("print shows the group table with details = TRUE", {
  a <- suppressWarnings(hosmerLemeshowTest(hfit))
  expect_output(print(a), "Number of groups: 10")
  expect_output(print(a, details = TRUE), "Observed vs Expected")
  # print() would otherwise write to the console during the test run
  expect_output(expect_invisible(print(a)))
})
