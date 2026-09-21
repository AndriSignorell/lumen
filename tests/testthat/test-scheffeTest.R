library(testthat)
library(lumen)

set.seed(1)
df <- data.frame(
  y = c(rnorm(20, 5), rnorm(20, 7), rnorm(20, 9)),
  g = factor(rep(c("A","B","C"), each = 20))
)
fit <- aov(y ~ g, data = df)

test_that("scheffeTest: aov method returns PostHocTest", {
  expect_s3_class(scheffeTest(fit), "PostHocTest")
})

test_that("scheffeTest: formula method works", {
  res <- scheffeTest(y ~ g, data = df)
  expect_s3_class(res, "PostHocTest")
})

test_that("scheffeTest: p-values in [0,1]", {
  res  <- scheffeTest(fit)
  pvals <- res[[1]][, "pval"]
  expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))
})

test_that("scheffeTest: A vs C significant", {
  res   <- scheffeTest(fit)
  pvals <- res[[1]][, "pval"]
  ac_p  <- pvals[grep("A-C|C-A", names(pvals))]
  expect_lt(min(ac_p), 0.05)
})

test_that("scheffeTest: CI contains diff for equal groups", {
  set.seed(42)
  df2 <- data.frame(y = rnorm(60), g = factor(rep(1:3, each=20)))
  res <- scheffeTest(aov(y ~ g, data = df2))
  # CI for all pairs should contain 0
  diffs <- res[[1]][, "diff"]
  lcis  <- res[[1]][, "lci"]
  ucis  <- res[[1]][, "uci"]
  expect_true(all(lcis <= 0 & 0 <= ucis, na.rm = TRUE))
})

test_that("scheffeTest: wider CI with higher conf.level", {
  res95 <- scheffeTest(fit, conf.level = 0.95)[[1]]
  res99 <- scheffeTest(fit, conf.level = 0.99)[[1]]
  w95   <- mean(res95[,"uci"] - res95[,"lci"], na.rm = TRUE)
  w99   <- mean(res99[,"uci"] - res99[,"lci"], na.rm = TRUE)
  expect_gt(w99, w95)
})


# -- added --------------------------------------------------------------------

wb <- aov(breaks ~ tension, data = warpbreaks)

test_that("scheffeTest: pairwise results equal postHocTest(method = 'scheffe')", {
  s <- scheffeTest(wb)$tension
  p <- postHocTest(wb, method = "scheffe")$tension
  expect_equal(unclass(s), unclass(p), ignore_attr = TRUE)
  expect_identical(rownames(s), rownames(p))
})

test_that("scheffeTest: result attributes", {
  r <- scheffeTest(wb)
  expect_identical(attr(r, "method"), "Scheff\u00e9 Test")
  expect_identical(attr(r, "conf.level"), 0.95)
  expect_false(attr(r, "ordered"))
  expect_equal(attr(r, "orig.call"), wb$call)
  expect_output(print(r), "Scheff")
})

test_that("scheffeTest: custom contrast by hand", {
  # L against the mean of M and H
  cc <- matrix(c(1, -0.5, -0.5), ncol = 1)
  r <- scheffeTest(wb, contrasts = cc)$tension
  m <- as.vector(tapply(warpbreaks$breaks, warpbreaks$tension, mean))
  n <- 18
  mse <- sum(wb$residuals^2) / wb$df.residual
  psi <- sum(cc * m)
  se <- sqrt(mse * sum(cc^2 / n))
  expect_identical(rownames(r), "L-M,H")
  expect_equal(unname(r[, "diff"]), psi)
  expect_equal(unname(r[, "pval"]),
               pf(psi^2 / (se^2 * 2), 2, wb$df.residual, lower.tail = FALSE))
  expect_equal(unname(r[, "uci"] - r[, "diff"]),
               sqrt(2 * qf(0.95, 2, wb$df.residual)) * se)
})

test_that("scheffeTest: several contrasts at once, and conf.level = NA", {
  cc <- cbind(c(1, -0.5, -0.5), c(0, 1, -1))
  r <- scheffeTest(wb, contrasts = cc)$tension
  expect_identical(rownames(r), c("L-M,H", "M-H"))
  expect_equal(r["M-H", "pval"], scheffeTest(wb)$tension["H-M", "pval"])

  p <- scheffeTest(wb, contrasts = cc, conf.level = NA)$tension
  expect_identical(colnames(p), c("diff", "pval"))
  expect_equal(p[, "pval"], r[, "pval"])
})

test_that("scheffeTest: conf.level = NA gives the lower triangle of p-values", {
  r <- scheffeTest(wb, conf.level = NA)$tension
  full <- scheffeTest(wb)$tension
  expect_identical(dimnames(r), list(c("M", "H"), c("L", "M")))
  expect_equal(c(r["M", "L"], r["H", "L"], r["H", "M"]), unname(full[, "pval"]))
  expect_true(is.na(r["M", "M"]))
})

test_that("scheffeTest: contrast checks", {
  expect_error(scheffeTest(wb, contrasts = matrix(c(1, 1, -1), ncol = 1)),
               "sum to zero")
  fit2 <- aov(breaks ~ wool + tension, data = warpbreaks)
  # regression: a 3-level contrast was recycled over the 2 wool levels
  expect_error(scheffeTest(fit2, contrasts = matrix(c(1, -0.5, -0.5), ncol = 1)),
               "select the term with 'which'")
  r <- scheffeTest(fit2, which = "tension",
                   contrasts = matrix(c(1, -0.5, -0.5), ncol = 1))
  expect_named(r, "tension")
})

test_that("scheffeTest: dfgrp is the term's own df in a two-factor model", {
  fit2 <- aov(breaks ~ wool + tension, data = warpbreaks)
  r <- scheffeTest(fit2)
  mse <- sum(fit2$residuals^2) / fit2$df.residual
  # wool: 2 levels -> df1 = 1, i.e. the plain t-test
  w <- r$wool
  se <- sqrt(mse * 2 / 27)
  expect_equal(unname(w[, "pval"]),
               2 * pt(abs(w[, "diff"]) / se, fit2$df.residual, lower.tail = FALSE),
               ignore_attr = TRUE)
})

test_that("scheffeTest: which", {
  fit2 <- aov(breaks ~ wool * tension, data = warpbreaks)
  expect_named(scheffeTest(fit2), c("wool", "tension", "wool:tension"))
  expect_named(scheffeTest(fit2, which = "tension"), "tension")
  expect_error(scheffeTest(fit2, which = "foo"), "specified no factors")
  expect_warning(r <- scheffeTest(fit2, which = c("wool", "foo")), "non-factors")
  expect_named(r, "wool")
})

test_that("scheffeTest: default method", {
  d <- with(warpbreaks, scheffeTest(breaks, tension))
  expect_equal(unclass(d$g), unclass(scheffeTest(wb)$tension), ignore_attr = TRUE)
  # regression: numeric group codes were treated as a covariate
  n <- with(warpbreaks, scheffeTest(breaks, as.integer(tension)))
  expect_equal(unname(n$g[, "pval"]), unname(d$g[, "pval"]))
  expect_error(scheffeTest(1:10), "'g' is missing")
  expect_error(scheffeTest(1:10, rep(1:2, 4)), "same length")
})

test_that("scheffeTest: formula method evaluates subset in the data", {
  s <- scheffeTest(breaks ~ tension, data = warpbreaks, subset = wool == "A")
  a <- scheffeTest(aov(breaks ~ tension, data = warpbreaks[warpbreaks$wool == "A", ]))
  expect_equal(s$tension, a$tension, ignore_attr = TRUE)
  # the stored call names the data instead of embedding it
  expect_lt(nchar(paste(deparse(attr(s, "orig.call")), collapse = "")), 200)
})

test_that("scheffeTest: formula method passes na.action and '...'", {
  d <- warpbreaks
  d$breaks[1:3] <- NA
  r <- scheffeTest(breaks ~ tension, data = d, na.action = na.omit,
                   conf.level = 0.9)
  expect_identical(attr(r, "conf.level"), 0.9)
  expect_equal(r$tension,
               scheffeTest(aov(breaks ~ tension, data = d[-(1:3), ]),
                           conf.level = 0.9)$tension, ignore_attr = TRUE)
})

test_that("scheffeTest: .contrasts builds all pairwise contrasts", {
  m <- lumen:::.contrasts(c("a", "b", "c"))
  expect_identical(names(m), c("b-a", "c-a", "c-b"))
  expect_identical(rownames(m), c("a", "b", "c"))
  expect_equal(colSums(m), c("b-a" = 0, "c-a" = 0, "c-b" = 0))
  expect_equal(m[["c-a"]], c(-1, 0, 1))
})
