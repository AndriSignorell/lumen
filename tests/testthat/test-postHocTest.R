library(testthat)
library(lumen)

set.seed(1)
df <- data.frame(
  y = c(rnorm(20, 5), rnorm(20, 7), rnorm(20, 9)),
  g = factor(rep(c("A","B","C"), each = 20))
)
fit <- aov(y ~ g, data = df)

test_that("postHocTest: returns PostHocTest", {
  expect_s3_class(postHocTest(fit), "PostHocTest")
})

test_that("postHocTest: result is a list", {
  expect_true(is.list(postHocTest(fit)))
})

test_that("postHocTest: hsd method works", {
  res <- postHocTest(fit, method = "hsd")
  expect_s3_class(res, "PostHocTest")
})

test_that("postHocTest: bonferroni method works", {
  res <- postHocTest(fit, method = "bonferroni")
  expect_s3_class(res, "PostHocTest")
})

test_that("postHocTest: lsd method works", {
  res <- postHocTest(fit, method = "lsd")
  expect_s3_class(res, "PostHocTest")
})

test_that("postHocTest: scheffe method works", {
  res <- postHocTest(fit, method = "scheffe")
  expect_s3_class(res, "PostHocTest")
})

test_that("postHocTest: p-values in [0,1]", {
  res <- postHocTest(fit, method = "hsd")
  pvals <- res[[1]][, "pval"]
  expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))
})

test_that("postHocTest: A vs C significant (large difference)", {
  res <- postHocTest(fit, method = "hsd")
  # A-C comparison should have small p-value given mean diff of ~4
  pvals <- res[[1]][, "pval"]
  ac_p  <- pvals[grep("A-C|C-A", names(pvals))]
  expect_lt(min(ac_p), 0.05)
})

test_that("postHocTest stops on covariates", {
  fit <- aov(temperature ~ driver + delivery_min, data = bedrock::Pizza)
  expect_error(postHocTest(fit), "covariate")
  expect_error(postHocTest(fit, which = "driver"), "covariate")
})


# postHocTest() --------------------------------------------------------------

fit1 <- aov(breaks ~ tension, data = warpbreaks)
fit2 <- aov(breaks ~ wool * tension, data = warpbreaks)

# the interval must touch 0 exactly at conf.level = 1 - pval; checks that
# width and p-value of a method carry the same multiplicity correction
dualityOK <- function(fit, method, term = 1L) {
  p <- postHocTest(fit, method = method)[[term]][, "pval"]
  ok <- p > 1e-6 & p < 0.99
  all(vapply(which(ok), function(i) {
    r <- postHocTest(fit, method = method, conf.level = 1 - p[i])[[term]]
    abs(abs(r[i, "diff"]) - (r[i, "uci"] - r[i, "lci"]) / 2) <
      1e-6 * max(1, abs(r[i, "diff"]))
  }, NA))
}

test_that("result structure", {
  r <- postHocTest(fit1)
  expect_s3_class(r, "PostHocTest")
  expect_named(r, "tension")
  expect_identical(colnames(r$tension), c("diff", "lci", "uci", "pval"))
  expect_identical(rownames(r$tension), c("M-L", "H-L", "H-M"))
  expect_identical(attr(r, "conf.level"), 0.95)
  expect_false(attr(r, "ordered"))
  expect_identical(attr(r, "method"), "Tukey HSD")
  expect_equal(attr(r, "orig.call"), fit1$call)
})

test_that("hsd equals TukeyHSD, including interaction terms", {
  r <- postHocTest(fit2, method = "hsd")
  t <- TukeyHSD(fit2)
  for (nm in names(t))
    expect_equal(unclass(r[[nm]]), unclass(t[[nm]]), ignore_attr = TRUE,
                 info = nm)
  expect_identical(rownames(r[["wool:tension"]]), rownames(t[["wool:tension"]]))
})

test_that("hsd ordered = TRUE equals TukeyHSD(ordered = TRUE)", {
  r <- postHocTest(fit1, method = "hsd", ordered = TRUE)
  t <- TukeyHSD(fit1, ordered = TRUE)
  expect_equal(unclass(r$tension), unclass(t$tension), ignore_attr = TRUE)
  expect_identical(rownames(r$tension), rownames(t$tension))
  expect_true(all(r$tension[, "diff"] > 0))
  expect_true(attr(r, "ordered"))
})

test_that("ordered is ignored for lsd / bonferroni / scheffe", {
  for (m in c("lsd", "bonferroni", "scheffe"))
    expect_identical(unclass(postHocTest(fit1, method = m, ordered = TRUE)$tension),
                     unclass(postHocTest(fit1, method = m)$tension), info = m)
})

test_that("lsd equals pooled pairwise t-tests without adjustment", {
  r <- postHocTest(fit1, method = "lsd")$tension
  pw <- pairwise.t.test(warpbreaks$breaks, warpbreaks$tension,
                        p.adjust.method = "none")$p.value
  expect_equal(unname(r[, "pval"]), pw[lower.tri(pw, diag = TRUE)])
  se <- sqrt(sum(fit1$residuals^2) / fit1$df.residual * 2 / 18)
  expect_equal(unname(r[, "uci"] - r[, "diff"]),
               rep(qt(0.975, fit1$df.residual) * se, 3))
  expect_identical(attr(postHocTest(fit1, method = "lsd"), "method"), "Fisher LSD")
})

test_that("bonferroni p-values equal pairwise.t.test(p.adjust = 'bonferroni')", {
  r <- postHocTest(fit1, method = "bonf")$tension
  pw <- pairwise.t.test(warpbreaks$breaks, warpbreaks$tension,
                        p.adjust.method = "bonferroni")$p.value
  expect_equal(unname(r[, "pval"]), pw[lower.tri(pw, diag = TRUE)])
})

test_that("bonferroni interval: two-sided, alpha / (2m)", {
  r <- postHocTest(fit1, method = "bonferroni")$tension
  se <- sqrt(sum(fit1$residuals^2) / fit1$df.residual * 2 / 18)
  expect_equal(unname(r[, "uci"] - r[, "diff"]),
               rep(qt(1 - 0.05 / (2 * 3), fit1$df.residual) * se, 3))
})

test_that("interval and p-value agree for every single-step method", {
  for (m in c("lsd", "bonferroni", "hsd", "scheffe", "newmankeuls", "duncan"))
    expect_true(dualityOK(fit1, m), info = m)
  expect_true(dualityOK(fit2, "bonferroni", "wool:tension"))
})

test_that("scheffe by hand", {
  r <- postHocTest(fit1, method = "scheffe")$tension
  mse <- sum(fit1$residuals^2) / fit1$df.residual
  se <- sqrt(mse * 2 / 18)
  est <- r[, "diff"] / se
  expect_equal(unname(r[, "pval"]), pf(unname(est)^2 / 2, 2, fit1$df.residual,
                                       lower.tail = FALSE))
  expect_equal(unname(r[, "uci"] - r[, "diff"]),
               rep(sqrt(2 * qf(0.95, 2, fit1$df.residual)) * se, 3))
})

test_that("newmankeuls / duncan reduce to lsd for adjacent means", {
  lsd <- postHocTest(fit1, method = "lsd")$tension
  for (m in c("newmankeuls", "duncan")) {
    r <- postHocTest(fit1, method = m)$tension
    # means L > M > H: M-L and H-M are adjacent, H-L spans 3 means
    expect_equal(r[c("M-L", "H-M"), "pval"], lsd[c("M-L", "H-M"), "pval"],
                 info = m)
    expect_gt(r["H-L", "pval"], lsd["H-L", "pval"])
  }
  # Newman-Keuls on the full span equals Tukey
  expect_equal(postHocTest(fit1, method = "newmankeuls")$tension["H-L", "pval"],
               postHocTest(fit1, method = "hsd")$tension["H-L", "pval"])
})

test_that("conf.level = NA returns the lower triangle of p-values", {
  r <- postHocTest(fit1, method = "hsd", conf.level = NA)
  ci <- postHocTest(fit1, method = "hsd")$tension
  m <- r$tension
  expect_identical(dimnames(m), list(c("M", "H"), c("L", "M")))
  expect_equal(c(m["M", "L"], m["H", "L"], m["H", "M"]), unname(ci[, "pval"]))
  expect_true(is.na(m["M", "M"]))
  expect_true(is.na(attr(r, "conf.level")))
})

test_that("which selects terms", {
  r <- postHocTest(fit2, which = "tension")
  expect_named(r, "tension")
  expect_equal(r$tension, postHocTest(fit2)$tension)
  expect_named(postHocTest(fit2, which = 2:3), c("tension", "wool:tension"))
  
  expect_error(postHocTest(fit2, which = "foo"), "specified no factors")
  expect_warning(r <- postHocTest(fit2, which = c("wool", "foo")),
                 "non-factors")
  expect_named(r, "wool")
})

test_that("method is matched", {
  expect_error(postHocTest(fit1, method = "foo"))
  expect_identical(attr(postHocTest(fit1, method = "sch"), "method"),
                   "Scheff\u00e9")
})

test_that("print methods run for both branches", {
  expect_output(print(postHocTest(fit1)), "95% family-wise confidence level")
  expect_output(print(postHocTest(fit1)), "Signif. codes", fixed = TRUE)
  expect_output(print(postHocTest(fit1, ordered = TRUE)), "have been ordered")
  expect_output(print(postHocTest(fit1, conf.level = NA)), "Tukey HSD")
  # wrapped in expect_output(): print() would write to the console during the run
  expect_output(expect_invisible(print(postHocTest(fit1, conf.level = NA))))
})

test_that("plot method runs", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(postHocTest(fit2)))
})

# -- tables ------------------------------------------------------------------

tab <- as.table(rbind(A = c(20, 30, 50),
                      B = c(35, 25, 40),
                      C = c(60, 20, 20)))

test_that("table: pairwise chi-square p-values on row pairs", {
  r <- postHocTest(tab)
  expect_s3_class(r, "PostHocTest")
  m <- r[[1L]]
  expect_identical(dimnames(m), list(c("B", "C"), c("A", "B")))
  ref <- function(i, j) chisq.test(tab[c(i, j), ])$p.value
  expect_equal(m["B", "A"], ref("A", "B"))
  expect_equal(m["C", "A"], ref("A", "C"))
  expect_equal(m["C", "B"], ref("B", "C"))
  expect_true(is.na(m["B", "B"]))
  expect_true(is.na(attr(r, "conf.level")))
})

test_that("table: p.adjust methods", {
  raw <- postHocTest(tab)[[1L]]
  p <- raw[lower.tri(raw, diag = TRUE)]
  for (m in c("holm", "bonferroni", "BH")) {
    adj <- postHocTest(tab, method = m)[[1L]]
    expect_equal(adj[lower.tri(adj, diag = TRUE)], p.adjust(p, m), info = m)
  }
})

test_that("table and matrix method agree; conf.level warns", {
  expect_equal(unclass(postHocTest(unclass(tab)))[[1L]],
               unclass(postHocTest(tab))[[1L]])
  expect_warning(postHocTest(tab, conf.level = 0.9), "not supported")
  # regression: the table method forwarded conf.level, so every call warned
  expect_no_warning(postHocTest(tab))
  expect_no_warning(postHocTest(tab, conf.level = NA))
  expect_output(print(postHocTest(tab, method = "holm")), "holm")
})


