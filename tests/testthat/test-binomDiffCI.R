library(testthat)
library(lumen)

# ===============================================================
# binomDiffCI TESTS
# ===============================================================

# the published reference values are rounded to four digits
tol <- 1e-3

methods_bdci <- c("wald", "wald-cc", "agresti-caffo", "newcombe-score",
                  "newcombe-score-cc", "mee-farrington-manning",
                  "miettinen-nurminen", "haldane", "jeffreys-perks",
                  "hauck-anderson")


# ---------------------------------------------------------------
# https://www.lexjansen.com/wuss/2016/127_Final_Paper_PDF.pdf, page 9

scenarios_sas <- list(
  s1 = list(x1 = 56, n1 = 70, x2 = 48, n2 = 80),
  s2 = list(x1 =  9, n1 = 10, x2 =  3, n2 = 10),
  s3 = list(x1 = 10, n1 = 10, x2 =  0, n2 = 20)
)

expected_sas <- list(

  wald = list(
    s1 = c(0.0575, 0.3425),
    s2 = c(0.2605, 0.9395),
    s3 = c(1.0000, 1.0000)
  ),

  `wald-cc` = list(
    s1 = c(0.0441, 0.3559),
    s2 = c(0.1605, 1.0000),
    s3 = c(0.9250, 1.0000)
  ),

  haldane = list(
    s1 = c(0.0535, 0.3351),
    s2 = c(0.1777, 0.8289),
    s3 = c(0.7482, 1.0000)
  ),

  `jeffreys-perks` = list(
    s1 = c(0.0531, 0.3355),
    s2 = c(0.1760, 0.8306),
    s3 = c(0.7431, 1.0000)
  ),

  `mee-farrington-manning` = list(
    s1 = c(0.0533, 0.3377),
    s2 = c(0.1821, 0.8370),
    s3 = c(0.7225, 1.0000)
  ),

  `miettinen-nurminen` = list(
    s1 = c(0.0528, 0.3382),
    s2 = c(0.1700, 0.8406),
    s3 = c(0.7156, 1.0000)
  ),

  # exact = list(
  #   s1 = c(0.0529, 0.3403),
  #   s2 = c(0.1393, 0.8836),
  #   s3 = c(0.6915, 1.0000)
  # ),

  `newcombe-score` = list(
    s1 = c(0.0524, 0.3339),
    s2 = c(0.1705, 0.8090),
    s3 = c(0.6791, 1.0000)
  ),

  `newcombe-score-cc` = list(
    s1 = c(0.0428, 0.3422),
    s2 = c(0.1013, 0.8387),
    s3 = c(0.6014, 1.0000)
  ),

  `hauck-anderson` = list(
    s1 = c(0.0494, 0.3506),
    s2 = c(0.1922, 1.0000),
    s3 = c(0.9500, 1.0000)
  ),

  `agresti-caffo` = list(
    s1 = c(0.0525, 0.3358),
    s2 = c(0.1600, 0.8400),
    s3 = c(0.6922, 1.0000)
  )
)


test_that("binomDiffCI matches the SAS reference values", {

  for (m in names(expected_sas)) {

    for (s in names(scenarios_sas)) {

      sc <- scenarios_sas[[s]]
      ci <- binomDiffCI(sc$x1, sc$n1, sc$x2, sc$n2, method = m)

      expect_equal(as.numeric(ci[c("lci", "uci")]),
                   expected_sas[[m]][[s]], tolerance = tol,
                   info = paste(m, s))
    }
  }
})


# ---------------------------------------------------------------
# same reference, HIV clinical trial

expected_hiv <- list(
  wald                     = c(-0.1162,  0.0843),
  `wald-cc`                = c(-0.1259,  0.0940),
  haldane                  = c(-0.1152,  0.0834),
  `jeffreys-perks`         = c(-0.1160,  0.0843),
  `mee-farrington-manning` = c(-0.1188,  0.0857),
  `miettinen-nurminen`     = c(-0.1191,  0.0860),
  `newcombe-score`         = c(-0.1177,  0.0851),
  `newcombe-score-cc`      = c(-0.1245,  0.0918),
  `hauck-anderson`         = c(-0.1216,  0.0898),
  `agresti-caffo`          = c(-0.1168,  0.0850)
)


test_that("binomDiffCI matches the HIV clinical trial reference values", {

  for (m in names(expected_hiv)) {

    ci <- binomDiffCI(84, 101, 89, 105, method = m)

    expect_equal(as.numeric(ci[c("lci", "uci")]),
                 expected_hiv[[m]], tolerance = tol, info = m)
  }
})


test_that("binomDiffCI: one-sided bounds follow design_rules 4.1", {

  for (m in methods_bdci) {

    left  <- binomDiffCI(56, 70, 48, 80, method = m, sides = "left")
    right <- binomDiffCI(56, 70, 48, 80, method = m, sides = "right")
    two   <- binomDiffCI(56, 70, 48, 80, method = m)

    # the free side is opened up to the end of the parameter range
    expect_equal(unname(left[["uci"]]),   1, info = m)
    expect_equal(unname(right[["lci"]]), -1, info = m)

    # the closed side is the two-sided bound at level 2 * conf.level - 1,
    # so it is tighter than the two-sided one
    two.90 <- binomDiffCI(56, 70, 48, 80, method = m, conf.level = 0.90)

    expect_equal(unname(left[["lci"]]),  unname(two.90[["lci"]]), info = m)
    expect_equal(unname(right[["uci"]]), unname(two.90[["uci"]]), info = m)

    expect_gte(left[["lci"]],  two[["lci"]])
    expect_lte(right[["uci"]], two[["uci"]])
  }
})


# ---------------------------------------------------------------
# Beal interval

# exact coverage over the full sample space of two binomials
.coverage_bdci <- function(n, p1, p2, method = "beal", conf.level = 0.95) {

  delta.true <- p1 - p2
  cov        <- 0

  for (x1 in 0:n) {
    for (x2 in 0:n) {

      pr <- dbinom(x1, n, p1) * dbinom(x2, n, p2)
      ci <- binomDiffCI(x1, n, x2, n, conf.level = conf.level, method = method)

      if (delta.true >= ci[["lci"]] && delta.true <= ci[["uci"]])
        cov <- cov + pr
    }
  }

  cov
}


test_that("Beal interval achieves near-nominal coverage, n=5, p1=p2=0.3", {

  cov <- .coverage_bdci(n = 5, p1 = 0.3, p2 = 0.3,
                        method = "beal", conf.level = 0.95)

  # a discrete interval cannot achieve the nominal level exactly, but it
  # should stay clearly above 0.90
  expect_gt(cov, 0.90)
  expect_lt(cov, 1.00)
})


test_that("Beal interval behaves reasonably in extreme case, n=5, p1=0.9, p2=0.05", {

  cov <- .coverage_bdci(n = 5, p1 = 0.9, p2 = 0.05,
                        method = "beal", conf.level = 0.95)

  # it must not collapse the way Wald does
  expect_gt(cov, 0.80)
  expect_lt(cov, 1.00)
})


test_that("Beal interval is symmetric when p1 = p2 and n1 = n2", {

  ci <- binomDiffCI(x1 = 8, n1 = 20, x2 = 8, n2 = 20, method = "beal")

  expect_equal(unname(ci[["lci"]]), -unname(ci[["uci"]]), tolerance = 1e-10)
  expect_equal(unname(ci[["est"]]), 0)
})


test_that("Beal interval respects the [-1, 1] range", {

  # extreme case: all successes against all failures
  ci <- binomDiffCI(x1 = 20, n1 = 20, x2 = 0, n2 = 20, method = "beal")

  expect_gte(ci[["lci"]], -1)
  expect_lte(ci[["uci"]],  1)
})


test_that("Beal interval point estimate equals the difference of proportions", {

  ci <- binomDiffCI(x1 = 15, n1 = 40, x2 = 10, n2 = 40, method = "beal")

  expect_equal(unname(ci[["est"]]), 15/40 - 10/40, tolerance = 1e-10)
})
