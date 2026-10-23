library(testthat)
library(lumen)

# ===============================================================
# binomPairedDiffCI TESTS
# ===============================================================

# Reference values computed with the implementations of Fagerland, Lydersen
# and Laake in the contingencytables package (functions *_paired_2x2), the
# companion code of their book and of Fagerland et al. (2014); the Bentur
# and Cavo tables are their worked examples. Rounded to four digits.
tol <- 1e-4

tabs <- list(
  bentur   = rbind(c(1, 1), c(7, 12)),
  cavo     = rbind(c(59, 6), c(16, 80)),
  zero     = rbind(c(10, 0), c(0, 5)),
  oneSided = rbind(c(3, 9), c(0, 4))
)

methods_bpdci <- c("tango", "newcombe-score", "bonett-price", "agresti-min",
                   "wald", "wald-cc")

expected <- list(
  "0.95" = list(
    bentur = rbind(tango = c(-0.5172, -0.0260), `newcombe-score` = c(-0.5069, -0.0256),
                   `bonett-price` = c(-0.5084, -0.0134), `agresti-min` = c(-0.4932, -0.0285),
                   wald = c(-0.5197, -0.0517), `wald-cc` = c(-0.5293, -0.0422)),
    cavo = rbind(tango = c(-0.1240, -0.0054), `newcombe-score` = c(-0.1186, -0.0046),
                 `bonett-price` = c(-0.1195, -0.0032), `agresti-min` = c(-0.1182, -0.0045),
                 wald = c(-0.1184, -0.0058), `wald-cc` = c(-0.1186, -0.0057)),
    zero = rbind(tango = c(-0.2039, 0.2039), `newcombe-score` = c(-0.1349, 0.1349),
                 `bonett-price` = c(-0.1630, 0.1630), `agresti-min` = c(-0.1153, 0.1153),
                 wald = c(0, 0), `wald-cc` = c(0, 0)),
    oneSided = rbind(tango = c(0.2600, 0.7690), `newcombe-score` = c(0.2341, 0.7453),
                     `bonett-price` = c(0.2224, 0.7776), `agresti-min` = c(0.2446, 0.7554),
                     wald = c(0.3194, 0.8056), `wald-cc` = c(0.2886, 0.8364))
  ),
  "0.9" = list(
    bentur = rbind(tango = c(-0.4816, -0.0746), `newcombe-score` = c(-0.4748, -0.0689),
                   `bonett-price` = c(-0.4686, -0.0532), `agresti-min` = c(-0.4559, -0.0659),
                   wald = c(-0.4821, -0.0893), `wald-cc` = c(-0.4901, -0.0813)),
    cavo = rbind(tango = c(-0.1131, -0.0151), `newcombe-score` = c(-0.1097, -0.0138),
                 `bonett-price` = c(-0.1101, -0.0125), `agresti-min` = c(-0.1091, -0.0136),
                 wald = c(-0.1094, -0.0149), `wald-cc` = c(-0.1095, -0.0147)),
    zero = rbind(tango = c(-0.1528, 0.1528), `newcombe-score` = c(-0.1131, 0.1131),
                 `bonett-price` = c(-0.1368, 0.1368), `agresti-min` = c(-0.0968, 0.0968),
                 wald = c(0, 0), `wald-cc` = c(0, 0)),
    oneSided = rbind(tango = c(0.3365, 0.7423), `newcombe-score` = c(0.2896, 0.7266),
                     `bonett-price` = c(0.2670, 0.7330), `agresti-min` = c(0.2857, 0.7143),
                     wald = c(0.3585, 0.7665), `wald-cc` = c(0.3326, 0.7924))
  )
)


test_that("all methods reproduce the reference implementation", {

  for (cl in names(expected))
    for (tb in names(tabs))
      for (m in methods_bpdci) {
        ci <- binomPairedDiffCI(tabs[[tb]], conf.level = as.numeric(cl),
                                method = m)
        # absolute: the references are rounded to four digits
        expect_lte(max(abs(unname(ci[c("lci", "uci")]) -
                           expected[[cl]][[tb]][m, ])), tol / 2 + 1e-10,
                   label = sprintf("%s, %s, conf.level %s", m, tb, cl))
        n <- tabs[[tb]]
        expect_equal(ci[["est"]], (n[1, 2] - n[2, 1]) / sum(n))
      }
})


test_that("the result has the documented structure", {

  ci <- binomPairedDiffCI(tabs$bentur)
  expect_named(ci, c("est", "lci", "uci"))
  expect_true(is.numeric(ci))

  # tango is the default
  expect_identical(ci, binomPairedDiffCI(tabs$bentur, method = "tango"))

  # several methods: one row per method
  df <- binomPairedDiffCI(tabs$bentur, method = methods_bpdci)
  expect_s3_class(df, "data.frame")
  expect_identical(nrow(df), length(methods_bpdci))
  expect_identical(names(df)[1:3], c("est", "lci", "uci"))
  expect_identical(df$method, methods_bpdci)
  expect_equal(df$lci[df$method == "newcombe-score"],
               binomPairedDiffCI(tabs$bentur, method = "newcombe-score")[["lci"]])

  # several levels
  df <- binomPairedDiffCI(tabs$cavo, conf.level = c(0.9, 0.95, 0.99))
  expect_identical(nrow(df), 3L)
  expect_true(all(diff(df$uci - df$lci) > 0))
})


test_that("two outcome vectors give the same result as their table", {

  set.seed(1)
  ok1 <- runif(150) < 0.8
  ok2 <- ifelse(runif(150) < 0.7, ok1, runif(150) < 0.75)
  tab <- table(factor(ok1, c(TRUE, FALSE)), factor(ok2, c(TRUE, FALSE)))

  for (m in methods_bpdci)
    expect_equal(binomPairedDiffCI(ok1, ok2, method = m),
                 binomPairedDiffCI(tab, method = m))

  # 0/1 coding as well
  expect_equal(binomPairedDiffCI(as.integer(ok1), as.integer(ok2)),
               binomPairedDiffCI(ok1, ok2))

  # est is the difference of the two proportions
  expect_equal(binomPairedDiffCI(ok1, ok2)[["est"]], mean(ok1) - mean(ok2))

  # swapping the measurements mirrors the interval
  a <- binomPairedDiffCI(ok1, ok2)
  b <- binomPairedDiffCI(ok2, ok1)
  expect_equal(unname(b), unname(-a[c("est", "uci", "lci")]), tolerance = 1e-8)
})


test_that("missing values follow the suite's NA rule", {

  x <- c(TRUE, TRUE, FALSE, NA, TRUE, FALSE, TRUE, TRUE)
  y <- c(TRUE, FALSE, FALSE, TRUE, NA, TRUE, TRUE, FALSE)

  # without na.rm: NA, no error, structure kept
  ci <- binomPairedDiffCI(x, y)
  expect_named(ci, c("est", "lci", "uci"))
  expect_true(all(is.na(ci)))

  # with na.rm: incomplete pairs removed as a whole
  ok <- !is.na(x) & !is.na(y)
  expect_equal(binomPairedDiffCI(x, y, na.rm = TRUE),
               binomPairedDiffCI(x[ok], y[ok]))

  expect_true(all(is.na(binomPairedDiffCI(matrix(c(1, NA, 2, 3), 2)))))
})


test_that("one-sided intervals are the ends of the two-sided ones", {

  two <- binomPairedDiffCI(tabs$bentur, conf.level = 0.90)
  left <- binomPairedDiffCI(tabs$bentur, conf.level = 0.95, sides = "left")
  right <- binomPairedDiffCI(tabs$bentur, conf.level = 0.95, sides = "right")

  expect_equal(left[["lci"]], two[["lci"]])
  expect_identical(left[["uci"]], 1)
  expect_equal(right[["uci"]], two[["uci"]])
  expect_identical(right[["lci"]], -1)
})


test_that("limits stay in [-1, 1] and bound the estimate", {

  for (n in list(rbind(c(0, 20), c(0, 0)), rbind(c(0, 0), c(20, 0)),
                 rbind(c(0, 3), c(0, 0)), rbind(c(5, 1), c(0, 0))))
    for (m in methods_bpdci) {
      ci <- binomPairedDiffCI(n, method = m)
      expect_true(ci[["lci"]] >= -1 && ci[["uci"]] <= 1, label = m)
      expect_true(ci[["lci"]] <= ci[["est"]] + 1e-8 || m == "bonett-price",
                  label = m)
    }

  # all pairs discordant in one direction: the estimate sits on the boundary
  expect_identical(binomPairedDiffCI(rbind(c(0, 20), c(0, 0)))[["uci"]], 1)
  expect_identical(binomPairedDiffCI(rbind(c(0, 0), c(20, 0)))[["lci"]], -1)
})


test_that("coverage of the recommended intervals is close to nominal", {

  skip_on_cran()

  set.seed(42)
  N <- 80
  # cell probabilities n11, n12, n21, n22, positively correlated measurements
  p <- c(0.55, 0.15, 0.08, 0.22)
  delta <- p[2] - p[3]
  tabsSim <- replicate(2000, matrix(rmultinom(1, N, p), 2, byrow = TRUE),
                       simplify = FALSE)

  for (m in c("tango", "newcombe-score", "bonett-price")) {
    cover <- vapply(tabsSim, function(n) {
      ci <- binomPairedDiffCI(n, method = m)
      ci[["lci"]] <= delta && delta <= ci[["uci"]]
    }, logical(1L))
    expect_gt(mean(cover), 0.93, label = m)
    expect_lt(mean(cover), 0.975, label = m)
  }
})


test_that("invalid arguments are rejected", {

  expect_error(binomPairedDiffCI(matrix(1:6, 2)), "2x2 table")
  expect_error(binomPairedDiffCI(c(1, 2)), "2x2 table")
  expect_error(binomPairedDiffCI(matrix(c(1, -1, 2, 3), 2)), "non-negative")
  expect_error(binomPairedDiffCI(matrix(c(1, 1.5, 2, 3), 2)), "integer")
  expect_error(binomPairedDiffCI(matrix(0, 2, 2)), "no pairs")
  expect_error(binomPairedDiffCI(c(TRUE, FALSE), c(TRUE)), "same length")
  expect_error(binomPairedDiffCI(c(1, 2), c(0, 1)), "logical or 0/1")
  expect_error(binomPairedDiffCI(tabs$bentur, method = "exact"))
  expect_error(binomPairedDiffCI(tabs$bentur, conf.level = 1.2), "conf.level")
  expect_error(binomPairedDiffCI(tabs$bentur, conf.level = 0.4, sides = "left"),
               "one-sided")
  expect_error(binomPairedDiffCI(tabs$bentur, na.rm = NA), "na.rm")
})
