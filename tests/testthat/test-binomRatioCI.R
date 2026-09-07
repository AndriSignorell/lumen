library(testthat)
library(lumen)

# ===============================================================
# binomRatioCI TESTS
# ===============================================================

methods_brc <- c("katz-log", "adj-log", "bailey", "koopman", "noether",
                 "sinh-1")

# helper: the structural invariants every result has to satisfy.
# Expectations raised inside a helper are counted by testthat, so this
# stays a function rather than being inlined everywhere.
.check_brc <- function(res) {

  expect_true(is.numeric(as.matrix(res)))

  nms <- if (is.null(dim(res))) names(res) else colnames(res)

  expect_true(all(c("est", "lci", "uci") %in% nms))
  expect_lte(res[["lci"]], res[["uci"]])

  invisible(TRUE)
}


test_that("binomRatioCI: every method returns est between the limits", {

  for (m in methods_brc) {

    res <- binomRatioCI(x1 = 12, n1 = 100, x2 = 20, n2 = 100, method = m)

    .check_brc(res)
    expect_lte(res[["lci"]], res[["est"]])
    expect_gte(res[["uci"]], res[["est"]])
  }
})


test_that("binomRatioCI: koopman reproduces Koopman (1984), Table 1", {

  # x1 = 36, n1 = 40, x2 = 16, n2 = 80, reported 95% CI (2.93, 7.15)
  res <- binomRatioCI(x1 = 36, n1 = 40, x2 = 16, n2 = 80, method = "koopman")

  expect_equal(unname(res[["est"]]), 4.5,    tolerance = 0.001)
  expect_equal(unname(res[["lci"]]), 2.9396, tolerance = 0.001)
  expect_equal(unname(res[["uci"]]), 7.1522, tolerance = 0.001)
})


test_that("binomRatioCI: a higher conf.level gives a wider interval", {

  for (m in methods_brc) {

    ci95 <- binomRatioCI(x1 = 12, n1 = 100, x2 = 20, n2 = 100,
                         method = m, conf.level = 0.95)
    ci99 <- binomRatioCI(x1 = 12, n1 = 100, x2 = 20, n2 = 100,
                         method = m, conf.level = 0.99)

    expect_lte(ci99[["lci"]], ci95[["lci"]])
    expect_gte(ci99[["uci"]], ci95[["uci"]])
  }
})


test_that("binomRatioCI: vectorized arguments give a data frame", {

  res <- binomRatioCI(
    x1 = c(5, 10, 20), n1 = c(50, 100, 200),
    x2 = c(4,  8, 25), n2 = c(50, 100, 200),
    method = c("katz-log", "koopman")
  )

  expect_s3_class(res, "data.frame")
})


test_that("binomRatioCI: both counts zero gives [0, Inf)", {

  for (m in methods_brc) {

    res <- binomRatioCI(x1 = 0, n1 = 100, x2 = 0, n2 = 100, method = m)

    expect_equal(unname(res[["est"]]), 0)
    expect_equal(unname(res[["lci"]]), 0)
    expect_true(is.infinite(res[["uci"]]))
  }
})


test_that("binomRatioCI: a zero numerator gives est = 0 and lci = 0", {

  for (m in methods_brc) {

    res <- binomRatioCI(x1 = 0, n1 = 100, x2 = 10, n2 = 100, method = m)

    expect_equal(unname(res[["est"]]), 0)
    expect_equal(unname(res[["lci"]]), 0)
  }
})


test_that("binomRatioCI: a zero denominator gives an infinite est and uci", {

  for (m in methods_brc) {

    res <- binomRatioCI(x1 = 10, n1 = 100, x2 = 0, n2 = 100, method = m)

    expect_true(is.infinite(res[["est"]]))
    expect_false(is.na(res[["lci"]]))
    expect_gte(res[["lci"]], 0)
    expect_true(is.infinite(res[["uci"]]))
  }
})


test_that("binomRatioCI: both proportions at one stays well defined", {

  for (m in methods_brc)
    .check_brc(binomRatioCI(x1 = 100, n1 = 100, x2 = 100, n2 = 100,
                            method = m))
})


test_that("binomRatioCI: one-sided intervals open the free side", {

  for (m in methods_brc) {

    res.left  <- binomRatioCI(x1 = 10, n1 = 100, x2 = 20, n2 = 100,
                              method = m, sides = "left")
    res.right <- binomRatioCI(x1 = 10, n1 = 100, x2 = 20, n2 = 100,
                              method = m, sides = "right")
    res.two   <- binomRatioCI(x1 = 10, n1 = 100, x2 = 20, n2 = 100,
                              method = m, sides = "two.sided")

    expect_true(is.infinite(res.left[["uci"]]))
    expect_equal(unname(res.right[["lci"]]), 0)

    # the one-sided bound is tighter on the closed side: it is the
    # two-sided bound at level 2 * conf.level - 1 (design_rules 4.1)
    expect_gte(res.left[["lci"]],  res.two[["lci"]])
    expect_lte(res.right[["uci"]], res.two[["uci"]])

    res.90 <- binomRatioCI(x1 = 10, n1 = 100, x2 = 20, n2 = 100,
                           method = m, conf.level = 0.90)

    expect_equal(unname(res.left[["lci"]]),  unname(res.90[["lci"]]))
    expect_equal(unname(res.right[["uci"]]), unname(res.90[["uci"]]))
  }
})


test_that("binomRatioCI: est increases with the numerator count", {

  for (m in methods_brc) {

    rr1 <- binomRatioCI(x1 =  5, n1 = 100, x2 = 10, n2 = 100, method = m)
    rr2 <- binomRatioCI(x1 = 20, n1 = 100, x2 = 10, n2 = 100, method = m)

    expect_gt(rr2[["est"]], rr1[["est"]])
  }
})


test_that("binomRatioCI: invalid input throws", {

  expect_error(binomRatioCI(x1 = 101, n1 = 100, x2 = 10, n2 = 100))

  for (bad.level in c(0, 1, -0.5, 1.5))
    expect_error(binomRatioCI(x1 = 10, n1 = 100, x2 = 10, n2 = 100,
                              conf.level = bad.level))
})


test_that("binomRatioCI: 500 random configurations stay well defined", {

  # the per-iteration checks are collected rather than asserted one by one:
  # 500 draws times six methods would otherwise add 9000 expectations to
  # the suite for a single property
  set.seed(123)
  bad <- character()

  for (i in 1:500) {

    n1 <- sample(10:500, 1)
    n2 <- sample(10:500, 1)
    x1 <- sample(0:n1, 1)
    x2 <- sample(0:n2, 1)

    for (m in methods_brc) {

      res <- try(binomRatioCI(x1 = x1, n1 = n1, x2 = x2, n2 = n2,
                              method = m),
                 silent = TRUE)

      label <- sprintf("x1=%d n1=%d x2=%d n2=%d method=%s", x1, n1, x2, n2, m)

      if (inherits(res, "try-error")) {
        bad <- c(bad, paste(label, "-> error"))
        next
      }

      if (anyNA(res[c("lci", "uci")]))
        bad <- c(bad, paste(label, "-> NA limit"))

      if (!(is.infinite(res[["lci"]]) || is.infinite(res[["uci"]]) ||
            res[["lci"]] <= res[["uci"]]))
        bad <- c(bad, paste(label, "-> lci > uci"))

      if (!isTRUE(res[["lci"]] >= 0))
        bad <- c(bad, paste(label, "-> negative lci"))
    }
  }

  expect_identical(bad, character())
})
