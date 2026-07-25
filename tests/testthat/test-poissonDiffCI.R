
# 15 events in 100 person-years against 6 events in 120 person-years.
# Reference values computed independently; the exact rate-ratio interval is
# additionally checked against stats::poisson.test, which uses the same
# conditional binomial construction.

test_that("poissonRatioCI reproduces the conditional exact interval", {

  r <- poissonRatioCI(15, 100, 6, 120)

  expect_named(r, c("est", "lci", "uci"))
  expect_equal(unname(r["est"]), 3)
  expect_equal(unname(r["lci"]), 1.099947155, tolerance = 1e-8)
  expect_equal(unname(r["uci"]), 9.437411052, tolerance = 1e-8)

  # same construction as the two-sample conditional test in stats
  expect_equal(as.vector(r[c("lci", "uci")]),
               as.vector(poisson.test(c(15, 6), c(100, 120))$conf.int),
               tolerance = 1e-8)
})


test_that("the other rate-ratio methods give their documented limits", {

  expect_equal(unname(poissonRatioCI(15, 100, 6, 120, method = "wald-log")[c("lci", "uci")]),
               c(1.163996961, 7.731978951), tolerance = 1e-8)
  # uniroot runs at tol = 1e-10, so the mid-p limits are pinned tightly
  expect_equal(unname(poissonRatioCI(15, 100, 6, 120, method = "midp")[c("lci", "uci")]),
               c(1.18888277839, 8.41548475954), tolerance = 1e-9)

  # mid-p is shorter than the exact interval it approximates
  ex <- poissonRatioCI(15, 100, 6, 120)
  mp <- poissonRatioCI(15, 100, 6, 120, method = "midp")
  expect_gt(mp[["lci"]], ex[["lci"]])
  expect_lt(mp[["uci"]], ex[["uci"]])
})


test_that("poissonDiffCI reproduces the MOVER and Wald intervals", {

  r <- poissonDiffCI(15, 100, 6, 120)
  expect_equal(unname(r["est"]), 0.1)
  expect_equal(unname(r[c("lci", "uci")]),
               c(0.01155262687, 0.2024156465), tolerance = 1e-8)

  w <- poissonDiffCI(15, 100, 6, 120, method = "wald")
  expect_equal(unname(w[c("lci", "uci")]),
               c(0.01419326324, 0.1858067368), tolerance = 1e-8)

  # MOVER is asymmetric about the estimate, Wald is symmetric by construction
  expect_equal(w[["uci"]] - w[["est"]], w[["est"]] - w[["lci"]])
  expect_false(isTRUE(all.equal(r[["uci"]] - r[["est"]], r[["est"]] - r[["lci"]])))
})


test_that("sides names the side carrying the finite limit", {

  # full alpha on one side, open end reported at the range limit
  l <- poissonRatioCI(15, 100, 6, 120, sides = "left")
  expect_equal(unname(l["lci"]), 1.262100115, tolerance = 1e-8)
  expect_equal(unname(l["uci"]), Inf)

  rr <- poissonRatioCI(15, 100, 6, 120, sides = "right")
  expect_equal(unname(rr["lci"]), 0)
  expect_true(is.finite(rr[["uci"]]))

  # the ratio is bounded below by zero, the difference is not
  expect_equal(unname(poissonDiffCI(15, 100, 6, 120, sides = "right")["lci"]), -Inf)
  expect_equal(unname(poissonDiffCI(15, 100, 6, 120, sides = "left")["uci"]), Inf)

  # one-sided limits are tighter than the corresponding two-sided ones
  expect_gt(l[["lci"]], poissonRatioCI(15, 100, 6, 120)[["lci"]])
})


test_that("arguments are recycled to a common length", {

  r <- poissonDiffCI(x1 = c(15, 20), n1 = 100, x2 = 6, n2 = 120)

  expect_s3_class(r, "data.frame")
  expect_equal(nrow(r), 2L)
  expect_equal(names(r)[1:3], c("est", "lci", "uci"))
  expect_equal(names(r)[4:8], c("x1", "n1", "x2", "n2", "conf.level"))
  expect_equal(r$x1, c(15, 20))
  expect_equal(r$n2, c(120, 120))
  expect_equal(r$est, c(0.1, 0.15))
  expect_equal(r$lci, c(0.01155262687, 0.05243411397), tolerance = 1e-8)
  expect_equal(r$uci, c(0.2024156465, 0.263390721), tolerance = 1e-8)

  q <- poissonRatioCI(x1 = c(15, 20), n1 = 100, x2 = 6, n2 = 120)
  expect_equal(q$lci, c(1.099947155, 1.549297746), tolerance = 1e-8)
  expect_equal(q$uci, c(9.437411052, 12.17194693), tolerance = 1e-8)

  # a single case still returns a plain named vector
  expect_true(is.numeric(poissonDiffCI(15, 100, 6, 120)))
  expect_null(dim(poissonDiffCI(15, 100, 6, 120)))
})


test_that("zero counts are handled at the boundary", {

  expect_equal(unname(poissonRatioCI(0, 100, 6, 120)[c("est", "lci", "uci")]),
               c(0, 0, 1.019173433), tolerance = 1e-8)
  expect_equal(unname(poissonRatioCI(15, 100, 0, 120)[c("lci", "uci")]),
               c(4.304098329, Inf), tolerance = 1e-8)

  # no events at all: the ratio is not estimable, the interval is uninformative
  z <- poissonRatioCI(0, 100, 0, 120)
  expect_equal(unname(z[c("lci", "uci")]), c(0, Inf))

  # the log scale is undefined with a zero count
  expect_equal(unname(poissonRatioCI(0, 100, 6, 120, method = "wald-log")[c("lci", "uci")]),
               c(0, Inf))

  expect_equal(unname(poissonDiffCI(0, 100, 6, 120)[c("est", "lci", "uci")]),
               c(-0.05, -0.1088289502, -0.001393812691), tolerance = 1e-7)
})


test_that("only an exhausted upper limit becomes infinite", {

  # x2 = 0 puts the whole conditional mass at x1, so the limit is genuinely Inf
  expect_equal(unname(poissonRatioCI(15, 100, 0, 120)["uci"]), Inf)

  # a limit merely close to the boundary must stay finite rather than be
  # rounded up to Inf: here the binomial upper limit is 0.99975
  r <- poissonRatioCI(100, 100, 1, 120)
  expect_equal(unname(r["est"]), 120)
  expect_true(is.finite(r[["uci"]]))
  expect_equal(unname(r[c("lci", "uci")]), c(21.05009771, 4786.544318),
               tolerance = 1e-7)
})


test_that("a wider confidence level gives a wider interval", {

  a <- poissonRatioCI(15, 100, 6, 120, conf.level = 0.90)
  b <- poissonRatioCI(15, 100, 6, 120, conf.level = 0.99)
  expect_gt(a[["lci"]], b[["lci"]])
  expect_lt(a[["uci"]], b[["uci"]])
  expect_equal(a[["est"]], b[["est"]])
})


test_that("invalid arguments are rejected", {

  expect_error(poissonRatioCI(15, 100, 6, 120, conf.level = 1), "'conf.level'")
  expect_error(poissonRatioCI(15, 100, 6, 120, conf.level = 0), "'conf.level'")
  expect_error(poissonRatioCI(-1, 100, 6, 120), "non-negative")
  expect_error(poissonRatioCI(15, 0, 6, 120), "positive")
  # logical NA fails the numeric check, NA_real_ the finiteness check
  expect_error(poissonDiffCI(15, 100, NA, 120), "must be numeric")
  expect_error(poissonDiffCI(15, 100, NA_real_, 120), "must be finite")
  expect_error(poissonRatioCI(15, 100, 6, 120, method = "katz"), "should be one of")
  expect_error(poissonDiffCI(15, 100, 6, 120, sides = "both"), "should be one of")
})
