
# ===============================================================
# binomCI TESTS
# ===============================================================
# Formerly a script with top-level stopifnot()/print(): its checks were not
# reported as expectations, a failure aborted the whole file, and every run
# printed the intervals to the console. Same checks, now as test_that blocks.

methods_bci <- c(
  "wald", "wald-cc", "wilson", "wilson-cc", "wilson-mod",
  "agresti-coull", "jeffreys", "jeffreys-mod",
  "clopper-pearson", "arcsine", "logit",
  "pratt", "mid-p", "blaker", "wang", "likelihood", "khouadji"
)
# witting excluded from deterministic tests (randomized)

# structure and range of a single interval
.expect_bci <- function(res, info = NULL) {
  nms <- if (is.null(dim(res))) names(res) else colnames(res)
  expect_true(is.numeric(as.matrix(res)), info = info)
  expect_true(all(c("est", "lci", "uci") %in% nms), info = info)
  expect_gte(res[["lci"]], 0, label = paste("lci", info))
  expect_lte(res[["uci"]], 1, label = paste("uci", info))
  expect_lte(res[["lci"]], res[["uci"]], label = paste("lci <= uci", info))
}


test_that("basic functionality: every method gives a valid interval around est", {
  for (m in methods_bci) {
    res <- binomCI(x = 37, n = 43, method = m)
    .expect_bci(res, info = m)
    expect_lte(res[["lci"]], res[["est"]], label = m)
    expect_gte(res[["uci"]], res[["est"]], label = m)
  }
})


test_that("reference values: prop.test and binom.test", {
  # Wilson == prop.test(correct = FALSE)
  expect_equal(unname(binomCI(37, 43, method = "wilson")[c("lci", "uci")]),
               prop.test(37, 43, correct = FALSE)$conf.int[1:2], tolerance = 1e-6)
  # Wilson-cc == prop.test(correct = TRUE)
  expect_equal(unname(binomCI(37, 43, method = "wilson-cc")[c("lci", "uci")]),
               prop.test(37, 43, correct = TRUE)$conf.int[1:2], tolerance = 1e-6)
  # Clopper-Pearson == binom.test
  expect_equal(unname(binomCI(42, 43, method = "clopper-pearson")[c("lci", "uci")]),
               binom.test(42, 43)$conf.int[1:2], tolerance = 1e-6)
})


test_that("reference values: Newcombe (1998) Table I, Wilson", {
  res <- binomCI(x = 81, n = 263, method = "wilson")
  # the table gives 4 decimals: absolute difference, not relative tolerance
  expect_lt(abs(res[["lci"]] - 0.2553), 0.001)
  expect_lt(abs(res[["uci"]] - 0.3662), 0.001)
})


test_that("point estimate is x/n for all standard methods", {
  for (m in methods_bci)
    expect_equal(binomCI(x = 15, n = 40, method = m)[["est"]], 15/40,
                 tolerance = 1e-10, label = m)
})


test_that("a higher conf.level gives a wider interval", {
  for (m in methods_bci) {
    ci95 <- binomCI(x = 15, n = 40, method = m, conf.level = 0.95)
    ci99 <- binomCI(x = 15, n = 40, method = m, conf.level = 0.99)
    expect_lte(ci99[["lci"]], ci95[["lci"]], label = m)
    expect_gte(ci99[["uci"]], ci95[["uci"]], label = m)
  }
})


test_that("bounds stay in [0, 1] near the edges", {
  for (m in methods_bci) {
    for (x in c(0, 1, 2, 38, 39, 40)) {
      # logit is undefined at x = 0 and x = n
      if (m == "logit" && x %in% c(0, 40)) next
      res <- binomCI(x = x, n = 40, method = m)
      info <- paste(m, "x =", x)
      expect_false(anyNA(res[c("lci", "uci")]), info = info)
      expect_gte(res[["lci"]], 0, label = info)
      expect_lte(res[["uci"]], 1, label = info)
    }
  }
})


test_that("edge case x = 0", {
  for (m in setdiff(methods_bci, "logit")) {
    res <- binomCI(x = 0, n = 20, method = m)
    expect_equal(res[["est"]], 0, label = m)
    # arcsine uses p.tilde > 0 internally, so lci is not exactly 0
    if (m != "arcsine") expect_equal(res[["lci"]], 0, label = m)
    else expect_true(res[["lci"]] >= 0 && res[["lci"]] < 0.05, label = m)
    # wald collapses to uci = 0 when x = 0 (known limitation)
    if (m != "wald") expect_gt(res[["uci"]], 0, label = m)
  }
})


test_that("edge case x = n", {
  for (m in setdiff(methods_bci, "logit")) {
    res <- binomCI(x = 20, n = 20, method = m)
    expect_equal(res[["est"]], 1, label = m)
    # arcsine uses p.tilde < 1 internally, so uci is not exactly 1
    if (m != "arcsine") expect_equal(res[["uci"]], 1, label = m)
    else expect_true(res[["uci"]] > 0.95 && res[["uci"]] <= 1, label = m)
    # wald collapses to lci = 1 when x = n (known limitation)
    if (m != "wald") expect_lt(res[["lci"]], 1, label = m)
  }
})


test_that("one-sided intervals", {
  for (m in setdiff(methods_bci, "khouadji")) {
    res.left  <- binomCI(x = 10, n = 40, method = m, sides = "left")
    res.right <- binomCI(x = 10, n = 40, method = m, sides = "right")
    res.two   <- binomCI(x = 10, n = 40, method = m, sides = "two.sided")

    expect_equal(res.left[["uci"]], 1, label = m)
    expect_equal(res.right[["lci"]], 0, label = m)

    # the one-sided bound is tighter on the closed side: it is the two-sided
    # bound at level 2 * conf.level - 1, the other side is opened up
    expect_gte(res.left[["lci"]],  res.two[["lci"]], label = m)
    expect_lte(res.right[["uci"]], res.two[["uci"]], label = m)
  }
})


test_that("stdEst = FALSE returns the adjusted estimator for agresti-coull", {
  res.std <- binomCI(x = 81, n = 263, method = "agresti-coull", stdEst = TRUE)
  res.adj <- binomCI(x = 81, n = 263, method = "agresti-coull", stdEst = FALSE)
  expect_equal(res.std[["est"]], 81/263)
  expect_false(isTRUE(all.equal(res.adj[["est"]], 81/263)))   # p.tilde != x/n
})


test_that("conf.level outside (0, 1) is an error", {
  for (bad.level in c(0, 1, -0.5, 1.5))
    expect_error(binomCI(x = 10, n = 40, conf.level = bad.level),
                 info = format(bad.level))
})


test_that("vectorization returns a data frame with one row per x", {
  res <- binomCI(x = c(42, 35, 23, 22), n = 43, method = "wilson")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 4L)

  res2 <- binomCI(x = c(42, 35, 23, 22), n = c(50, 60, 70, 80),
                  method = "jeffreys")
  expect_s3_class(res2, "data.frame")
  expect_equal(nrow(res2), 4L)
})


test_that("witting: bounds in [0, 1]", {
  set.seed(42)
  .expect_bci(binomCI(x = 10, n = 30, method = "witting"), info = "witting")
})


test_that("witting works at x = 0 and x = n", {
  # uniroot() failed in about 95 % of the calls here
  set.seed(1)
  for (i in 1:50) {
    r0 <- binomCI(0, 20, method = "witting")
    rn <- binomCI(20, 20, method = "witting")
    expect_equal(r0[["lci"]], 0)
    expect_equal(rn[["uci"]], 1)
    .expect_bci(r0, info = "x = 0"); .expect_bci(rn, info = "x = n")
  }
})


test_that("witting bounds solve their defining equation with alpha/2", {
  # P(X + U <= x.tilde | p) at the bound equals 1 - alpha/2 (lci) and
  # alpha/2 (uci); x.tilde is recovered from the p.tilde attribute
  pAbs <- function(p, t, n)
    pbinom(floor(t) - 1, n, p) + (t - floor(t)) * dbinom(floor(t), n, p)
  set.seed(7)
  for (i in 1:20) {
    n <- 25; x <- sample(1:(n - 1), 1)
    r <- lumen:::.binomCI.witting(x, n, 0.05)
    t <- attr(r, "p.tilde") * n
    expect_equal(pAbs(r[["lci"]], t, n), 0.975, tolerance = 1e-7)
    expect_equal(pAbs(r[["uci"]], t, n), 0.025, tolerance = 1e-7)
  }
})


test_that("witting keeps the level exactly (Monte Carlo)", {
  # randomized exactness: coverage 1 - alpha for every p; with alpha per
  # tail it was about 0.885 two-sided and 0.90 one-sided at 95 %
  set.seed(11)
  n <- 20; p <- 0.3; R <- 4000
  x  <- rbinom(R, n, p)
  ci <- t(vapply(x, function(xx)
    binomCI(xx, n, method = "witting")[c("lci", "uci")], numeric(2)))
  lo <- vapply(x, function(xx)
    binomCI(xx, n, sides = "left", method = "witting")[["lci"]], numeric(1))
  se <- sqrt(0.95 * 0.05 / R)
  expect_lt(abs(mean(ci[, 1] <= p & p <= ci[, 2]) - 0.95), 4 * se)
  expect_lt(abs(mean(lo <= p) - 0.95), 4 * se)
})


# --- regression test ------------------------------------------------

test_that("blaker limits are nested in clopper-pearson for all scales", {
  
  for (nn in c(19, 43, 263, 1e4, 1e5)) {
    for (xx in unique(round(c(1, 3, 0.1, 0.5, 0.9) * c(1, 1, nn, nn, nn)))) {
      
      cp <- binomCI(xx, nn, method = "clopper-pearson")
      bl <- binomCI(xx, nn, method = "blaker")
      
      expect_gte(bl[["lci"]], cp[["lci"]])
      expect_lte(bl[["uci"]], cp[["uci"]])
      expect_lte(bl[["lci"]], xx / nn)
      expect_gte(bl[["uci"]], xx / nn)
    }
  }
  
  # the case that used to return lci = 1, uci = 0
  expect_equal(unname(binomCI(50000, 100000, method = "blaker")[c("lci", "uci")]),
               c(0.4968999707, 0.5031000293), tolerance = 1e-8)
})


test_that("stress test: random x, n for all methods", {

  # collect failures instead of one expectation per call (300 x 16), so a
  # failure reports every offending case at once
  set.seed(123)
  bad <- character()

  for (i in 1:300) {
    n <- sample(5:200, 1)
    x <- sample(0:n, 1)

    for (m in methods_bci) {
      # logit can be NA at x = 0 or x = n
      if (m == "logit" && x %in% c(0, n)) next

      case <- sprintf("%s x=%d n=%d", m, x, n)
      res  <- tryCatch(binomCI(x = x, n = n, method = m),
                       error = function(e) e)

      if (inherits(res, "error")) {
        bad <- c(bad, paste(case, "error:", conditionMessage(res)))
        next
      }
      lu <- res[c("lci", "uci")]
      if (anyNA(lu) || lu[1] < 0 || lu[2] > 1 || lu[1] > lu[2])
        bad <- c(bad, sprintf("%s lci=%g uci=%g", case, lu[1], lu[2]))
    }
  }

  expect_length(bad, 0)
  if (length(bad)) message(paste(head(bad, 20), collapse = "\n"))
})
test_that("wang reproduces Wang (2014) / ExactCIone::WbinoCI", {
  # ExactCIone 1.0.5, WbinoCI(x, 5, 0.95, details = TRUE)$CIM
  ref <- rbind(c(0,          0.5000000), c(0.01020614, 0.6574084),
               c(0.07644030, 0.8107447), c(0.18925530, 0.9235597),
               c(0.34259163, 0.9897939), c(0.49999997, 1))
  got <- t(sapply(0:5, function(x)
    binomCI(x, 5, method = "wang")[c("lci", "uci")]))
  expect_equal(unname(got), ref, tolerance = 1e-6)
})

test_that("wang is nested in clopper-pearson and keeps exact coverage", {
  for (n in c(7, 20, 45)) {
    w  <- binomCI(0:n, n, method = "wang")
    cp <- binomCI(0:n, n, method = "clopper-pearson")
    expect_true(all(w$lci >= cp$lci - 1e-12 & w$uci <= cp$uci + 1e-12))
    pts <- unique(c(w$lci, w$uci)); pts <- pts[pts > 0 & pts < 1]
    pp  <- c(pts - 1e-9, pts + 1e-9, seq(0.001, 0.999, length.out = 500))
    cov <- vapply(pp, function(p)
      sum(dbinom(0:n, n, p)[w$lci <= p & p <= w$uci]), numeric(1))
    expect_gte(min(cov), 0.95 - 1e-9)
  }
})

test_that("wang keeps the level also between two nearly equal limits", {
  # the one-sided coverage limits at every breakpoint; shifting p by a
  # relative 1e-12 instead missed a piece of width 2e-13 at n = 11 between
  # U(0) and L(6) with coverage 0.919
  infcov <- function(n, conf.level) {
    w <- binomCI(0:n, n, conf.level = conf.level, method = "wang")
    p <- unique(c(w$lci, w$uci)); p <- p[p > 0 & p < 1]
    min(vapply(p, function(pp) {
      d <- dbinom(0:n, n, pp)
      min(sum(d[w$lci < pp & w$uci >= pp]), sum(d[w$lci <= pp & w$uci > pp]))
    }, numeric(1)))
  }
  expect_gte(infcov(11, 0.95), 0.95 - 1e-12)
  for (n in c(30, 41, 47, 56)) for (cl in c(0.95, 0.99, 0.999))
    expect_gte(infcov(n, cl), cl - 1e-12, label = sprintf("n=%d cl=%g", n, cl))
})

test_that("wang rejects non-integer counts instead of truncating", {
  expect_error(.wangBinomCI(2.5, 5, 0.05), "integer")
  expect_error(.wangBinomCI(2, 5, 1.2), "alpha")
})

test_that("blaker and wang one-sided equal one-sided clopper-pearson", {
  # an end of the two-sided blaker interval at the doubled alpha is no
  # one-sided bound: x = 1, n = 1 gave the lower bound 0.1 at 95 %, which
  # covers p just below 0.1 with probability 0.9 only
  for (m in c("blaker", "wang")) for (s in c("left", "right"))
    expect_equal(binomCI(3, 17, method = m, sides = s),
                 binomCI(3, 17, method = "clopper-pearson", sides = s),
                 label = paste(m, s))
  expect_equal(binomCI(1, 1, sides = "left", method = "blaker")[["lci"]],
               0.05)
})


test_that("wang stays valid at low levels", {
  for (n in c(2, 7, 13, 30)) for (cl in c(0.05, 0.2, 0.3, 0.5)) {
    w <- binomCI(0:n, n, conf.level = cl, method = "wang")
    lab <- sprintf("n=%d cl=%g", n, cl)
    expect_true(all(w$lci <= w$uci) && !is.unsorted(w$lci), label = lab)
    p <- unique(c(w$lci, w$uci)); p <- p[p > 0 & p < 1]
    cov <- vapply(p, function(pp) {
      d <- dbinom(0:n, n, pp)
      min(sum(d[w$lci < pp & w$uci >= pp]), sum(d[w$lci <= pp & w$uci > pp]))
    }, numeric(1))
    expect_gte(min(cov), cl - 1e-12, label = lab)
  }
})

test_that("wang is symmetric in x <-> n - x", {
  a <- binomCI(4, 31, method = "wang"); b <- binomCI(27, 31, method = "wang")
  expect_equal(a[["lci"]], 1 - b[["uci"]])
  expect_equal(a[["uci"]], 1 - b[["lci"]])
})
