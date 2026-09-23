
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
  "pratt", "mid-p", "blaker", "likelihood", "khouadji"
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
