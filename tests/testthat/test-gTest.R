
# Agresti (2007) p.39 - gender x party
M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
dimnames(M) <- list(gender = c("M","F"), party = c("Democrat","Independent","Republican"))

test_that("gTest: returns htest", {
  expect_s3_class(gTest(M), "htest")
})

test_that("gTest: statistic named G", {
  expect_named(gTest(M)$statistic, "G")
})

test_that("gTest: p.value in [0,1]", {
  res <- gTest(M)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("gTest: independence test df = (r-1)(c-1)", {
  res <- gTest(M)
  expect_equal(unname(res$parameter), (2-1) * (3-1))
})

test_that("gTest: Agresti example significant", {
  expect_lt(gTest(M)$p.value, 0.001)
})

test_that("gTest: GOF uniform p agrees with chisq.test for large n", {
  x <- c(100, 100, 100, 100)
  expect_gt(gTest(x)$p.value, 0.99)
})

test_that("gTest: GOF non-uniform gives small p", {
  x <- c(90, 10, 50, 50)
  expect_lt(gTest(x)$p.value, 0.05)
})

test_that("gTest: GOF df = k-1", {
  x <- c(20, 15, 25)
  expect_equal(unname(gTest(x)$parameter), 2L)
})

test_that("gTest: G statistic >= 0", {
  expect_gte(unname(gTest(M)$statistic), 0)
  expect_gte(unname(gTest(c(20, 15, 25))$statistic), 0)
})

test_that("gTest: observed and expected in result", {
  res <- gTest(M)
  expect_false(is.null(res$observed))
  expect_false(is.null(res$expected))
  expect_equal(sum(res$expected), sum(res$observed), tolerance = 1e-10)
})

test_that("gTest: rescaleP=TRUE allows non-summing p", {
  x <- c(89, 37, 30, 28, 2)
  p <- c(40, 20, 20, 15, 5)
  res <- gTest(x, p = p, rescaleP = TRUE)
  expect_s3_class(res, "htest")
})

test_that("gTest: williams correction gives larger p than none", {
  res_none <- gTest(M, correct = "none")
  res_will <- gTest(M, correct = "williams")
  expect_gte(res_will$p.value, res_none$p.value)
})


test_that("gTest: observed keeps original counts under Yates correction", {
  # regression test: the continuity-corrected (+-0.5) counts used to be
  # returned as 'observed'
  tab <- as.table(matrix(c(10, 4, 3, 11), 2))

  res <- gTest(tab, correct = "yates")

  expect_equal(unname(as.vector(res$observed)), c(10, 4, 3, 11))
  expect_equal(sum(res$expected), sum(res$observed), tolerance = 1e-10)
})


test_that("gTest: input validation", {
  expect_error(gTest(c(10, -1, 5)), "nonnegative")
  expect_error(gTest(c(10, 20), p = c(0.5, 0.4)), "sum to 1")
  expect_error(gTest(c(A = 5), y = NULL), "2 elements")
  expect_no_error(gTest(matrix(1:4, 2), correct = "yates"))
})


M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
dimnames(M) <- list(gender = c("M", "F"),
                    party = c("Democrat", "Independent", "Republican"))

Gref <- function(O, E) 2 * sum(O[O > 0] * log(O[O > 0] / E[O > 0]))

test_that("independence: Agresti (2007), p. 39", {
  r <- gTest(M)
  E <- outer(rowSums(M), colSums(M)) / sum(M)
  
  expect_s3_class(r, "htest")
  expect_equal(r$statistic, c(G = Gref(M, E)))
  expect_equal(round(unname(r$statistic), 1), 30.0)   # G^2 = 30.0 in Agresti
  expect_equal(r$parameter, c(df = 2L))
  expect_equal(r$p.value, pchisq(Gref(M, E), 2, lower.tail = FALSE))
  expect_equal(r$expected, E)
  expect_identical(r$observed, M)
  expect_match(r$method, "without correction")
})

test_that("independence: factor input equals the table", {
  df <- as.data.frame(M)
  g <- rep(df$gender, df$Freq)
  p <- rep(df$party, df$Freq)
  r <- gTest(g, p)
  expect_equal(unname(r$statistic), unname(gTest(M)$statistic))
  expect_identical(r$data.name, "g and p")
  
  # NA pairs are dropped
  expect_equal(gTest(c(g, NA), c(p, "Democrat"))$statistic, r$statistic)
})

test_that("independence: data frame input", {
  expect_equal(gTest(as.data.frame.matrix(M))$statistic, gTest(M)$statistic)
})

test_that("Williams' correction", {
  n <- sum(M)
  q <- 1 + ((n * sum(1 / rowSums(M)) - 1) * (n * sum(1 / colSums(M)) - 1)) /
    (6 * n * 2 * 1)
  r <- gTest(M, correct = "williams")
  expect_equal(unname(r$statistic), unname(gTest(M)$statistic) / q)
  expect_match(r$method, "Williams")
})

test_that("Yates' correction shifts by min(0.5, |O - E|)", {
  for (tab in list(matrix(c(12, 5, 7, 14), 2),     # ad > bc
                   matrix(c(3, 11, 9, 4), 2),      # ad < bc
                   matrix(c(10, 10, 10, 10.2), 2))) {
    E <- outer(rowSums(tab), colSums(tab)) / sum(tab)
    s <- min(0.5, abs(tab[1, 1] - E[1, 1]))
    O <- tab - sign(tab - E) * s
    r <- gTest(tab, correct = "yates")
    expect_equal(unname(r$statistic), Gref(O, E))
    # observed stays uncorrected
    expect_identical(r$observed, tab)
    expect_lte(unname(r$statistic), unname(gTest(tab)$statistic))
  }
})

test_that("Yates' correction never overshoots", {
  # regression: perfect independence gave G = 0.1
  expect_equal(unname(gTest(matrix(10, 2, 2), correct = "yates")$statistic), 0)
  # regression: an empty row produced negative cells and G = Inf
  r <- gTest(matrix(c(3, 0, 6, 0), 2), correct = "yates")
  expect_equal(unname(r$statistic), 0)
})

test_that("Yates' correction requires 2 x 2", {
  expect_error(gTest(M, correct = "yates"), "2 x 2")
})

test_that("goodness of fit: equal and given probabilities", {
  x <- c(A = 20, B = 15, C = 25)
  r <- gTest(x)
  E <- rep(sum(x) / 3, 3)
  expect_equal(unname(r$statistic), Gref(x, E))
  expect_equal(r$parameter, c(df = 2L))
  expect_named(r$expected, names(x))
  expect_match(r$method, "goodness of fit")
  expect_equal(gTest(as.table(x))$statistic, r$statistic)
  
  p <- c(0.5, 0.2, 0.3)
  expect_equal(unname(gTest(x, p = p)$statistic), Gref(x, sum(x) * p))
})

test_that("goodness of fit: one-row / one-column matrix", {
  x <- c(20, 15, 25)
  expect_equal(gTest(matrix(x, nrow = 1))$statistic, gTest(x)$statistic)
  expect_equal(gTest(matrix(x, ncol = 1))$statistic, gTest(x)$statistic)
})

test_that("goodness of fit: rescaleP", {
  x <- c(89, 37, 30, 28, 2)
  p <- c(40, 20, 20, 15, 5)
  expect_error(gTest(x, p = p), "sum to 1")
  expect_equal(gTest(x, p = p, rescaleP = TRUE)$statistic,
               gTest(x, p = p / sum(p))$statistic)
})

test_that("goodness of fit: Williams", {
  x <- c(20, 15, 25)
  q <- 1 + (3 + 1) / (6 * sum(x))
  expect_equal(unname(gTest(x, correct = "williams")$statistic),
               unname(gTest(x)$statistic) / q)
})

test_that("goodness of fit: Yates", {
  E <- c(20, 20)
  # x1 above E1 by more than 0.25 -> moved down by 0.5
  expect_equal(unname(gTest(c(25, 15), correct = "yates")$statistic),
               Gref(c(24.5, 15.5), E))
  # below -> moved up
  expect_equal(unname(gTest(c(15, 25), correct = "yates")$statistic),
               Gref(c(15.5, 24.5), E))
  # within 0.25 -> unchanged
  x <- c(20.2, 19.8)
  expect_equal(gTest(x, correct = "yates")$statistic, gTest(x)$statistic)
  
  expect_error(gTest(c(1, 2, 3), correct = "yates"), "2 data values")
})

test_that("input checks", {
  expect_error(gTest(1:3, 1:2), "same length")
  expect_error(gTest(c("a", "a", "a"), c("x", "y", "x")), "at least 2 levels")
  expect_error(gTest(c(3, -1, 2)), "nonnegative")
  expect_error(gTest(c(3, NA, 2)), "nonnegative")
  expect_error(gTest(c(0, 0, 0)), "at least one entry")
  expect_error(gTest(5), "at least have 2 elements")
  expect_error(gTest(c(1, 2, 3), p = c(0.5, 0.5)), "same number of elements")
  expect_error(gTest(c(1, 2, 3), p = c(1.2, -0.1, -0.1)), "non-negative")
  expect_error(gTest(array(1:8, c(2, 2, 2))), "invalid 'x'")
  expect_error(gTest(c(1, 2), correct = "foo"))
})
