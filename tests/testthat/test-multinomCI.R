library(testthat)
library(lumen)

x <- c(A = 20, B = 15, C = 25)

test_that("multinomCI: returns matrix with 3 columns", {
  res <- multinomCI(x)
  expect_true(is.matrix(res))
  expect_equal(ncol(res), 3L)
})

test_that("multinomCI: column names est/lci/uci", {
  expect_equal(colnames(multinomCI(x)), c("est", "lci", "uci"))
})

test_that("multinomCI: est sums to 1", {
  res <- multinomCI(x)
  expect_equal(sum(res[,"est"]), 1, tolerance = 1e-10)
})

test_that("multinomCI: lci <= est <= uci for all categories", {
  res <- multinomCI(x)
  expect_true(all(res[,"lci"] <= res[,"est"]))
  expect_true(all(res[,"est"] <= res[,"uci"]))
})

test_that("multinomCI: lci >= 0 and uci <= 1", {
  res <- multinomCI(x)
  expect_true(all(res[,"lci"] >= 0))
  expect_true(all(res[,"uci"] <= 1))
})

test_that("multinomCI: est = proportions", {
  res <- multinomCI(x)
  expect_equal(unname(res[,"est"]), unname(x / sum(x)), tolerance = 1e-10)
})

test_that("multinomCI: wider CI with higher conf.level", {
  w95 <- mean(multinomCI(x, conf.level = 0.95)[,"uci"] -
                multinomCI(x, conf.level = 0.95)[,"lci"])
  w99 <- mean(multinomCI(x, conf.level = 0.99)[,"uci"] -
                multinomCI(x, conf.level = 0.99)[,"lci"])
  expect_gt(w99, w95)
})

test_that("multinomCI: all methods return valid result", {
  methods <- c("sison-glaz","goodman","wald","waldcc","wilson")
  for (m in methods) {
    res <- multinomCI(x, method = m)
    expect_true(all(res[,"lci"] >= 0), label = paste(m, "lci>=0"))
    expect_true(all(res[,"uci"] <= 1), label = paste(m, "uci<=1"))
  }
})

test_that("multinomCI: sides='left' gives uci=1", {
  res <- multinomCI(x, sides = "left")
  expect_true(all(res[,"uci"] == 1))
})

test_that("multinomCI: sides='right' gives lci=0", {
  res <- multinomCI(x, sides = "right")
  expect_true(all(res[,"lci"] == 0))
})


# -- added --------------------------------------------------------------------

xs <- list(c(56, 72, 73, 59, 62, 87, 58), c(A = 20, B = 15, C = 25),
           c(1, 10, 3, 0), c(5, 5))

test_that("multinomCI sison-glaz equals MultinomialCI::multinomialCI", {
  skip_if_not_installed("MultinomialCI")
  for (x in xs) for (cl in c(0.9, 0.95, 0.99)) {
    ref <- MultinomialCI::multinomialCI(x, alpha = 1 - cl)
    expect_equal(unname(multinomCI(x, conf.level = cl)[, c("lci", "uci")]),
                 unname(ref), info = paste(c(x, cl), collapse = " "))
  }
})

# score interval of prop.test() at a given level
wilsonCI <- function(x, n, cl)
  unname(t(vapply(x, function(xi) prop.test(xi, n, conf.level = cl,
                                     correct = FALSE)$conf.int, numeric(2))))

test_that("goodman = Wilson score interval at the Bonferroni level 1 - alpha/k", {
  for (x in xs) {
    k <- length(x); n <- sum(x)
    r <- multinomCI(x, method = "goodman", conf.level = 0.95)
    expect_equal(unname(r[, c("lci", "uci")]),
                 wilsonCI(x, n, 1 - 0.05 / k), tolerance = 1e-10)
  }
})

test_that("wilson and qh are score intervals at level cl resp. chi2(k-1)", {
  x <- c(A = 20, B = 15, C = 25); n <- 60; k <- 3
  r <- multinomCI(x, method = "wilson")
  expect_equal(unname(r[, c("lci", "uci")]), wilsonCI(x, n, 0.95),
               tolerance = 1e-10)
  # Quesenberry-Hurst: z^2 = qchisq(cl, k - 1)
  clq <- 2 * pnorm(sqrt(qchisq(0.95, k - 1))) - 1
  r <- multinomCI(x, method = "qh")
  expect_equal(unname(r[, c("lci", "uci")]), wilsonCI(x, n, clq),
               tolerance = 1e-10)
})

test_that("wald, waldcc and fs by hand", {
  x <- c(A = 20, B = 15, C = 25); n <- 60; p <- x / n
  z <- qnorm(0.975)
  w <- z * sqrt(p * (1 - p) / n)
  expect_equal(unname(multinomCI(x, method = "wald")[, 2:3]),
               unname(cbind(p - w, p + w)))
  expect_equal(unname(multinomCI(x, method = "waldcc")[, 2:3]),
               unname(cbind(p - w - 1 / (2 * n), p + w + 1 / (2 * n))))
  expect_equal(unname(multinomCI(x, method = "fs")[, 2:3]),
               unname(cbind(p - z / (2 * sqrt(n)), p + z / (2 * sqrt(n)))))
})

test_that("cplus1 is the symmetric Sison-Glaz interval widened by 1/n", {
  x <- c(56, 72, 73, 59, 62, 87, 58); n <- sum(x); p <- x / n
  sg <- multinomCI(x)
  cp <- multinomCI(x, method = "cplus1")
  # lower bound: p - c/n  ->  p - c/n - 1/n
  expect_equal(unname(cp[, "lci"]), unname(pmax(0, sg[, "lci"] - 1 / n)))
  expect_equal(unname(cp[, "uci"] - p), unname(p - cp[, "lci"]))
})

test_that("bounds are clipped to [0, 1] and names are kept", {
  x <- c(A = 1, B = 0, C = 29)
  for (m in c("wald", "waldcc", "fs", "goodman", "wilson", "qh",
              "sison-glaz", "cplus1")) {
    r <- multinomCI(x, method = m)
    expect_identical(rownames(r), c("A", "B", "C"), info = m)
    expect_true(all(r[, "lci"] >= 0 & r[, "uci"] <= 1), info = m)
  }
  expect_equal(unname(multinomCI(x, method = "wald")["B", ]), c(0, 0, 0))
})

test_that("one-sided = end of the two-sided interval at 2*cl - 1", {
  x <- c(A = 20, B = 15, C = 25)
  for (m in c("goodman", "wilson", "sison-glaz")) {
    two <- multinomCI(x, method = m, conf.level = 0.9)
    lft <- multinomCI(x, method = m, conf.level = 0.95, sides = "left")
    rgt <- multinomCI(x, method = m, conf.level = 0.95, sides = "r")
    expect_equal(lft[, "lci"], two[, "lci"], info = m)
    expect_equal(rgt[, "uci"], two[, "uci"], info = m)
  }
  expect_error(multinomCI(x, sides = "left", conf.level = 0.5), "above 0.5")
})

test_that("input checks", {
  expect_error(multinomCI(c(5, -1, 3)), "non-negative")
  expect_error(multinomCI(c(5, NA, 3)), "non-negative")
  expect_error(multinomCI(c(5, Inf, 3)), "non-negative")
  expect_error(multinomCI(c("a", "b")), "non-negative")
  expect_error(multinomCI(5), "at least 2")
  expect_error(multinomCI(c(0, 0, 0)), "positive")
  expect_error(multinomCI(c(5, 3), conf.level = 1), "'conf.level'")
  expect_error(multinomCI(c(5, 3), conf.level = c(0.9, 0.95)), "'conf.level'")
  expect_error(multinomCI(c(5, 3), method = "foo"))
  expect_error(multinomCI(c(5, 3), sides = "foo"))
})
