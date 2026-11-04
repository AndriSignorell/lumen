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
  methods <- c("sison-glaz","goodman","wald","wald-cc","wilson")
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
  # reference values from MultinomialCI 1.2, multinomialCI(x, alpha = 1 - cl);
  # columns lci, uci; one matrix per x (outer) and cl in 0.90/0.95/0.99 (inner)
  ref <- list(
    matrix(c(0.0835117773, 0.1177730193, 0.1199143469, 0.0899357602, 0.0963597430, 0.1498929336, 0.0877944325,
             0.1601075553, 0.1943687972, 0.1965101249, 0.1665315381, 0.1729555210, 0.2264887116, 0.1643902105), ncol = 2),
    matrix(c(0.0792291221, 0.1134903640, 0.1156316916, 0.0856531049, 0.0920770878, 0.1456102784, 0.0835117773,
             0.1643557340, 0.1986169760, 0.2007583036, 0.1707797169, 0.1772036998, 0.2307368904, 0.1686383893), ncol = 2),
    matrix(c(0.0706638116, 0.1049250535, 0.1070663812, 0.0770877944, 0.0835117773, 0.1370449679, 0.0749464668,
             0.1733812016, 0.2076424435, 0.2097837712, 0.1798051844, 0.1862291673, 0.2397623579, 0.1776638568), ncol = 2),
    matrix(c(0.2333333333, 0.1500000000, 0.3166666667,
             0.4630515591, 0.3797182258, 0.5463848924), ncol = 2),
    matrix(c(0.2166666667, 0.1333333333, 0.3000000000,
             0.4809808206, 0.3976474873, 0.5643141539), ncol = 2),
    matrix(c(0.1833333333, 0.1000000000, 0.2666666667,
             0.5152597643, 0.4319264310, 0.5985930977), ncol = 2),
    matrix(c(0.0000000000, 0.5714285714, 0.0714285714, 0.0000000000,
             0.2896647163, 0.9325218592, 0.4325218592, 0.2182361449), ncol = 2),
    matrix(c(0.0000000000, 0.5714285714, 0.0714285714, 0.0000000000,
             0.3406616714, 0.9835188142, 0.4835188142, 0.2692330999), ncol = 2),
    matrix(c(0.0000000000, 0.5000000000, 0.0000000000, 0.0000000000,
             0.3679571860, 1.0000000000, 0.5108143289, 0.2965286146), ncol = 2),
    matrix(c(0.3000000000, 0.3000000000,
             0.7685378195, 0.7685378195), ncol = 2),
    matrix(c(0.3000000000, 0.3000000000,
             0.8607856120, 0.8607856120), ncol = 2),
    matrix(c(0.2000000000, 0.2000000000,
             0.9267463609, 0.9267463609), ncol = 2) 
  )
  i <- 0L
  for (x in xs) for (cl in c(0.90, 0.95, 0.99)) {
    i <- i + 1L
    expect_equal(unname(multinomCI(x, conf.level = cl)[, c("lci", "uci")]),
                 ref[[i]], tolerance = 1e-9,
                 info = paste(c(x, cl), collapse = " "))
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
  r <- multinomCI(x, method = "quesenberry-hurst")
  expect_equal(unname(r[, c("lci", "uci")]), wilsonCI(x, n, clq),
               tolerance = 1e-10)
})

test_that("wald, wald-cc and fs by hand", {
  x <- c(A = 20, B = 15, C = 25); n <- 60; p <- x / n
  z <- qnorm(0.975)
  w <- z * sqrt(p * (1 - p) / n)
  expect_equal(unname(multinomCI(x, method = "wald")[, 2:3]),
               unname(cbind(p - w, p + w)))
  expect_equal(unname(multinomCI(x, method = "wald-cc")[, 2:3]),
               unname(cbind(p - w - 1 / (2 * n), p + w + 1 / (2 * n))))
  expect_equal(unname(multinomCI(x, method = "fitzpatrick-scott")[, 2:3]),
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
  for (m in c("wald", "wald-cc", "fitzpatrick-scott", "goodman", "wilson", "quesenberry-hurst",
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
