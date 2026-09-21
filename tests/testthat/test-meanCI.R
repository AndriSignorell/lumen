set.seed(1)
x <- rnorm(50, mean = 5, sd = 2)

# fixed data with two outliers; reference values for the trimmed intervals
# from scipy.stats.mstats.trimmed_mean_ci (Tukey-McLaughlin, df = n - 2g - 1)
xr <- c(2.1, 3.4, 1.9, 5.6, 4.4, 3.3, 12.8, 2.7, 3.9, 4.1,
        0.4, 3.6, 2.2, 4.8, 3.0, 9.5, 3.1, 2.9, 4.0, 3.7)

# winsorized variance of the trimmed mean, written out independently
seTrim2 <- function(x, trim) {
  n <- length(x)
  g <- floor(trim * n)
  s <- sort(x)
  var(pmin(pmax(x, s[g + 1]), s[n - g])) / ((1 - 2 * trim)^2 * n)
}

bootLimits <- function(ci, type) {
  slot <- c(norm = "normal", basic = "basic", perc = "percent",
            bca = "bca", stud = "student")[[type]]
  m <- ci[[slot]]
  if (type == "norm") m[2:3] else m[4:5]
}


# -- classic ------------------------------------------------------------------

test_that("returns est/lci/uci with est = mean(x)", {
  res <- meanCI(x)
  expect_named(res, c("est", "lci", "uci"))
  expect_equal(res[["est"]], mean(x))
  expect_true(res[["lci"]] <= res[["est"]] && res[["est"]] <= res[["uci"]])
})


test_that("classic interval equals t.test() for all sides and levels", {
  for (cl in c(0.8, 0.95, 0.99)) {
    expect_equal(unname(meanCI(x, conf.level = cl)[c("lci", "uci")]),
                 as.vector(t.test(x, conf.level = cl)$conf.int))
    expect_equal(meanCI(x, conf.level = cl, sides = "left")[["lci"]],
                 t.test(x, conf.level = cl, alternative = "greater")$conf.int[1])
    expect_equal(meanCI(x, conf.level = cl, sides = "right")[["uci"]],
                 t.test(x, conf.level = cl, alternative = "less")$conf.int[2])
  }
  expect_identical(meanCI(x, sides = "left")[["uci"]], Inf)
  expect_identical(meanCI(x, sides = "r")[["lci"]], -Inf)   # partial matching
})


test_that("width shrinks with lower level and larger n", {
  w <- function(r) r[["uci"]] - r[["lci"]]
  expect_gt(w(meanCI(x, conf.level = 0.99)), w(meanCI(x)))
  set.seed(1)
  expect_gt(w(meanCI(rnorm(50))), w(meanCI(rnorm(500))))
})


test_that("known sd gives the z-interval", {
  res <- meanCI(x, sd = 2, conf.level = 0.9)
  expect_equal(unname(res[c("lci", "uci")]),
               mean(x) + c(-1, 1) * qnorm(0.95) * 2 / sqrt(length(x)))
  expect_equal(meanCI(x, sd = 2, sides = "left")[["lci"]],
               mean(x) - qnorm(0.95) * 2 / sqrt(length(x)))
})


test_that("trimmed mean: Tukey-McLaughlin interval", {
  r10 <- meanCI(xr, trim = 0.1)
  expect_equal(unname(r10), c(3.55, 2.862014258876056, 4.237985741123944),
               tolerance = 1e-10)
  r20 <- meanCI(xr, trim = 0.2)
  expect_equal(unname(r20), c(3.5083333333333333, 2.947647728364524,
                              4.0690189383021425), tolerance = 1e-10)
  expect_equal(unname(meanCI(xr, trim = 0.2, conf.level = 0.9)[c("lci", "uci")]),
               c(3.050844212421926, 3.9658224542447407), tolerance = 1e-10)

  # same thing written out, unequal n and a trim that does not divide n
  set.seed(4)
  z <- rt(37, df = 3)
  g <- floor(0.15 * 37)
  expect_equal(unname(meanCI(z, trim = 0.15)[c("lci", "uci")]),
               mean(z, trim = 0.15) + c(-1, 1) * qt(0.975, 37 - 2 * g - 1) *
                 sqrt(seTrim2(z, 0.15)))

  # one-sided
  expect_equal(meanCI(xr, trim = 0.2, sides = "left")[["lci"]],
               meanCI(xr, trim = 0.2, conf.level = 0.9)[["lci"]])
})


test_that(".winvar(): winsorized variance and df", {
  wv <- .winvar(xr, 0.2)
  expect_named(wv, c("var", "DF"))
  expect_equal(wv[["var"]], seTrim2(xr, 0.2) * (1 - 2 * 0.2)^2 * length(xr))
  expect_equal(wv[["DF"]], 20 - 2 * 4 - 1)
  expect_equal(.winvar(xr, 0)[["var"]], var(xr))
})


test_that("missing values", {
  expect_equal(meanCI(c(x, NA), na.rm = TRUE), meanCI(x))
  expect_true(all(is.na(meanCI(c(x, NA)))))
})


test_that("input validation", {
  expect_error(meanCI(letters), "'x' must be numeric")
  expect_error(meanCI(1), "at least two")
  expect_error(meanCI(c(1, NA), na.rm = TRUE), "at least two")
  expect_error(meanCI(x, trim = -0.1), "'trim' must be")
  expect_error(meanCI(x, trim = 0.5), "'trim' must be")
  expect_error(meanCI(x, trim = c(0.1, 0.2)), "'trim' must be")
  expect_error(meanCI(x, trim = "a"), "'trim' must be")
  expect_error(meanCI(x, sides = "both"))
  expect_error(meanCI(x, method = "exact"))
})


# -- bootstrap ----------------------------------------------------------------
# wiring test: meanCI(method = "boot") must equal boot::boot + boot::boot.ci
# with the same statistic and the same RNG stream

test_that("boot, untrimmed: all five types equal boot.ci()", {
  stat <- function(d, i) c(mean(d[i]), var(d[i]) / length(i))
  for (type in c("norm", "basic", "perc", "bca", "stud")) {
    set.seed(1)
    res <- meanCI(x, method = "boot", type = type, R = 999, parallel = "no")
    set.seed(1)
    b  <- boot::boot(x, stat, R = 999, parallel = "no")
    ci <- boot::boot.ci(b, conf = 0.95, type = type)
    expect_equal(unname(res), c(mean(x), bootLimits(ci, type)),
                 label = paste("type =", type))
  }
})


test_that("boot, trimmed: statistic uses the resampled winsorized variance", {
  stat <- function(d, i) c(mean(d[i], trim = 0.1), seTrim2(d[i], 0.1))
  for (type in c("perc", "stud")) {
    set.seed(2)
    res <- meanCI(xr, trim = 0.1, method = "boot", type = type, R = 999,
                  parallel = "no")
    set.seed(2)
    b  <- boot::boot(xr, stat, R = 999, parallel = "no")
    ci <- boot::boot.ci(b, conf = 0.95, type = type)
    expect_equal(unname(res), c(mean(xr, trim = 0.1), bootLimits(ci, type)),
                 label = paste("type =", type))
  }
})


test_that("boot: sides and conf.level", {
  set.seed(3)
  two <- meanCI(x, method = "boot", type = "perc", R = 999, parallel = "no",
                conf.level = 0.9)
  set.seed(3)
  left <- meanCI(x, method = "boot", type = "perc", R = 999, parallel = "no",
                 sides = "left")
  expect_equal(left[c("est", "lci")], two[c("est", "lci")])
  expect_identical(left[["uci"]], Inf)
})
