
tol <- 1e-4
set.seed(1); x <- rnorm(50, mean = 5, sd = 2)

test_that("zTest: one-sample returns htest", {
  expect_s3_class(zTest(x, sd_pop = 2), "htest")
})

test_that("zTest: p.value in [0,1]", {
  res <- zTest(x, sd_pop = 2)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("zTest: H0 true (mu=5) gives large p", {
  set.seed(42); xh <- rnorm(200, mean = 5, sd = 2)
  expect_gt(zTest(xh, mu = 5, sd_pop = 2)$p.value, 0.05)
})

test_that("zTest: H0 false gives small p", {
  set.seed(42); xh <- rnorm(200, mean = 5, sd = 2)
  expect_lt(zTest(xh, mu = 0, sd_pop = 2)$p.value, 0.001)
})

test_that("zTest: statistic matches manual z calculation", {
  res <- zTest(x, mu = 0, sd_pop = 2)
  z_manual <- (mean(x) - 0) / (2 / sqrt(length(x)))
  expect_equal(unname(res$statistic), z_manual, tolerance = tol)
})

test_that("zTest: CI contains true mean for large n", {
  set.seed(1); xh <- rnorm(10000, mean = 5, sd = 2)
  res <- zTest(xh, mu = 5, sd_pop = 2)
  expect_true(res$conf.int[1] <= 5 && 5 <= res$conf.int[2])
})

test_that("zTest: alternative='less' p <= two.sided p", {
  # use mu > mean(x) so zstat < 0 and the left-tail p-value is small
  res2 <- zTest(x, mu = 10, sd_pop = 2, alternative = "two.sided")
  resL <- zTest(x, mu = 10, sd_pop = 2, alternative = "less")
  expect_lte(resL$p.value, res2$p.value + tol)
})

test_that("zTest: two-sample works", {
  set.seed(2); y <- rnorm(40, mean = 3, sd = 1.5)
  res <- zTest(x, y, sd_pop = 2)
  expect_s3_class(res, "htest")
})

# zTest() --------------------------------------------------------------------

x <- c(102.1, 98.4, 105.3, 99.8, 101.2, 97.6, 103.9, 100.5)
y <- c( 96.2, 99.1,  94.8, 98.3,  97.7, 95.5)

test_that("one-sample: statistic, p-value and interval by hand", {
  r <- zTest(x, mu = 99, sd_pop = 3)
  se <- 3 / sqrt(length(x))
  z <- (mean(x) - 99) / se
  
  expect_s3_class(r, "htest")
  expect_equal(r$statistic, c(z = z))
  expect_equal(r$p.value, 2 * pnorm(-abs(z)))
  expect_equal(as.vector(r$conf.int), mean(x) + c(-1, 1) * qnorm(0.975) * se)
  expect_equal(attr(r$conf.int, "conf.level"), 0.95)
  expect_equal(r$estimate, c("mean of x" = mean(x)))
  expect_equal(r$null.value, c(mean = 99))
  expect_equal(r$stderr, se)
  expect_identical(r$method, "One Sample z-test")
  expect_identical(r$data.name, "x")
})

test_that("one-sided alternatives", {
  se <- 3 / sqrt(length(x))
  z <- (mean(x) - 99) / se
  
  g <- zTest(x, mu = 99, sd_pop = 3, alternative = "greater", conf.level = 0.9)
  expect_equal(g$p.value, pnorm(z, lower.tail = FALSE))
  expect_equal(as.vector(g$conf.int), c(mean(x) - qnorm(0.9) * se, Inf))
  
  l <- zTest(x, mu = 99, sd_pop = 3, alternative = "l")
  expect_equal(l$p.value, pnorm(z))
  expect_equal(as.vector(l$conf.int), c(-Inf, mean(x) + qnorm(0.95) * se))
  
  expect_equal(g$p.value + l$p.value, 1)
})

test_that("two-sample", {
  r <- zTest(x, y, sd_pop = 2, mu = 1)
  se <- 2 * sqrt(1 / length(x) + 1 / length(y))
  z <- (mean(x) - mean(y) - 1) / se
  
  expect_equal(r$statistic, c(z = z))
  expect_equal(r$p.value, 2 * pnorm(-abs(z)))
  expect_equal(as.vector(r$conf.int),
               mean(x) - mean(y) + c(-1, 1) * qnorm(0.975) * se)
  expect_equal(unname(r$estimate), c(mean(x), mean(y)))
  expect_named(r$estimate, c("mean of x", "mean of y"))
  expect_equal(r$null.value, c("difference in means" = 1))
  expect_identical(r$method, "Two Sample z-test")
  expect_identical(r$data.name, "x and y")
})

test_that("paired test equals one-sample test on the differences", {
  a <- c(44.5, 55, 52.5, 50.2, 45.3, 46.1, 52.1, 50.5, 50.6, 49.2)
  b <- c(44.9, 54.8, 55.6, 55.2, 55.6, 47.7, 53, 49.1, 52.3, 50.7)
  p <- zTest(a, b, sd_pop = 3, paired = TRUE)
  o <- zTest(a - b, sd_pop = 3)
  
  expect_equal(p$statistic, o$statistic)
  expect_equal(p$p.value, o$p.value)
  expect_equal(p$conf.int, o$conf.int)
  expect_identical(p$method, "Paired z-test")
  expect_named(p$estimate, "mean of the differences")
  expect_named(p$null.value, "difference in means")
})

test_that("missing values: pairwise for paired, separately otherwise", {
  a <- c(1.2, NA, 3.1, 4.0, 5.2, 2.2)
  b <- c(2.0, 2.5, NA, 4.4, 5.9, 2.0)
  
  ok <- complete.cases(a, b)
  expect_equal(zTest(a, b, sd_pop = 1, paired = TRUE)$statistic,
               zTest(a[ok] - b[ok], sd_pop = 1)$statistic)
  expect_equal(zTest(a, b, sd_pop = 1)$statistic,
               zTest(na.omit(a), na.omit(b), sd_pop = 1)$statistic)
  expect_equal(zTest(a, sd_pop = 1)$estimate[[1]], mean(a, na.rm = TRUE))
})

test_that("argument checks", {
  expect_error(zTest(x, mu = c(1, 2), sd_pop = 1), "'mu' must be a single number")
  expect_error(zTest(x, mu = NA, sd_pop = 1), "'mu' must be a single number")
  expect_error(zTest(x, sd_pop = 1, conf.level = 1.5), "'conf.level'")
  expect_error(zTest(x, sd_pop = 1, conf.level = NA), "'conf.level'")
  expect_error(zTest(x, sd_pop = 1, paired = TRUE), "'y' is missing")
  expect_error(zTest(c(NA, NA), 1:3, sd_pop = 1), "not enough 'x' observations")
  expect_error(zTest(1:3, NA, sd_pop = 1), "not enough 'y' observations")
  expect_error(zTest(c(5, 5, 5), sd_pop = 0), "'sd_pop'")
  expect_error(zTest(c(5, 5), c(5, 5), sd_pop = 0), "'sd_pop'")
  expect_error(zTest(x, sd_pop = 1, alternative = "foo"))
})

test_that("formula interface equals the default method", {
  f <- zTest(extra ~ group, data = sleep, sd_pop = 2)
  d <- with(sleep, zTest(extra[group == 1], extra[group == 2], sd_pop = 2))
  
  expect_equal(f$statistic, d$statistic)
  expect_equal(f$p.value, d$p.value)
  expect_equal(f$conf.int, d$conf.int)
  expect_equal(f$estimate, d$estimate)
  
  # regression: the name was read from a component resolveFormula() does
  # not return, which deleted data.name from the result
  expect_type(f$data.name, "character")
  expect_match(f$data.name, "extra")
})

test_that("formula interface passes '...' on", {
  f <- zTest(extra ~ group, data = sleep, sd_pop = 2,
             alternative = "less", mu = -1)
  d <- with(sleep, zTest(extra[group == 1], extra[group == 2], sd_pop = 2,
                         alternative = "less", mu = -1))
  expect_equal(f$p.value, d$p.value)
  expect_identical(f$alternative, "less")
})

test_that("formula interface honours subset", {
  # regression: do.call() evaluated the subset expression outside the data
  f <- zTest(extra ~ group, data = sleep, subset = ID != "1", sd_pop = 2)
  s <- sleep[sleep$ID != "1", ]
  d <- zTest(s$extra[s$group == 1], s$extra[s$group == 2], sd_pop = 2)
  expect_equal(f$statistic, d$statistic)
})

test_that("formula interface rejects one-sided formulas", {
  expect_error(zTest(~ extra, data = sleep, sd_pop = 2),
               "'formula' missing or incorrect")
})


test_that("formula interface rejects paired, also abbreviated", {

  msg <- "formula interface"
  expect_error(zTest(extra ~ group, data = sleep, sd_pop = 2, paired = TRUE), msg)
  # abbreviated: would be matched partially by zTest.default()
  expect_error(zTest(extra ~ group, data = sleep, sd_pop = 2, pair = TRUE), msg)
  expect_error(zTest(extra ~ group, data = sleep, sd_pop = 2, paired = NA), msg)

  # paired = FALSE stays allowed
  expect_s3_class(zTest(extra ~ group, data = sleep, sd_pop = 2,
                        paired = FALSE), "htest")
})


test_that("guards on sd_pop, mu, paired and the data", {

  expect_error(zTest(x), "'sd_pop'.*required")
  expect_error(zTest(x, sd_pop = -1), "'sd_pop'")
  expect_error(zTest(x, sd_pop = NA), "'sd_pop'")
  expect_error(zTest(x, sd_pop = Inf), "'sd_pop'")
  expect_error(zTest(x, sd_pop = c(1, 2)), "'sd_pop'")
  expect_error(zTest(x, sd_pop = "2"), "'sd_pop'")

  expect_error(zTest(x, mu = Inf, sd_pop = 1), "'mu'")
  expect_error(zTest(x, mu = "1", sd_pop = 1), "'mu'")

  expect_error(zTest(x, x, sd_pop = 1, paired = NA), "'paired'")
  expect_error(zTest(x, x, sd_pop = 1, paired = c(TRUE, FALSE)), "'paired'")
  expect_error(zTest(x, x[-1], sd_pop = 1, paired = TRUE), "same length")

  expect_error(zTest(letters, sd_pop = 1), "numeric")
})


test_that("constant data are fine: the standard error does not depend on them", {
  r <- zTest(c(5, 5, 5, 5), mu = 4, sd_pop = 2)
  expect_equal(unname(r$statistic), (5 - 4) / (2 / 2))
})


test_that("non-finite values are removed, as documented", {
  expect_equal(zTest(c(x, Inf, -Inf), sd_pop = 3)$statistic,
               zTest(x, sd_pop = 3)$statistic)
  expect_equal(zTest(c(x, Inf), c(y, NA), sd_pop = 3)$statistic,
               zTest(x, y, sd_pop = 3)$statistic)
  expect_equal(zTest(c(x[1:6], Inf), c(y, 1), sd_pop = 3, paired = TRUE)$statistic,
               zTest(x[1:6], y, sd_pop = 3, paired = TRUE)$statistic)
})


test_that("a single observation suffices when sd_pop is known", {
  r <- zTest(3, mu = 1, sd_pop = 2)
  expect_equal(unname(r$statistic), (3 - 1) / 2)
  expect_equal(unname(r$conf.int[1:2]), 3 + c(-1, 1) * qnorm(0.975) * 2)

  r2 <- zTest(3, 1, sd_pop = 2)
  expect_equal(unname(r2$statistic), (3 - 1) / (2 * sqrt(2)))

  expect_error(zTest(numeric(0), sd_pop = 1), "not enough 'x' observations")
})
