
# Reference values computed independently from the 2 x k tables; the pooled
# median of warpbreaks$breaks is 26, which several observations attain, so the
# three 'ties' conventions genuinely differ here.

test_that("moodMedianTest reproduces the warpbreaks table and statistic", {

  r <- moodMedianTest(breaks ~ tension, data = warpbreaks)

  expect_equal(r$grand.median, 26)
  expect_equal(unname(r$observed["above", ]), c(12, 9, 4))
  expect_equal(unname(r$observed["below", ]), c(6, 9, 14))
  expect_equal(colnames(r$observed), c("L", "M", "H"))

  expect_equal(unname(r$statistic), 7.299310345, tolerance = 1e-9)
  expect_equal(unname(r$parameter), 2)
  expect_equal(r$p.value, 0.02600009278, tolerance = 1e-9)
  expect_equal(unname(r$estimate), c(29.5, 27, 20.5))

  expect_s3_class(r, "htest")
  expect_identical(names(r$statistic), "X-squared")
  expect_equal(sum(r$expected), sum(r$observed))
})


test_that("the ties conventions differ as documented", {

  below <- moodMedianTest(breaks ~ tension, data = warpbreaks, ties = "below")
  above <- moodMedianTest(breaks ~ tension, data = warpbreaks, ties = "above")
  drop  <- moodMedianTest(breaks ~ tension, data = warpbreaks, ties = "drop")

  expect_equal(unname(above$observed["above", ]), c(14, 10, 5))
  expect_equal(unname(drop$observed["above", ]), c(12, 9, 4))
  expect_equal(unname(drop$observed["below", ]), c(4, 8, 13))

  expect_equal(unname(above$statistic), 9.086896552, tolerance = 1e-9)
  expect_equal(unname(drop$statistic),  8.823529412, tolerance = 1e-9)

  # "drop" discards the observations sitting exactly on the pooled median
  expect_lt(sum(drop$observed), sum(below$observed))
  # the grand median itself does not depend on the convention
  expect_equal(above$grand.median, drop$grand.median)
})


test_that("the two-sample case matches a 2 x 2 chi-squared test", {

  x <- c(12, 15, 11, 18, 22, 14, 17, 13, 25, 16)
  y <- c(21, 24, 19, 28, 23, 30, 20, 27, 26, 22)
  v <- c(x, y)
  g <- rep(c("a", "b"), each = 10)

  r <- moodMedianTest(v, g)
  expect_equal(r$grand.median, 20.5)
  expect_equal(as.vector(r$observed), c(2, 8, 8, 2))

  # Yates' correction is on by default for a 2 x 2 table
  expect_equal(unname(r$statistic), 5, tolerance = 1e-12)
  expect_equal(r$p.value, 0.02534731868, tolerance = 1e-9)
  expect_match(r$method, "continuity correction", fixed = TRUE)

  r0 <- moodMedianTest(v, g, correct = FALSE)
  expect_equal(unname(r0$statistic), 7.2, tolerance = 1e-12)
  expect_equal(r0$p.value, 0.007290358092, tolerance = 1e-9)
  expect_false(grepl("continuity", r0$method))

  expect_equal(unname(r$estimate), c(15.5, 23.5))
  expect_equal(unname(chisq.test(r$observed)$statistic), unname(r$statistic))
})


test_that("correct = TRUE is ignored for more than two groups", {

  a <- moodMedianTest(breaks ~ tension, data = warpbreaks, correct = TRUE)
  b <- moodMedianTest(breaks ~ tension, data = warpbreaks, correct = FALSE)
  expect_equal(a$statistic, b$statistic)
  expect_false(grepl("continuity", a$method))
})


test_that("method = 'exact' uses Fisher's test on the same table", {

  r <- moodMedianTest(breaks ~ tension, data = warpbreaks, method = "exact")

  expect_equal(r$p.value, 0.03153901764, tolerance = 1e-8)
  expect_null(r$statistic)
  expect_null(r$parameter)
  expect_match(r$method, "Fisher's exact test", fixed = TRUE)
  # the table itself is unchanged
  expect_equal(r$observed,
               moodMedianTest(breaks ~ tension, data = warpbreaks)$observed)
})


test_that("the test is invariant under monotone transformation", {

  a <- moodMedianTest(breaks ~ tension, data = warpbreaks)
  b <- moodMedianTest(log(breaks) ~ tension, data = warpbreaks)
  expect_equal(unname(a$statistic), unname(b$statistic))
  expect_equal(a$observed, b$observed)
})


test_that("the formula method matches the default method", {

  f <- moodMedianTest(breaks ~ tension, data = warpbreaks)
  g <- moodMedianTest(warpbreaks$breaks, warpbreaks$tension)

  expect_equal(unname(f$statistic), unname(g$statistic))
  expect_equal(f$data.name, "breaks ~ tension")

  # subset is evaluated in 'data'
  s <- moodMedianTest(breaks ~ tension, data = warpbreaks, subset = wool == "A")
  expect_equal(sum(s$observed), 27)
})


test_that("missing values are removed casewise", {

  d <- warpbreaks
  d$breaks[c(1, 5, 50)] <- NA
  r <- moodMedianTest(breaks ~ tension, data = d)
  expect_equal(sum(r$observed), 51)
})


test_that("invalid arguments are rejected", {

  expect_error(moodMedianTest(1:10, rep("a", 10)), "at least two levels")
  expect_error(moodMedianTest(1:10, rep(c("a", "b"), 3)), "same length")
  expect_error(moodMedianTest(letters[1:10], rep(c("a", "b"), 5)), "numeric")
  expect_error(moodMedianTest(1:10, rep(c("a", "b"), 5), correct = NA), "'correct'")
  expect_error(moodMedianTest(1:10, rep(c("a", "b"), 5), ties = "middle"),
               "should be one of")

  # every observation identical: "drop" leaves nothing behind
  expect_error(moodMedianTest(rep(1, 10), rep(c("a", "b"), 5), ties = "drop"),
               "no observations left")
})
