
test_that("coxStuartTest counts signs and pairs correctly", {

  r <- coxStuartTest(1:12)

  expect_s3_class(r, "htest")
  expect_equal(unname(r$statistic), 6)
  expect_equal(unname(r$parameter), 6)
  expect_identical(names(r$statistic), "S")
  expect_equal(r$p.value, 0.03125)
  expect_equal(unname(r$estimate), 1)
  expect_equal(r$method, "Cox-Stuart trend test")
  expect_equal(r$data.name, "1:12")
})


test_that("the alternatives point in the documented direction", {

  expect_equal(coxStuartTest(1:12, alternative = "increasing")$p.value, 0.015625)
  expect_equal(coxStuartTest(1:12, alternative = "decreasing")$p.value, 1)

  # a decreasing series mirrors the result
  expect_equal(coxStuartTest(12:1, alternative = "decreasing")$p.value, 0.015625)
  expect_equal(coxStuartTest(12:1, alternative = "increasing")$p.value, 1)
  expect_equal(coxStuartTest(12:1)$p.value, coxStuartTest(1:12)$p.value)
  expect_equal(unname(coxStuartTest(12:1)$statistic), 0)
})


test_that("an odd-length series drops its central observation", {

  r <- coxStuartTest(1:9)
  expect_equal(unname(r$statistic), 4)
  expect_equal(unname(r$parameter), 4)
  expect_equal(r$p.value, 0.125)

  # the central value itself does not influence the result
  x <- c(1, 2, 3, 4, 999, 6, 7, 8, 9)
  expect_equal(coxStuartTest(x)$p.value, r$p.value)
  expect_equal(unname(coxStuartTest(x)$statistic), unname(r$statistic))
})


test_that("tied pairs are dropped rather than counted", {

  # pairs are (1,1), (2,5), (3,3), (4,9): two ties, two increases
  r <- coxStuartTest(c(1, 2, 3, 4, 1, 5, 3, 9))
  expect_equal(unname(r$statistic), 2)
  expect_equal(unname(r$parameter), 2)
  expect_equal(r$p.value, 0.5)
})


test_that("the test is invariant under monotone transformation", {

  set.seed(3)
  x <- cumsum(abs(rnorm(30))) + 1
  expect_equal(coxStuartTest(x)$p.value, coxStuartTest(log(x))$p.value)
  expect_equal(unname(coxStuartTest(x)$statistic),
               unname(coxStuartTest(exp(x / 100))$statistic))
})


test_that("missing values are removed before the series is split", {

  x <- c(1:6, NA, 7:12)
  expect_equal(coxStuartTest(x)$p.value, coxStuartTest(1:12)$p.value)
})


test_that("invalid arguments are rejected", {

  expect_error(coxStuartTest(letters[1:10]), "numeric")
  expect_error(coxStuartTest(1:3), "at least 4")
  expect_error(coxStuartTest(c(1, 2, NA, NA, 5)), "at least 4")
  expect_error(coxStuartTest(1:12, alternative = "up"), "should be one of")

  # a perfectly periodic series pairs every value with itself
  expect_error(coxStuartTest(rep(1:3, 4)), "no directional information")
})
