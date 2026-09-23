

# tests/testthat/test-siegelTukeyTest.R

test_that("siegelTukeyTest.default returns htest object", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 3, 5, 7, 9)
  
  res <- siegelTukeyTest(x, y)
  
  expect_s3_class(res, "htest")
  expect_named(res, c(
    "statistic", "parameter", "p.value", "null.value",
    "alternative", "method", "data.name", "exact", "ties"
  ))
  expect_named(res$statistic, "W")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
  expect_equal(res$method, "Siegel-Tukey test for scale differences")
})


test_that("siegelTukeyTest requires y", {
  x <- 1:5
  
  expect_error(
    siegelTukeyTest(x),
    "'y' is missing",
    fixed = TRUE
  )
})


test_that("alternative argument is validated", {
  x <- 1:5
  y <- 6:10
  
  expect_error(
    siegelTukeyTest(x, y, alternative = "invalid")
  )
  
  expect_equal(
    siegelTukeyTest(x, y, alternative = "less")$alternative,
    "less"
  )
  
  expect_equal(
    siegelTukeyTest(x, y, alternative = "greater")$alternative,
    "greater"
  )
})


test_that("mu must be a single finite number", {
  x <- 1:5
  y <- 6:10
  
  expect_error(
    siegelTukeyTest(x, y, mu = c(0, 1)),
    "'mu' must be a single number",
    fixed = TRUE
  )
  
  expect_error(
    siegelTukeyTest(x, y, mu = Inf),
    "'mu' must be a single number",
    fixed = TRUE
  )
  
  expect_equal(
    siegelTukeyTest(x, y, mu = 2)$null.value,
    c(mu = 2)
  )
})


test_that("formula interface works for two independent groups", {
  dat <- data.frame(
    value = c(1, 2, 3, 4, 1, 3, 5, 7),
    group = rep(c("a", "b"), each = 4)
  )
  
  res_formula <- siegelTukeyTest(value ~ group, data = dat)
  res_default <- siegelTukeyTest(
    x = dat$value[dat$group == "a"],
    y = dat$value[dat$group == "b"]
  )
  
  expect_s3_class(res_formula, "htest")
  expect_equal(res_formula$p.value, res_default$p.value)
  expect_equal(res_formula$statistic, res_default$statistic)
})


test_that("formula interface rejects incorrect formula", {
  dat <- data.frame(
    value = 1:5,
    group = rep(c("a", "b"), length.out = 5)
  )
  
  expect_error(
    siegelTukeyTest(~ value, data = dat),
    "'formula' missing or incorrect",
    fixed = TRUE
  )
})


test_that("adjustMedian changes result when medians differ", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(10, 12, 14, 16, 18)
  
  res_unadjusted <- siegelTukeyTest(x, y, adjustMedian = FALSE)
  res_adjusted   <- siegelTukeyTest(x, y, adjustMedian = TRUE)
  
  expect_s3_class(res_adjusted, "htest")
  expect_false(isTRUE(all.equal(
    res_unadjusted$statistic,
    res_adjusted$statistic
  )))
})


test_that("ties are detected", {
  x <- c(1, 2, 2, 3, 4)
  y <- c(1, 2, 3, 3, 5)
  
  res <- siegelTukeyTest(x, y)
  
  expect_true(res$ties)
  expect_false(res$exact)
})


test_that("exact p-value is used for small samples without ties", {
  x <- c(1, 2, 3, 4)
  y <- c(5, 6, 7, 8)
  
  res <- siegelTukeyTest(x, y)
  
  expect_true(res$exact)
  expect_false(res$ties)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})


test_that("normal approximation is used for large samples", {
  x <- 1:60
  y <- 101:160
  
  res <- siegelTukeyTest(x, y)
  
  expect_false(res$exact)
  expect_false(res$ties)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})


test_that(".siegelTukeyRank returns expected columns and length for even n", {
  x <- c(1, 2, 3, 4, 5, 6)
  g <- c(0, 0, 0, 1, 1, 1)
  
  res <- .siegelTukeyRank(x, g)
  
  expect_s3_class(res, "data.frame")
  expect_named(res, c("sort.x", "sort.id", "unique.ranks", "raw.ranks"))
  expect_equal(nrow(res), length(x))
  expect_equal(res$sort.x, sort(x))
})


test_that(".siegelTukeyRank drops one median observation for odd n", {
  x <- c(1, 2, 3, 4, 5)
  g <- c(0, 0, 1, 1, 1)
  
  res <- .siegelTukeyRank(x, g, dropMedian = TRUE)
  
  expect_equal(nrow(res), length(x) - 1)
  expect_false(3 %in% res$sort.x)
})


test_that(".siegelTukeyRank averages ranks for tied x values", {
  x <- c(1, 2, 2, 3)
  g <- c(0, 0, 1, 1)
  
  res <- .siegelTukeyRank(x, g, dropMedian = FALSE)
  
  tied_ranks <- res$unique.ranks[res$sort.x == 2]
  
  expect_equal(length(unique(tied_ranks)), 1)
})



# reference: two-sided normal approximation with the permutation variance
# of the rank sum, computed independently from the Siegel-Tukey ranks
.refStP <- function(x, y, correct = TRUE) {
  st <- .siegelTukeyRank(c(x, y), g = rep(0:1, c(length(x), length(y))))
  r  <- st$unique.ranks
  m  <- sum(st$sort.id == 0)
  n  <- sum(st$sort.id == 1)
  N  <- m + n
  U  <- sum(r[st$sort.id == 1]) - n * (n + 1) / 2
  z  <- U - m * n / 2
  if (correct) z <- z - sign(z) * 0.5
  V  <- m * n / (N * (N - 1)) * sum((r - mean(r))^2)
  2 * pnorm(-abs(z) / sqrt(V))
}


test_that("without ties the approximation equals wilcox.test on the ST ranks", {
  
  x <- c(23, 18, 17, 25, 22, 19, 31, 26, 29, 33)
  y <- c(21, 28, 32, 30, 41, 24, 35, 34, 27, 39, 36)
  
  st <- .siegelTukeyRank(c(x, y), g = rep(0:1, c(length(x), length(y))))
  expect_false(anyDuplicated(st$sort.x) > 0)
  
  # without ties the ST ranks are a permutation of 1..N, so ranking them
  # again changes nothing
  ref <- wilcox.test(st$unique.ranks[st$sort.id == 1],
                     st$unique.ranks[st$sort.id == 0],
                     exact = FALSE)$p.value
  
  expect_equal(siegelTukeyTest(x, y, exact = FALSE)$p.value, ref)
})


test_that("with ties the variance is the permutation variance of the ST ranks", {
  
  # example from the DescTools issue: aligning the medians creates ties
  x <- c(23, 18, 17, 25, 22, 19, 31, 26, 29, 33)
  y <- c(21, 28, 32, 30, 41, 24, 35, 34, 27, 39, 36)
  
  res <- suppressWarnings(
    siegelTukeyTest(x, y, adjustMedian = TRUE, exact = FALSE))
  
  expect_true(res$ties)
  expect_equal(res$p.value,
               .refStP(x - (median(x) - median(y)), y))
})


test_that("tie groups sharing a mean rank are not merged", {
  
  # the two smallest (ranks 1, 4) and the two largest (ranks 3, 2) both
  # average to 2.5; table(unique.ranks) counted them as one group of four
  x <- c(1, 3, 5, 8)
  y <- c(1, 4, 6, 8)
  
  st <- .siegelTukeyRank(c(x, y), g = rep(0:1, each = 4))
  expect_equal(sum(st$unique.ranks == 2.5), 4L)
  
  for (cc in c(TRUE, FALSE))
    expect_equal(siegelTukeyTest(x, y, exact = FALSE, correct = cc)$p.value,
                 .refStP(x, y, correct = cc))
})
