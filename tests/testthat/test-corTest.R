

set.seed(42)
X <- matrix(rnorm(200), nrow = 50, ncol = 4)

test_that("corTest: returns list with cor/pValue/n", {
  res <- corTest(X)
  expect_named(res, c("cor", "pValue", "n"))
})

test_that("corTest: cor matrix dimensions", {
  res <- corTest(X)
  expect_equal(dim(res$cor), c(4L, 4L))
})

test_that("corTest: diagonal of cor is 1", {
  res <- corTest(X)
  expect_equal(diag(res$cor), rep(1, 4), tolerance = 1e-10)
})

test_that("corTest: diagonal of pValue is NA", {
  res <- corTest(X)
  expect_true(all(is.na(diag(res$pValue))))
})

test_that("corTest: pValues >= 0 (off-diagonal)", {
  # Note: the 2*pbeta() formula can yield p > 1 for weak correlations;
  # the source caps via pmin(..., 1) after the fix in corTest.R
  res <- corTest(X)
  p_off <- res$pValue[!is.na(res$pValue)]
  expect_true(all(p_off >= 0))
})

test_that("corTest: cor matrix is symmetric", {
  res <- corTest(X)
  expect_equal(res$cor, t(res$cor))
})

test_that("corTest: pValue matrix is symmetric", {
  res <- corTest(X)
  p <- res$pValue
  p[is.na(p)] <- 0
  expect_equal(p, t(p))
})

test_that("corTest: n matrix has correct sample size", {
  res <- corTest(X)
  expect_true(all(res$n == 50L))
})

test_that("corTest: cor agrees with cor.test on single pair", {
  res <- corTest(X)
  ref <- cor.test(X[,1], X[,2])
  expect_equal(res$cor[1,2],    ref$estimate[[1]], tolerance = 1e-10)
  expect_equal(res$pValue[1,2], ref$p.value,       tolerance = 1e-6)
})

test_that("corTest: pValues in [0,1] (off-diagonal)", {
  res   <- corTest(X)
  p_off <- res$pValue[!is.na(res$pValue)]
  expect_true(all(p_off >= 0 & p_off <= 1))
})

test_that("corTest: triangle='upper' NAs lower triangle", {
  res <- corTest(X, triangle = "upper")
  expect_true(all(is.na(res$cor[lower.tri(res$cor)])))
})

test_that("corTest: triangle='lower' NAs upper triangle", {
  res <- corTest(X, triangle = "lower")
  expect_true(all(is.na(res$cor[upper.tri(res$cor)])))
})

test_that("corTest: maxPValue filters non-significant", {
  res <- corTest(X, maxPValue = 0.05)
  # NAs where p > 0.05 (ignoring diagonal NAs)
  p <- corTest(X)$pValue
  should_be_na <- !is.na(p) & p > 0.05
  expect_true(all(is.na(res$cor[should_be_na])))
})

test_that("corTest: NA values handled pairwise", {
  Xna <- X
  Xna[1:5, 1] <- NA
  res <- corTest(Xna)
  expect_equal(res$n[1,2], 45L)
  expect_equal(res$n[2,3], 50L)
})


