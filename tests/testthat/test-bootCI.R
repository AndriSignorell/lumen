# bootCI() -------------------------------------------------------------------

test_that("bootCI returns est/lci/uci with est = FUN(x)", {
  set.seed(1)
  r <- bootCI(mtcars$mpg, FUN = mean, bci.method = "perc", R = 199)

  expect_named(r, c("est", "lci", "uci"))
  expect_equal(r[["est"]], mean(mtcars$mpg))
  expect_lt(r[["lci"]], r[["est"]])
  expect_gt(r[["uci"]], r[["est"]])
})

test_that("bootCI bounds equal boot::boot.ci() for every usable type", {
  x <- mtcars$mpg
  for (type in c("norm", "basic", "perc", "bca")) {
    set.seed(7)
    r <- bootCI(x, FUN = mean, bci.method = type, R = 499)
    set.seed(7)
    b <- boot::boot(x, function(x, d) mean(x[d]), R = 499)
    ref <- lumen:::.bootCIBounds(boot::boot.ci(b, type = type), type)
    expect_equal(unname(r[c("lci", "uci")]), ref, info = type)
  }
})

test_that("bootCI evaluates the dots in the caller's frame", {
  # regression: substitute() made 'tr' unresolvable inside boot()
  f <- function(tr) bootCI(mtcars$mpg, FUN = mean, trim = tr,
                           bci.method = "perc", R = 99)
  set.seed(1)
  expect_equal(f(0.1)[["est"]], mean(mtcars$mpg, trim = 0.1))
})

test_that("bootCI strips names returned by FUN", {
  set.seed(1)
  r <- bootCI(mtcars$mpg, FUN = function(z) quantile(z, 0.5),
              bci.method = "perc", R = 99)
  expect_named(r, c("est", "lci", "uci"))
  expect_equal(r[["est"]], median(mtcars$mpg))
})

test_that("bootCI resamples rows of a matrix / data frame", {
  d <- mtcars[, c("mpg", "hp")]
  set.seed(3)
  r <- bootCI(d, FUN = function(z) cor(z[, 1], z[, 2]),
              bci.method = "perc", R = 199)
  expect_equal(r[["est"]], cor(d$mpg, d$hp))
  expect_true(r[["lci"]] >= -1 && r[["uci"]] <= 1)

  set.seed(3)
  r2 <- bootCI(as.matrix(d), FUN = function(z) cor(z[, 1], z[, 2]),
               bci.method = "perc", R = 199)
  expect_equal(r2, r)
})

test_that("bootCI bivariate: y is resampled jointly with x", {
  sp <- function(x, y) cor(x, y, method = "spearman")
  set.seed(5)
  r <- bootCI(mtcars$mpg, mtcars$hp, FUN = sp, bci.method = "perc", R = 199)
  set.seed(5)
  r2 <- bootCI(mtcars[, c("mpg", "hp")],
               FUN = function(z) sp(z[, 1], z[, 2]),
               bci.method = "perc", R = 199)
  expect_equal(r, r2)
})

test_that("bootCI one-sided = end of the two-sided interval at 2*cl - 1", {
  x <- mtcars$mpg
  set.seed(11)
  two <- bootCI(x, FUN = mean, bci.method = "norm", conf.level = 0.90, R = 199)
  set.seed(11)
  lft <- bootCI(x, FUN = mean, bci.method = "norm", conf.level = 0.95,
                sides = "left", R = 199)
  set.seed(11)
  rgt <- bootCI(x, FUN = mean, bci.method = "norm", conf.level = 0.95,
                sides = "r", R = 199)

  expect_equal(lft[["lci"]], two[["lci"]])
  expect_identical(lft[["uci"]], Inf)
  expect_equal(rgt[["uci"]], two[["uci"]])
  expect_identical(rgt[["lci"]], -Inf)
})

test_that("bootCI rejects one-sided conf.level <= 0.5", {
  expect_error(bootCI(1:10, FUN = mean, sides = "left", conf.level = 0.5),
               "above 0.5")
})

test_that("bootCI 'stud' fails with a message naming the interval", {
  set.seed(1)
  expect_error(
    suppressWarnings(bootCI(mtcars$mpg, FUN = mean, bci.method = "stud", R = 99)),
    "stud")
})

test_that("bootCI rejects unknown bci.method", {
  expect_error(bootCI(1:10, FUN = mean, bci.method = "all"))
})
