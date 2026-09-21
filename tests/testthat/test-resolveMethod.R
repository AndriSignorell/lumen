# .resolveMethod() -----------------------------------------------------------

rm_fun <- function(method = c("alpha", "beta", "gamma"), several.ok = FALSE)
  lumen:::.resolveMethod(method, several.ok = several.ok)

test_that(".resolveMethod reads the choices from the caller's formals", {
  expect_identical(rm_fun("beta"), "beta")
  expect_identical(rm_fun("g"), "gamma")          # abbreviation
  # the full default vector resolves to the first entry, as match.arg()
  expect_identical(rm_fun(), "alpha")
})

test_that(".resolveMethod: NULL gives the first choice", {
  expect_identical(rm_fun(NULL), "alpha")
})

test_that(".resolveMethod: '.all' gives every choice", {
  expect_identical(rm_fun(".all"), c("alpha", "beta", "gamma"))
  expect_identical(rm_fun(".all", several.ok = TRUE), c("alpha", "beta", "gamma"))
})

test_that(".resolveMethod: several.ok", {
  expect_identical(rm_fun(c("beta", "alpha"), several.ok = TRUE),
                   c("beta", "alpha"))
  expect_error(rm_fun(c("beta", "alpha")))
  expect_error(rm_fun("delta"))
})

test_that(".resolveMethod: explicit fn", {
  f <- function(method = c("x1", "x2")) NULL
  expect_identical(lumen:::.resolveMethod("x2", fn = f), "x2")
  expect_identical(lumen:::.resolveMethod(".all", fn = f), c("x1", "x2"))
  expect_identical(lumen:::.resolveMethod(NULL, fn = f), "x1")
})
