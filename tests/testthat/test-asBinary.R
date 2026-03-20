
# tests/testthat/test-asBinary.R

test_that(".asBinary works for logical", {
  x <- c(TRUE, FALSE, TRUE)
  expect_equal(.asBinary(x), c(1, 0, 1))
})

test_that(".asBinary works for numeric 0/1", {
  x <- c(0, 1, 1, 0)
  expect_equal(.asBinary(x), x)
})

test_that(".asBinary rejects non-binary numeric", {
  x <- c(0, 1, 2)
  expect_error(.asBinary(x), "binary")
})

test_that(".asBinary works for factor (default)", {
  x <- factor(c("F","U","F"))
  
  expect_warning(
    res <- .asBinary(x),
    "coercing factor"
  )
  
  expect_equal(res, c(0,1,0))
})

test_that(".asBinary works for factor with ref", {
  x <- factor(c("F","U","F"))
  
  expect_no_warning(
    res <- .asBinary(x, ref = "F")
  )
  
  expect_equal(res, c(1,0,1))
})

test_that(".asBinary fails if ref not in levels", {
  x <- factor(c("F","U"))
  
  expect_error(.asBinary(x, ref = "X"))
})

test_that(".asBinary fails for factor with >2 levels", {
  x <- factor(c("A","B","C"))
  
  expect_error(.asBinary(x), "exactly 2 levels")
})

test_that(".asBinary works for character", {
  x <- c("yes","no","yes", NA)
  
  expect_warning(
    res <- .asBinary(x),
    "coercing"
  )
  
  expect_true(all(res[!is.na(res)] %in% c(0,1)))
  expect_equal(is.na(res), is.na(x))
})

test_that(".asBinary works for character with ref", {
  x <- c("yes","no","yes")
  
  expect_no_warning(
    res <- .asBinary(x, ref = "yes")
  )
  
  expect_equal(res, c(1,0,1))
})

test_that(".asBinary fails for character with !=2 values", {
  x <- c("a","b","c")
  
  expect_error(.asBinary(x), "exactly 2")
})

test_that(".asBinary handles NA correctly", {
  x <- c(0, 1, NA)
  
  expect_equal(.asBinary(x), c(0,1,NA))
})

test_that(".asBinary removes names", {
  x <- c(a=0, b=1)
  
  res <- .asBinary(x)
  
  expect_null(names(res))
})

test_that(".asBinary stable under level order", {
  x <- factor(c("A","B","A"), levels=c("B","A"))
  
  res <- .asBinary(x)
  
  expect_equal(res, c(1,0,1))
})