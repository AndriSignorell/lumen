# .rowsFromNames() / .orderIndex() -------------------------------------------

test_that(".rowsFromNames", {
  rfn <- lumen:::.rowsFromNames
  expect_null(rfn(NULL))

  d <- data.frame(a = 1:5, row.names = letters[1:5])
  expect_identical(rfn(c("b", "d"), d), c(2L, 4L))
  expect_null(rfn(c("b", "z"), d))

  # default row names, matched as names
  expect_identical(rfn(c("2", "5"), data.frame(a = 1:5)), c(2L, 5L))

  # no data: row names read as numbers
  expect_identical(rfn(c("3", "1")), c(3L, 1L))
  expect_null(rfn(c("3", "x")))

  # consistent with what lm() keeps after subset and na.action
  d <- data.frame(y = c(1, NA, 3, 4, 5, 6), x = 1:6)
  fit <- lm(y ~ x, data = d, subset = x != 5)
  expect_identical(rfn(rownames(model.frame(fit)), d), c(1L, 3L, 4L, 6L))
})

test_that(".orderIndex: NULL and plain vectors", {
  oi <- lumen:::.orderIndex
  expect_null(oi(NULL, 5))
  z <- c(3, 1, 2)
  expect_identical(oi(z, 3), order(z))
  # NA goes last, as with order()
  expect_identical(oi(c(2, NA, 1), 3), c(3L, 1L, 2L))
})

test_that(".orderIndex: data frame and formula give successive keys", {
  oi <- lumen:::.orderIndex
  d <- data.frame(a = c(2, 1, 2, 1), b = c(1, 2, 0, 1))
  ref <- order(d$a, d$b)
  expect_identical(oi(d, 4), ref)
  expect_identical(oi(~ a + b, 4, data = d), ref)
  # evaluated in the formula environment when data is missing
  a <- d$a; b <- d$b
  expect_identical(oi(~ a + b, 4), ref)
})

test_that(".orderIndex aligns via rows", {
  oi <- lumen:::.orderIndex
  z <- c(50, 40, 30, 20, 10)
  rows <- c(1L, 3L, 5L)
  expect_identical(oi(z, 3, rows = rows), order(z[rows]))
  expect_identical(oi(~ t, 3, data = data.frame(t = z), rows = rows),
                   order(z[rows]))
})

test_that(".orderIndex errors", {
  oi <- lumen:::.orderIndex
  expect_error(oi(numeric(0), 3), "is empty")
  expect_error(oi(data.frame(), 3), "is empty")
  expect_error(oi(~ poly(t, 2), 5, data = data.frame(t = 1:5)),
               "one value per observation")
  expect_error(oi(1:5, 3), "cannot be aligned")
  expect_error(oi(1:5, 3, rows = 1:2), "cannot be aligned")
  expect_error(oi(1:5, 3, rows = c(1L, 2L, 9L)), "cannot be aligned")
})
