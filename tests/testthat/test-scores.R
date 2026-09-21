# scores() -------------------------------------------------------------------

tab <- matrix(c(10, 20, 30,
                 5,  5, 10), nrow = 3,
              dimnames = list(dose = c("0", "2.5", "10"), resp = c("a", "b")))

test_that("scores 'table' uses numeric dimnames, else 1..k", {
  expect_identical(scores(tab), c(0, 2.5, 10))
  expect_identical(scores(tab, MARGIN = 2), 1:2)
  expect_identical(scores(unname(tab)), 1:3)
  expect_identical(scores(unname(tab), MARGIN = 2), 1:2)
  # partly numeric dimnames fall back to 1..k
  t2 <- tab; rownames(t2) <- c("0", "x", "10")
  expect_identical(scores(t2), 1:3)
})

test_that("scores ranks / ridit / modridit are midranks of the margin", {
  # row totals 15, 25, 40 -> midranks 8, 28, 60.5
  mr <- c(8, 28, 60.5)
  expect_equal(scores(tab, method = "ranks"), mr, ignore_attr = TRUE)
  expect_equal(scores(tab, method = "ridit"), mr / 80, ignore_attr = TRUE)
  expect_equal(scores(tab, method = "modridit"), mr / 81, ignore_attr = TRUE)

  # column totals 60, 20 -> 30.5, 70.5
  expect_equal(scores(tab, MARGIN = 2, method = "ranks"), c(30.5, 70.5),
               ignore_attr = TRUE)
})

test_that("scores ranks agree with rank() on the expanded data", {
  lv <- rep(seq_len(nrow(tab)), rowSums(tab))
  expect_equal(unname(scores(tab, method = "ranks")),
               as.vector(tapply(rank(lv), lv, mean)))
})

test_that("scores default method and abbreviation", {
  expect_identical(scores(tab), scores(tab, method = "table"))
  expect_identical(scores(tab, method = "mod"), scores(tab, method = "modridit"))
  expect_error(scores(tab, method = "foo"))
})
