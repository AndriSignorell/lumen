# .extractBootArgs() / .bootCIBounds() ------------------------------------------

xba <- function(...) lumen:::.extractBootArgs(...)

test_that(".extractBootArgs defaults", {
  a <- xba(list())
  expect_identical(a$type, "bca")
  expect_identical(a$R, 999)
  expect_identical(a$parallel, "no")
  expect_true(!is.null(a$ncpus))
})

test_that(".extractBootArgs honours user values and caller default", {
  a <- xba(list(type = "perc", R = 199, parallel = "snow", ncpus = 2))
  expect_identical(a[c("type", "R", "parallel", "ncpus")],
                   list(type = "perc", R = 199, parallel = "snow", ncpus = 2))
  expect_identical(xba(list(), default = "basic")$type, "basic")
})

test_that(".extractBootArgs internal contract on 'default'", {
  expect_error(xba(list(), types = c("perc", "bca"), default = "norm"),
               "internal")
  expect_error(xba(list(), default = c("perc", "bca")), "internal")
})

test_that(".extractBootArgs validates 'type'", {
  expect_error(xba(list(type = "all")), "'type' must be one of")
  expect_error(xba(list(type = c("perc", "bca"))), "'type' must be one of")
  expect_error(xba(list(type = 1)), "'type' must be one of")
  # narrower set offered by a compiled caller
  expect_error(xba(list(type = "norm"), types = c("perc", "bca")),
               "'type' must be one of")
  expect_identical(xba(list(type = "perc"), types = c("perc", "bca"))$type,
                   "perc")
})

test_that(".extractBootArgs validates 'R'", {
  for (bad in list(0, -1, 999.5, Inf, NA_real_, c(99, 199), "999",
                   .Machine$integer.max + 1))
    expect_error(xba(list(type = "perc", R = bad)),
                 "'R' must be a single positive whole number",
                 info = format(bad))
})

test_that(".extractBootArgs: bca replicate floor and soft warning", {
  expect_error(xba(list(type = "bca", R = 20)), "at least 49")
  expect_warning(xba(list(type = "bca", R = 49)), "unstable tails")
  expect_warning(xba(list(type = "bca", R = 100)), "unstable tails")
  expect_no_warning(xba(list(type = "bca", R = 199)))
  # the floor applies to bca only
  expect_no_warning(xba(list(type = "perc", R = 20)))
})

test_that(".extractBootArgs validates 'ncpus' and 'parallel'", {
  for (bad in list(0, 1.5, NA_real_, c(1, 2), "2"))
    expect_error(xba(list(type = "perc", ncpus = bad)), "'ncpus'",
                 info = format(bad))
  for (bad in list("yes", c("no", "snow"), 1))
    expect_error(xba(list(type = "perc", parallel = bad)), "'parallel'",
                 info = format(bad))
})

test_that(".extractBootArgs parallel = FALSE rejects parallel/ncpus", {
  expect_error(xba(list(parallel = "snow"), parallel = FALSE),
               "compiled bootstrap")
  expect_error(xba(list(ncpus = 2), parallel = FALSE), "ncpus")
  a <- xba(list(type = "perc"), parallel = FALSE)
  expect_null(a$parallel)
  expect_null(a$ncpus)
})

test_that(".bootCIBounds maps every type onto its boot.ci component", {
  set.seed(1)
  b <- boot::boot(mtcars$mpg, function(x, d) mean(x[d]), R = 199)
  ci <- boot::boot.ci(b, type = c("norm", "basic", "perc", "bca"))

  expect_identical(lumen:::.bootCIBounds(ci, "norm"),  unname(ci$normal[2:3]))
  expect_identical(lumen:::.bootCIBounds(ci, "basic"), unname(ci$basic[4:5]))
  expect_identical(lumen:::.bootCIBounds(ci, "perc"),  unname(ci$percent[4:5]))
  expect_identical(lumen:::.bootCIBounds(ci, "bca"),   unname(ci$bca[4:5]))

  expect_error(lumen:::.bootCIBounds(ci, "stud"), "returned no")
  expect_error(lumen:::.bootCIBounds(ci, "all"), "no boot.ci\\(\\) component")
})
