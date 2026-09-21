# .checkVerbose() / .printSignifCodes() --------------------------------------

test_that(".checkVerbose: argument > option > default", {
  op <- options(DescTools.verbose = NULL)
  on.exit(options(op))

  expect_identical(lumen:::.checkVerbose(), 2L)
  expect_identical(lumen:::.checkVerbose(3), 3L)

  options(DescTools.verbose = 1L)
  expect_identical(lumen:::.checkVerbose(), 1L)
  expect_identical(lumen:::.checkVerbose(3L), 3L)
})

test_that(".checkVerbose rejects everything but a single 1, 2 or 3", {
  for (bad in list(0, 4, 2.5, NA, "2", c(1, 2), numeric(0)))
    expect_error(lumen:::.checkVerbose(bad), "verbose must be",
                 info = format(bad))

  op <- options(DescTools.verbose = 7)
  on.exit(options(op))
  expect_error(lumen:::.checkVerbose(), "verbose must be")
})

test_that(".printSignifCodes prints the legend", {
  expect_output(lumen:::.printSignifCodes(),
                "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
                fixed = TRUE)
  expect_output(lumen:::.printSignifCodes(), "^---")
})
