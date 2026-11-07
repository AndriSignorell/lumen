tol <- 1e-4

test_that("poissonCI: default returns named vector with est/lci/uci", {
  res <- poissonCI(10, 1)
  expect_type(res, "double")
  expect_true(all(c("est", "lci", "uci") %in% names(res)))
})

test_that("poissonCI: est = x/n", {
  expect_equal(unname(poissonCI(10, 2)["est"]), 5)
  expect_equal(unname(poissonCI(0,  1)["est"]), 0)
})

test_that("poissonCI: lci <= est <= uci", {
  for (m in c("exact", "score", "wald", "byar")) {
    res <- poissonCI(10, 1, method = m)
    expect_false(is.na(res["lci"]), label = paste(m, "not NA"))
    expect_lte(unname(res["lci"]), unname(res["est"]), label = paste(m, "lci<=est"))
    expect_lte(unname(res["est"]), unname(res["uci"]), label = paste(m, "est<=uci"))
  }
})

test_that("poissonCI: lci >= 0 for x=0", {
  for (m in c("exact", "score", "byar")) {
    expect_gte(unname(poissonCI(0, 1, method = m)["lci"]), 0, label = m)
  }
})

test_that("poissonCI: exact matches poisson.test", {
  ref <- poisson.test(10, 1)$conf.int
  res <- poissonCI(10, 1, method = "exact")
  expect_equal(unname(res["lci"]), ref[1], tolerance = tol)
  expect_equal(unname(res["uci"]), ref[2], tolerance = tol)
})

test_that("poissonCI: wider CI with higher conf.level", {
  w95 <- diff(unname(poissonCI(10, 1, conf.level = 0.95)[c("lci","uci")]))
  w80 <- diff(unname(poissonCI(10, 1, conf.level = 0.80)[c("lci","uci")]))
  expect_gt(w95, w80)
})

test_that("poissonCI: sides=left gives uci=Inf", {
  expect_equal(unname(poissonCI(10, 1, sides = "left")["uci"]), Inf)
})

test_that("poissonCI: sides=right gives lci=0", {
  expect_equal(unname(poissonCI(10, 1, sides = "right")["lci"]), 0)
})

test_that("poissonCI: n>1 scales est correctly", {
  expect_equal(unname(poissonCI(10, 2)["est"]),
               unname(poissonCI(10, 1)["est"]) / 2, tolerance = tol)
})

test_that("poissonCI: multiple methods returns data.frame with est/lci/uci", {
  res <- poissonCI(10, 1, method = c("exact", "score", "wald", "byar"))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("est", "lci", "uci") %in% names(res)))
  expect_equal(nrow(res), 4L)
  expect_true(all(res[["lci"]] <= res[["est"]]))
  expect_true(all(res[["est"]] <= res[["uci"]]))
})


test_that("poissonCI: byar follows Rothman & Boice (x for lci, x + 1 for uci)", {
  # reference values computed independently from the closed-form formula
  res <- poissonCI(10, 1, method = "byar")
  expect_equal(unname(res[c("lci", "uci")]),
               c(4.787449913967563, 18.391458605756473), tolerance = 1e-10)

  # x = 0: lower limit 0, upper limit from x + 1 = 1, scaled by n
  res0 <- poissonCI(0, 2, method = "byar")
  expect_identical(unname(res0["lci"]), 0)
  expect_equal(unname(res0["uci"]), 1.8340059160868754, tolerance = 1e-10)

  # close to the exact interval for moderate counts
  expect_equal(unname(res[c("lci", "uci")]),
               unname(poissonCI(10, 1, method = "exact")[c("lci", "uci")]),
               tolerance = 5e-3)
})


test_that("poissonCI: invalid x and n are rejected for every method", {
  for (m in c("exact", "score", "wald", "byar")) {
    expect_error(poissonCI(2.5, 1, method = m), label = m)
    expect_error(poissonCI(-1, 1, method = m), label = m)
    expect_error(poissonCI(5, 0, method = m), "'n'", label = m)
    expect_error(poissonCI(5, -2, method = m), "'n'", label = m)
    expect_error(poissonCI(5, Inf, method = m), "'n'", label = m)
  }
  # non-integer exposure is legitimate
  expect_no_error(poissonCI(137, 24.19893))
})
