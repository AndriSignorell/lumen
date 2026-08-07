test_that("madDiff/madRatio resample the two samples independently", {

  # Der eigentliche Regressionstest. Bei ungleichen Laengen zog die alte
  # Fassung ueber die auf n_max aufgefuellte Matrix, wodurch die
  # effektive Stichprobengroesse pro Replikation binomialverteilt war.
  # Sichtbar wird das an der Streuung: sie war zu gross und haing an
  # max(nx, ny) statt an nx und ny.

  set.seed(1)
  x <- rnorm(50,  sd = 1)
  y <- rnorm(400, sd = 1)

  set.seed(2)
  ci <- mad_diff_boot_cpp(x, y, R = 2000, seed = 7L)

  expect_named(ci, c("est", "lci", "uci"))
  expect_equal(unname(ci["est"]), mad(x) - mad(y), tolerance = 1e-12)
  expect_true(ci[["lci"]] < ci[["est"]] && ci[["est"]] < ci[["uci"]])

  # Die Breite muss der Theorie folgen: die Varianz der Differenz wird
  # von der KLEINEREN Stichprobe dominiert. Vervierfacht man nur die
  # grosse, aendert sich fast nichts.
  set.seed(1)
  y2 <- rnorm(1600, sd = 1)
  set.seed(2)
  ci2 <- mad_diff_boot_cpp(x, y2, R = 2000, seed = 7L)

  w1 <- diff(ci[c("lci", "uci")])
  w2 <- diff(ci2[c("lci", "uci")])
  expect_lt(abs(w1 - w2) / w1, 0.25)
})


test_that("the interval covers the true difference at roughly the right rate", {

  skip_on_cran()

  hits <- replicate(300, {
    x  <- rnorm(60,  sd = 1)
    y  <- rnorm(90,  sd = 2)
    ci <- mad_diff_boot_cpp(x, y, R = 399)
    ci[["lci"]] <= (1 - 2) && (1 - 2) <= ci[["uci"]]
  })

  # 95% nominal; der MAD-Bootstrap deckt etwas schlechter, aber nicht
  # dramatisch. Weit unter 0.88 deutet auf ein falsches Ziehschema hin -
  # genau das war der Fehler.
  expect_gt(mean(hits), 0.88)
})


test_that("the ratio statistic is the squared MAD ratio", {

  set.seed(3)
  x <- rnorm(80, sd = 3)
  y <- rnorm(80, sd = 1)

  ci <- mad_ratio_boot_cpp(x, y, R = 999, seed = 11L)

  expect_equal(unname(ci["est"]), (mad(x) / mad(y))^2, tolerance = 1e-12)
  expect_true(ci[["lci"]] > 0)
  expect_true(ci[["lci"]] < ci[["est"]] && ci[["est"]] < ci[["uci"]])
})


test_that("unequal lengths need no padding any more", {

  set.seed(4)
  x <- rnorm(7)
  y <- rnorm(300)

  expect_silent(mad_diff_boot_cpp(x, y, R = 500, seed = 1L))
  expect_silent(mad_diff_boot_cpp(y, x, R = 500, seed = 1L))
})


test_that("set.seed governs the result, and the seed argument too", {

  set.seed(5)
  x <- rnorm(50); y <- rnorm(60)

  # expliziter Seed
  a <- mad_diff_boot_cpp(x, y, R = 500, seed = 99L)
  b <- mad_diff_boot_cpp(x, y, R = 500, seed = 99L)
  expect_identical(a, b)

  # ohne expliziten Seed: aus R's RNG, nicht aus random_device
  set.seed(6); p <- mad_diff_boot_cpp(x, y, R = 500)
  set.seed(6); q <- mad_diff_boot_cpp(x, y, R = 500)
  expect_identical(p, q)

  set.seed(7); r <- mad_diff_boot_cpp(x, y, R = 500)
  expect_false(isTRUE(all.equal(p[["lci"]], r[["lci"]])))
})


test_that("the result does not depend on the number of threads", {

  set.seed(8)
  x <- rnorm(60); y <- rnorm(90)

  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = 2), add = TRUE)

  RcppParallel::setThreadOptions(numThreads = 1)
  a <- mad_diff_boot_cpp(x, y, R = 999, seed = 3L)
  RcppParallel::setThreadOptions(numThreads = min(4, old))
  b <- mad_diff_boot_cpp(x, y, R = 999, seed = 3L)

  expect_identical(a, b)
})


test_that("missing values are rejected rather than silently dropped", {

  set.seed(9)
  x <- rnorm(40); y <- rnorm(40)
  x[3] <- NA

  expect_error(mad_diff_boot_cpp(x, y, R = 200, seed = 1L), "missing")
  expect_error(mad_ratio_boot_cpp(y, x, R = 200, seed = 1L), "missing")
})


test_that("bca and perc agree on the estimate and stay ordered", {

  set.seed(10)
  x <- rnorm(70, sd = 2); y <- rnorm(70)

  p <- mad_diff_boot_cpp(x, y, R = 999, seed = 4L, method = "perc")
  b <- mad_diff_boot_cpp(x, y, R = 999, seed = 4L, method = "bca")

  expect_equal(p[["est"]], b[["est"]])
  expect_true(b[["lci"]] < b[["uci"]])

  # gleiche Replikationen, also nahe beieinander - aber nicht gleich
  expect_equal(b[["lci"]], p[["lci"]], tolerance = 0.3)
})


test_that("a degenerate sample does not produce a reversed interval", {

  # konstantes y: MAD(y) = 0, jeder Jackknife-Wert gleich. Fruher gab
  # das eine erfundene Beschleunigung und konnte adjl > adju liefern.
  set.seed(11)
  x <- rnorm(50)
  y <- rep(2, 50)

  ci <- suppressWarnings(
    mad_diff_boot_cpp(x, y, R = 999, seed = 5L, method = "bca"))

  expect_true(ci[["lci"]] <= ci[["uci"]])
  expect_equal(unname(ci["est"]), mad(x), tolerance = 1e-12)
})
