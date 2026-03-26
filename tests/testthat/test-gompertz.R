
test_that("dgompertz matches flexsurv", {
  skip_if_not_installed("flexsurv")
  
  x     <- c(0, 0.5, 1, 2, 5)
  shape <- c(0.1, 0.5, 1)
  rate  <- c(0.5, 1, 2)
  
  # standard cases
  expect_equal(
    dgompertz(x, shape = 0.5, rate = 1),
    flexsurv::dgompertz(x, shape = 0.5, rate = 1)
  )
  
  # log density
  expect_equal(
    dgompertz(x, shape = 0.5, rate = 1, log = TRUE),
    flexsurv::dgompertz(x, shape = 0.5, rate = 1, log = TRUE)
  )
  
  # shape = 0 (exponential special case)
  expect_equal(
    dgompertz(x, shape = 0, rate = 1),
    flexsurv::dgompertz(x, shape = 0, rate = 1)
  )
  
  # negative shape
  expect_equal(
    dgompertz(x, shape = -0.5, rate = 1),
    flexsurv::dgompertz(x, shape = -0.5, rate = 1)
  )
  
  # vectorised parameters
  expect_equal(
    dgompertz(x, shape = shape[1:5 %% 3 + 1], rate = 1),
    flexsurv::dgompertz(x, shape = shape[1:5 %% 3 + 1], rate = 1)
  )
})

test_that("pgompertz matches flexsurv", {
  skip_if_not_installed("flexsurv")
  
  x <- c(0, 0.5, 1, 2, 5)
  
  # standard
  expect_equal(
    pgompertz(x, shape = 0.5, rate = 1),
    flexsurv::pgompertz(x, shape = 0.5, rate = 1)
  )
  
  # upper tail
  expect_equal(
    pgompertz(x, shape = 0.5, rate = 1, lower.tail = FALSE),
    flexsurv::pgompertz(x, shape = 0.5, rate = 1, lower.tail = FALSE)
  )
  
  # log.p
  expect_equal(
    pgompertz(x, shape = 0.5, rate = 1, log.p = TRUE),
    flexsurv::pgompertz(x, shape = 0.5, rate = 1, log.p = TRUE)
  )
  
  # shape = 0
  expect_equal(
    pgompertz(x, shape = 0, rate = 1),
    flexsurv::pgompertz(x, shape = 0, rate = 1)
  )
  
  # negative shape (immortal case)
  expect_equal(
    pgompertz(x, shape = -0.5, rate = 1),
    flexsurv::pgompertz(x, shape = -0.5, rate = 1)
  )
})

test_that("qgompertz matches flexsurv", {
  skip_if_not_installed("flexsurv")
  
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  
  # standard
  expect_equal(
    qgompertz(p, shape = 0.5, rate = 1),
    flexsurv::qgompertz(p, shape = 0.5, rate = 1)
  )
  
  # upper tail
  expect_equal(
    qgompertz(p, shape = 0.5, rate = 1, lower.tail = FALSE),
    flexsurv::qgompertz(p, shape = 0.5, rate = 1, lower.tail = FALSE)
  )
  
  # shape = 0 (exponential)
  expect_equal(
    qgompertz(p, shape = 0, rate = 1),
    flexsurv::qgompertz(p, shape = 0, rate = 1)
  )
  
  # negative shape — immortal case returns Inf
  expect_equal(
    qgompertz(0.99, shape = -0.5, rate = 0.1),
    flexsurv::qgompertz(0.99, shape = -0.5, rate = 0.1)
  )
  
  # roundtrip p -> q -> p
  expect_equal(
    pgompertz(qgompertz(p, shape = 0.5, rate = 1), shape = 0.5, rate = 1),
    p
  )
})

test_that("rgompertz matches flexsurv distribution", {
  skip_if_not_installed("flexsurv")
  
  # rgompertz is stochastic so we test distributional properties
  # rather than exact equality
  set.seed(42)
  x <- rgompertz(10000, shape = 0.5, rate = 1)
  
  # mean should be close to theoretical value
  expect_true(all(x >= 0))
  expect_equal(mean(x), 
               unname(mgompertz(shape = 0.5, rate = 1)["mean"]), 
               tolerance = 0.05)
  
  # KS test against flexsurv
  set.seed(42)
  x_flex <- flexsurv::rgompertz(10000, shape = 0.5, rate = 1)
  ks <- ks.test(x, x_flex)
  expect_gt(ks$p.value, 0.01)
})

test_that("gompertz boundary and edge cases", {
  
  # negative x returns 0 density
  expect_equal(dgompertz(-1, shape = 0.5, rate = 1), 0)
  
  # NA input returns NA
  expect_true(is.na(dgompertz(NA, shape = 0.5, rate = 1)))
  
  # invalid rate returns NA
  expect_true(is.na(dgompertz(1, shape = 0.5, rate = -1)))
  
  # p = 0 and p = 1
  expect_equal(qgompertz(0, shape = 0.5, rate = 1), 0)
  expect_equal(qgompertz(1, shape = 0.5, rate = 1), Inf)
})