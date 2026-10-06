# madCI(), madDiffCI(), madRatioCI() -- R layer. The bootstrap kernels
# (mad_boot_cpp & co.) have their own tests from the BCa round.

set.seed(20260911)
x <- rlnorm(80)
y <- rlnorm(120, meanlog = 0.5)
z975 <- qnorm(0.975)


# -- asymptotic variance ------------------------------------------------------

test_that(".asv.mad() matches the closed-form asymptotic variance of mad()", {
  # noise-free samples, so only the GLD approximation error remains
  # Normal: n Var(mad(x)) -> 1.4826^2 / (16 phi(q75)^2) = 1.3605
  expect_equal(.asv.mad(qnorm(ppoints(1000))),
               1.4826^2 / (16 * dnorm(qnorm(0.75))^2), tolerance = 0.05)

  # Exp(1): m = log 2, zeta = asinh(1/2); f(m -/+ zeta) = exp(+/-zeta)/2,
  # f(m) = 1/2, so A = cosh(zeta), C = sinh(zeta) = 1/2 and
  # 1 - F(m + zeta) - F(m - zeta) = cosh(zeta) - 1
  zeta <- asinh(0.5)
  A <- cosh(zeta)
  C <- 0.5
  B <- C^2 + 4 * C * 0.5 * (cosh(zeta) - 1)
  expect_equal(.asv.mad(qexp(ppoints(1000))),
               1.4826^2 * (1 + B / 0.25) / (4 * A^2), tolerance = 0.05)
})


# -- input checks -------------------------------------------------------------

test_that("input validation", {
  expect_error(madCI("a"), "'x' must be a non-empty numeric")
  expect_error(madCI(numeric(0)), "'x' must be a non-empty numeric")
  expect_error(madCI(x, conf.level = 0), "conf.level")
  expect_error(madCI(x, conf.level = 1), "conf.level")
  expect_error(madCI(x, conf.level = c(0.9, 0.95)), "conf.level")
  expect_error(madCI(x, conf.level = "0.95"), "conf.level")
  expect_error(madCI(c(NA_real_, NA_real_), na.rm = TRUE), "No non-missing values in 'x'")
  expect_error(madCI(x, sides = "both"))
  expect_error(madCI(x, method = "exact"))

  for (f in list(madDiffCI, madRatioCI)) {
    expect_error(f("a", y), "'x' must be a non-empty numeric")
    expect_error(f(x, "a"), "'y' must be a non-empty numeric")
    expect_error(f(x, numeric(0)), "'y' must be a non-empty numeric")
    expect_error(f(x, y, conf.level = 1.5), "conf.level")
    expect_error(f(c(NA_real_, NA_real_), y, na.rm = TRUE), "No non-missing values in 'x'")
    expect_error(f(x, c(NA_real_, NA_real_), na.rm = TRUE), "No non-missing values in 'y'")
  }

  expect_error(madRatioCI(x, c(1, 1, 1, 2, 3)), "MAD of 'y' is zero")
  expect_error(madCI(x, gldMethod = "XYZ"))
})


# -- classic ------------------------------------------------------------------

test_that("classic: estimates and Wald form", {
  r <- madCI(x)
  expect_named(r, c("est", "lci", "uci"))
  expect_equal(r[["est"]], mad(x))
  expect_equal(unname(r[c("lci", "uci")]),
               mad(x) + c(-1, 1) * z975 * sqrt(.asv.mad(x) / length(x)))

  d <- madDiffCI(x, y)
  expect_equal(d[["est"]], mad(x) - mad(y))
  expect_equal(unname(d[c("lci", "uci")]),
               d[["est"]] + c(-1, 1) * z975 *
                 sqrt(.asv.mad(x) / length(x) + .asv.mad(y) / length(y)))
})


test_that("classic: ratio interval is the delta-method interval on the log scale", {
  q <- madRatioCI(x, y)
  expect_equal(q[["est"]], (mad(x) / mad(y))^2)
  expect_gt(q[["lci"]], 0)

  # log((a/b)^2) = 2 (log a - log b)  =>  se = 2 sqrt(Va/a^2 + Vb/b^2)
  seLog <- 2 * sqrt(.asv.mad(x) / length(x) / mad(x)^2 +
                    .asv.mad(y) / length(y) / mad(y)^2)
  expect_equal(unname(log(q[c("lci", "uci")] / q[["est"]])),
               c(-1, 1) * z975 * seLog)
})


test_that("classic: conf.level, gldMethod and na.rm are honoured", {
  w <- function(r) r[["uci"]] - r[["lci"]]
  expect_gt(w(madCI(x, conf.level = 0.99)), w(madCI(x)))
  expect_equal(w(madCI(x, conf.level = 0.90)) / w(madCI(x)),
               qnorm(0.95) / z975)

  expect_equal(madCI(x, gldMethod = "Lmom")[["uci"]],
               mad(x) + z975 * sqrt(.asv.mad(x, method = "Lmom") / length(x)))

  expect_equal(madCI(c(x, NA), na.rm = TRUE), madCI(x))
  expect_equal(madDiffCI(c(NA, x), c(y, NA), na.rm = TRUE), madDiffCI(x, y))
  expect_equal(madRatioCI(c(NA, x), c(y, NA), na.rm = TRUE), madRatioCI(x, y))
})


test_that("sides: the finite bound equals the two-sided bound at 1 - 2 alpha", {
  for (f in list(function(...) madCI(x, ...),
                 function(...) madDiffCI(x, y, ...),
                 function(...) madRatioCI(x, y, ...))) {
    two   <- f(conf.level = 0.90)
    left  <- f(sides = "left")
    right <- f(sides = "r")                               # partial matching
    expect_equal(left[c("est", "lci")], two[c("est", "lci")])
    expect_identical(left[["uci"]], Inf)
    expect_equal(right[c("est", "uci")], two[c("est", "uci")])
    expect_identical(right[["lci"]], -Inf)
  }
})


# -- bootstrap ----------------------------------------------------------------

test_that("boot: estimates, ordering, sides", {
  for (type in c("perc", "bca")) {
    r <- madCI(x, method = "boot", R = 499, type = type)
    expect_named(r, c("est", "lci", "uci"))
    expect_equal(r[["est"]], mad(x))
    expect_true(r[["lci"]] < r[["est"]] && r[["est"]] < r[["uci"]])

    d <- madDiffCI(x, y, method = "boot", R = 499, type = type)
    expect_equal(d[["est"]], mad(x) - mad(y))
    expect_true(d[["lci"]] < d[["est"]] && d[["est"]] < d[["uci"]])

    q <- madRatioCI(x, y, method = "boot", R = 499, type = type)
    expect_equal(q[["est"]], (mad(x) / mad(y))^2)
    expect_true(0 < q[["lci"]] && q[["lci"]] < q[["est"]] && q[["est"]] < q[["uci"]])
  }

  expect_identical(madCI(x, method = "boot", R = 199, type = "perc",
                         sides = "left")[["uci"]], Inf)
  expect_identical(madDiffCI(x, y, method = "boot", R = 199, type = "perc",
                             sides = "right")[["lci"]], -Inf)
})


test_that("boot: set.seed() determines the result", {
  run <- function(f) { set.seed(7); f() }
  expect_identical(run(function() madCI(x, method = "boot", R = 299, type = "bca")),
                   run(function() madCI(x, method = "boot", R = 299, type = "bca")))
  expect_identical(run(function() madDiffCI(x, y, method = "boot", R = 299, type = "perc")),
                   run(function() madDiffCI(x, y, method = "boot", R = 299, type = "perc")))
  expect_identical(run(function() madRatioCI(x, y, method = "boot", R = 299, type = "perc")),
                   run(function() madRatioCI(x, y, method = "boot", R = 299, type = "perc")))
})


test_that("boot and classic agree roughly for a large normal sample", {
  set.seed(5)
  xn <- rnorm(400)
  cl <- madCI(xn)
  bt <- madCI(xn, method = "boot", R = 999, type = "perc")
  expect_equal(bt[["uci"]] - bt[["lci"]], cl[["uci"]] - cl[["lci"]],
               tolerance = 0.25)
})
