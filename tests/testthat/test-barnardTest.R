
# Mehta and Senchaudhuri (2003), vaccine example
tab_mehta <- matrix(c(7, 12, 8, 3), nrow = 2,
                    dimnames = list(treat = c("vaccine", "placebo"),
                                    infection = c("yes", "no")))
tab_small <- matrix(c(8, 14, 1, 3), nrow = 2)


# -------------------------------------------------------------------------
# htest contract
# -------------------------------------------------------------------------
test_that("barnardTest: htest contract", {
  r <- barnardTest(tab_mehta, method = "z-pooled")
  expect_s3_class(r, "htest")
  expect_named(r$statistic, "z")
  expect_named(r$parameter, c("n1", "n2"))
  expect_equal(unname(r$parameter), c(15, 15))
  expect_named(r$estimate, "difference in proportion")
  expect_equal(unname(r$estimate), 7 / 15 - 12 / 15)
  expect_identical(r$model, "binomial")
  expect_identical(r$tsmethod, "square")
  expect_identical(r$data.name, "tab_mehta")
  expect_match(r$method, "z-pooled ordering")
  expect_true(r$p.value >= 0 && r$p.value <= 1)
  expect_output(print(r), "Barnard")
  # CSM has no ordering statistic
  expect_null(barnardTest(tab_mehta)$statistic)
})


# -------------------------------------------------------------------------
# reference values: Exact 3.3, exact.test(..., useStoredCSM = FALSE)
# -------------------------------------------------------------------------
test_that("barnardTest: p-values of the Mehta example match Exact", {
  ref <- rbind(
    "csm"           = c(0.0829997385869, 0.0366894992456, 0.963501305316),
    "z-pooled"      = c(0.0682183092995, 0.0341091546498, 1),
    "z-unpooled"    = c(0.0682183092995, 0.0341091546498, 1),
    "boschloo"      = c(0.0682183092995, 0.0341091546498, 0.977121534756),
    "santner-snell" = c(0.0987371467054, 0.0493685733527, 1))
  alts <- c("two.sided", "less", "greater")
  for (m in rownames(ref)) for (k in 1:3)
    expect_equal(barnardTest(tab_mehta, method = m,
                             alternative = alts[k])$p.value,
                 unname(ref[m, k]), tolerance = 1e-10,
                 info = paste(m, alts[k]))

  expect_equal(unname(barnardTest(tab_mehta, method = "z-pooled")$statistic),
               -1.89433807606, tolerance = 1e-10)
})


# p-values are refined with optimize() and agree to about 1e-9 relative
# (the exact value depends on the search interval and the platform);
# statistics are deterministic and are compared at 1e-10

test_that("barnardTest: further reference values", {
  expect_equal(barnardTest(tab_small, fixed = 2)$p.value,
               0.721796847403, tolerance = 1e-10)
  expect_equal(barnardTest(tab_mehta, tsMethod = "central")$p.value,
               0.0733789984911, tolerance = 1e-10)
  # Exact's grid stops at 0.00001 / 0.99999 and does not refine the
  # multinomial maximum: here the supremum is slightly larger
  expect_gte(barnardTest(tab_mehta, method = "z-pooled",
                         alternative = "greater", delta = -0.1)$p.value,
             0.944435232161)
  expect_equal(barnardTest(tab_mehta, delta = 0.1)$p.value,
               0.0149546414721, tolerance = 1e-8)
  expect_equal(barnardTest(tab_mehta, method = "z-pooled",
                           npInterval = TRUE)$p.value,
               0.0692183092515, tolerance = 1e-8)
  # multinomial model: at least Exact's grid maximum, and close to it
  for (r in list(list(tab_small, "z-pooled", 0.756858245484),
                 list(tab_mehta, "boschloo", 0.0836644715438))) {
    p <- barnardTest(r[[1]], method = r[[2]], fixed = NA)$p.value
    expect_gte(p, r[[3]])
    expect_lt(p - r[[3]], 1e-4)
  }
  expect_equal(barnardTest(tab_mehta, method = "z-unpooled", fixed = NA,
                           alternative = "greater")$p.value, 1)
})


test_that("barnardTest: confidence intervals match Exact", {
  r <- barnardTest(tab_mehta, method = "z-pooled", conf.level = 0.95)
  expect_equal(as.vector(r$conf.int), c(-0.61172895191, 0.02723180294),
               tolerance = 1e-8)
  expect_identical(attr(r$conf.int, "conf.level"), 0.95)
  expect_true(r$conf.int[1] < r$estimate && r$estimate < r$conf.int[2])

  r <- barnardTest(tab_mehta, conf.level = 0.95)
  expect_equal(as.vector(r$conf.int), c(-0.59112655605, 0.02713747241),
               tolerance = 1e-8)

  # one-sided intervals are open on the other side
  expect_equal(barnardTest(tab_mehta, method = "z-pooled", alternative = "less",
                           conf.level = 0.95)$conf.int[1], -1)
  # conf.level = NA (default): no interval
  expect_null(barnardTest(tab_mehta, method = "z-pooled")$conf.int)
  ci90 <- barnardTest(tab_mehta, method = "z-pooled", conf.level = 0.9)$conf.int
  ci95 <- barnardTest(tab_mehta, method = "z-pooled", conf.level = 0.95)$conf.int
  expect_lt(diff(as.vector(ci90)), diff(as.vector(ci95)))
})


# -------------------------------------------------------------------------
# internal consistency
# -------------------------------------------------------------------------
test_that("barnardTest: fixed = 2 equals fixed = 1 on the transposed table", {
  for (m in c("csm", "z-pooled", "boschloo"))
    expect_equal(barnardTest(tab_mehta, fixed = 2, method = m)$p.value,
                 barnardTest(t(tab_mehta), fixed = 1, method = m)$p.value,
                 info = m)
})

test_that("barnardTest: swapping the groups swaps less and greater", {
  sw <- tab_mehta[2:1, ]
  for (m in c("csm", "z-pooled", "santner-snell")) {
    expect_equal(barnardTest(tab_mehta, method = m, alternative = "less")$p.value,
                 barnardTest(sw, method = m, alternative = "greater")$p.value,
                 info = m)
    expect_equal(barnardTest(tab_mehta, method = m)$p.value,
                 barnardTest(sw, method = m)$p.value, info = m)
  }
  expect_equal(unname(barnardTest(sw, method = "z-pooled")$statistic),
               -unname(barnardTest(tab_mehta, method = "z-pooled")$statistic))
})

test_that("barnardTest: two-sided >= one-sided; central = twice one-sided", {
  for (m in c("csm", "z-pooled", "boschloo")) {
    two  <- barnardTest(tab_mehta, method = m)$p.value
    less <- barnardTest(tab_mehta, method = m, alternative = "less")$p.value
    expect_gte(two, less, label = m)
  }
  expect_equal(barnardTest(tab_mehta, method = "z-pooled",
                           tsMethod = "central")$p.value,
               2 * barnardTest(tab_mehta, method = "z-pooled",
                               alternative = "less")$p.value)
})

test_that("barnardTest: symmetric problems report both maximizing nuisance values", {
  np <- barnardTest(tab_mehta, method = "z-pooled")$np[[1]]
  expect_length(np, 2L)
  expect_equal(sum(np), 1, tolerance = 1e-6)
})

test_that("barnardTest: two factors give the same result as their table", {
  df  <- as.data.frame(as.table(tab_mehta))
  tr  <- rep(df$treat, df$Freq)
  inf <- rep(df$infection, df$Freq)
  r <- barnardTest(tr, inf, method = "z-pooled")
  expect_equal(r$p.value, barnardTest(tab_mehta, method = "z-pooled")$p.value)
  expect_identical(r$data.name, "tr and inf")
})

test_that("barnardTest: data.name is a single string", {
  r <- barnardTest(matrix(c(7, 12, 8, 3), nrow = 2, dimnames = list(
    treatment_group = c("vaccine", "placebo"),
    infection_status = c("yes", "no"))), method = "z-pooled")
  expect_length(r$data.name, 1L)
})

test_that("barnardTest: multinomial result structure", {
  r <- barnardTest(tab_mehta, method = "boschloo", fixed = NA)
  expect_identical(r$model, "multinomial")
  expect_named(r$parameter, "N")
  expect_named(r$np, c("first nuisance parameter", "second nuisance parameter"))
  expect_equal(unname(r$estimate), (7 * 3 - 8 * 12) / 30^2)
})

test_that("barnardTest: CSM is fast without stored orderings", {
  skip_on_cran()
  tab <- matrix(c(30, 70, 50, 50), 2, byrow = TRUE)
  t <- system.time(r <- barnardTest(tab))[["elapsed"]]
  expect_lt(t, 10)
  expect_lt(r$p.value, 0.01)
})


# -------------------------------------------------------------------------
# input validation
# -------------------------------------------------------------------------
test_that("barnardTest: input validation", {
  expect_error(barnardTest(matrix(1:6, nrow = 3)), "2 x 2")
  expect_error(barnardTest(tab_small, fixed = c(1, 2)), "fisher.test")
  expect_error(barnardTest(tab_small, fixed = c(2, 1)), "fisher.test")
  for (bad in list(3, 0, "1", c(1, 1)))
    expect_error(barnardTest(tab_small, fixed = bad), "'fixed' must be",
                 info = format(bad))
  expect_error(barnardTest(tab_small, method = "fisher"))
  expect_error(barnardTest(tab_small, alternative = "foo"))
  expect_error(barnardTest(tab_small, fixed = NA), "multinomial")
  expect_error(barnardTest(tab_small, fixed = NA, method = "z-pooled",
                           delta = 0.1), "binomial model")
  expect_error(barnardTest(tab_small, npInterval = TRUE), "CSM")
  expect_s3_class(barnardTest(tab_small, method = "z-unpooled", delta = 0.1),
                  "htest")
  expect_error(barnardTest(tab_small, method = "boschloo", delta = 0.1),
               "central")
  expect_error(barnardTest(tab_small, delta = 1), "delta")
  expect_error(barnardTest(tab_small, method = "z-pooled", conf.level = 1.5))
  expect_error(barnardTest(tab_small, method = "z-pooled", fixed = NA,
                           conf.level = 0.95), "binomial model")
  expect_error(barnardTest(tab_small, npNumbers = 0), "npNumbers")
  expect_error(barnardTest(matrix(c(1.5, 2, 3, 4), 2)), "integer")
  expect_error(barnardTest(matrix(c(0, 0, 3, 4), 2), fixed = 2),
               "at least one observation")
})


# -------------------------------------------------------------------------
# independent reference: brute force in plain R
# -------------------------------------------------------------------------

# restricted ML estimate of the second group's probability under
# pi1 = pi2 + delta, found numerically (independent of the closed-form cubic
# root of Miettinen and Nurminen used in the package)
mlRestricted <- function(a, n1, b, n2, delta) {
  lo <- max(0, -delta); hi <- min(1, 1 - delta)
  ll <- function(p) dbinom(a, n1, p + delta, log = TRUE) +
    dbinom(b, n2, p, log = TRUE)
  o <- optimize(ll, c(lo, hi), maximum = TRUE, tol = 1e-12)
  # the maximum may sit on a bound, where optimize() does not evaluate
  cand <- c(lo, o$maximum, hi)
  cand[which.max(vapply(cand, ll, numeric(1)))]
}

# Barnard's p-value by definition, groups as columns of 'tab': all tables at
# least as extreme as the observed one, maximised over a fine grid of the
# nuisance parameter (a lower bound of the supremum, accurate to ~1e-6)
bfBarnard <- function(tab, method, alternative, delta = 0) {
  n1 <- sum(tab[, 1]); n2 <- sum(tab[, 2])
  stat1 <- function(a, b) {
    p1 <- a / n1; p2 <- b / n2
    num <- p1 - p2 - delta
    se <- switch(method,
      "z-pooled"   = {
        q2 <- if (delta == 0) (a + b) / (n1 + n2)
              else mlRestricted(a, n1, b, n2, delta)
        q1 <- q2 + delta
        sqrt(q1 * (1 - q1) / n1 + q2 * (1 - q2) / n2) },
      "z-unpooled" = sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2),
      "santner-snell" = 1)
    if (abs(num) < 1e-12) 0 else num / se
  }
  stat <- function(a, b) mapply(stat1, a, b)
  g <- expand.grid(a = 0:n1, b = 0:n2)
  t <- stat(g$a, g$b)
  t0 <- stat(tab[1, 1], tab[1, 2])
  eps <- if (is.finite(t0)) 1e-7 * max(1, abs(t0)) else 0
  keep <- switch(alternative,
                 two.sided = abs(t) >= abs(t0) - eps,
                 less      = t <= t0 + eps,
                 greater   = t >= t0 - eps)
  pp <- seq(max(0, -delta), min(1, 1 - delta), length.out = 20001)
  max(vapply(pp, function(p)
    sum(dbinom(g$a[keep], n1, min(1, p + delta)) * dbinom(g$b[keep], n2, p)),
    numeric(1)))
}

test_that("barnardTest: two-sided region contains only tables as extreme as observed", {
  # Exact's boundary walk admits (1, 1) and (3, 1) with |Z| = 0.612 < 0.775
  # and reports 0.7834; scipy.stats.barnard_exact() gives 0.5681192439
  tab <- matrix(c(3, 1, 2, 0), nrow = 2, byrow = TRUE)
  # scipy maximises with shgo, we with optimize(): agreement to ~1e-9
  expect_equal(barnardTest(tab, method = "z-pooled")$p.value, 0.5681192439,
               tolerance = 1e-7)
  expect_equal(barnardTest(tab, method = "z-pooled")$p.value,
               bfBarnard(t(tab), "z-pooled", "two.sided"), tolerance = 1e-6)
})

test_that("barnardTest: one-sided p-values agree with scipy", {
  # scipy.stats.barnard_exact / boschloo_exact (1.18), groups as columns;
  # two-sided values are not compared: scipy drops tables whose |Z| equals
  # the observed one only up to rounding
  tab <- matrix(c(7, 12, 8, 3), 2)
  ref <- c("z-pooled" = 0.0480281572, "z-unpooled" = 0.0553392293,
           "boschloo" = 0.0334806620)
  for (m in names(ref))
    expect_equal(barnardTest(tab, method = m, alternative = "less",
                             fixed = 2)$p.value, unname(ref[m]),
                 tolerance = 1e-8, info = m)
  tab <- matrix(c(8, 14, 1, 3), 2)
  ref <- c("z-pooled" = 0.4599281499, "z-unpooled" = 0.6021256654,
           "boschloo" = 0.3825281901)
  for (m in names(ref))
    expect_equal(barnardTest(tab, method = m, alternative = "greater",
                             fixed = 2)$p.value, unname(ref[m]),
                 tolerance = 1e-7, info = m)
})

test_that("barnardTest: agrees with the brute-force definition", {
  set.seed(20261106)
  for (r in 1:60) {
    n <- sample(1:10, 2)
    x <- rbinom(2, n, runif(2))
    tab <- matrix(c(x[1], n[1] - x[1], x[2], n[2] - x[2]), 2)
    method <- sample(c("z-pooled", "z-unpooled", "santner-snell"), 1)
    alt <- sample(c("two.sided", "less", "greater"), 1)
    delta <- if (runif(1) < 0.4) sample(c(-0.3, 0.2), 1) else 0

    info <- paste(paste(tab, collapse = ","), method, alt, delta)
    ours <- barnardTest(tab, method = method, alternative = alt, fixed = 2,
                        delta = delta)$p.value
    ref <- bfBarnard(tab, method, alt, delta)
    # ours refines the grid maximum numerically; ref is a grid lower bound
    expect_gte(ours, ref - 1e-8, label = info)
    expect_lt(ours - ref, 1e-5, label = info)
  }
})

test_that("bfBarnard() handles infinite statistics", {
  # z-unpooled of diag(2) is +Inf (zero variance); region = both diagonals
  expect_equal(bfBarnard(diag(2), "z-unpooled", "two.sided"), 0.5)
})

test_that("barnardTest: z-pooled with delta != 0 matches the numerical restricted ML", {
  tab <- matrix(c(9, 1, 5, 5), 2)            # groups as columns
  for (d in c(-0.2, 0.15)) {
    r <- barnardTest(tab, method = "z-pooled", fixed = 2, delta = d,
                     alternative = "greater")
    q2 <- mlRestricted(9, 10, 5, 10, d); q1 <- q2 + d
    expect_equal(unname(r$statistic),
                 (0.9 - 0.5 - d) / sqrt(q1 * (1 - q1) / 10 + q2 * (1 - q2) / 10),
                 tolerance = 1e-8, info = d)
    expect_equal(r$p.value, bfBarnard(tab, "z-pooled", "greater", d),
                 tolerance = 1e-5, info = d)
  }
})

test_that("barnardTest: z-unpooled includes delta in the statistic", {
  tab <- matrix(c(9, 1, 5, 5), 2)            # groups as columns
  r <- barnardTest(tab, method = "z-unpooled", fixed = 2, delta = 0.2,
                   alternative = "greater")
  p1 <- 0.9; p2 <- 0.5
  expect_equal(unname(r$statistic),
               (p1 - p2 - 0.2) / sqrt(p1 * (1 - p1) / 10 + p2 * (1 - p2) / 10),
               tolerance = 1e-10)
  expect_equal(r$p.value, bfBarnard(tab, "z-unpooled", "greater", 0.2),
               tolerance = 1e-5)
})


# -------------------------------------------------------------------------
# nuisance grid
# -------------------------------------------------------------------------
test_that("barnardTest: the grid includes the bounds of the parameter space", {
  # one observation per group, no successes: the one-sided score p-value is
  # attained at pi = 0 and is exactly 1 (Exact: 0.9999900001)
  tab <- matrix(c(0, 0, 1, 1), 2)
  expect_identical(barnardTest(tab, method = "z-pooled",
                               alternative = "less")$p.value, 1)
  # delta close to +-1: the grid stays within the admissible range
  tab <- matrix(c(9, 1, 1, 9), 2, byrow = TRUE)
  for (d in c(-0.999, 0.999)) {
    r <- barnardTest(tab, method = "z-pooled", alternative = "less", delta = d)
    expect_true(r$p.value >= 0 && r$p.value <= 1)
    expect_true(all(r$np[[1]] >= max(0, -d) & r$np[[1]] <= min(1, 1 - d)))
    expect_equal(r$np.range[[1]], c(max(0, -d), min(1, 1 - d)))
  }
})

test_that("barnardTest: multinomial maximum is found off the grid", {
  # region = the two diagonal tables, P = 4 p1 (1 - p1) p2 (1 - p2): zero
  # with zero gradient at all corners of the 2-point grid, 0.25 at the centre
  r <- barnardTest(diag(2), method = "z-pooled", fixed = NA, npNumbers = 2)
  expect_equal(r$p.value, 0.25, tolerance = 1e-8)
  expect_equal(c(r$np[[1]], r$np[[2]]), c(0.5, 0.5), tolerance = 1e-4)
})

test_that("barnardTest: np refers to the second group also after the internal swap", {
  # observed p1 - p2 > delta: the groups are swapped internally
  tab <- matrix(c(9, 1, 2, 8), 2)            # groups as columns
  d <- 0.1
  r <- barnardTest(tab, method = "z-pooled", fixed = 2, delta = d)
  expect_equal(r$np.range[[1]], c(max(0, -d), min(1, 1 - d)))
  np <- r$np[[1]]
  expect_true(all(np >= 0 & np <= 1 - d))
  # the p-value function, in the original parametrisation, peaks there
  M <- outer(0:10, 0:10, function(a, b)
    abs(.barnardStat(10, 10, "z-pooled", "two.sided", d)[cbind(a + 1, b + 1)]) >=
      abs(r$statistic))
  pf <- function(p) sum(outer(dbinom(0:10, 10, p + d), dbinom(0:10, 10, p)) * M)
  expect_equal(max(vapply(np, pf, numeric(1))), r$p.value, tolerance = 1e-8)
})

test_that("barnardTest: npNumbers and npInterval are validated", {
  expect_error(barnardTest(tab_small, npNumbers = 1), "at least 2")
  expect_true(is.finite(barnardTest(tab_small, method = "z-pooled",
                                    npNumbers = 2)$p.value))
  # the pooled successes are not binomial for delta != 0
  tab <- matrix(c(90, 10, 10, 90), 2, byrow = TRUE)
  expect_error(barnardTest(tab, method = "z-pooled", alternative = "less",
                           delta = 0.8, npInterval = TRUE), "delta = 0")
  expect_error(barnardTest(tab, method = "z-pooled", npInterval = TRUE,
                           conf.level = 0.95), "delta = 0")
})


# -------------------------------------------------------------------------
# CSM against Exact::exact.test()
# -------------------------------------------------------------------------
test_that("barnardTest: CSM agrees with Exact::exact.test", {
  skip_on_cran()
  skip_if_not_installed("Exact")

  set.seed(20261106)
  for (r in 1:25) {
    n <- sample(2:12, 2)
    x <- rbinom(2, n, runif(2))
    tab <- matrix(c(x[1], n[1] - x[1], x[2], n[2] - x[2]), 2)
    alt <- sample(c("two.sided", "less", "greater"), 1)
    delta <- if (runif(1) < 0.3) sample(c(-0.2, 0.1), 1) else 0

    info <- paste(paste(tab, collapse = ","), alt, delta)
    ours <- barnardTest(tab, alternative = alt, fixed = 2, delta = delta)$p.value
    ref  <- Exact::exact.test(tab, alternative = alt, method = "csm",
                              cond.row = FALSE, delta = delta, to.plot = FALSE,
                              useStoredCSM = FALSE)$p.value
    # identical up to the grid bounds (Exact stops at 0.00001 / 0.99999)
    expect_gte(ours, ref - 1e-8, label = info)
    expect_lt(ours - ref, 1e-4, label = info)
  }
})
