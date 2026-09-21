# Regressions from the pattern sweep ------------------------------------------

# formula + subset must equal the formula on the pre-filtered data, and the
# data name must come from resolveFormula()'s 'dataName'

test_that("two-sample formula methods: subset is evaluated in the data", {
  wb <- warpbreaks
  s  <- wb[wb$tension == "L", ]
  for (fn in list(yuenTTest, mosesTest, siegelTukeyTest)) {
    a <- suppressWarnings(fn(breaks ~ wool, data = wb, subset = tension == "L"))
    b <- suppressWarnings(fn(breaks ~ wool, data = s))
    expect_equal(a$statistic, b$statistic)
    expect_equal(a$p.value, b$p.value)
    expect_type(a$data.name, "character")
    expect_match(a$data.name, "breaks")
  }
})

test_that("block-design formula methods: subset is evaluated in the data", {
  set.seed(1)
  d <- expand.grid(id = factor(1:12), time = factor(1:3))
  d$resp <- rbinom(nrow(d), 1, c(0.3, 0.5, 0.7)[d$time])
  d$x <- rnorm(nrow(d), as.integer(d$time))
  s <- d[d$id != "1", ]
  s$id <- droplevels(s$id)

  a <- cochranQTest(resp ~ time | id, data = d, subset = id != "1")
  b <- cochranQTest(resp ~ time | id, data = s)
  expect_equal(a$statistic, b$statistic)
  expect_match(a$data.name, "resp")

  a <- pageTest(x ~ time | id, data = d, subset = id != "1")
  b <- pageTest(x ~ time | id, data = s)
  expect_equal(a$statistic, b$statistic)
  expect_match(a$data.name, "x")
})

test_that("block designs: unused factor levels are no empty blocks", {
  set.seed(1)
  d <- expand.grid(id = factor(1:12), time = factor(1:3))
  d$resp <- rbinom(nrow(d), 1, c(0.3, 0.5, 0.7)[d$time])
  d$x <- rnorm(nrow(d), as.integer(d$time))
  s <- d[d$id != "1", ]          # level "1" of id stays, but is empty
  expect_equal(cochranQTest(s$resp, s$time, s$id)$statistic,
               cochranQTest(s$resp, s$time, droplevels(s$id))$statistic)
  expect_equal(pageTest(s$x, s$time, s$id)$statistic,
               pageTest(s$x, s$time, droplevels(s$id))$statistic)
})

test_that("bpTest is Koenker's studentized BP test: identical to lmtest::bptest", {
  skip_if_not_installed("lmtest")
  set.seed(1)
  d <- data.frame(x = rnorm(60), z = rnorm(60))
  d$y <- d$x + rnorm(60) * (1 + abs(d$x))
  for (f in list(y ~ x, y ~ x + z, y ~ poly(x, 2) + z)) {
    a <- bpTest(lm(f, data = d))
    b <- lmtest::bptest(lm(f, data = d))
    expect_equal(unname(a$statistic), unname(b$statistic))
    expect_equal(unname(a$parameter), unname(b$parameter))
    expect_equal(a$p.value, unname(b$p.value))
  }
})

test_that("bpTest: na.exclude does not inflate n", {
  set.seed(2)
  d <- data.frame(x = rnorm(50)); d$y <- d$x + rnorm(50) * (1 + abs(d$x))
  d$x[1:5] <- NA
  expect_equal(bpTest(lm(y ~ x, d, na.action = na.exclude))$statistic,
               bpTest(lm(y ~ x, d))$statistic)
})

test_that("bpTest: rejects weights and glm", {
  d <- data.frame(x = rnorm(30), y = rnorm(30), w = runif(30))
  expect_error(bpTest(lm(y ~ x, d, weights = w)), "weighted")
  expect_error(bpTest(glm(y ~ x, data = d)), "lm object")
  expect_identical(bpTest(lm(y ~ x, d))$data.name, "y ~ x")
})


# -- one-sided intervals need conf.level > 0.5 ---------------------------------

test_that("CI functions reject one-sided conf.level <= 0.5", {
  x <- c(4.1, 5.3, 2.2, 6.7, 3.9, 5.0, 4.4, 6.1, 3.3, 5.8)
  y <- c(3.1, 4.0, 2.9, 5.2, 3.6, 4.8)
  calls <- list(
    meanCI     = function(cl, s) meanCI(x, conf.level = cl, sides = s),
    meanDiffCI = function(cl, s) meanDiffCI(x, y, conf.level = cl, sides = s),
    medianCI   = function(cl, s) medianCI(x, conf.level = cl, sides = s),
    quantileCI = function(cl, s) quantileCI(x, probs = 0.5, conf.level = cl,
                                            sides = s),
    varCI      = function(cl, s) varCI(x, conf.level = cl, sides = s),
    poissonCI  = function(cl, s) poissonCI(7, 2, conf.level = cl, sides = s)
  )
  for (nm in names(calls)) {
    f <- calls[[nm]]
    expect_error(f(0.5, "left"), "above 0.5", info = nm)
    expect_error(f(0.3, "right"), "above 0.5", info = nm)
    expect_error(f(1.2, "two.sided"), "'conf.level'", info = nm)
    expect_error(f(NA, "two.sided"), "'conf.level'", info = nm)
    # two-sided below 0.5 stays legal
    expect_no_error(f(0.4, "two.sided"))
    expect_no_error(f(0.9, "left"))
  }
})

test_that("madCI family rejects one-sided conf.level <= 0.5", {
  x <- c(4.1, 5.3, 2.2, 6.7, 3.9, 5.0, 4.4, 6.1, 3.3, 5.8)
  y <- c(3.1, 4.0, 2.9, 5.2, 3.6, 4.8, 3.9)
  expect_error(madCI(x, conf.level = 0.5, sides = "left"), "above 0.5")
  expect_error(madDiffCI(x, y, conf.level = 0.4, sides = "right"), "above 0.5")
  expect_error(madRatioCI(x, y, conf.level = 0.5, sides = "left"), "above 0.5")
})

# -- bootstrap bounds by name ---------------------------------------------------

test_that("bootstrap CIs read their bounds from the named boot.ci component", {
  x <- mtcars$mpg
  for (type in c("norm", "basic", "perc")) {
    set.seed(3)
    r <- varCI(x, method = "boot", type = type, R = 199)
    set.seed(3)
    b <- boot::boot(x, function(x, d) var(x[d]), R = 199)
    expect_equal(unname(r[c("lci", "uci")]),
                 lumen:::.bootCIBounds(boot::boot.ci(b, type = type), type),
                 info = type)
  }
  # no variances: a clear message instead of "subscript out of bounds"
  set.seed(3)
  expect_error(suppressWarnings(varCI(x, method = "boot", type = "stud", R = 99)),
               "stud")
  set.seed(3)
  expect_error(suppressWarnings(quantileCI(x, probs = 0.5, method = "boot",
                                           type = "stud", R = 99)),
               "stud")
})

# -- normality tests: degenerate input -----------------------------------------

test_that("pearsonTest and shapiroFranciaTest reject degenerate input", {
  expect_error(pearsonTest(rep(3, 20)), "identical")
  expect_error(pearsonTest(c(1:20, Inf)), "infinite")
  expect_error(pearsonTest(letters), "numeric")
  expect_error(pearsonTest(c(1.2, 3.4)), "degree of freedom")
  expect_error(shapiroFranciaTest(rep(3, 20)), "identical")
  expect_error(shapiroFranciaTest(c(1:20, -Inf)), "infinite")
  expect_error(shapiroFranciaTest(letters), "numeric")
})

test_that("pearsonTest and shapiroFranciaTest still equal nortest", {
  skip_if_not_installed("nortest")
  set.seed(4)
  for (x in list(rnorm(30), rexp(80), rt(200, 3))) {
    expect_equal(pearsonTest(x)$p.value, nortest::pearson.test(x)$p.value)
    expect_equal(shapiroFranciaTest(x)$p.value, nortest::sf.test(x)$p.value)
  }
})

# -- weighted aov --------------------------------------------------------------

test_that("postHocTest and scheffeTest reject weighted aov", {
  w <- rep(c(1, 2), length.out = nrow(warpbreaks))
  fit <- aov(breaks ~ tension, data = warpbreaks, weights = w)
  expect_error(postHocTest(fit), "weighted")
  expect_error(scheffeTest(fit), "weighted")
  one <- aov(breaks ~ tension, data = warpbreaks,
             weights = rep(1, nrow(warpbreaks)))
  expect_equal(postHocTest(one)$tension,
               postHocTest(aov(breaks ~ tension, data = warpbreaks))$tension)
})
