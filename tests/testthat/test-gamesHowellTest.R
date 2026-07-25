
# Reference values computed independently from the group means, variances and
# sizes of warpbreaks, using the studentized range distribution.

test_that("gamesHowellTest reproduces the warpbreaks comparisons", {

  r <- gamesHowellTest(breaks ~ tension, data = warpbreaks)

  expect_s3_class(r, "PostHocTest")
  expect_named(r, "tension")
  expect_equal(colnames(r$tension), c("diff", "lci", "uci", "pval"))
  expect_equal(rownames(r$tension), c("M-L", "H-L", "H-M"))

  expect_equal(unname(r$tension[, "diff"]),
               c(-10, -14.72222222, -4.722222222), tolerance = 1e-8)
  expect_equal(unname(r$tension[, "lci"]),
               c(-21.00112832, -25.54578586, -11.86790457), tolerance = 1e-8)
  expect_equal(unname(r$tension[, "uci"]),
               c(1.001128322, -3.898658584, 2.423460122), tolerance = 1e-8)
  expect_equal(unname(r$tension[, "pval"]),
               c(0.080082272, 0.006355163064, 0.2512982436), tolerance = 1e-8)

  expect_equal(attr(r, "conf.level"), 0.95)
  expect_equal(attr(r, "method"), "Games-Howell")
})


test_that("diff is the second group minus the first, as in TukeyHSD", {

  r <- gamesHowellTest(breaks ~ tension, data = warpbreaks)
  m <- tapply(warpbreaks$breaks, warpbreaks$tension, mean)
  expect_equal(unname(r$tension["M-L", "diff"]), unname(m["M"] - m["L"]))
  expect_equal(unname(r$tension["H-M", "diff"]), unname(m["H"] - m["M"]))
})


test_that("the interval and the p-value agree at the confidence level", {

  r <- gamesHowellTest(breaks ~ tension, data = warpbreaks)
  # an interval excluding zero must correspond to pval < 0.05 and vice versa
  excludes <- r$tension[, "lci"] > 0 | r$tension[, "uci"] < 0
  expect_equal(excludes, r$tension[, "pval"] < 0.05)

  # a wider level gives wider intervals
  w <- gamesHowellTest(breaks ~ tension, data = warpbreaks, conf.level = 0.99)
  expect_true(all(w$tension[, "uci"] - w$tension[, "lci"] >
                  r$tension[, "uci"] - r$tension[, "lci"]))
  expect_equal(unname(w$tension[, "diff"]), unname(r$tension[, "diff"]))
})


test_that("with equal sizes and equal variances the SE reduces to Tukey's", {

  # three shifted copies of one vector: identical variances and identical n,
  # so the Games-Howell and Tukey standard errors coincide exactly and the
  # procedures can differ only through the degrees of freedom
  base <- c(-3, -1, 0, 1, 2, 4, 5, -2, 3, -4)
  v <- c(base, base + 1, base + 2)
  g <- rep(letters[1:3], each = 10)

  r <- gamesHowellTest(v, g)
  tk <- TukeyHSD(aov(v ~ factor(g)))[[1L]]

  expect_equal(unname(r$g[, "diff"]), c(1, 2, 1))
  expect_equal(unname(r$g[, "diff"]), unname(tk[, "diff"]), tolerance = 1e-12)

  # pairwise df 2(n-1) = 18 against pooled k(n-1) = 27 makes Games-Howell the
  # slightly more conservative of the two, by the ratio of the two quantiles
  hwGH <- (r$g[, "uci"] - r$g[, "lci"]) / 2
  hwTK <- (tk[, "upr"] - tk[, "lwr"]) / 2
  expect_true(all(hwGH > hwTK))
  expect_equal(unname(hwGH / hwTK),
               rep(qtukey(0.95, 3, 18) / qtukey(0.95, 3, 27), 3),
               tolerance = 1e-10)
})


test_that("the aov and formula methods agree with the default method", {

  a <- gamesHowellTest(breaks ~ tension, data = warpbreaks)
  b <- gamesHowellTest(warpbreaks$breaks, warpbreaks$tension)
  d <- gamesHowellTest(aov(breaks ~ tension, data = warpbreaks))

  expect_equal(unname(a$tension), unname(b[[1L]]))
  expect_equal(unname(a$tension), unname(d$tension))
  expect_named(d, "tension")

  s <- gamesHowellTest(breaks ~ tension, data = warpbreaks, subset = wool == "A")
  expect_false(isTRUE(all.equal(unname(s$tension), unname(a$tension))))

  expect_error(gamesHowellTest(aov(breaks ~ tension + wool, data = warpbreaks)),
               "one-way model")

  # an offset need not show up among the term labels
  d <- warpbreaks
  d$off <- 1
  expect_error(gamesHowellTest(aov(breaks ~ tension + offset(off), data = d)),
               "offset")
})


test_that("invalid arguments are rejected", {

  v <- c(1, 2, 3, 4, 5, 6)
  g <- rep(c("a", "b"), each = 3)

  expect_error(gamesHowellTest(letters[1:6], g), "numeric")
  expect_error(gamesHowellTest(v, g[1:4]), "same length")
  expect_error(gamesHowellTest(v, rep("a", 6)), "at least two levels")
  expect_error(gamesHowellTest(v, g, conf.level = 1), "'conf.level'")

  # each group contributes its own variance, so singletons cannot be used
  expect_error(gamesHowellTest(c(1, 2, 3), c("a", "a", "b")),
               "at least 2 non-missing observations")

  # an infinite value would make var() NaN and break the degeneracy test
  expect_error(gamesHowellTest(c(1, 2, Inf, 4, 5, 6), g), "only finite values")
  expect_error(gamesHowellTest(c(1, 2, -Inf, 4, 5, 6), g), "only finite values")
})


test_that("missing values are removed casewise", {

  d <- warpbreaks
  d$breaks[1:3] <- NA
  r <- gamesHowellTest(breaks ~ tension, data = d)
  s <- gamesHowellTest(breaks ~ tension, data = d[-(1:3), ])
  expect_equal(unname(r$tension), unname(s$tension))
})


test_that("pairs that are constant in both groups are reported as NA", {

  # a and b are both constant: no standard error and no degrees of freedom
  expect_warning(r <- gamesHowellTest(c(1, 1, 2, 2), rep(c("a", "b"), each = 2)),
                 "both constant")
  expect_equal(unname(r[[1L]][, "diff"]), 1)
  expect_true(is.na(r[[1L]][, "pval"]))
  expect_true(all(is.na(r[[1L]][, c("lci", "uci")])))

  # one constant group among several is harmless, provided the remaining pairs
  # clear df >= 2; with five per group the Welch df is 4 or 8
  v <- c(rep(1, 5), rep(2, 5), 3:7, 4:8)
  g <- rep(c("a", "b", "c", "d"), each = 5)
  expect_warning(r <- gamesHowellTest(v, g), "1 of 6 comparisons")
  expect_true(is.na(r$g["b-a", "pval"]))
  expect_false(any(is.na(r$g[setdiff(rownames(r$g), "b-a"), "pval"])))
})


test_that("pairs below two Welch degrees of freedom are reported as NA", {

  # ptukey and qtukey are defined only for df >= 2. Two observations in each of
  # two groups give df = (a+b)^2/(a^2+b^2), which lies in [1, 2] and reaches 2
  # only when the two variances are exactly equal.
  v <- c(1, 1, 3, 5)
  g <- rep(c("a", "b"), each = 2)
  expect_warning(r <- gamesHowellTest(v, g), "fewer than 2")
  expect_true(is.na(r$g[, "pval"]))
  expect_equal(unname(r$g[, "diff"]), 3)

  # the obstacle is the tiny group size, not the procedure: three observations
  # per group with equal variances already give df = 4. Note that the list is
  # named after the deparsed 'g' argument, so this has to go through a variable
  # for r$g to resolve.
  v <- c(1, 3, 5, 5, 7, 9)
  g <- rep(c("a", "b"), each = 3)
  expect_silent(r <- gamesHowellTest(v, g))
  expect_false(is.na(r$g[, "pval"]))
  expect_equal(unname(r$g[, "diff"]), 4)
})
