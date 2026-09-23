# signifDiff() ---------------------------------------------------------------

fit <- aov(breaks ~ tension, data = warpbreaks)
ph  <- postHocTest(fit, method = "hsd")
# HSD p-values: M-L 0.038, H-L 0.0014, H-M 0.46; means L > M > H

test_that("PostHocTest: levels, labels and signed differences", {
  s <- signifDiff(ph)
  expect_s3_class(s, "signifDiff")
  expect_named(s, "tension")
  expect_equal(attr(s, "alpha"), 0.05)
  expect_true(attr(s, "signed"))

  d <- s$tension
  expect_identical(rownames(d), c("L", "M", "H"))
  expect_identical(d$label, c("1", "2", "3"))
  expect_identical(d$diff, c("2+, 3+", "1-", "1-"))
})

test_that("alpha: default from conf.level, explicit value wins", {
  s <- signifDiff(postHocTest(fit, conf.level = 0.99))
  expect_equal(attr(s, "alpha"), 0.01)
  expect_identical(s$tension$diff, c("3+", "", "1-"))

  s <- signifDiff(ph, alpha = 0.01)
  expect_identical(s$tension$diff, c("3+", "", "1-"))
  expect_identical(signifDiff(ph, alpha = 1e-6)$tension$diff, c("", "", ""))
})

test_that("direction = FALSE drops the signs", {
  s <- signifDiff(ph, direction = FALSE)
  expect_false(attr(s, "signed"))
  expect_identical(s$tension$diff, c("2, 3", "1", "1"))
})

test_that("labels and sep", {
  expect_identical(signifDiff(ph, labels = "letters")$tension$diff,
                   c("b+, c+", "a-", "a-"))
  expect_identical(signifDiff(ph, labels = "LETTERS")$tension$label,
                   c("A", "B", "C"))
  expect_identical(signifDiff(ph, labels = "names", sep = "/")$tension$diff,
                   c("M+/H+", "L-", "L-"))
  expect_identical(signifDiff(ph, labels = c("lo", "mid", "hi"))$tension$diff,
                   c("mid+, hi+", "lo-", "lo-"))
  lv <- c("control", "treatment", "placebo")
  expect_identical(lumen:::.makeLabels(lv, "abbreviate", minlength = 4),
                   unname(abbreviate(lv, minlength = 4)))
  expect_error(signifDiff(ph, labels = c("a", "b")), "one entry per level")
  expect_error(signifDiff(ph, labels = "foo"))
})

test_that("works on every term of a multi-factor model", {
  fit2 <- aov(breaks ~ wool + tension, data = warpbreaks)
  s <- signifDiff(postHocTest(fit2))
  expect_named(s, c("wool", "tension"))
  expect_identical(rownames(s$wool), c("A", "B"))
})

test_that("ordered = TRUE: same verdicts, other level order", {
  s <- signifDiff(postHocTest(fit, ordered = TRUE), labels = "names")
  expect_identical(rownames(s$tension), c("H", "M", "L"))
  expect_identical(s$tension["L", "diff"], "H+, M+")
  expect_identical(s$tension["H", "diff"], "L-")
})

test_that("p-value branch: no signs, warning only if asked for", {
  pv <- postHocTest(fit, conf.level = NA)
  expect_no_warning(s <- signifDiff(pv))
  expect_equal(attr(s, "alpha"), 0.05)
  expect_false(attr(s, "signed"))
  expect_identical(s$tension$diff, c("2, 3", "1", "1"))
  expect_warning(signifDiff(pv, direction = TRUE), "no differences")
})

test_that("pairwise.htest", {
  pw <- pairwise.t.test(warpbreaks$breaks, warpbreaks$tension)
  s <- signifDiff(pw)
  expect_s3_class(s, "signifDiff")
  expect_named(s, pw$data.name)
  expect_false(attr(s, "signed"))

  # reference from the object's own p-values
  p <- pw$p.value
  sig <- c(L = paste(c("2", "3")[c(p["M", "L"], p["H", "L"]) < 0.05], collapse = ", "))
  expect_identical(s[[1L]]["L", "diff"], sig[["L"]])

  expect_warning(signifDiff(pw, direction = TRUE), "no differences")
  expect_identical(signifDiff(pw, alpha = 1e-8)[[1L]]$diff, c("", "", ""))
})

test_that("print: legend only for signed results", {
  expect_output(print(signifDiff(ph)), "Sign codes")
  expect_output(print(signifDiff(ph)), "alpha = 0.05", fixed = TRUE)
  out <- capture.output(print(signifDiff(ph), legend = FALSE))
  expect_false(any(grepl("Sign codes", out)))
  out <- capture.output(print(signifDiff(ph, direction = FALSE)))
  expect_false(any(grepl("Sign codes", out)))
  # wrapped in expect_output(): print() would write to the console during the run
  expect_output(expect_invisible(print(signifDiff(ph))))
})

test_that(".pairMatrices: symmetric p, antisymmetric d", {
  m <- lumen:::.pairMatrices(ph$tension)
  expect_equal(m$p, t(m$p))
  expect_equal(m$d, -t(m$d))
  expect_equal(m$d["M", "L"], ph$tension["M-L", "diff"])
  expect_equal(m$p["L", "H"], ph$tension["H-L", "pval"])
  expect_true(all(is.na(diag(m$p))))

  # levels attribute is preferred over the pair labels
  z <- ph$tension
  attr(z, "levels") <- c("L", "M", "H")
  expect_equal(lumen:::.pairMatrices(z), m)

  attr(z, "levels") <- c("L", "M", "H", "X")
  expect_error(lumen:::.pairMatrices(z), "does not match")
})

test_that(".deduceLevels", {
  expect_identical(lumen:::.deduceLevels(c("M-L", "H-L", "H-M")),
                   c("L", "M", "H"))
  expect_identical(lumen:::.deduceLevels(c("b-a", "c-a", "d-a", "c-b", "d-b", "d-c")),
                   c("a", "b", "c", "d"))
  # a separator inside a level name cannot be resolved
  expect_error(lumen:::.deduceLevels(c("b-1-a-1", "c-1-a-1", "c-1-b-1")),
               "cannot be recovered")
})

test_that(".letterSeq continues beyond the alphabet", {
  ls <- lumen:::.letterSeq
  expect_identical(ls(3), c("a", "b", "c"))
  expect_identical(ls(26), letters)
  expect_identical(ls(28)[27:28], c("aa", "ab"))
  expect_identical(ls(28, LETTERS)[27:28], c("AA", "AB"))
  expect_length(unique(ls(700)), 700)
})
