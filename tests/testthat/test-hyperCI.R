
# ===============================================================
# hyperCI TESTS
# ===============================================================

methods_hci <- c("wilson", "wald", "clopper-pearson", "mid-p", "blaker", "wang")
exact_hci   <- c("clopper-pearson", "blaker", "wang")

# all intervals of a method over the sample space, one row per x = 0..n
.hciFamily <- function(n, N, method, conf.level = 0.95)
  t(vapply(0:n, function(x)
    hyperCI(x, n, N, conf.level, method = method)[c("lci", "uci")],
    numeric(2)))

# coverage minimum over all M = 0..N
.hciInfCov <- function(ci, n, N)
  min(vapply(0:N, function(M)
    sum(dhyper(0:n, M, N - M, n)[ci[, 1] <= M & M <= ci[, 2]]), numeric(1)))


test_that("every method gives a valid integer interval around est", {
  for (m in methods_hci) {
    res <- hyperCI(x = 10, n = 50, N = 2000, method = m)
    expect_equal(res[["est"]], 400, label = m)
    expect_true(res[["lci"]] <= res[["est"]] && res[["est"]] <= res[["uci"]],
                label = m)
    expect_equal(res[c("lci", "uci")], round(res[c("lci", "uci")]), label = m)
  }
})


test_that("wang reproduces Wang (2015) / ExactCIone::WhyperCI_M", {
  # ExactCIone 1.0.5, WhyperCI_M(x, n, N, 0.95)$CI
  expect_equal(unname(hyperCI(10, 50, 2000, method = "wang")[c("lci", "uci")]),
               c(211, 661))
  expect_equal(unname(hyperCI(30, 500, 1e5, method = "wang")[c("lci", "uci")]),
               c(4151, 8448))
  # WhyperCI_M(0, 5, 20, 0.95, details = TRUE)$CIM
  expect_equal(unname(.hciFamily(5, 20, "wang")),
               cbind(c(0, 1, 2, 5, 8, 12), c(8, 12, 15, 18, 19, 20)))
})


test_that("clopper-pearson limits satisfy the defining tail conditions", {
  x <- 7; n <- 30; N <- 200; a2 <- 0.025
  ci <- hyperCI(x, n, N, method = "clopper-pearson")
  L <- ci[["lci"]]; U <- ci[["uci"]]
  expect_gt(phyper(x - 1, L, N - L, n, lower.tail = FALSE), a2)
  expect_lte(phyper(x - 1, L - 1, N - L + 1, n, lower.tail = FALSE), a2)
  expect_gt(phyper(x, U, N - U, n), a2)
  expect_lte(phyper(x, U + 1, N - U - 1, n), a2)
})


test_that("exact methods keep the level, blaker and wang lie within cp", {
  for (d in list(c(10, 40), c(7, 23), c(12, 12 * 5))) {
    n <- d[1]; N <- d[2]
    cp <- .hciFamily(n, N, "clopper-pearson")
    for (m in exact_hci) {
      ci <- .hciFamily(n, N, m)
      expect_gte(.hciInfCov(ci, n, N), 0.95 - 1e-12,
                 label = sprintf("%s n=%d N=%d", m, n, N))
      expect_true(all(ci[, 1] >= cp[, 1] & ci[, 2] <= cp[, 2]), label = m)
    }
  }
})


test_that("blaker equals the brute force inversion over all M", {
  x <- 4; n <- 15; N <- 60
  acc <- vapply(0:N, function(M)
    if (x < max(0, n - N + M) || x > min(n, M)) 0
    else lumen:::.hyperBlakerAccept(x, n, N, M), numeric(1))
  expect_equal(unname(hyperCI(x, n, N, method = "blaker")[c("lci", "uci")]),
               range(which(acc > 0.05) - 1))
})


test_that("limits are symmetric in x <-> n - x", {
  for (m in methods_hci) {
    a <- hyperCI(4, 31, 150, method = m); b <- hyperCI(27, 31, 150, method = m)
    expect_equal(a[["lci"]], 150 - b[["uci"]], label = m)
    expect_equal(a[["uci"]], 150 - b[["lci"]], label = m)
  }
})


test_that("limits stay within x <= M <= N - n + x", {
  for (m in methods_hci) for (x in c(0, 3, 20)) {
    ci <- hyperCI(x, 20, 45, method = m)
    expect_gte(ci[["lci"]], x, label = m)
    expect_lte(ci[["uci"]], 45 - 20 + x, label = m)
  }
})


test_that("a census returns M = x for every method", {
  for (m in methods_hci)
    expect_equal(unname(hyperCI(3, 20, 20, method = m)), c(3, 3, 3), label = m)
})


test_that("large N approaches binomCI", {
  N <- 1e7
  for (m in c("wilson", "wald", "clopper-pearson")) {
    h <- hyperCI(10, 50, N, method = m)[c("lci", "uci")] / N
    b <- binomCI(10, 50, method = m)[c("lci", "uci")]
    expect_equal(unname(h), unname(b), tolerance = 1e-5, label = m)
  }
})


test_that("one-sided intervals", {
  x <- 10; n <- 50; N <- 2000
  for (m in methods_hci) {
    l <- hyperCI(x, n, N, method = m, sides = "left")
    r <- hyperCI(x, n, N, method = m, sides = "right")
    t <- hyperCI(x, n, N, method = m)
    expect_equal(l[["uci"]], N - n + x, label = m)
    expect_equal(r[["lci"]], x, label = m)
    expect_gte(l[["lci"]], t[["lci"]], label = m)
    expect_lte(r[["uci"]], t[["uci"]], label = m)
  }
  # blaker and wang one-sided are the one-sided clopper-pearson bound
  for (m in c("blaker", "wang")) for (s in c("left", "right"))
    expect_equal(hyperCI(x, n, N, method = m, sides = s),
                 hyperCI(x, n, N, method = "clopper-pearson", sides = s),
                 label = paste(m, s))
})


test_that("one-sided exact bounds keep the level", {
  # an end of the two-sided blaker interval at the doubled alpha gave the
  # lower bound 2 here, which misses M = 1 with probability 0.1
  expect_equal(hyperCI(1, 1, 10, sides = "left", method = "blaker")[["lci"]], 1)
  n <- 12; N <- 40
  for (m in exact_hci) {
    lo <- vapply(0:n, function(x)
      hyperCI(x, n, N, sides = "left", method = m)[["lci"]], numeric(1))
    cov <- vapply(0:N, function(M)
      sum(dhyper(0:n, M, N - M, n)[lo <= M]), numeric(1))
    expect_gte(min(cov), 0.95 - 1e-12, label = m)
  }
})


test_that("wang stays valid at low levels", {
  # crossed limits (lci 4, uci 3) and a middle limit passing its neighbour
  # (coverage 0 at M = 18) before the steps were checked for validity
  ci <- hyperCI(3, 6, 7, conf.level = 0.2, method = "wang")
  expect_lte(ci[["lci"]], ci[["uci"]])
  for (d in list(c(6, 7), c(13, 36), c(23, 31), c(15, 40)))
    for (cl in c(0.05, 0.2, 0.3, 0.5)) {
      n <- d[1]; N <- d[2]
      ci <- .hciFamily(n, N, "wang", cl)
      lab <- sprintf("n=%d N=%d cl=%g", n, N, cl)
      expect_true(all(ci[, 1] <= ci[, 2]) && !is.unsorted(ci[, 1]), label = lab)
      expect_gte(.hciInfCov(ci, n, N), cl - 1e-12, label = lab)
    }
})


test_that("recycling and conf.level = NA", {
  res <- hyperCI(x = c(0, 5, 10), n = 50, N = c(100, 1000, 2000))
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 3)
  expect_equal(hyperCI(5, 50, 1000, conf.level = NA), 100)
})


test_that("input validation", {
  expect_error(hyperCI(6, 5, 20), "larger than 'n'")
  expect_error(hyperCI(2, 30, 20), "larger than 'N'")
  expect_error(hyperCI(2, 5, 20, conf.level = 0.4, sides = "left"), "above 0.5")
  expect_error(lumen:::.wangHyperCI(2.5, 5, 20, 0.05), "integer")
})
