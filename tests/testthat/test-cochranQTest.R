
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("cochranQTest returns an htest object", {
  
  d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
    rep(1:8, c(6,2,2,6,16,4,4,6)), ]
  row.names(d.frm) <- NULL
  d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[1:3],
                    v.names="resp", direction="long")
  
  res <- cochranQTest(resp ~ time | id, data = d.long)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "Cochran's Q")
  expect_named(res$parameter, "df")
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# SAS reference values
# -------------------------------------------------------------------------
test_that("SAS drugs example gives correct Q and p-value", {
  
  d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
    rep(1:8, c(6,2,2,6,16,4,4,6)), ]
  row.names(d.frm) <- NULL
  d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[1:3],
                    v.names="resp", direction="long")
  
  res <- cochranQTest(resp ~ time | id, data = d.long)
  
  # SAS PROC FREQ: Q = 8.4706, p = 0.0144
  expect_equal(unname(res$statistic["Cochran's Q"]), 8.4706, tolerance = 1e-3)
  expect_equal(res$p.value, 0.01447555, tolerance = 1e-4)
  expect_equal(unname(res$parameter["df"]), 2L)
})
# -------------------------------------------------------------------------
# Matrix interface
# -------------------------------------------------------------------------
test_that("matrix interface gives same result as formula interface", {
  
  d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
    rep(1:8, c(6,2,2,6,16,4,4,6)), ]
  row.names(d.frm) <- NULL
  d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[1:3],
                    v.names="resp", direction="long")
  
  res_formula <- cochranQTest(resp ~ time | id, data = d.long)
  
  # build binary matrix manually
  d_ord <- d.long[order(d.long$id, d.long$time), ]
  mat   <- matrix(as.integer(d_ord$resp == "U"),
                  ncol = 3, byrow = TRUE)
  
  res_mat <- cochranQTest(mat)
  
  expect_equal(unname(res_formula$statistic), unname(res_mat$statistic),
               tolerance = 1e-6)
  expect_equal(res_formula$p.value, res_mat$p.value, tolerance = 1e-6)
})

test_that("vector interface gives same result as formula interface", {
  
  d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
    rep(1:8, c(6,2,2,6,16,4,4,6)), ]
  row.names(d.frm) <- NULL
  d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[1:3],
                    v.names="resp", direction="long")
  
  res_formula <- cochranQTest(resp ~ time | id, data = d.long)
  res_default <- cochranQTest(d.long$resp, d.long$time, d.long$id)
  
  expect_equal(unname(res_formula$statistic), unname(res_default$statistic),
               tolerance = 1e-6)
  expect_equal(res_formula$p.value, res_default$p.value, tolerance = 1e-6)
})
# -------------------------------------------------------------------------
# Known Q = 0 case
# -------------------------------------------------------------------------
test_that("identical responses give Q = 0", {
  
  # all blocks respond the same in every group
  mat <- matrix(c(1,1,1,
                  0,0,0,
                  1,1,1,
                  0,0,0), nrow = 4, byrow = TRUE)
  
  res <- cochranQTest(mat)
  
  expect_equal(unname(res$statistic["Cochran's Q"]), 0)
  expect_equal(res$p.value, 1)
})
# -------------------------------------------------------------------------
# degrees of freedom
# -------------------------------------------------------------------------
test_that("df equals k - 1", {
  
  mat <- matrix(sample(0:1, 20, replace = TRUE), nrow = 5, ncol = 4)
  
  res <- cochranQTest(mat)
  
  expect_equal(unname(res$parameter["df"]), 3L)
})
# -------------------------------------------------------------------------
# Approximate method
# -------------------------------------------------------------------------
test_that("approximate method gives plausible p-value", {
  
  skip_if_not_installed("coin")
  
  set.seed(1)
  
  d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
    rep(1:8, c(6,2,2,6,16,4,4,6)), ]
  row.names(d.frm) <- NULL
  d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[1:3],
                    v.names="resp", direction="long")
  
  res <- cochranQTest(resp ~ time | id, data = d.long,
                      method = "approximate", nresample = 5000)
  
  expect_match(res$method, "approximate")
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
  # should be in ballpark of asymptotic p = 0.0144
  expect_lt(res$p.value, 0.05)
})

test_that("approximate method uses nresample argument", {
  
  skip_if_not_installed("coin")
  
  set.seed(1)
  
  d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
    rep(1:8, c(6,2,2,6,16,4,4,6)), ]
  row.names(d.frm) <- NULL
  d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[1:3],
                    v.names="resp", direction="long")
  
  res <- cochranQTest(resp ~ time | id, data = d.long,
                      method = "approximate", nresample = 999)
  
  expect_match(res$method, "999")
})
# -------------------------------------------------------------------------
# NA handling
# -------------------------------------------------------------------------
test_that("incomplete blocks are removed", {
  
  mat    <- matrix(c(1,0,1, 0,1,0, 1,1,0, 0,0,1), nrow = 4, byrow = TRUE)
  mat_na <- mat
  mat_na[2, 2] <- NA
  
  res_full <- cochranQTest(mat[-2, ])
  res_na   <- cochranQTest(mat_na)
  
  expect_equal(unname(res_full$statistic), unname(res_na$statistic),
               tolerance = 1e-10)
  expect_equal(res_full$p.value, res_na$p.value, tolerance = 1e-10)
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("non-binary matrix throws error", {
  
  mat <- matrix(c(0,1,2,0,1,0), nrow = 2)
  
  expect_error(cochranQTest(mat), "0/1")
})

test_that("single group throws error", {
  
  mat <- matrix(sample(0:1, 5, replace = TRUE), nrow = 5, ncol = 1)
  
  expect_error(cochranQTest(mat), "2 groups")
})

test_that("single block throws error", {
  
  mat <- matrix(c(0,1,1), nrow = 1)
  
  expect_error(cochranQTest(mat), "2 blocks")
})

test_that("replicated design throws error", {
  
  expect_error(
    cochranQTest(
      y      = c(0,1,0,1,0,1,0,1),
      groups = rep(1:2, 4),
      blocks = rep(1:2, each = 4)
    ),
    "unreplicated"
  )
})

test_that("invalid method throws error", {
  
  mat <- matrix(sample(0:1, 12, replace = TRUE), nrow = 4)
  
  expect_error(cochranQTest(mat, method = "invalid"))
})

test_that("no variation gives Q = 0 and p = 1", {
  
  mat <- matrix(rep(1, 12), nrow = 4)
  
  res <- cochranQTest(mat)
  
  expect_equal(unname(res$statistic["Cochran's Q"]), 0)
  expect_equal(res$p.value, 1)
})

# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
    rep(1:8, c(6,2,2,6,16,4,4,6)), ]
  row.names(d.frm) <- NULL
  d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[1:3],
                    v.names="resp", direction="long")
  
  res <- cochranQTest(resp ~ time | id, data = d.long)
  
  expect_output(print(res), "Cochran")
})

