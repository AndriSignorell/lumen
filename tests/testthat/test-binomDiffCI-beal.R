
library(testthat)

bdci <- binomDiffCI

# helper: exact coverage calculation
.coverage_bdci <- function(n, p1, p2,
                           method     = "beal",
                           conf.level = 0.95) {
  
  delta_true <- p1 - p2
  cov        <- 0
  
  for (x1 in 0:n) {
    for (x2 in 0:n) {
      
      pr <- dbinom(x1, n, p1) * dbinom(x2, n, p2)
      
      ci <- bdci(x1, n, x2, n,
                 conf.level = conf.level,
                 method     = method)
      
      if (delta_true >= ci["lci"] && delta_true <= ci["uci"])
        cov <- cov + pr
    }
  }
  
  cov
}

test_that("Beal interval achieves near-nominal coverage, n=5, p1=p2=0.3", {
  
  cov <- .coverage_bdci(n  = 5,
                        p1 = 0.3,
                        p2 = 0.3,
                        method     = "beal",
                        conf.level = 0.95)
  
  # discrete CIs cannot achieve exact nominal coverage;
  # a reasonable interval should stay clearly above 0.90
  expect_gt(cov, 0.90)
  expect_lt(cov, 1.00)
})

test_that("Beal interval behaves reasonably in extreme case, n=5, p1=0.9, p2=0.05", {
  
  cov <- .coverage_bdci(n  = 5,
                        p1 = 0.9,
                        p2 = 0.05,
                        method     = "beal",
                        conf.level = 0.95)
  
  # should not collapse like Wald; liberal lower bound for stress test
  expect_gt(cov, 0.80)
  expect_lt(cov, 1.00)
})


test_that("Beal interval is symmetric when p1 = p2 and n1 = n2", {
  
  ci <- bdci(x1 = 8, n1 = 20,
             x2 = 8, n2 = 20,
             method = "beal")
  
  expect_equal(unname(ci["lci"]), -unname(ci["uci"]), tolerance = 1e-10)
  expect_equal(unname(ci["est"]), 0)
})

test_that("Beal interval respects [-1, 1] bounds", {
  
  # extreme case: all successes vs all failures
  ci <- bdci(x1 = 20, n1 = 20,
             x2 =  0, n2 = 20,
             method = "beal")
  
  expect_gte(ci["lci"], -1)
  expect_lte(ci["uci"],  1)
})

test_that("Beal interval point estimate equals difference of proportions", {
  
  ci <- bdci(x1 = 15, n1 = 40,
             x2 = 10, n2 = 40,
             method = "beal")
  
  expect_equal(unname(ci["est"]), 15/40 - 10/40, tolerance = 1e-10)
})

