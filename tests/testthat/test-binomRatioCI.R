
# ===============================================================
# binomRatioCI TESTS
# ===============================================================
tol <- 1e-7
# ===============================================================

# helper
.check <- function(res){
  
  stopifnot(is.numeric(as.matrix(res)))
  
  if(is.null(dim(res))){
    nms <- names(res)
  } else {
    nms <- colnames(res)
  }
  
  stopifnot(all(c("est","lci","uci") %in% nms))
  stopifnot(res["lci"] <= res["uci"])
  
  invisible(TRUE)
  
}
# ===============================================================
# methods
methods <- c(
  "katz-log",
  "adj-log",
  "bailey",
  "koopman",
  "noether",
  "sinh-1"
)
# ===============================================================
# basic functionality
for(m in methods){
  
  cat("\n----------------------------------\n")
  cat("Testing:", m, "\n")
  
  res <- binomRatioCI(
    x1 = 12, n1 = 100,
    x2 = 20, n2 = 100,
    method = m
  )
  
  print(res)
  .check(res)
  
  stopifnot(res["lci"] <= res["est"])
  stopifnot(res["uci"] >= res["est"])
  
}
# ===============================================================
# Koopman (1984) reference value
# Example from Koopman (1984): x1=36, n1=40, x2=16, n2=80
# Reported 95% CI: (2.93, 7.15) – Table 1
res.koop <- binomRatioCI(
  x1 = 36, n1 = 40,
  x2 = 16, n2 = 80,
  method = "koopman"
)
cat("\nKoopman (1984) reference:\n")
print(round(res.koop, 4))

stopifnot(abs(res.koop["est"] - 4.5)    < 0.001)
stopifnot(abs(res.koop["lci"] - 2.9396) < 0.001)
stopifnot(abs(res.koop["uci"] - 7.1522) < 0.001)


# ===============================================================
# conf.level effect: wider interval at higher level
for(m in methods){
  
  ci95 <- binomRatioCI(
    x1 = 12, n1 = 100,
    x2 = 20, n2 = 100,
    method     = m,
    conf.level = 0.95
  )
  
  ci99 <- binomRatioCI(
    x1 = 12, n1 = 100,
    x2 = 20, n2 = 100,
    method     = m,
    conf.level = 0.99
  )
  
  stopifnot(ci99["lci"] <= ci95["lci"])
  stopifnot(ci99["uci"] >= ci95["uci"])
  
}


# ===============================================================
# vectorization
res <- binomRatioCI(
  x1 = c(5, 10, 20),
  n1 = c(50, 100, 200),
  x2 = c(4,   8,  25),
  n2 = c(50, 100, 200),
  method = c("katz-log", "koopman")
)
print(res)
stopifnot(is.data.frame(res))


# ===============================================================
# edge cases
# ---------------------------------------------------------------
# both zero
for(m in methods){
  
  res <- binomRatioCI(
    x1 = 0, n1 = 100,
    x2 = 0, n2 = 100,
    method = m
  )
  
  print(res)
  stopifnot(res["est"] == 0)
  stopifnot(res["lci"] == 0)
  stopifnot(is.infinite(res["uci"]))
  
}
# ---------------------------------------------------------------
# x1 = 0, x2 > 0  →  est = 0, lci = 0
for(m in methods){
  
  res <- binomRatioCI(
    x1 =  0, n1 = 100,
    x2 = 10, n2 = 100,
    method = m
  )
  
  print(res)
  stopifnot(res["est"] == 0)
  stopifnot(res["lci"] == 0)
  
}
# ---------------------------------------------------------------
# x2 = 0, x1 > 0  →  est = Inf, uci = Inf
for(m in methods){
  
  res <- binomRatioCI(
    x1 = 10, n1 = 100,
    x2 =  0, n2 = 100,
    method = m
  )
  
  print(res)
  stopifnot(is.infinite(res["est"]))
  stopifnot(!is.na(res["lci"]))
  stopifnot(res["lci"] >= 0)
  stopifnot(is.infinite(res["uci"]))
  
}


# ---------------------------------------------------------------
# both one
for(m in methods){
  
  res <- binomRatioCI(
    x1 = 100, n1 = 100,
    x2 = 100, n2 = 100,
    method = m
  )
  
  print(res)
  .check(res)
  
}


# ===============================================================
# one-sided intervals
for(m in methods){
  
  # left: uci must be Inf
  res.left <- binomRatioCI(
    x1 = 10, n1 = 100,
    x2 = 20, n2 = 100,
    method = m,
    sides  = "left"
  )
  
  print(res.left)
  stopifnot(is.infinite(res.left["uci"]))
  
  # right: lci must be 0
  res.right <- binomRatioCI(
    x1 = 10, n1 = 100,
    x2 = 20, n2 = 100,
    method = m,
    sides  = "right"
  )
  
  print(res.right)
  stopifnot(res.right["lci"] == 0)
  
  # one-sided must be wider on the open side than two-sided
  res.two <- binomRatioCI(
    x1 = 10, n1 = 100,
    x2 = 20, n2 = 100,
    method = m,
    sides  = "two.sided"
  )
  
  stopifnot(res.left["lci"]  <= res.two["lci"])
  stopifnot(res.right["uci"] >= res.two["uci"])
  
}

# ===============================================================
# monotonicity sanity check
for(m in methods){
  
  rr1 <- binomRatioCI(
    x1 =  5, n1 = 100,
    x2 = 10, n2 = 100,
    method = m
  )["est"]
  
  rr2 <- binomRatioCI(
    x1 = 20, n1 = 100,
    x2 = 10, n2 = 100,
    method = m
  )["est"]
  
  stopifnot(rr2 > rr1)
  
}


# ===============================================================
# compare methods
tab <- sapply(methods, function(m){
  binomRatioCI(
    x1 = 25, n1 = 100,
    x2 = 10, n2 = 100,
    method = m
  )
})
print(round(tab, 5))


# ===============================================================
# invalid input: x > n
ok <- FALSE
tryCatch({
  binomRatioCI(
    x1 = 101, n1 = 100,
    x2 =  10, n2 = 100
  )
}, error = function(e){
  ok <<- TRUE
})
stopifnot(ok)


# ===============================================================
# invalid input: conf.level out of range
for(bad.level in c(0, 1, -0.5, 1.5)){
  
  ok <- FALSE
  tryCatch({
    binomRatioCI(
      x1 = 10, n1 = 100,
      x2 = 10, n2 = 100,
      conf.level = bad.level
    )
  }, error = function(e){
    ok <<- TRUE
  })
  stopifnot(ok)
  
}


# ===============================================================
# stress test
set.seed(123)
for(i in 1:500){
  
  n1 <- sample(10:500, 1)
  n2 <- sample(10:500, 1)
  x1 <- sample(0:n1, 1)
  x2 <- sample(0:n2, 1)
  
  for(m in methods){
    
    res <- try(
      binomRatioCI(
        x1 = x1, n1 = n1,
        x2 = x2, n2 = n2,
        method = m
      ),
      silent = TRUE
    )
    
    if(inherits(res, "try-error")){
      cat("\nFAILED:\n")
      print(list(x1=x1, n1=n1, x2=x2, n2=n2, method=m))
      stop("Stress test failed.")
    }
    
    # no NA/NaN in output
    stopifnot(!anyNA(res[c("lci","uci")]))
    
    # lci <= uci (allow Inf)
    stopifnot(
      is.infinite(res["lci"]) ||
        is.infinite(res["uci"]) ||
        res["lci"] <= res["uci"]
    )
    
    # lci must be non-negative for ratio
    stopifnot(res["lci"] >= 0)
    
  }
  
}


# ===============================================================
# coverage simulation (optional / slow)
# set.seed(1)
# ...
cat("\n====================================\n")
cat("ALL TESTS PASSED\n")
cat("====================================\n")


