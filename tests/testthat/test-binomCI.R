
# ===============================================================
# binomCI TESTS
# ===============================================================

# ===============================================================
# helper
.check_bci <- function(res) {
  
  stopifnot(is.numeric(as.matrix(res)))
  
  if (is.null(dim(res))) {
    nms <- names(res)
  } else {
    nms <- colnames(res)
  }
  
  stopifnot(all(c("est", "lci", "uci") %in% nms))
  stopifnot(res["lci"] >= 0)
  stopifnot(res["uci"] <= 1)
  stopifnot(res["lci"] <= res["uci"])
  
  invisible(TRUE)
  
}

# ===============================================================
# methods
methods_bci <- c(
  "wald", "wald-cc", "wilson", "wilson-cc", "wilson-mod",
  "agresti-coull", "jeffreys", "jeffreys-mod",
  "clopper-pearson", "arcsine", "logit",
  "pratt", "mid-p", "blaker", "likelihood", "khouadji"
)
# witting excluded from deterministic tests (randomized)

# ===============================================================
# basic functionality
for (m in methods_bci) {
  
  cat("\n----------------------------------\n")
  cat("Testing:", m, "\n")
  
  res <- binomCI(x = 37, n = 43, method = m)
  
  print(res)
  .check_bci(res)
  
  stopifnot(res["lci"] <= res["est"])
  stopifnot(res["uci"] >= res["est"])
  
}

# ===============================================================
# reference values
# ---------------------------------------------------------------
# Wilson == prop.test(correct=FALSE)
res.wilson <- binomCI(x = 37, n = 43, method = "wilson")
res.prop   <- prop.test(x = 37, n = 43, correct = FALSE)

cat("\nWilson vs prop.test:\n")
stopifnot(abs(res.wilson["lci"] - res.prop$conf.int[1]) < 1e-6)
stopifnot(abs(res.wilson["uci"] - res.prop$conf.int[2]) < 1e-6)

# Wilson-cc == prop.test(correct=TRUE)
res.wilsoncc <- binomCI(x = 37, n = 43, method = "wilson-cc")
res.prop.cc  <- prop.test(x = 37, n = 43, correct = TRUE)

cat("\nWilson-cc vs prop.test(correct=TRUE):\n")
stopifnot(abs(res.wilsoncc["lci"] - res.prop.cc$conf.int[1]) < 1e-6)
stopifnot(abs(res.wilsoncc["uci"] - res.prop.cc$conf.int[2]) < 1e-6)

# Clopper-Pearson == binom.test
res.cp   <- binomCI(x = 42, n = 43, method = "clopper-pearson")
res.binom <- binom.test(x = 42, n = 43)$conf.int

cat("\nClopper-Pearson vs binom.test:\n")
stopifnot(abs(res.cp["lci"] - res.binom[1]) < 1e-6)
stopifnot(abs(res.cp["uci"] - res.binom[2]) < 1e-6)

# ===============================================================
# point estimate = x/n for all standard methods
for (m in methods_bci) {
  
  res <- binomCI(x = 15, n = 40, method = m)
  stopifnot(abs(res["est"] - 15/40) < 1e-10)
  
}

# ===============================================================
# conf.level effect: wider interval at higher level
for (m in methods_bci) {
  
  ci95 <- binomCI(x = 15, n = 40, method = m, conf.level = 0.95)
  ci99 <- binomCI(x = 15, n = 40, method = m, conf.level = 0.99)
  
  stopifnot(ci99["lci"] <= ci95["lci"])
  stopifnot(ci99["uci"] >= ci95["uci"])
  
}

# ===============================================================
# bounds always in [0, 1]
for (m in methods_bci) {
  for (x in c(0, 1, 2, 38, 39, 40)) {
    
    # logit ist bei x=0 und x=n nicht definiert
    if (m == "logit" && x %in% c(0, 40)) next
    
    res <- binomCI(x = x, n = 40, method = m)
    stopifnot(!anyNA(res[c("lci","uci")]))
    stopifnot(res["lci"] >= 0)
    stopifnot(res["uci"] <= 1)
  }
}




# ===============================================================
# edge cases: x = 0
for (m in methods_bci) {
  
  if (m == "logit") next
  
  res <- binomCI(x = 0, n = 20, method = m)
  
  stopifnot(res["est"] == 0)
  
  # arcsine uses p.tilde > 0 internally, so lci is not exactly 0
  if (m != "arcsine")
    stopifnot(res["lci"] == 0)
  else
    stopifnot(res["lci"] >= 0 && res["lci"] < 0.05)
  
  # wald collapses to uci=0 when x=0 (known limitation)
  if (m != "wald")
    stopifnot(res["uci"] > 0)
  
}


# ===============================================================
# edge cases: x = n
for (m in methods_bci) {
  
  if (m == "logit") next
  
  res <- binomCI(x = 20, n = 20, method = m)
  
  stopifnot(res["est"] == 1)
  
  # arcsine uses p.tilde < 1 internally, so uci is not exactly 1
  if (m != "arcsine")
    stopifnot(res["uci"] == 1)
  else
    stopifnot(res["uci"] > 0.95 && res["uci"] <= 1)
  
  # wald collapses to lci=1 when x=n (known limitation)
  if (m != "wald")
    stopifnot(res["lci"] < 1)
  
}

# ===============================================================
# one-sided intervals
for (m in methods_bci) {
  
  res.left  <- binomCI(x = 10, n = 40, method = m, sides = "left")
  res.right <- binomCI(x = 10, n = 40, method = m, sides = "right")
  res.two   <- binomCI(x = 10, n = 40, method = m, sides = "two.sided")
  
  # left: uci = 1
  stopifnot(res.left["uci"] == 1)
  
  # right: lci = 0
  stopifnot(res.right["lci"] == 0)
  
  # one-sided wider on the open side
  stopifnot(res.left["lci"]  <= res.two["lci"])
  stopifnot(res.right["uci"] >= res.two["uci"])
  
}

# ===============================================================
# vectorization
res <- binomCI(
  x = c(42, 35, 23, 22),
  n = 43,
  method = "wilson"
)
print(res)
stopifnot(is.data.frame(res))
stopifnot(nrow(res) == 4)

res2 <- binomCI(
  x = c(42, 35, 23, 22),
  n = c(50, 60, 70, 80),
  method = "jeffreys"
)
print(res2)
stopifnot(is.data.frame(res2))
stopifnot(nrow(res2) == 4)

# ===============================================================
# witting: reproducible with set.seed, bounds in [0,1]
set.seed(42)
res.wit <- binomCI(x = 10, n = 30, method = "witting")
print(res.wit)
stopifnot(res.wit["lci"] >= 0)
stopifnot(res.wit["uci"] <= 1)
stopifnot(res.wit["lci"] <= res.wit["uci"])

# ===============================================================
# Newcombe (1998) Table I reference values
# x=81, n=263 – Wilson interval
res.newc <- binomCI(x = 81, n = 263, method = "wilson")
cat("\nNewcombe (1998) Wilson, x=81, n=263:\n")
print(round(res.newc, 4))
stopifnot(abs(res.newc["lci"] - 0.2553) < 0.001)
stopifnot(abs(res.newc["uci"] - 0.3662) < 0.001)

# ===============================================================
# stdEst = FALSE returns adjusted estimator for agresti-coull
res.std  <- binomCI(x = 81, n = 263, method = "agresti-coull", stdEst = TRUE)
res.adj  <- binomCI(x = 81, n = 263, method = "agresti-coull", stdEst = FALSE)

stopifnot(res.std["est"] == 81/263)
stopifnot(res.adj["est"] != 81/263)   # p.tilde != x/n

# ===============================================================
# invalid input: conf.level out of range
for (bad.level in c(0, 1, -0.5, 1.5)) {
  
  ok <- FALSE
  tryCatch({
    binomCI(x = 10, n = 40, conf.level = bad.level)
  }, error = function(e) {
    ok <<- TRUE
  })
  stopifnot(ok)
  
}

# ===============================================================
# stress test
set.seed(123)
for (i in 1:300) {
  
  n  <- sample(5:200, 1)
  x  <- sample(0:n, 1)
  
  for (m in methods_bci) {
    
    res <- try(
      binomCI(x = x, n = n, method = m),
      silent = TRUE
    )
    
    if (inherits(res, "try-error")) {
      cat("\nFAILED:\n")
      print(list(x = x, n = n, method = m))
      stop("Stress test failed.")
    }

    # logit can be NA at x=0 or x=n
    if (m == "logit" && x %in% c(0, n)) next

    if (!is.na(res["lci"]) && !is.na(res["uci"]) && res["lci"] > res["uci"])
      cat("lci > uci:", m, "x=", x, "n=", n, 
          "lci=", res["lci"], "uci=", res["uci"], "\n")
    
    stopifnot(!anyNA(res[c("lci", "uci")]))
    stopifnot(res["lci"] >= 0)
    stopifnot(res["uci"] <= 1)
    stopifnot(res["lci"] <= res["uci"])
    
  }
  
}

cat("\n====================================\n")
cat("ALL binomCI TESTS PASSED\n")
cat("====================================\n")


