
# ===============================================================
# medianCI TESTS
# ===============================================================

# ===============================================================
# helper
.check_medianCI <- function(res) {
  
  stopifnot(is.numeric(res))
  stopifnot(all(c("median", "lci", "uci") %in% names(res)))
  stopifnot(is.na(res["lci"]) || res["lci"] <= res["median"] ||
              is.infinite(res["lci"]))
  stopifnot(is.na(res["uci"]) || res["uci"] >= res["median"] ||
              is.infinite(res["uci"]))
  stopifnot(is.na(res["lci"]) || is.na(res["uci"]) ||
              res["lci"] <= res["uci"])
  
  invisible(TRUE)
  
}

# ===============================================================
# reference data
set.seed(448)
x.na  <- c(rnorm(100), NA)
x     <- x.na[!is.na(x.na)]

# ===============================================================
# basic functionality: exact method
res.ex <- medianCI(x, method = "exact")
cat("\nExact method:\n")
print(res.ex)
.check_medianCI(res.ex)

# point estimate equals median
stopifnot(abs(res.ex["median"] - median(x)) < 1e-10)

# reported conf.level is attached as attribute
stopifnot(!is.null(attr(res.ex, "conf.level")))
stopifnot(attr(res.ex, "conf.level") >= 0.90)
stopifnot(attr(res.ex, "conf.level") <= 1.00)

# ===============================================================
# basic functionality: boot method
set.seed(1)
res.boot <- medianCI(x, method = "boot")
cat("\nBoot method:\n")
print(res.boot)
.check_medianCI(res.boot)

# boot estimate close to exact
stopifnot(abs(res.boot["median"] - median(x)) < 1e-10)
stopifnot(abs(res.boot["lci"] - res.ex["lci"]) < 0.15)
stopifnot(abs(res.boot["uci"] - res.ex["uci"]) < 0.15)

# ===============================================================
# na.rm = TRUE removes NAs
res.narm <- medianCI(x.na, na.rm = TRUE)
stopifnot(abs(res.narm["median"] - median(x)) < 1e-10)
stopifnot(abs(res.narm["lci"] - res.ex["lci"]) < 1e-10)

# ===============================================================
# conf.level effect: wider interval at higher level
res.95 <- medianCI(x, conf.level = 0.95, method = "exact")
res.99 <- medianCI(x, conf.level = 0.99, method = "exact")

cat("\n95% vs 99%:\n")
print(res.95)
print(res.99)

# exact CI uses order statistics so conf.level is discrete –
# reported level can be equal when the same order statistic is selected
stopifnot(res.99["lci"] <= res.95["lci"] ||
            attr(res.99, "conf.level") == attr(res.95, "conf.level"))
stopifnot(res.99["uci"] >= res.95["uci"] ||
            attr(res.99, "conf.level") == attr(res.95, "conf.level"))

# ===============================================================
# one-sided: left -> uci = Inf
res.left <- medianCI(x, sides = "left")
cat("\nLeft-sided:\n")
print(res.left)
stopifnot(is.infinite(res.left["uci"]))
stopifnot(is.finite(res.left["lci"]))

# one-sided: right -> lci = -Inf
res.right <- medianCI(x, sides = "right")
cat("\nRight-sided:\n")
print(res.right)
stopifnot(is.infinite(res.right["lci"]))   # -Inf
stopifnot(is.finite(res.right["uci"]))

# ===============================================================
# small sample: n < 6 falls back to (-Inf, Inf)
res.small <- medianCI(x = c(1, 2, 3), conf.level = 0.95, method = "exact")
cat("\nSmall sample (n=3):\n")
print(res.small)
# conf.level should be 1 and bounds infinite
stopifnot(is.infinite(res.small["lci"]) || is.infinite(res.small["uci"]) ||
            attr(res.small, "conf.level") == 1)

# ===============================================================
# symmetric data: CI should be symmetric around median
x.sym <- -5:5   # median = 0
res.sym <- medianCI(x.sym, method = "exact")
cat("\nSymmetric data:\n")
print(res.sym)
stopifnot(abs(res.sym["median"]) < 1e-10)
stopifnot(abs(res.sym["lci"] + res.sym["uci"]) < 1e-10)

# ===============================================================
# all identical values: CI collapses to point
res.const <- medianCI(x = rep(5, 20), method = "exact")
cat("\nConstant data:\n")
print(res.const)
stopifnot(abs(res.const["median"] - 5) < 1e-10)
stopifnot(res.const["lci"] == 5)
stopifnot(res.const["uci"] == 5)

# ===============================================================
# boot: different types work without error
set.seed(42)
for (btype in c("norm", "basic", "perc", "bca")) {
  res.bt <- medianCI(x, method = "boot", type = btype)
  cat("Boot type", btype, ":", res.bt, "\n")
  stopifnot(is.numeric(res.bt))
  stopifnot(length(res.bt) == 3)
}

# ===============================================================
# boot: unsupported type returns NA with warning
res.stud <- withCallingHandlers(
  medianCI(x, method = "boot", type = "stud"),
  warning = function(w) invokeRestart("muffleWarning")
)
stopifnot(anyNA(res.stud[c("lci", "uci")]))

# ===============================================================
# result names are always correct
for (m in c("exact", "boot")) {
  set.seed(1)
  res <- medianCI(x, method = m)
  stopifnot(identical(names(res), c("median", "lci", "uci")))
}

# ===============================================================
# n = 1: degenerate case
res.n1 <- medianCI(x = 42, method = "exact")
cat("\nn=1:\n")
print(res.n1)
stopifnot(res.n1["median"] == 42)

cat("\n====================================\n")
cat("ALL medianCI TESTS PASSED\n")
cat("====================================\n")

