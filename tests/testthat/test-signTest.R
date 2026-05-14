
# ===============================================================
# signTest TESTS
# ===============================================================

# ===============================================================
# helper: extract fields from htest object
.check_signtest <- function(res) {
  
  stopifnot(inherits(res, "htest"))
  stopifnot(!is.null(res$statistic))
  stopifnot(!is.null(res$p.value))
  stopifnot(!is.null(res$conf.int))
  stopifnot(!is.null(res$estimate))
  stopifnot(res$p.value >= 0 && res$p.value <= 1)
  stopifnot(res$conf.int[1] <= res$conf.int[2])
  
  invisible(TRUE)
  
}

# ===============================================================
# reference data from documentation examples
x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

d.light <- data.frame(
  black = c(25.85, 28.84, 32.05, 25.74, 20.89, 41.05, 25.01, 24.96, 27.47),
  white = c(18.23, 20.84, 22.96, 19.68, 19.5,  24.98, 16.61, 16.07, 24.59),
  d     = c(7.62, 8.00,  9.09,  6.06,  1.39,  16.07,  8.40,  8.89,  2.88)
)

# ===============================================================
# basic structure: two-sample
res.two <- signTest(x, y)
cat("\nTwo-sample sign test:\n")
print(res.two)
.check_signtest(res.two)
stopifnot(res.two$method == "Dependent-samples Sign-Test")
stopifnot(names(res.two$statistic) == "S")
stopifnot(names(res.two$parameter) == "number of differences")

# ===============================================================
# basic structure: one-sample
res.one <- signTest(x = d.light$d, mu = 4)
cat("\nOne-sample sign test:\n")
print(res.one)
.check_signtest(res.one)
stopifnot(res.one$method == "One-sample Sign-Test")

# ===============================================================
# two-sample == one-sample on differences
res.diff <- signTest(x = d.light$black - d.light$white)
res.pair <- signTest(x = d.light$black, y = d.light$white)

cat("\nPaired vs difference:\n")
stopifnot(abs(res.diff$p.value   - res.pair$p.value)   < 1e-10)
stopifnot(abs(res.diff$statistic - res.pair$statistic) < 1e-10)

# ===============================================================
# S statistic == number of positive differences
d   <- d.light$d - 4
s   <- sum(d > 0)
res <- signTest(x = d.light$d, mu = 4)
stopifnot(unname(res$statistic) == s)

# ===============================================================
# p-value consistent with binom.test
d       <- x - y
n.valid <- sum(d != 0)
s       <- sum(d > 0)
ref     <- binom.test(x = s, n = n.valid, p = 0.5)
stopifnot(abs(res.two$p.value - ref$p.value) < 1e-10)

# ===============================================================
# alternative hypotheses
res.ts <- signTest(x = d.light$d, mu = 4, alternative = "two.sided")
res.gt <- signTest(x = d.light$d, mu = 4, alternative = "greater")
res.lt <- signTest(x = d.light$d, mu = 4, alternative = "less")

cat("\nAlternatives:\n")
print(res.ts$p.value)
print(res.gt$p.value)
print(res.lt$p.value)

.check_signtest(res.ts)
.check_signtest(res.gt)
.check_signtest(res.lt)

# p-values of one-sided tests sum to approximately 1 + p(two-sided)/2
# (standard relationship for discrete distributions)
stopifnot(res.gt$p.value + res.lt$p.value > 1)

# one-sided p <= two-sided p
stopifnot(res.gt$p.value <= res.ts$p.value ||
            res.lt$p.value <= res.ts$p.value)

# ===============================================================
# conf.level effect: wider interval at higher level
res.95 <- signTest(x = d.light$d, mu = 4, conf.level = 0.95)
res.99 <- signTest(x = d.light$d, mu = 4, conf.level = 0.99)

stopifnot(res.99$conf.int[1] <= res.95$conf.int[1])
stopifnot(res.99$conf.int[2] >= res.95$conf.int[2])

# ===============================================================
# mu = 0 default: median of x with no shift
res.mu0 <- signTest(x = c(-2, -1, 1, 2, 3))
# 3 positive, 2 negative -> S = 3
stopifnot(unname(res.mu0$statistic) == 3)

# ===============================================================
# all positive differences -> S = n, p-value large for two-sided
res.pos <- signTest(x = c(1, 2, 3, 4, 5))
stopifnot(unname(res.pos$statistic) == 5)

# ===============================================================
# all negative differences -> S = 0
res.neg <- signTest(x = c(-1, -2, -3, -4, -5))
stopifnot(unname(res.neg$statistic) == 0)

# ===============================================================
# ties at mu are excluded from n
# x = c(1, 2, 4) with mu=2: d = c(-1, 0, 2), valid = 2
res.tie <- signTest(x = c(1, 2, 4), mu = 2)
stopifnot(unname(res.tie$parameter) == 2)
stopifnot(unname(res.tie$statistic) == 1)

# ===============================================================
# NA handling: NAs are removed
res.na <- signTest(x = c(1, 2, NA, 3, 4))
res.clean <- signTest(x = c(1, 2, 3, 4))
stopifnot(abs(res.na$p.value - res.clean$p.value) < 1e-10)

# ===============================================================
# invalid inputs
ok <- FALSE
tryCatch(signTest(x = c(1, 2, 3), mu = c(1, 2)),
         error = function(e) ok <<- TRUE)
stopifnot(ok)

ok <- FALSE
tryCatch(signTest(x = c(1, 2, 3), conf.level = 1.5),
         error = function(e) ok <<- TRUE)
stopifnot(ok)

ok <- FALSE
tryCatch(signTest(x = c(1, 2, 3), y = c(1, 2)),
         error = function(e) ok <<- TRUE)
stopifnot(ok)

ok <- FALSE
tryCatch(signTest(x = "a"),
         error = function(e) ok <<- TRUE)
stopifnot(ok)

# ===============================================================
# null.value is mu
res <- signTest(x = d.light$d, mu = 4)
stopifnot(res$null.value == 4)
stopifnot(names(res$null.value) == "median")

res.two2 <- signTest(x, y)
stopifnot(names(res.two2$null.value) == "median difference")

# ===============================================================
# estimate is median of (x - mu) + mu = median(x) for mu=0
res.med <- signTest(x = c(1, 2, 3, 4, 5))
stopifnot(abs(unname(res.med$estimate) - median(c(1, 2, 3, 4, 5))) < 1e-10)

cat("\n====================================\n")
cat("ALL signTest TESTS PASSED\n")
cat("====================================\n")

