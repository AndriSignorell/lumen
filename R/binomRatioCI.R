
#' Confidence Intervals for the Ratio of Two Binomial Proportions
#'
#' Computes confidence intervals for the ratio of two independent binomial
#' proportions (relative risk). Several classical and modern methods are
#' available, which may behave quite differently for small samples or extreme
#' proportions.
#'
#' All arguments are vectorized and recycled according to standard R rules.
#'
#' The ratio of proportions is estimated by
#' \deqn{
#' \hat{\theta} =
#' \frac{\hat{p}_1}{\hat{p}_2} =
#' \frac{x_1 / n_1}{x_2 / n_2}.
#' }
#'
#' \strong{Katz-log}:
#' The classical large-sample log-transformed Wald interval described by
#' Katz et al. (1978). This method is simple and widely used but may perform
#' poorly for small sample sizes or proportions close to 0 or 1.
#'
#' \strong{Adjusted-log}:
#' A continuity-adjusted modification of the Katz interval using additive
#' corrections to reduce bias and improve performance in sparse data settings
#' (Walter, 1975; Pettigrew et al., 1986).
#'
#' \strong{Bailey}:
#' A skewness-corrected interval proposed by Bailey (1987) based on a cubic
#' transformation. Often provides improved coverage probabilities relative to
#' standard Wald-type intervals.
#'
#' \strong{Koopman}:
#' An asymptotic score interval obtained by inverting Pearson's chi-square
#' statistic (Koopman, 1984). Confidence limits are determined iteratively via
#' root-finding and the method is generally considered among the most reliable
#' asymptotic procedures.
#'
#' \strong{Noether}:
#' A large-sample interval based directly on the asymptotic variance of the
#' estimated ratio (Noether, 1957).
#'
#' \strong{Inverse hyperbolic sine}:
#' Based on the variance-stabilizing inverse hyperbolic sine transformation
#' proposed by Newcombe (2001).
#'
#' Several methods require special handling in cases where
#' \eqn{x_1 = 0},
#' \eqn{x_2 = 0},
#' \eqn{x_1 = n_1}, or
#' \eqn{x_2 = n_2}.
#' Small continuity corrections are therefore applied internally where needed.
#'
#' Some methods may produce infinite limits when one observed proportion is
#' zero. This is expected behavior for ratio parameters.
#'
#' The choice of method remains an active topic of discussion. The Koopman
#' interval is often recommended due to its comparatively good coverage
#' properties across a broad range of scenarios.
#'
#' @param x1 Number of successes in the first group (numerator).
#' @param n1 Number of trials in the first group.
#' @param x2 Number of successes in the second group (denominator).
#' @param n2 Number of trials in the second group.
#' @param conf.level Confidence level, default is 0.95.
#' @param sides A character string specifying the type of confidence interval:
#'   \code{"two.sided"} (default), \code{"left"}, or \code{"right"}.
#'   Partial matching is allowed.
#' @param method One of:
#'   \code{"koopman"},
#'   \code{"bailey"},
#'   \code{"adj-log"},
#'   \code{"katz-log"},
#'   \code{"sinh-1"},
#'   \code{"noether"}.
#' @param tol Desired accuracy (convergence tolerance) for the iterative
#'   root-finding procedure used by the Koopman interval.
#'
#' @return A matrix with three columns containing the estimated ratio and the
#' lower and upper confidence limits.
#'
#' @references
#' Bailey BJR (1987).
#' Confidence limits to the risk ratio.
#' \emph{Biometrics}, 43(1), 201-205.
#'
#' Katz D, Baptista J, Azen SP, Pike MC (1978).
#' Obtaining confidence intervals for the risk ratio in cohort studies.
#' \emph{Biometrics}, 34, 469-474.
#'
#' Koopman PAR (1984).
#' Confidence intervals for the ratio of two binomial proportions.
#' \emph{Biometrics}, 40, 513-517.
#'
#' Newcombe RG (2001).
#' Logit confidence intervals and the inverse sinh transformation.
#' \emph{The American Statistician}, 55, 200-202.
#'
#' Noether GE (1957).
#' Sample size determination for some common nonparametric tests.
#' \emph{Journal of the American Statistical Association}, 52, 645-647.
#'
#' Pettigrew HM, Gart JJ, Thomas DG (1986).
#' The bias and higher cumulants of the logarithm of a binomial variate.
#' \emph{Biometrika}, 73(2), 425-435.
#'
#' Walter SD (1975).
#' The distribution of Levin's measure of attributable risk.
#' \emph{Biometrika}, 62(2), 371-374.
#'
#' @seealso \code{\link{binom.test}}, \code{\link{prop.test}}
#'
#' @family topic.categoricalData
#' @concept categorical data
#' @concept confidence intervals
#'
#' @examples
#'
#' # Example from Koopman (1984)
#'
#' binomRatioCI(
#'   x1 = 36, n1 = 40,
#'   x2 = 16, n2 = 80,
#'   method = "katz-log"
#' )
#'
#' binomRatioCI(
#'   x1 = 36, n1 = 40,
#'   x2 = 16, n2 = 80,
#'   method = "koopman"
#' )
#'
#'
#' # Compare several methods
#'
#' meths <- c(
#'   "katz-log",
#'   "adj-log",
#'   "bailey",
#'   "koopman",
#'   "noether",
#'   "sinh-1"
#' )
#'
#' binomRatioCI(
#'   x1 = 25, n1 = 100,
#'   x2 = 10, n2 = 100,
#'   method = meths
#' )
#'
#'
#' # Sparse data
#'
#' binomRatioCI(
#'   x1 = 1, n1 = 20,
#'   x2 = 0, n2 = 20,
#'   method = meths
#' )
#'
#'
#' @family ci.proportion
#' @concept confidence-intervals
#' @concept descriptive-statistics
#'


#' @export
binomRatioCI <- function(
    x1, n1,
    x2, n2,
    conf.level = 0.95,
    sides = c("two.sided", "left", "right"),
    method = c(
      "koopman",
      "bailey",
      "adj-log",
      "katz-log",
      "sinh-1",
      "noether"
    ),
    tol = .Machine$double.eps^0.25
) {
  
  
  if(conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be between 0 and 1 (exclusive).")
  
  sides <- match.arg(sides)
  
  if (missing(method)) {
    method <- eval(formals(sys.function())$method)[1]
  } else {
    method <- .resolveMethod(method, several.ok = TRUE)
  }
  
  res <- .recycleApply(
    .binomRatioCI_engine,
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    conf.level = conf.level,
    sides = sides,
    method = method,
    tol = tol
  )
  
  if(length(res) == 1){
    
    out <- res[[1]]
    
  } else {
    
    out <- as.data.frame(attr(res, "recycle"))
    out <- data.frame(do.call(rbind, res), out)
    
  }
  
  return(out)
  
}


# == internal helper functions ===============================================


.binomRatioCI_engine <- function(
    x1, n1,
    x2, n2,
    conf.level,
    sides,
    method,
    tol
) {
  
  if((x1 > n1) || (x2 > n2))
    stop("x cannot be larger than n.")
  
  alpha <- 1 - conf.level
  
  if (sides != "two.sided")
    alpha <- alpha / 2
  
  CI <- switch( method,
    
    "katz-log" = { .brci.katzlog(x1, n1, x2, n2, alpha) },
    "adj-log"  = { .brci.adjlog(x1, n1, x2, n2, alpha) },
    "bailey"   = { .brci.bailey(x1, n1, x2, n2, alpha) },
    "koopman"  = { .brci.koopman(x1, n1, x2, n2, alpha, tol) },
    "noether"  = { .brci.noether(x1, n1, x2, n2, alpha) },
    "sinh-1"   = { .brci.sinh1(x1, n1, x2, n2, alpha) },

    stop(gettextf("Unknown method '%s'.", method))
  )
  
  est <- (x1/n1)/(x2/n2)
  
  if(x1 == 0 && x2 == 0)
    est <- 0
  
  if(x1 > 0 && x2 == 0)
    est <- Inf
  
  ci <- c(
    est = unname(est),
    lci = unname(CI["lci"]),
    uci = unname(CI["uci"])
  )
  
  if(sides == "left")
    ci["uci"] <- Inf
  else if(sides == "right")
    ci["lci"] <- 0
  
  return(ci)
  
}



# ===============================================================
# helper


.brci.z <- function(alpha) {
  qnorm(1 - alpha/2)
}


.brci.rr <- function(x1, n1, x2, n2) {
  (x1/n1)/(x2/n2)
}



# ===============================================================
# katz-log


.brci.katzlog <- function(x1, n1, x2, n2, alpha){
  
  z <- .brci.z(alpha)
  
  if(x1 == 0 && x2 == 0)
    return(c(lci = 0, uci = Inf))
  
  x1a <- ifelse(x1 == 0, 0.5, x1)
  x2a <- ifelse(x2 == 0, 0.5, x2)
  
  if(x1 == n1)
    x1a <- n1 - 0.5
  
  if(x2 == n2)
    x2a <- n2 - 0.5
  
  rr <- .brci.rr(x1a, n1, x2a, n2)
  
  v <- (1/x1a) - (1/n1) +
    (1/x2a) - (1/n2)
  
  lci <- rr * exp(-z * sqrt(v))
  uci <- rr * exp(+z * sqrt(v))
  
  if(x1 == 0 && x2 > 0)
    lci <- 0
  
  if(x1 > 0 && x2 == 0)
    uci <- Inf
  
  c(lci = lci, uci = uci)
  
}




.brci.adjlog <- function(x1, n1, x2, n2, alpha){
  
  z <- .brci.z(alpha)
  
  if(x1 == 0 && x2 == 0)
    return(c(lci = 0, uci = Inf))
  
  x1a <- ifelse(x1 == n1, n1 - 0.5, x1)
  x2a <- ifelse(x2 == n2, n2 - 0.5, x2)
  
  rr <- ((x1a + 0.5)/(n1 + 0.5)) /
    ((x2a + 0.5)/(n2 + 0.5))
  
  v <- (1/(x1a + 0.5)) - (1/(n1 + 0.5)) +
    (1/(x2a + 0.5)) - (1/(n2 + 0.5))
  
  lci <- rr * exp(-z * sqrt(v))
  uci <- rr * exp(+z * sqrt(v))
  
  # NEU: Randfall x1 = 0
  if(x1 == 0 && x2 > 0)
    lci <- 0
  
  # NEU: Randfall x2 = 0
  if(x1 > 0 && x2 == 0)
    uci <- Inf
  
  c(lci = lci, uci = uci)
  
}


.brci.noether <- function(x1, n1, x2, n2, alpha){
  
  z <- .brci.z(alpha)
  
  if(x1 == 0 && x2 == 0)
    return(c(lci = 0, uci = Inf))
  
  x1a <- ifelse(x1 == 0, 0.5, x1)
  x2a <- ifelse(x2 == 0, 0.5, x2)
  
  if(x1 == n1)
    x1a <- n1 - 0.5
  
  if(x2 == n2)
    x2a <- n2 - 0.5
  
  rr <- .brci.rr(x1a, n1, x2a, n2)
  
  se <- rr * sqrt(
    (1/x1a) - (1/n1) +
      (1/x2a) - (1/n2)
  )
  
  lci <- max(0, rr - z * se)
  uci <- rr + z * se
  
  if(x1 > 0 && x2 == 0)
    uci <- Inf
  
  c(lci = lci, uci = uci)
  
}



.brci.sinh1 <- function(x1, n1, x2, n2, alpha){
  
  z <- .brci.z(alpha)
  
  if(x1 == 0 && x2 == 0)
    return(c(lci = 0, uci = Inf))
  
  x1a <- ifelse(x1 == 0, z, x1)
  x2a <- ifelse(x2 == 0, z, x2)
  
  if(x1 == n1)
    x1a <- n1 - 0.5
  
  if(x2 == n2)
    x2a <- n2 - 0.5
  
  rr <- .brci.rr(x1a, n1, x2a, n2)
  
  v <- 2 * asinh(
    (z/2) * sqrt(
      1/x1a + 1/x2a -
        1/n1  - 1/n2
    )
  )
  
  lci <- exp(log(rr) - v)
  uci <- exp(log(rr) + v)
  
  if(x1 == 0 && x2 > 0)
    lci <- 0
  
  if(x1 > 0 && x2 == 0)
    uci <- Inf
  
  c(lci = lci, uci = uci)
  
}



.brci.bailey <- function(x1, n1, x2, n2, alpha){
  
  z <- .brci.z(alpha)
  
  rr <- .brci.rr(x1, n1, x2, n2)
  
  if(x1 == 0 && x2 == 0)
    return(c(lci = 0, uci = Inf))
  
  x1a <- x1
  x2a <- x2
  
  if(x1a == 0)
    x1a <- 0.5
  
  if(x2a == 0)
    x2a <- 0.5
  
  if(x1a == n1)
    x1a <- n1 - 0.5
  
  if(x2a == n2)
    x2a <- n2 - 0.5
  
  rr.adj <- .brci.rr(x1a, n1, x2a, n2)
  
  p1 <- x1a / n1
  p2 <- x2a / n2
  
  q1 <- 1 - p1
  q2 <- 1 - p2
  
  term <- sqrt(
    (q1 / x1a) +
      (q2 / x2a) -
      (z^2 * q1 * q2) /
      (9 * x1a * x2a)
  )
  
  denom <- 1 - (z^2 * q2) / (9 * x2a)
  
  lci <- rr.adj *
    ((1 - z * term / 3) / denom)^3
  
  uci <- rr.adj *
    ((1 + z * term / 3) / denom)^3
  
  if(x1 == 0 && x2 > 0)
    lci <- 0
  
  if(x1 > 0 && x2 == 0)
    uci <- Inf
  
  c(lci = lci, uci = uci)
  
}



.brci.koopman <- function(x1, n1, x2, n2, alpha, tol){
  
  z <- .brci.z(alpha)
  x2q <- qchisq(1 - alpha, 1)
  
  if(x1 == 0 && x2 == 0)
    return(c(lci = 0, uci = Inf))
  
  rr <- .brci.rr(x1, n1, x2, n2)
  
  rr.adj <- rr
  
  if(x1 == 0 && x2 > 0)
    rr.adj <- .brci.rr(x1 + 0.5, n1, x2, n2)
  
  if(x2 == 0 && x1 > 0)
    rr.adj <- .brci.rr(x1, n1, x2 + 0.5, n2)
  
  if(x1 == n1 && x2 == n2)
    rr.adj <- 1
  
  varhat <- (1/max(x1, 0.5)) - (1/n1) +
    (1/max(x2, 0.5)) - (1/n2)
  
  lower.start <- rr.adj * exp(-z * sqrt(varhat)) / 4
  upper.start <- rr.adj * exp(+z * sqrt(varhat)) * 4
  
  score.fun <- function(r){
    
    p.hat <- (
      r * (n1 + x2) + x1 + n2 -
        sqrt(
          (r * (n1 + x2) + x1 + n2)^2 -
            4 * r * (n1 + n2) * (x1 + x2)
        )
    ) / (2 * (n1 + n2))
    
    (
      ((x1 - n1 * p.hat)^2) /
        (n1 * p.hat * (1 - p.hat))
    ) *
      (
        1 + (n1 * (r - p.hat)) /
          (n2 * (1 - p.hat))
      ) -
      x2q
    
  }
  
  # ============================================================
  # x1 = 0
  
  if(x1 == 0 && x2 != n2){
    
    upper <- tryCatch(
      uniroot(
        score.fun,
        c(rr.adj, upper.start),
        tol = tol
      )$root,
      error = function(e) Inf
    )
    
    return(c(
      lci = 0,
      uci = upper
    ))
    
  }
  
  # ============================================================
  # x2 = 0
  
  if(x2 == 0 && x1 != 0){
    
    lower <- tryCatch(
      uniroot(
        score.fun,
        c(lower.start, rr.adj),
        tol = tol
      )$root,
      error = function(e) 0
    )
    
    return(c(
      lci = lower,
      uci = Inf
    ))
    
  }
  
  # ============================================================
  # interior solution
  
  lower <- tryCatch(
    uniroot(
      score.fun,
      c(lower.start, rr.adj),
      tol = tol
    )$root,
    error = function(e) 0
  )
  
  upper <- tryCatch(
    uniroot(
      score.fun,
      c(rr.adj, upper.start),
      tol = tol
    )$root,
    error = function(e) Inf
  )
  
  c(
    lci = lower,
    uci = upper
  )
  
}

