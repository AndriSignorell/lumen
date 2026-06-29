
# -------------------------------------------------------------------------
# Jarque-Bera Test
# -------------------------------------------------------------------------

#' (Robust) Jarque-Bera Test
#'
#' A goodness-of-fit test for normality based on sample skewness and kurtosis.
#'
#' This function performs either the classical Jarque-Bera test or the robust
#' Jarque-Bera test proposed by Gel and Gastwirth (2008).
#'
#' The robust version replaces the classical standard deviation by the mean
#' absolute deviation from the median.
#'
#' @param x A numeric vector of observations.
#' @param robust Logical. Should the robust version be used?
#' @param method Character string specifying how p-values are obtained:
#'   \code{"chisq"} (asymptotic chi-square approximation) or
#'   \code{"mc"} (Monte Carlo simulation).
#' @param R Number of Monte Carlo simulations.
#' @param na.rm Logical. Should missing values be removed?
#'
#' @return An object of class \code{"htest"}.
#'
#' @note
#' Based on the implementation of \code{rjb.test} from the
#' \pkg{lawstat} package by W. Wallace Hui, Yulia R. Gel,
#' Joseph L. Gastwirth, and Weiwen Miao.
#' The code was refactored and adapted to conform to package standards.
#'  
#' @references
#' Gel, Y. R. and Gastwirth, J. L. (2008).
#' \emph{A robust modification of the Jarque-Bera test of normality}.
#' Economics Letters 99, 30-32.
#'
#' Jarque, C. and Bera, A. (1980).
#' \emph{Efficient tests for normality, homoscedasticity and serial independence
#' of regression residuals}.
#' Economics Letters 6, 255-259.


#'

#' @family test.normality  
#' @concept normality-test  
#' @concept goodness-of-fit
#'
#'
#' @export

jarqueBeraTest <- function(
    x,
    robust = TRUE,
    method = c("chisq", "mc"),
    R = 1000,
    na.rm = FALSE
) {
  
  dname <- paste(deparse(substitute(x)), collapse = " ")
  
  method <- match.arg(method)
  
  if(!is.null(dim(x))) {
    stop("x must be a numeric vector or univariate time series")
  }
  
  x <- as.numeric(x)
  
  if(anyNA(x)) {
    
    if(na.rm) {
      x <- x[!is.na(x)]
    } else {
      stop("missing values in x")
    }
  }
  
  n <- length(x)
  
  if(n < 3L) {
    stop("sample size must be at least 3")
  }
  
  if(method == "mc") {
    
    R <- as.integer(R)
    
    if(is.na(R) || R < 1L) {
      stop("R must be a positive integer")
    }
  }
  
  statistic <- .jb_statistic(x, robust = robust)
  
  p.value <- switch(
    
    method,
    
    chisq = {
      1 - pchisq(statistic, df = 2)
    },
    
    mc = {
      .jb_pvalue_mc(
        statistic = statistic,
        n = n,
        robust = robust,
        R = R
      )
    }
  )
  
  method_name <- paste0(
    if(robust) "Robust " else "",
    "Jarque-Bera Test (",
    if(method == "mc") {
      "Monte Carlo"
    } else {
      "Chi-square approximation"
    },
    ")"
  )
  
  parameter_out <- switch(
    method,
    chisq = c(df = 2),
    mc = c(R = R)
  )
  
  statistic_out <- c("X-squared" = statistic)
  
  structure(
    list(
      statistic = statistic_out,
      parameter = parameter_out,
      p.value = p.value,
      method = method_name,
      data.name = dname
    ),
    class = "htest"
  )
}


# == internal helper functions ============================================


# -------------------------------------------------------------------------
# Jarque-Bera statistic

.jb_statistic <- function(x, robust = TRUE) {
  
  if(diff(range(x)) == 0) {
    stop("all values are identical")
  }
  
  n <- length(x)
  
  zc <- x - mean(x)
  
  m2 <- mean(zc^2)
  m3 <- mean(zc^3)
  m4 <- mean(zc^4)
  
  if(robust) {
    
    # Normalising constant such that:
    # E[sqrt(pi / 2) * |Z - median(Z)|] = sigma
    # for Z ~ N(0, sigma^2)
    
    J <- sqrt(pi / 2) * mean(abs(x - median(x)))
    
    if(J <= 0) {
      stop("robust scale estimate is zero")
    }
    
    sq_skew <- (m3 / J^3)^2
    kurtosis <- m4 / J^4
    
    vs <- 6 / n
    
    # Gel & Gastwirth (2008)
    vk <- 384 / n
    
  } else {
    
    sq_skew <- (m3 / m2^(3 / 2))^2
    kurtosis <- m4 / m2^2
    
    vs <- 6 / n
    vk <- 24 / n
  }
  
  sq_skew / vs + (kurtosis - 3)^2 / vk
}


# -------------------------------------------------------------------------
# Monte Carlo p-value

.jb_pvalue_mc <- function(statistic, n, robust, R) {
  
  sim <- replicate(
    R,
    .jb_statistic(rnorm(n), robust = robust)
  )
  
  # Monte Carlo p-value with finite-sample correction
  (sum(sim >= statistic) + 1) / (R + 1)
}


