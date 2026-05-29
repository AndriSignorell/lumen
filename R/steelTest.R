
#' Steel Test for Comparing Several Treatments With a Control
#'
#' A nonparametric multiple-comparison procedure for comparing several
#' treatment groups with a single control group based on Wilcoxon
#' rank-sum statistics.
#'
#' Performs Steel's many-to-one rank test for independent samples.
#' The procedure is the nonparametric analogue of Dunnett's test and
#' controls the family-wise error rate when comparing multiple treatment
#' groups against a common control.
#'
#' The test statistic is based on pairwise Wilcoxon rank-sum statistics
#' comparing each treatment group with the control group. Ties are
#' handled using the asymptotic variance and covariance corrections
#' described by Scholz (2023). Adjusted p-values are obtained from the
#' asymptotic multivariate normal distribution of the standardized
#' statistics.
#'
#' If \code{x} is a list, its elements are taken as the samples to be
#' compared and hence have to be numeric data vectors. In this case,
#' \code{g} is ignored and one can simply use \code{steelTest(x)}.
#'
#' Otherwise, \code{x} must be a numeric vector and \code{g} a grouping
#' factor (or vector coercible to a factor) of the same length.
#'
#' @name steelTest
#' @aliases steelTest steelTest.default steelTest.formula
#'
#' @param x a numeric vector of observations, or a list of numeric
#'   vectors.
#' @param g a grouping factor corresponding to \code{x}. Ignored if
#'   \code{x} is a list.
#' @param control the level of the control group against which all
#'   treatment groups are compared. Defaults to the first group.
#' @param alternative character string specifying the alternative
#'   hypothesis. One of \code{"two.sided"}, \code{"greater"} or
#'   \code{"less"}.
#' @param output character string specifying the output format. One of
#'   \code{"list"} (default) or \code{"matrix"}.
#' @param formula a formula of the form \code{response ~ group}.
#' @param data an optional data frame containing the variables in
#'   \code{formula}.
#' @param subset an optional expression specifying a subset of
#'   observations to be used.
#' @param na.action a function specifying how missing values should be
#'   handled.
#' @param \dots further arguments passed to methods.
#'
#' @return An object of class \code{c("rankTest", "htest")} containing:
#' \describe{
#'   \item{res}{
#'     comparison results. For \code{output="list"} a matrix with
#'     columns \code{W}, \code{z} and \code{pval}; for
#'     \code{output="matrix"} a many-to-one matrix of adjusted
#'     p-values.
#'   }
#'   \item{pmat}{
#'     matrix of adjusted p-values.
#'   }
#'   \item{statistic}{
#'     observed Steel test statistic.
#'   }
#'   \item{p.value}{
#'     asymptotic p-value for the global Steel test.
#'   }
#' }
#'
#' Additional information is stored in attributes:
#' \code{method}, \code{alternative}, \code{output},
#' \code{main}, and \code{data.name}.
#'
#' @details
#' Steel's test is the nonparametric analogue of Dunnett's test.
#' Unlike Dunn's or Conover's procedures, only comparisons against
#' the designated control group are performed.
#'
#' The asymptotic covariance structure follows Scholz (2023), including
#' tie corrections for the Wilcoxon rank-sum statistics. P-values are
#' obtained from the multivariate normal distribution corresponding to
#' the joint asymptotic distribution of the standardized Steel
#' statistics.
#'
#' @references
#'
#' Steel, R. G. D. (1959).
#' A multiple comparison rank sum test: Treatments versus control.
#' \emph{Biometrics}, \strong{15}, 560--572.
#'
#' Scholz, F. W. (2023).
#' Improved tie correction methods for Steel-type rank tests.
#' \emph{Journal of Nonparametric Statistics}, \strong{35}, 541--563.
#'
#' @seealso
#' \code{\link{wilcox.test}}
#'
#' @examples
#' ## Hollander & Wolfe example
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
#' y <- c(3.8, 2.7, 4.0, 2.4)
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
#'
#' steelTest(list(x, y, z))
#'
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)))
#'
#' steelTest(x, g)
#'
#' ## Specify control group
#' steelTest(x, g, control = "1")
#'
#' ## Formula interface
#' steelTest(Ozone ~ Month, data = airquality)
#'
#' @rdname steelTest
#' @family test.posthoc
#' @concept multiple-comparisons
#' @concept nonparametric
#' @concept hypothesis-testing


#' @export
steelTest <- function(x, ...)
  UseMethod("steelTest")



#' @rdname steelTest
#' @export
steelTest.formula <- function(
    formula,
    data,
    subset,
    na.action,
    ...
) {
  
  if (missing(formula) ||
      (length(formula) != 3L) ||
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  
  subset_expr <- if (!missing(subset))
    substitute(subset)
  else
    NULL
  
  na_expr <- if (!missing(na.action))
    substitute(na.action)
  else
    NULL
  
  pf <- resolveFormula(
    formula   = formula,
    data      = data,
    subset    = subset_expr,
    na.action = na_expr,
    allowed   = "n-sample-independent"
  )
  
  out <- steelTest(
    x = pf$x,
    g = pf$group,
    ...
  )
  
  out$data.name <- pf$data.name
  
  out
}


#' @rdname steelTest
#' @export
steelTest.default <- function(
    x,
    g,
    control = NULL,
    alternative = c(
      "two.sided",
      "greater",
      "less"
    ),
    output = c(
      "list",
      "matrix"
    ),
    ...
) {
  
  alternative <- match.arg(alternative)
  output      <- match.arg(output)
  
  DG <- resolveGroups(x, g)
  
  x <- DG$x
  g <- DG$g
  
  gn <- DG$group.names
  
  if (is.null(control))
    control <- gn[1L]
  
  if (!control %in% gn)
    stop(gettextf(
      "control level '%s' not found in grouping variable",
      control
    ))
  
  N <- DG$n
  
  ## ------------------------------------------------------------
  ## Tie corrections (Scholz 2023)
  ## ------------------------------------------------------------
  
  rx <- rank(x, ties.method = "average")
  d  <- as.vector(table(rx))
  
  d2 <- sum(d * (d - 1))
  d3 <- sum(d * (d - 1) * (d - 2))
  
  n2 <- N * (N - 1)
  n3 <- n2 * (N - 2)
  
  ## ------------------------------------------------------------
  ## Control group
  ## ------------------------------------------------------------
  
  ctrl <- x[g == control]
  n0   <- length(ctrl)
  
  trt.names <- setdiff(gn, control)
  m         <- length(trt.names)
  
  if (m < 1L)
    stop("at least one treatment group required")
  
  ## ------------------------------------------------------------
  ## Steel pairwise statistics
  ## ------------------------------------------------------------
  
  W      <- numeric(m)
  mu     <- numeric(m)
  tau    <- numeric(m)
  ni.vec <- numeric(m)
  
  sig02 <- (n0 / 12) * (1 - d3 / n3)
  
  for (i in seq_len(m)) {
    
    trt <- x[g == trt.names[i]]
    ni  <- length(trt)
    
    ni.vec[i] <- ni
    
    W[i] <-
      sum(rank(c(ctrl, trt),
               ties.method = "average")[(n0 + 1L):(n0 + ni)]) -
      ni * (ni + 1L) / 2
    
    mu[i] <- n0 * ni / 2
    
    sig2 <-
      (n0 * ni / 12) *
      (n0 + 1 - 3 * d2 / n2 - (n0 - 2) * d3 / n3)
    
    tau[i] <- sqrt(ni^2 * sig02 + sig2)
  }
  
  ## ------------------------------------------------------------
  ## Standardized statistics
  ## ------------------------------------------------------------
  
  z <- (W - mu) / tau
  
  ## ------------------------------------------------------------
  ## m = 1 special case: pmvnorm degenerates, use pnorm directly
  ## ------------------------------------------------------------
  
  if (m == 1L) {
    
    S <- switch(
      alternative,
      greater   =  z[1L],
      less      =  z[1L],
      two.sided = abs(z[1L])
    )
    
    p.global <- switch(
      alternative,
      greater   = pnorm(z[1L], lower.tail = FALSE),
      less      = pnorm(z[1L]),
      two.sided = 2 * pnorm(abs(z[1L]), lower.tail = FALSE)
    )
    
    pval <- p.global
    R    <- matrix(1, nrow = 1L, ncol = 1L)
    
  } else {
    
    ## ----------------------------------------------------------
    ## Correlation matrix
    ## ----------------------------------------------------------
    
    R <- diag(m)
    
    for (i in seq_len(m - 1L)) {
      for (j in (i + 1L):m) {
        covij    <- n0 * ni.vec[i] * ni.vec[j] / 12 * (1 - d3 / n3)
        R[i, j]  <- covij / (tau[i] * tau[j])
        R[j, i]  <- R[i, j]
      }
    }
    
    diag(R) <- 1
    
    if (anyNA(R))
      stop("correlation matrix contains missing values")
    
    eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    
    if (min(eig) < -1e-8)
      warning(
        "correlation matrix is not positive semidefinite; ",
        "p-values may be unreliable",
        call. = FALSE
      )
    
    ## ----------------------------------------------------------
    ## Steel statistic
    ## ----------------------------------------------------------
    
    S <- switch(
      alternative,
      greater   = max(z),
      less      = min(z),
      two.sided = max(abs(z))
    )
    
    ## ----------------------------------------------------------
    ## Global p-value
    ## ----------------------------------------------------------
    
    p.global <- switch(
      alternative,
      greater = {
        1 - mvtnorm::pmvnorm(upper = rep(S, m), corr = R)[1L]
      },
      less = {
        1 - mvtnorm::pmvnorm(lower = rep(S, m), corr = R)[1L]
      },
      two.sided = {
        1 - mvtnorm::pmvnorm(
          lower = rep(-S, m),
          upper = rep( S, m),
          corr  = R
        )[1L]
      }
    )
    
    ## ----------------------------------------------------------
    ## Steel-adjusted pairwise p-values
    ## ----------------------------------------------------------
    
    pval <- numeric(m)
    
    for (i in seq_len(m)) {
      pval[i] <- switch(
        alternative,
        greater = {
          1 - mvtnorm::pmvnorm(upper = rep(z[i], m), corr = R)[1L]
        },
        less = {
          1 - mvtnorm::pmvnorm(lower = rep(z[i], m), corr = R)[1L]
        },
        two.sided = {
          zi <- abs(z[i])
          1 - mvtnorm::pmvnorm(
            lower = rep(-zi, m),
            upper = rep( zi, m),
            corr  = R
          )[1L]
        }
      )
    }
    
  } # end m > 1
  
  ## ------------------------------------------------------------
  ## Result table
  ## ------------------------------------------------------------
  
  res.mat <- cbind(W = W, z = z, pval = pval)
  rownames(res.mat) <- paste(trt.names, control, sep = "-")
  
  if (output == "list") {
    
    res <- res.mat
    
  } else {
    
    res <- matrix(
      NA_real_,
      nrow = length(gn),
      ncol = length(gn),
      dimnames = list(gn, gn)
    )
    res[trt.names, control] <- pval
  }
  
  ## ------------------------------------------------------------
  ## p-value matrix
  ## ------------------------------------------------------------
  
  pmat <- matrix(
    NA_real_,
    nrow = length(gn),
    ncol = length(gn),
    dimnames = list(gn, gn)
  )
  
  pmat[trt.names, control] <- pval
  diag(pmat) <- 1
  
  attr(pmat, "lbl") <- apply(
    pmat,
    1,
    function(x)
      paste(rownames(pmat)[!is.na(x) & x < 0.05], collapse = ",")
  )
  
  ## ------------------------------------------------------------
  ## Return
  ## ------------------------------------------------------------
  
  out <- list(
    statistic = S,
    p.value   = p.global,
    res       = res,
    pmat      = pmat,
    corr      = R
  )
  
  class(out) <- c("rankTest", "htest")
  
  attr(out, "main")        <- "Steel test for multiple comparisons with a control"
  attr(out, "method")      <- "Steel"
  attr(out, "alternative") <- alternative
  attr(out, "output")      <- output
  attr(out, "data.name")   <- DG$data.name
  
  out
}

