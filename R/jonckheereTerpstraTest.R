
#' Jonckheere-Terpstra Test for Ordered Alternatives
#'
#' A nonparametric test for monotonic trends across ordered independent groups.
#'
#' @description
#' Performs the Jonckheere-Terpstra test for ordered alternatives across
#' \eqn{k} independent groups.
#'
#' The test assesses whether observations tend to increase (or decrease)
#' systematically with group order.
#'
#' Exact p-values are computed whenever feasible and valid
#' (small samples without ties). For larger samples or tied data,
#' permutation or asymptotic approximations are used.
#'
#' @details
#' The Jonckheere-Terpstra statistic is:
#'
#' \deqn{
#' JT =
#' \sum_{k<l}
#' \sum_{ij}
#' \left[
#' I(X_{ik} < X_{jl})
#' +
#' \frac12 I(X_{ik} = X_{jl})
#' \right]
#' }
#'
#' where \eqn{i,j} index observations from ordered groups
#' \eqn{k,l}.
#'
#' Large values of the statistic indicate increasing trends across groups.
#'
#' ## Exact inference
#'
#' Exact p-values are computed from the exact permutation distribution
#' using a dynamic programming recursion implemented in C++.
#'
#' Exact inference is only available:
#' \itemize{
#'   \item without ties,
#'   \item for total sample sizes \eqn{n \le 100}.
#' }
#'
#' ## Permutation inference
#'
#' When ties are present or sample sizes are large, permutation p-values
#' can be computed by permuting group labels under the null hypothesis.
#' The number of permutations is controlled by \code{R}.
#'
#' This approach remains valid in the presence of ties.
#'
#' ## Asymptotic approximation
#'
#' If \code{method = "asymptotic"} or no other method applies,
#' a normal approximation with tie-corrected variance is applied.
#'
#' @name jonckheereTerpstraTest
#' @aliases jonckheereTerpstraTest jonckheereTerpstraTest.default jonckheereTerpstraTest.formula
#'
#' @param x Numeric vector of observations,
#'   or a list of numeric vectors.
#'
#' @param g Grouping variable.
#'   Ignored if \code{x} is a list.
#'
#' @param alternative Character string specifying the alternative:
#'   \code{"increasing"},
#'   \code{"decreasing"},
#'   or \code{"two.sided"}.
#'
#' @param method Character string specifying the inference method:
#'   \code{"auto"} (default), \code{"exact"}, \code{"permutation"},
#'   or \code{"asymptotic"}.
#'   \code{"auto"} uses exact inference when possible
#'   (no ties, \eqn{n \le 100}), otherwise asymptotic.
#'
#' @param R Number of permutations. Required when \code{method = "permutation"}.
#'
#' @param formula Formula interface of the form \code{lhs ~ rhs}.
#'
#' @param data Optional data frame for the formula interface.
#'
#' @param subset Optional subset expression.
#'
#' @param na.action NA handling function.
#'
#' @param \dots Further arguments passed to methods.
#'
#' @return
#' An object of class \code{"htest"} containing:
#'
#' \itemize{
#'   \item statistic,
#'   \item p-value,
#'   \item method,
#'   \item alternative,
#'   \item sample size information.
#' }
#'
#' @references
#'
#' Jonckheere, A. R. (1954).
#' A distribution-free k-sample test against ordered alternatives.
#' \emph{Biometrika}, 41, 133--145.
#'
#' Terpstra, T. J. (1952).
#' The asymptotic normality and consistency of Kendall's test
#' against trend.
#' \emph{Indagationes Mathematicae}, 14, 327--333.
#'
#' @examples
#'
#' set.seed(1)
#'
#' g <- ordered(rep(1:4, each = 10))
#'
#' x <- rnorm(40) + 0.5 * as.numeric(g)
#'
#' jonckheereTerpstraTest(x, g)
#'
#' x[1:2] <- mean(x[1:2])
#'
#' jonckheereTerpstraTest(x, g, method = "permutation", R = 2000)
#'
#' coffee <- list(
#'   c_4 = c(447,396,383,410),
#'   c_2 = c(438,521,468,391,504,472),
#'   c_0 = c(513,543,506,489,407)
#' )
#'
#' jonckheereTerpstraTest(coffee)
#'
#' # Hollander-Wolfe-Chicken Example 6.2 
#' #   Motivational Effect of Knowledge of Performance
#' 
#' motiv <- list(
#'   no    = c(40,35,38,43,44,41),
#'   rough = c(38,40,47,44,40,42),
#'   acc   = c(48,40,45,43,46,44))
#' 
#' jonckheereTerpstraTest(motiv, method = "auto")    # here: exact
#' jonckheereTerpstraTest(motiv, method = "exact")   # exact:       p ~ 0.038
#' jonckheereTerpstraTest(motiv, method = "asymp")   # asymptotic:  p ~ 0.032
#' 
#' set.seed(42)
#' jonckheereTerpstraTest(motiv, method = "perm", R= 10000)   # permutation: p ~ 0.038
#' 
#'
#' @rdname jonckheereTerpstraTest
#' @family test.trend
#' @concept hypothesis-testing
#' @concept nonparametric


#' @export
jonckheereTerpstraTest <- function(x, ...)
  UseMethod("jonckheereTerpstraTest")



#' @rdname jonckheereTerpstraTest
#' @export
jonckheereTerpstraTest.formula <- function(formula,
                                           data,
                                           subset,
                                           na.action = na.pass,
                                           ...) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")
  
  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "n.sample.independent"
  )
  
  if (!missing(data))
    args$data <- data
  
  if (!missing(subset))
    args$subset <- substitute(subset)
  
  d <- do.call(bedrock::resolveFormula, args)
  
  res <- jonckheereTerpstraTest.default(
    x = d$x,
    g = d$g,
    ...
  )
  
  res$data.name <- paste(
    deparse(formula[[2L]]),
    "by",
    deparse(formula[[3L]])
  )
  
  res
}



#' @rdname jonckheereTerpstraTest
#' @export
jonckheereTerpstraTest.default <- function(
    x,
    g,
    alternative = c("increasing", "decreasing", "two.sided"),
    method      = c("auto", "exact", "permutation", "asymptotic"),
    R           = NULL,
    ...
) {
  
  alternative <- match.arg(alternative)
  method      <- match.arg(method)
  
  ## ---------------------------------------------------------------------
  ## List interface
  ## ---------------------------------------------------------------------
  
  if (is.list(x)) {
    
    if (length(x) < 2L)
      stop("'x' must contain at least two groups")
    
    if (!missing(g))
      warning("'g' ignored because 'x' is a list")
    
    DNAME <- deparse1(substitute(x))
    
    x <- lapply(x, function(z) z[is.finite(z)])
    
    if (!all(vapply(x, is.numeric, logical(1))))
      stop("all groups must be numeric")
    
    sizes <- lengths(x)
    
    if (any(sizes == 0L))
      stop("all groups must contain observations")
    
    g <- ordered(rep.int(seq_along(x), sizes))
    x <- unlist(x)
    
  } else {
    
    ## -------------------------------------------------------------------
    ## Standard interface
    ## -------------------------------------------------------------------
    
    if (missing(g))
      stop("'g' is missing")
    
    if (length(x) != length(g))
      stop("'x' and 'g' must have same length")
    
    DNAME <- paste(
      deparse1(substitute(x)),
      "by",
      deparse1(substitute(g))
    )
    
    ok <- complete.cases(x, g)
    x  <- x[ok]
    g  <- g[ok]
    
    if (!(is.numeric(g) || is.factor(g) || is.ordered(g)))
      stop("'g' must be numeric or a factor")
    
    g <- ordered(g)
    
  }
  
  n <- length(x)
  
  if (n < 2L)
    stop("not enough observations")
  
  k <- nlevels(g)
  
  if (k < 2L)
    stop("all observations belong to same group")
  
  ## ---------------------------------------------------------------------
  ## Ordering
  ## ---------------------------------------------------------------------
  
  ord <- order(g)
  x   <- x[ord]
  g   <- g[ord]
  
  gsize  <- as.integer(table(g))
  cgsize <- c(0L, cumsum(gsize))
  
  TIES <- anyDuplicated(x) > 0L
  
  ## ---------------------------------------------------------------------
  ## JT statistic
  ## ---------------------------------------------------------------------
  
  JT <- 0
  
  for (i in seq_len(k - 1L)) {
    
    idx1 <- (cgsize[i] + 1L):cgsize[i + 1L]
    idx2 <- (cgsize[i + 1L] + 1L):n
    
    JT <- JT +
      sum(
        outer(
          x[idx1],
          x[idx2],
          function(a, b) (a < b) + 0.5 * (a == b)
        )
      )
  }
  
  STATISTIC <- c(JT = JT)
  JT_int    <- as.integer(round(JT))
  
  ## ---------------------------------------------------------------------
  ## Resolve method
  ## ---------------------------------------------------------------------
  
  if (method == "auto")
    method <- if (!TIES && n <= 100L) "exact" else "asymptotic"
  
  if (method == "exact" && TIES) {
    warning(
      "exact inference not available with ties; ",
      "falling back to asymptotic approximation"
    )
    method <- "asymptotic"
  }
  
  if (method == "exact" && n > 100L) {
    warning(
      "exact inference requested for n = ", n, " > 100; ",
      "this may be slow or fail"
    )
  }
  
  if (method == "permutation" && is.null(R))
    stop("'R' must be specified when method = \"permutation\"")
  
  if (!is.null(R) && method != "permutation")
    warning("'R' is ignored when method != \"permutation\"")
  
  
  ## ---------------------------------------------------------------------
  ## P-value calculation
  ## ---------------------------------------------------------------------
  
  METHOD <- "Jonckheere-Terpstra Test for Ordered Alternatives"
  
  if (method == "permutation") {
    
    PVAL <- .permPvalue(
      x           = x,
      g           = g,
      observed    = JT,
      R           = R,
      alternative = alternative
    )
    
    METHOD <- paste0(METHOD, " (permutation, R = ", R, ")")
    
  } else if (method == "exact") {
    
    pdf <- .jtpdf(gsize)
    
    lower_tail <- sum(pdf[seq_len(JT_int + 1L)])
    upper_tail <- if (JT_int == 0L) 1 else 1 - sum(pdf[seq_len(JT_int)])
    
    PVAL <- switch(
      alternative,
      "increasing" = upper_tail,
      "decreasing" = lower_tail,
      "two.sided"  = min(2 * min(lower_tail, upper_tail), 1)
    )
    
    METHOD <- paste(METHOD, "(exact)")
    
  } else {
    
    ## -------------------------------------------------------------------
    ## Tie-corrected asymptotic variance (Hollander & Wolfe, eq. 6.19)
    ##
    ## Var(JT) = (A - B - C) / 72
    ##         + (D * E) / (36 * N * (N-1) * (N-2))
    ##         + (F * G) / (8  * N * (N-1))
    ##
    ## A = N^2 (2N + 3)
    ## B = sum n_i^2 (2 n_i + 3)
    ## C = sum t_j (t_j - 1)(2 t_j + 5)       [tie correction]
    ## D = sum n_i (n_i - 1)(n_i - 2)
    ## E = sum t_j (t_j - 1)(t_j - 2)         [tie correction]
    ## F = sum n_i (n_i - 1)
    ## G = sum t_j (t_j - 1)                   [tie correction]
    ## -------------------------------------------------------------------
    
    tie_tab <- as.numeric(table(x))
    
    A <- n^2 * (2 * n + 3)
    B <- sum(gsize^2 * (2 * gsize + 3))
    C <- sum(tie_tab * (tie_tab - 1) * (2 * tie_tab + 5))
    D <- sum(gsize * (gsize - 1) * (gsize - 2))
    E <- sum(tie_tab * (tie_tab - 1) * (tie_tab - 2))
    F <- sum(gsize * (gsize - 1))
    G <- sum(tie_tab * (tie_tab - 1))
    
    sigma2 <- (A - B - C) / 72 +
      (D * E) / (36 * n * (n - 1) * (n - 2)) +
      (F * G) / (8  * n * (n - 1))
    
    muJT <- (n^2 - sum(gsize^2)) / 4
    
    z <- (JT - muJT) / sqrt(sigma2)
    
    PVAL <- switch(
      alternative,
      "increasing" = stats::pnorm(z, lower.tail = FALSE),
      "decreasing" = stats::pnorm(z),
      "two.sided"  = min(2 * stats::pnorm(abs(z), lower.tail = FALSE), 1)
    )
    
    METHOD <- paste(METHOD, "(asymptotic)")
  }
  
  ## ---------------------------------------------------------------------
  ## Return object
  ## ---------------------------------------------------------------------
  
  rval <- list(
    statistic   = STATISTIC,
    parameter   = c(k = k, n = n),
    p.value     = as.numeric(PVAL),
    alternative = alternative,
    method      = METHOD,
    data.name   = DNAME
  )
  
  class(rval) <- "htest"
  
  rval
}



# == internal helper functions ================================================


# Permutation test

.permPvalue <- function(x, g, observed, R, alternative) {
  
  perm_stats <- numeric(R)
  
  for (b in seq_len(R)) {
    
    gperm <- sample(g)
    ordp  <- order(gperm)
    
    xp <- x[ordp]
    gp <- gperm[ordp]
    
    gsizep  <- as.integer(table(gp))
    cgsizep <- c(0L, cumsum(gsizep))
    
    stat <- 0
    
    for (i in seq_len(length(gsizep) - 1L)) {
      
      idx1 <- (cgsizep[i] + 1L):cgsizep[i + 1L]
      idx2 <- (cgsizep[i + 1L] + 1L):length(xp)
      
      stat <- stat +
        sum(
          outer(
            xp[idx1],
            xp[idx2],
            function(a, b) (a < b) + 0.5 * (a == b)
          )
        )
    }
    
    perm_stats[b] <- stat
  }
  
  if (alternative == "increasing") {
    mean(perm_stats >= observed)
  } else if (alternative == "decreasing") {
    mean(perm_stats <= observed)
  } else {
    mu <- mean(perm_stats)
    mean(abs(perm_stats - mu) >= abs(observed - mu))
  }
}



# Exact null distribution of JT via DP recursion

.jtpdf <- function(gsize) {
  jtpdf_cpp(as.integer(gsize))
}

