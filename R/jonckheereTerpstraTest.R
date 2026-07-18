
#' Jonckheere-Terpstra Test for Ordered Alternatives
#'
#' A nonparametric test for monotonic trends across ordered independent
#' groups, assessing whether observations tend to increase (or decrease)
#' systematically with group order.
#'
#' The Jonckheere-Terpstra statistic is
#' \deqn{JT = \sum_{k<l} \sum_{ij} \left[ I(X_{ik} < X_{jl}) +
#' \frac{1}{2} I(X_{ik} = X_{jl}) \right]}{JT = sum_{k<l} sum_{ij}
#' [I(X_ik < X_jl) + 1/2 I(X_ik = X_jl)]}
#' where \eqn{i, j} index observations from ordered groups \eqn{k, l}.
#' Large values of the statistic indicate increasing trends across groups.
#'
#' Exact p-values are computed from the exact permutation distribution
#' using a dynamic programming recursion implemented in C++. Exact
#' inference is only valid without ties and, for practical reasons, is
#' offered for total sample sizes \eqn{n \le 100}.
#'
#' When ties are present or sample sizes are large, permutation p-values
#' can be computed by permuting group labels under the null hypothesis
#' (\code{method = "permutation"}); the number of permutations is
#' controlled by \code{R}, and the reported p-value uses the finite-sample
#' correction \eqn{(m + 1)/(R + 1)}. This approach remains valid in the
#' presence of ties.
#'
#' With \code{method = "asymptotic"} (the fallback of \code{"auto"} when
#' exact inference does not apply), a normal approximation with the
#' tie-corrected variance of Hollander and Wolfe (1999, Eq. 6.19) is used.
#'
#' @name jonckheereTerpstraTest
#' @aliases jonckheereTerpstraTest jonckheereTerpstraTest.default jonckheereTerpstraTest.formula
#'
#' @param x a numeric vector of observations, or a list of numeric vectors.
#' @param g a grouping variable corresponding to \code{x}, whose (factor)
#' level order defines the hypothesised ordering; ignored when \code{x} is
#' a list.
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of \code{"two.sided"} (default),
#' \code{"increasing"} or \code{"decreasing"}.
#' @param method a character string specifying the inference method, one of
#' \code{"auto"} (default), \code{"exact"}, \code{"permutation"} or
#' \code{"asymptotic"}. \code{"auto"} uses exact inference when possible
#' (no ties, \eqn{n \le 100}), otherwise the asymptotic approximation.
#' @param R the number of permutations, required when
#' \code{method = "permutation"}.
#' @param formula a formula of the form \code{response ~ group}.
#' @param data an optional data frame containing the variables in
#' \code{formula}.
#' @param subset an optional expression specifying a subset of observations.
#' @param na.action a function indicating how missing values should be
#' handled.
#' @param \dots further arguments passed to methods.
#' @return A list with class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the value of the Jonckheere-Terpstra statistic.}
#' \item{parameter}{the number of groups \code{k} and the total sample
#' size \code{n}.}
#' \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.}
#' \item{method}{a character string indicating the test performed and the
#' inference method used.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @references
#' Jonckheere, A. R. (1954) A distribution-free k-sample test against
#' ordered alternatives. \emph{Biometrika}, 41, 133--145.
#'
#' Terpstra, T. J. (1952) The asymptotic normality and consistency of
#' Kendall's test against trend, when ties are present in one ranking.
#' \emph{Indagationes Mathematicae}, 14, 327--333.
#'
#' Hollander, M. and Wolfe, D. A. (1999) \emph{Nonparametric Statistical
#' Methods}, 2nd ed., New York: Wiley.
#'
#' @seealso \code{\link{kruskal.test}}, \code{\link{cochranArmitageTest}}
#'
#' @examples
#' set.seed(1)
#' g <- ordered(rep(1:4, each = 10))
#' x <- rnorm(40) + 0.5 * as.numeric(g)
#'
#' jonckheereTerpstraTest(x, g)
#'
#' # with ties: permutation inference
#' x[1:2] <- mean(x[1:2])
#' jonckheereTerpstraTest(x, g, method = "permutation", R = 2000)
#'
#' coffee <- list(
#'   c_4 = c(447, 396, 383, 410),
#'   c_2 = c(438, 521, 468, 391, 504, 472),
#'   c_0 = c(513, 543, 506, 489, 407)
#' )
#' jonckheereTerpstraTest(coffee)
#'
#' # Hollander & Wolfe, Example 6.2:
#' # motivational effect of knowledge of performance
#' motiv <- list(
#'   no    = c(40, 35, 38, 43, 44, 41),
#'   rough = c(38, 40, 47, 44, 40, 42),
#'   acc   = c(48, 40, 45, 43, 46, 44))
#'
#' jonckheereTerpstraTest(motiv, alternative = "increasing")
#' ## exact one-sided p-value 0.0379 as in Hollander & Wolfe
#'
#' jonckheereTerpstraTest(motiv, alternative = "increasing",
#'                        method = "asymptotic")
#'
#' set.seed(42)
#' jonckheereTerpstraTest(motiv, alternative = "increasing",
#'                        method = "permutation", R = 10000)
#'
#' @family test.trend
#' @concept trend-test
#' @concept nonparametric
#' @concept k-sample
#'
#' @export
jonckheereTerpstraTest <- function(x, ...)
  UseMethod("jonckheereTerpstraTest")



#' @rdname jonckheereTerpstraTest
#' @export
jonckheereTerpstraTest.formula <- function(formula,
                                           data,
                                           subset,
                                           na.action,
                                           ...) {

  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")

  # capture subset / na.action here, before they are evaluated
  subset_expr <- if (!missing(subset)) substitute(subset) else NULL
  na_expr     <- if (!missing(na.action)) substitute(na.action) else NULL

  pf <- resolveFormula(
    formula   = formula,
    data      = data,
    subset    = subset_expr,
    na.action = na_expr,
    allowed   = "n-sample-independent"
  )

  y <- jonckheereTerpstraTest(x = pf$x, g = pf$group, ...)

  y$data.name <- pf$data.name

  y
}



#' @rdname jonckheereTerpstraTest
#' @export
jonckheereTerpstraTest.default <- function(
    x,
    g,
    alternative = c("two.sided", "increasing", "decreasing"),
    method      = c("auto", "exact", "permutation", "asymptotic"),
    R           = NULL,
    ...
) {

  alternative <- match.arg(alternative)
  method      <- match.arg(method)

  DG <- resolveGroups(x, g)

  x <- DG$x
  g <- ordered(DG$groups)

  n <- DG$n
  k <- DG$k

  DNAME <- DG$data.name

  # order by group so that observations form contiguous blocks
  ord <- order(g)
  x   <- x[ord]
  g   <- g[ord]

  gsize  <- as.integer(table(g))
  cgsize <- c(0L, cumsum(gsize))

  TIES <- anyDuplicated(x) > 0L

  JT        <- .jtStatistic(x, cgsize)
  STATISTIC <- c(JT = JT)
  JT_int    <- as.integer(round(JT))

  ## resolve method -------------------------------------------------------

  if (method == "auto")
    method <- if (!TIES && n <= 100L) "exact" else "asymptotic"

  if (method == "exact" && TIES) {
    warning("exact inference not available with ties; ",
            "falling back to asymptotic approximation")
    method <- "asymptotic"
  }

  if (method == "exact" && n > 100L)
    warning("exact inference requested for n = ", n, " > 100; ",
            "this may be slow or fail")

  if (method == "permutation" && is.null(R))
    stop("'R' must be specified when method = \"permutation\"")

  if (!is.null(R) && method != "permutation")
    warning("'R' is ignored when method != \"permutation\"")

  ## p-value --------------------------------------------------------------

  METHOD <- "Jonckheere-Terpstra test for ordered alternatives"

  # exact mean of JT under the null hypothesis
  muJT <- (n^2 - sum(gsize^2)) / 4

  if (method == "permutation") {

    PVAL <- .jtPvaluePerm(x = x, g = g, observed = JT, mu = muJT,
                          R = R, alternative = alternative)

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

    # tie-corrected asymptotic variance,
    # Hollander and Wolfe (1999), Eq. 6.19
    tie_tab <- as.numeric(table(x))

    a1 <- n^2 * (2 * n + 3)
    b1 <- sum(gsize^2 * (2 * gsize + 3))
    c1 <- sum(tie_tab * (tie_tab - 1) * (2 * tie_tab + 5))

    a2 <- sum(gsize * (gsize - 1) * (gsize - 2))
    b2 <- sum(tie_tab * (tie_tab - 1) * (tie_tab - 2))

    a3 <- sum(gsize * (gsize - 1))
    b3 <- sum(tie_tab * (tie_tab - 1))

    sigma2 <- (a1 - b1 - c1) / 72 +
      (a2 * b2) / (36 * n * (n - 1) * (n - 2)) +
      (a3 * b3) / (8 * n * (n - 1))

    z <- (JT - muJT) / sqrt(sigma2)

    PVAL <- switch(
      alternative,
      "increasing" = pnorm(z, lower.tail = FALSE),
      "decreasing" = pnorm(z),
      "two.sided"  = min(2 * pnorm(abs(z), lower.tail = FALSE), 1)
    )

    METHOD <- paste(METHOD, "(asymptotic)")
  }

  structure(
    list(
      statistic   = STATISTIC,
      parameter   = c(k = k, n = n),
      p.value     = as.numeric(PVAL),
      alternative = alternative,
      method      = METHOD,
      data.name   = DNAME
    ),
    class = "htest"
  )
}



# == internal helper functions ============================================


# JT statistic for observations ordered by group;
# cgsize is c(0, cumsum(group sizes))

.jtStatistic <- function(x, cgsize) {

  n <- length(x)
  k <- length(cgsize) - 1L

  JT <- 0

  for (i in seq_len(k - 1L)) {

    idx1 <- (cgsize[i] + 1L):cgsize[i + 1L]
    idx2 <- (cgsize[i + 1L] + 1L):n

    JT <- JT + sum(outer(x[idx1], x[idx2],
                         function(a, b) (a < b) + 0.5 * (a == b)))
  }

  JT
}


# permutation p-value with finite-sample correction

.jtPvaluePerm <- function(x, g, observed, mu, R, alternative) {

  gsize  <- as.integer(table(g))
  cgsize <- c(0L, cumsum(gsize))

  perm_stats <- vapply(seq_len(R), function(b) {
    xp <- x[order(sample(g))]
    .jtStatistic(xp, cgsize)
  }, numeric(1))

  m <- switch(
    alternative,
    "increasing" = sum(perm_stats >= observed),
    "decreasing" = sum(perm_stats <= observed),
    # center at the exact null mean, not the empirical one
    "two.sided"  = sum(abs(perm_stats - mu) >= abs(observed - mu))
  )

  (m + 1) / (R + 1)
}


# exact null distribution of JT via DP recursion (C++)

.jtpdf <- function(gsize) {
  jtpdf_cpp(as.integer(gsize))
}
