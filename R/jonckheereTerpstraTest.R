
#' Jonckheere-Terpstra Test for Detecting Ordered Differences Across Independent Groups
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
#' using dynamic programming recursions implemented in C++, with and
#' without ties.
#'
#' Without ties the null distribution depends on the group sizes alone and
#' is obtained from the classical recursion; it is offered for total sample
#' sizes \eqn{n \le 100}.
#'
#' With ties the statistic depends on the data through the table of counts
#' of group by distinct value, which under the null hypothesis follows the
#' multiple hypergeometric distribution with the group sizes and the tie
#' counts as its margins. The distribution is built by splitting off one
#' row of that table at a time, the state being the counts not yet
#' assigned. Its cost grows with \eqn{\prod_a (c_a + 1)}, the \eqn{c_a}
#' being the tie counts, so it is cheapest where ties are heaviest, and
#' prohibitive where they are few and the sample is large. `"auto"`
#' therefore turns to the asymptotic approximation beyond a cost of about
#' \eqn{10^7} table cells, and `method = "exact"` warns above
#' \eqn{2 \times 10^8} and falls back to the approximation as well. Note
#' that with ties the statistic and the support of its distribution are
#' half-integral.
#'
#' For large samples permutation p-values can be computed by permuting
#' group labels under the null hypothesis (`method = "permutation"`);
#' the number of permutations is controlled by `R`, and the reported
#' p-value uses the finite-sample correction \eqn{(m + 1)/(R + 1)}.
#'
#' Two-sided p-values are the smaller one-sided p-value doubled, in each of
#' the three methods. The null distribution of the statistic is symmetric
#' only for equal group sizes, so this is not the same as counting the
#' values lying at least as far from the null mean as the observed one.
#'
#' With `method = "asymptotic"` (the fallback of `"auto"` when
#' exact inference does not apply), a normal approximation with the
#' tie-corrected variance of Hollander and Wolfe (1999, Eq. 6.19) is used.
#'
#' @name jonckheereTerpstraTest
#' @aliases jonckheereTerpstraTest jonckheereTerpstraTest.default jonckheereTerpstraTest.formula
#'
#' @param x a numeric vector of observations, or a list of numeric vectors.
#' @param g a grouping variable corresponding to `x`, whose (factor)
#' level order defines the hypothesised ordering; ignored when `x` is
#' a list.
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of `"two.sided"` (default),
#' `"increasing"` or `"decreasing"`.
#' @param method a character string specifying the inference method, one of
#' `"auto"` (default), `"exact"`, `"permutation"` or
#' `"asymptotic"`. `"auto"` uses exact inference where it is
#' affordable (without ties \eqn{n \le 100}, with ties a cost below
#' \eqn{10^7} table cells), otherwise the asymptotic approximation.
#' @param R the number of permutations, a single positive integer,
#' required when `method = "permutation"`.
#' @param formula a formula of the form `response ~ group`.
#' @param data an optional data frame containing the variables in
#' `formula`.
#' @param subset an optional expression specifying a subset of observations,
#' evaluated in `data` (`subset = dose > 0.5`), as in [kruskal.test()].
#' @param na.action a function indicating how missing values should be
#' handled. Defaults to [na.omit()].
#' @param \dots further arguments passed to methods.
#' @return A list with class `"htest"` containing the following
#' components:
#' \item{statistic}{the value of the Jonckheere-Terpstra statistic.}
#' \item{parameter}{the number of groups `k` and the total sample
#' size `n`.}
#' \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.}
#' \item{method}{a character string indicating the test performed and the
#' inference method used.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @references
#' Jonckheere, A. R. (1954) A distribution-free k-sample test against
#' ordered alternatives. *Biometrika*, 41, 133--145.
#'
#' Terpstra, T. J. (1952) The asymptotic normality and consistency of
#' Kendall's test against trend, when ties are present in one ranking.
#' *Indagationes Mathematicae*, 14, 327--333.
#'
#' Hollander, M. and Wolfe, D. A. (1999) *Nonparametric Statistical
#' Methods*, 2nd ed., New York: Wiley.
#'
#' @seealso [kruskal.test()]
#'
#' @examples
#' set.seed(1)
#' g <- ordered(rep(1:4, each = 10))
#' x <- rnorm(40) + 0.5 * as.numeric(g)
#'
#' jonckheereTerpstraTest(x, g)
#'
#' # with ties: exact inference as long as it is affordable,
#' # permutation inference otherwise
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
#' ## exact one-sided p-value 0.0210, the data being tied. Hollander and
#' ## Wolfe report 0.0231 from the tie-free null distribution and 0.0207
#' ## from the tie-corrected normal approximation
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
                                           na.action = na.omit,
                                           ...) {

  # formula, data and subset are forwarded unevaluated, so that 'subset' is
  # evaluated in 'data' as in kruskal.test(); y ~ a:b compares the cells
  pf <- resolveFormulaFromCall(
    allowed   = "n-sample-independent",
    na.action = na.action
  )

  y <- jonckheereTerpstraTest(x = pf$x, g = pf$group, ...)

  y$data.name <- pf$dataName

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

  DNAME <- DG$dataName

  # order by group so that observations form contiguous blocks
  ord <- order(g)
  x   <- x[ord]
  g   <- g[ord]

  gsize  <- as.integer(table(g))
  cgsize <- c(0L, cumsum(gsize))

  tieTab <- as.integer(table(x))
  TIES   <- any(tieTab > 1L)

  JT        <- .jtStatistic(x, cgsize)
  STATISTIC <- c(JT = JT)
  JT_int    <- as.integer(round(JT))

  ## resolve method -------------------------------------------------------

  # cost of the exact distribution with ties, in table cells
  cells <- if (TIES) .jtTiesCells(gsize, tieTab) else 0

  # without ties the cost is governed by the sample size alone, with ties
  # by the table the recursion has to walk
  affordable <- if (TIES) cells <= .jtTiesAutoCells else n <= 100L

  if (method == "auto")
    method <- if (affordable) "exact" else "asymptotic"

  if (method == "exact" && cells > .jtTiesMaxCells) {

    warning(gettextf(
      paste("exact inference with ties would need %.1e table cells,",
            "the limit is %.1e; falling back to the asymptotic",
            "approximation, method = \"permutation\" is the",
            "distribution-free alternative"),
      cells, .jtTiesMaxCells), call. = FALSE)

    method <- "asymptotic"
  }

  if (method == "exact" && !TIES && n > 100L)
    warning("exact inference requested for n = ", n, " > 100; ",
            "this may be slow or fail")

  if (method == "permutation") {

    if (is.null(R))
      stop("'R' must be specified when method = \"permutation\"")

    if (!is.numeric(R) || length(R) != 1L || !is.finite(R) ||
        R < 1 || R != round(R))
      stop("'R' must be a single positive integer")

    R <- as.integer(R)
  }

  if (!is.null(R) && method != "permutation")
    warning("'R' is ignored when method != \"permutation\"")

  ## p-value --------------------------------------------------------------

  METHOD <- "Jonckheere-Terpstra test for ordered alternatives"

  # exact mean of JT under the null hypothesis
  muJT <- (n^2 - sum(gsize^2)) / 4

  if (method == "permutation") {

    PVAL <- .jtPvaluePerm(x = x, g = g, observed = JT,
                          R = R, alternative = alternative)

    METHOD <- paste0(METHOD, " (permutation, R = ", R, ")")

  } else if (method == "exact") {

    # with ties the distribution runs over 2 * JT, the half weights of the
    # ties making the statistic half-integral
    if (TIES) {
      pdf <- .jtpdfTies(gsize, tieTab)
      at  <- as.integer(round(2 * JT)) + 1L
    } else {
      pdf <- .jtpdf(gsize)
      at  <- JT_int + 1L
    }

    lower_tail <- sum(pdf[seq_len(at)])
    upper_tail <- sum(pdf[at:length(pdf)])

    PVAL <- switch(
      alternative,
      "increasing" = upper_tail,
      "decreasing" = lower_tail,
      "two.sided"  = min(2 * min(lower_tail, upper_tail), 1)
    )

    METHOD <- paste(METHOD, if (TIES) "(exact, ties)" else "(exact)")

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

    # constant data leave nothing to compare: the statistic equals its null
    # expectation with probability one, and z would be 0/0
    if (sigma2 <= 0) {

      PVAL <- 1

    } else {

      z <- (JT - muJT) / sqrt(sigma2)

      PVAL <- switch(
        alternative,
        "increasing" = pnorm(z, lower.tail = FALSE),
        "decreasing" = pnorm(z),
        "two.sided"  = min(2 * pnorm(abs(z), lower.tail = FALSE), 1)
      )
    }

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


# permutation p-value with finite-sample correction; two-sided as the
# doubled smaller tail, as in the exact and the asymptotic branch

.jtPvaluePerm <- function(x, g, observed, R, alternative) {

  gsize  <- as.integer(table(g))
  cgsize <- c(0L, cumsum(gsize))

  perm_stats <- vapply(seq_len(R), function(b) {
    xp <- x[order(sample(g))]
    .jtStatistic(xp, cgsize)
  }, numeric(1))

  upper_tail <- (sum(perm_stats >= observed) + 1) / (R + 1)
  lower_tail <- (sum(perm_stats <= observed) + 1) / (R + 1)

  switch(
    alternative,
    "increasing" = upper_tail,
    "decreasing" = lower_tail,
    "two.sided"  = min(2 * min(lower_tail, upper_tail), 1)
  )
}


# exact null distribution of JT via DP recursion (C++), no ties

.jtpdf <- function(gsize) {
  jtpdf_cpp(as.integer(gsize))
}


# exact null distribution of 2 * JT via DP over the tables of group by
# distinct value (C++), ties included

.jtpdfTies <- function(gsize, tieTab) {

  # the recursion holds the states it can reach in memory, so the cost is
  # checked here as well: a direct call that skipped the test above would
  # otherwise take the session down with it
  if (.jtTiesCells(gsize, tieTab) > .jtTiesMaxCells)
    stop("the table of group by distinct value is too large to enumerate",
         call. = FALSE)

  jtpdfTies_cpp(as.integer(gsize), as.integer(tieTab))
}


# cost of that recursion: states times the support of 2 * JT.
# Roughly 1.5e7 cells per second, and the memory follows the same measure

.jtTiesCells <- function(gsize, tieTab) {

  gsize <- as.numeric(gsize)
  ahead <- sum(gsize) - cumsum(gsize)

  prod(as.numeric(tieTab) + 1) * (2 * sum(gsize * ahead) + 1)
}


# beyond this many cells "auto" prefers the asymptotic approximation,
# beyond the second constant "exact" refuses the computation altogether

.jtTiesAutoCells <- 1e7

.jtTiesMaxCells <- 2e8
