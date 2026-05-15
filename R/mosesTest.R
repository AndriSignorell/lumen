
#' Moses Test of Extreme Reactions
#'
#' A nonparametric test comparing the spread (range) of two independent
#' groups, assessing whether the control group shows greater variability
#' than the experiment group.
#'
#' @description
#' Performs the Moses test of extreme reactions, a distribution-free
#' nonparametric test for the difference in extremity of scores between two
#' independent groups. Scores from both groups are pooled and converted to
#' ranks. The test statistic is the span of the control group's ranks
#' (range + 1). An exact one-tailed probability is computed for the raw
#' span, then recomputed after trimming \code{h} extreme ranks from each end
#' of the control group.
#'
#' @details
#' **Distributional derivation.**
#' Under the null hypothesis, any assignment of \eqn{n_k + n_e} pooled ranks
#' to the two groups is equally likely. The exact distribution of the span
#' \eqn{S} of the control group (a subset of size \eqn{n_k} drawn without
#' replacement from \eqn{\{1, \ldots, N\}}) is given by Moses (1952) as:
#'
#' \deqn{P(S \le s) = \frac{1}{\binom{N}{n_k}}
#'   \sum_{i=0}^{\min(s - n_k,\, n_e)}
#'   \binom{i + n_k - 2}{i}\,
#'   \binom{n_e + 1 - i}{n_e - i}}
#'
#' This is a cumulative distribution function (CDF):
#' \eqn{P(\mathrm{span} \le s)}.
#'
#' The one-tailed p-value for the observed span \eqn{S_{obs}} is:
#'
#' \deqn{
#' p = P(S \ge S_{obs})
#'   = 1 - P(S \le S_{obs} - 1)
#' }
#'
#' This distinction matters because some textbook presentations write the
#' summation formula as \eqn{P(S = s)}, whereas the upper-tail probability
#' required for hypothesis testing is obtained from the cumulative form.
#'
#' **Generalisation with trimming (\eqn{h > 0}).**
#' After removing \eqn{h} extreme ranks from each end of the sorted control
#' group ranks, the effective control-group size becomes
#' \eqn{n_k^\prime = n_k - 2h}, yielding:
#'
#' \deqn{
#' P(S_h \le s) =
#' \frac{1}{\binom{N}{n_k}}
#' \sum_{i=0}^{\min(s - n_k^\prime,\, n_e)}
#' \binom{i + n_k^\prime - 2}{i}\,
#' \binom{n_e + 2h + 1 - i}{n_e - i}
#' }
#'
#' **Tie handling.**
#' The exact combinatorial distribution assumes a continuous underlying
#' distribution (no ties). By default,
#' \code{ties.method = "first"} is passed to \code{\link{rank}},
#' producing integer-valued ranks that remain compatible with the exact
#' formula. Tied observations are ordered according to their occurrence in
#' the pooled sample, matching SPSS behaviour.
#'
#' Alternative tie methods such as \code{"average"} may produce fractional
#' mid-ranks. In such cases the span is rounded upward via
#' \code{ceiling()} to retain an integer-valued test statistic, but the
#' resulting p-values should be regarded as approximate rather than exact.
#'
#' **Numerical stability.**
#' The exact distribution is evaluated entirely in log-space using
#' \code{\link{lchoose}} and a log-sum-exp transformation, avoiding
#' overflow for large sample sizes.
#'
#' @name mosesTest
#' @aliases mosesTest mosesTest.default mosesTest.formula
#'
#' @param x Numeric vector of observations from the control group.
#'   Non-finite values are removed.
#' @param y Numeric vector of observations from the experiment group.
#'   Non-finite values are removed.
#' @param formula A formula of the form \code{lhs ~ rhs}, where
#'   \code{lhs} gives the observations and \code{rhs} the grouping
#'   variable with exactly two levels.
#' @param data Optional data frame containing the variables in
#'   \code{formula}.
#' @param subset Optional vector specifying a subset of observations.
#' @param na.action Function specifying how missing values are handled.
#' @param extreme Non-negative integer \eqn{h}. Number of extreme ranks
#'   trimmed from each end of the control group before recomputing the
#'   span. If \code{NULL}, defaults to
#'   \code{max(floor(0.05 * length(x)), 1)}.
#' @param ties.method Character string passed to \code{\link{rank}}.
#'   Default \code{"first"} preserves integer-valued ranks and exact
#'   combinatorial validity.
#' @param \dots Further arguments passed to or from methods.
#'
#' @return An object of class \code{"mosesTestResult"} inheriting from
#'   \code{"htest"} with components:
#'   \describe{
#'     \item{\code{statistic}}{Named vector containing
#'       \code{S_raw} and \code{S_trimmed}.}
#'     \item{\code{p.value}}{Named vector containing
#'       \code{p_raw} and \code{p_trimmed}.}
#'     \item{\code{extreme}}{Effective trimming parameter \eqn{h}.}
#'     \item{\code{parameter}}{Vector containing
#'       \code{n_control} and \code{n_experiment}.}
#'     \item{\code{null.value}}{Null hypothesis of equal extremity.}
#'     \item{\code{alternative}}{Character string:
#'       \code{"greater"}.}
#'     \item{\code{method}}{Character string describing the test.}
#'     \item{\code{data.name}}{Character string describing the data.}
#'   }
#'
#' @seealso
#' \code{\link{wilcox.test}},
#' \code{\link{ks.test}},
#' \code{\link{ansari.test}}
#'
#' @references
#' Moses, L.E. (1952).
#' A two-sample test.
#' \emph{Psychometrika},
#' \strong{17}, 239--247.
#' \doi{10.1007/BF02288735}
#'
#' @examples
#' x <- c(0.80, 0.83, 1.89, 1.04, 1.45,
#'        1.38, 1.91, 1.64, 0.73, 1.46)
#'
#' y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
#'
#' mosesTest(x, y)
#'
#' set.seed(1479)
#'
#' x2 <- sample(1:20, 10, replace = TRUE)
#' y2 <- sample(5:25, 6, replace = TRUE)
#'
#' mosesTest(x2, y2, extreme = 2)
#'
#' df <- data.frame(
#'   score = c(x, y),
#'   group = factor(rep(
#'     c("control", "experiment"),
#'     c(length(x), length(y))
#'   ))
#' )
#'
#' mosesTest(score ~ group, data = df)
#'
#' @rdname mosesTest
#' @family test.scale
#' @concept hypothesis-testing
#' @concept nonparametric

#' @export
mosesTest <- function(x, ...)
  UseMethod("mosesTest")


# -------------------------------------------------------------------------
# Formula method
# -------------------------------------------------------------------------

#' @rdname mosesTest
#' @export
mosesTest.formula <- function(formula,
                              data,
                              subset,
                              na.action = na.pass,
                              ...) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")
  
  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "two.sample.independent"
  )
  
  if (!missing(data))
    args$data <- data
  
  if (!missing(subset))
    args$subset <- substitute(subset)
  
  d <- do.call(bedrock::resolveFormula, args)
  
  mosesTest.default(
    x = d$x,
    y = d$y,
    ...
  )
}


# -------------------------------------------------------------------------
# Default method
# -------------------------------------------------------------------------

#' @rdname mosesTest
#' @export
mosesTest.default <- function(
    x,
    y,
    extreme = NULL,
    ties.method = c("first", "average", "random", "max", "min"),
    ...
) {
  
  DNAME <- paste(
    deparse(substitute(x)),
    "and",
    deparse(substitute(y))
  )
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 4L)
    stop("'x' must contain at least 4 finite observations")
  
  if (length(y) < 1L)
    stop("'y' must contain at least 1 finite observation")
  
  nk <- length(x)
  ne <- length(y)
  N  <- nk + ne
  
  ## --- trimming parameter -----------------------------------------------
  
  if (is.null(extreme)) {
    
    h <- max(floor(0.05 * nk), 1L)
    
  } else {
    
    if (length(extreme) != 1L ||
        !is.numeric(extreme) ||
        is.na(extreme) ||
        extreme < 0 ||
        extreme != floor(extreme)) {
      
      stop("'extreme' must be a single non-negative integer")
    }
    
    h <- as.integer(extreme)
  }
  
  h_max <- floor((nk - 2L) / 2L)
  
  if (h > h_max) {
    
    warning(sprintf(
      "'extreme' reduced from %d to %d (maximum for n_control = %d)",
      h,
      h_max,
      nk
    ))
    
    h <- h_max
  }
  
  ## --- ties handling ----------------------------------------------------
  
  ties.method <- match.arg(ties.method)
  
  if (ties.method != "first" &&
      anyDuplicated(c(x, y))) {
    
    warning(
      "non-integer tie handling produces approximate p-values only"
    )
  }
  
  ## --- pooled ranks -----------------------------------------------------
  
  R_all <- rank(
    -c(x, y),
    ties.method = ties.method
  )
  
  ## --- numerically stable log-sum-exp -----------------------------------
  
  log_sum_exp <- function(lx) {
    
    m <- max(lx)
    
    m + log(sum(exp(lx - m)))
  }
  
  ## --- exact CDF of Moses span distribution -----------------------------
  
  moses_cdf <- function(s,
                        n_ctrl,
                        n_exp,
                        trim = 0L) {
    
    n_eff <- n_ctrl - 2L * trim
    
    i_max <- min(s - n_eff, n_exp)
    
    if (i_max < 0L)
      return(0)
    
    log_terms <- vapply(
      seq.int(0L, i_max),
      function(i) {
        
        lchoose(i + n_eff - 2L, i) +
          lchoose(
            n_exp + 2L * trim + 1L - i,
            n_exp - i
          )
      },
      numeric(1L)
    )
    
    exp(
      log_sum_exp(log_terms) -
        lchoose(N, n_ctrl)
    )
  }
  
  ## --- one-tailed p-value -----------------------------------------------
  
  moses_pval <- function(trim) {
    
    R_ctrl <- sort(R_all[seq_len(nk)])
    
    if (trim > 0L) {
      
      R_ctrl <- R_ctrl[
        seq(trim + 1L,
            length(R_ctrl) - trim)
      ]
    }
    
    S <- as.integer(
      ceiling(
        max(R_ctrl) -
          min(R_ctrl) +
          1L
      )
    )
    
    p <- 1 - moses_cdf(
      S - 1L,
      nk,
      ne,
      trim
    )
    
    p <- min(max(p, 0), 1)
    
    list(
      S = S,
      p = p
    )
  }
  
  ## --- results ----------------------------------------------------------
  
  res_raw <- moses_pval(trim = 0L)
  
  res_trimmed <- moses_pval(trim = h)
  
  structure(
    list(
      statistic = c(
        S_raw     = res_raw$S,
        S_trimmed = res_trimmed$S
      ),
      
      p.value = c(
        p_raw     = res_raw$p,
        p_trimmed = res_trimmed$p
      ),
      
      extreme = h,
      
      parameter = c(
        n_control    = nk,
        n_experiment = ne
      ),
      
      null.value = c(
        "equal extremity" = 0
      ),
      
      alternative = "greater",
      
      method = "Moses Test of Extreme Reactions",
      
      data.name = DNAME
    ),
    
    class = c(
      "mosesTestResult",
      "htest"
    )
  )
}


# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------

#' @export
print.mosesTestResult <- function(
    x,
    digits = getOption("digits"),
    ...
) {
  
  cat("\n\t", x$method, "\n\n")
  
  cat("data: ", x$data.name, "\n")
  
  cat(sprintf(
    "n_control = %d,  n_experiment = %d\n\n",
    x$parameter["n_control"],
    x$parameter["n_experiment"]
  ))
  
  cat(sprintf(
    "  %-32s  span = %4d,  p-value = %s\n",
    "Without trimming:",
    x$statistic["S_raw"],
    format.pval(
      x$p.value["p_raw"],
      digits = digits
    )
  ))
  
  cat(sprintf(
    "  %-32s  span = %4d,  p-value = %s\n",
    sprintf(
      "After trimming %d from each end:",
      x$extreme
    ),
    x$statistic["S_trimmed"],
    format.pval(
      x$p.value["p_trimmed"],
      digits = digits
    )
  ))
  
  cat(
    "\nalternative hypothesis:",
    "extreme values are more likely in the control group\n\n"
  )
  
  invisible(x)
}

