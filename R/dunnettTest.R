#' Dunnett's Test for Comparing Treatments With a Control
#'
#' Performs Dunnett's parametric post hoc test for comparing multiple
#' treatment groups with a control group while controlling the familywise
#' error rate.
#'
#' If \code{x} is a list, its elements are taken as the samples to be
#' compared and must therefore be numeric vectors. In this case, \code{g}
#' is ignored. The test can be performed with \code{dunnettTest(x)} or,
#' if the samples are stored in separate objects, with
#' \code{dunnettTest(list(x, ...))}.
#'
#' Otherwise, \code{x} must be a numeric vector and \code{g} a vector or
#' factor of the same length that identifies the group of each observation.
#'
#' The adjusted p-values and simultaneous confidence intervals are
#' two-sided. The confidence limits use the equicoordinate quantile of
#' \eqn{\max_j |T_j|} (\code{tail = "both.tails"} in
#' \code{mvtnorm::qmvt()}), in accordance with the adjusted p-values
#' \eqn{1 - P(\max_j |T_j| \le |t_i|)}. Multivariate t probabilities are
#' evaluated by randomized quasi-Monte Carlo integration using
#' \pkg{mvtnorm}. A fixed seed ensures that repeated calls produce identical
#' results up to the numerical integration accuracy.
#'
#' @name dunnettTest
#' @aliases dunnettTest dunnettTest.default dunnettTest.formula
#' 
#' @param x a numeric vector of data values or a list of numeric vectors.
#' @param g a vector or factor identifying the group of each element of
#'   \code{x}; ignored if \code{x} is a list.
#' @param control a character vector identifying one or more control levels.
#'   Each specified control is compared separately with all remaining groups.
#'   Defaults to the first group.
#' @param conf.level the confidence level for the simultaneous confidence
#'   intervals. Defaults to \code{0.95}.
#' @param formula a formula of the form \code{lhs ~ rhs}, where \code{lhs}
#'   contains the data values and \code{rhs} defines the corresponding groups.
#' @param data an optional matrix or data frame (or a similar object; see
#'   \code{\link{model.frame}}) containing the variables in \code{formula}.
#'   By default, the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to
#'   be used.
#' @param na.action a function indicating how missing values should be
#'   handled. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments passed to or from methods.
#' 
#' @return An object of class \code{"PostHocTest"}: a list containing one
#'   matrix for each control level. Each matrix has columns \code{diff} for
#'   the observed mean difference (treatment minus control), \code{lci}
#'   and \code{uci} for the simultaneous confidence limits, and
#'   \code{pval} for the multiplicity-adjusted p-value.
#'
#' Print and plot methods are available for class \code{"PostHocTest"}.
#' The plot method supplies its own axis labels and title and therefore does
#' not accept \code{xlab}, \code{ylab}, or \code{main}.
#'
#' @seealso
#' \code{\link{print.PostHocTest}},
#' \code{\link{plot.PostHocTest}},
#' \code{\link[mvtnorm]{pmvt}},
#' \code{\link[mvtnorm]{qmvt}}
#'
#' @references
#' Dunnett, C. W. (1955). A multiple comparison procedure for comparing
#' several treatments with a control. \emph{Journal of the American
#' Statistical Association}, \bold{50}, 1096--1121.
#'
#' @examples
#' ## Hollander and Wolfe (1973, p. 116)
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ##  subjects, subjects with obstructive airway disease, and subjects
#' ##  with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#'
#' dunnettTest(list(x, y, z))
#'
#' ## Equivalent vector-and-group interface
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)),
#'             labels = c("Normal subjects",
#'                        "Subjects with obstructive airway disease",
#'                        "Subjects with asbestosis"))
#'
#' dunnettTest(x, g)
#'
#' ## Formula interface
#' boxplot(Ozone ~ factor(Month), data = airquality)
#' dunnettTest(Ozone ~ factor(Month), data = airquality)
#'
#' ## Single control level with a 90% simultaneous confidence level
#' dunnettTest(Ozone ~ factor(Month), data = airquality,
#'             control = "8", conf.level = 0.9)
#'
#' ## Multiple control levels
#' dunnettTest(Ozone ~ factor(Month), data = airquality,
#'             control = c("5", "8"))
#'
#' @family test.posthoc
#' @concept k-sample
#' @concept parametric
#'
#' @export
dunnettTest <- function(x, ...)
  UseMethod("dunnettTest")


# ======================================================================
# Formula method
# ======================================================================

#' @rdname dunnettTest
#' @export
dunnettTest.formula <- function(formula,
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

  rval <- dunnettTest(x = pf$x, g = pf$group, ...)

  attr(rval, "data.name") <- pf$data.name

  rval
}


# ======================================================================
# Default method
# ======================================================================

#' @rdname dunnettTest
#' @export
dunnettTest.default <- function(x,
                                g,
                                control = NULL,
                                conf.level = 0.95,
                                ...) {

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number between 0 and 1")

  DG <- resolveGroups(x, g)

  x <- DG$x
  g <- DG$groups

  N <- DG$n
  k <- DG$k

  gn <- DG$group.names

  ni <- as.numeric(DG$group.sizes)
  names(ni) <- gn

  DNAME <- DG$data.name

  if (is.null(control))
    control <- gn[1]

  ctrls <- control

  if (!all(ctrls %in% gn)) {
    stop(gettextf("control level '%s' not found in grouping variable",
                  ctrls[!ctrls %in% gn][1]))
  }

  means <- tapply(x, g, mean)

  rss <- tapply(x, g, function(z) sum((z - mean(z))^2))

  s <- sqrt(sum(rss) / (N - k))

  out <- vector("list", length(ctrls))

  for (ii in seq_along(ctrls)) {

    control <- ctrls[ii]

    meandiffs <- means[names(means) != control] - means[control]

    fittedn  <- ni[names(ni) != control]
    controln <- ni[control]

    se <- s * sqrt(1 / fittedn + 1 / controln)
    Dj <- meandiffs / se

    Rij <- sqrt(fittedn / (fittedn + controln))

    R <- outer(Rij, Rij, "*")
    diag(R) <- 1

    # two-sided equicoordinate quantile of max |T_j|, consistent with
    # the two-sided adjusted p-values below
    qvt <- mvtnorm::qmvt(
      p     = conf.level,
      df    = N - k,
      sigma = R,
      tail  = "both.tails",
      seed  = 5L
    )$quantile

    lower <- meandiffs - se * qvt
    upper <- meandiffs + se * qvt

    pval <- numeric(length(meandiffs))

    for (i in seq_along(meandiffs)) {
      pval[i] <- 1 - mvtnorm::pmvt(
        lower = -abs(Dj[i]),
        upper =  abs(Dj[i]),
        corr  = R,
        delta = rep(0, length(meandiffs)),
        df    = N - k,
        seed  = 5L
      )[1]
    }

    out[[ii]] <- cbind(
      diff   = meandiffs,
      lci = lower,
      uci = upper,
      pval   = pval
    )

    rownames(out[[ii]]) <- paste(names(meandiffs), control, sep = "-")
  }

  names(out) <- ctrls

  class(out) <- "PostHocTest"

  attr(out, "conf.level") <- conf.level
  attr(out, "ordered")    <- FALSE
  attr(out, "method")     <- "Dunnett"
  attr(out, "data.name")  <- DNAME

  attr(out, "method.str") <- gettextf(
    "\n  Dunnett's test for comparing several treatments with a control : %s\n",
    attr(out, "method")
  )

  out
}
