
#' Dunn's Test for Pairwise Rank Comparisons After a Kruskal–Wallis Test
#'
#' A nonparametric post hoc test for multiple pairwise comparisons
#' following a significant Kruskal-Wallis test, based on rank sums with
#' adjustment for multiple testing.
#'
#' \code{dunnTest} performs the post hoc pairwise multiple-comparison
#' procedure appropriate after rejection of the Kruskal-Wallis null
#' hypothesis. In contrast to performing separate Wilcoxon rank-sum tests,
#' Dunn's procedure preserves the pooled ranking and variance estimate
#' underlying the Kruskal-Wallis test. It is intended as a post hoc
#' procedure following a significant Kruskal-Wallis test, i.e. typically
#' for three or more groups.
#'
#' If \code{x} is a list, its elements are taken as the samples to be
#' compared and must be numeric vectors. In this case \code{g} is ignored.
#' Otherwise, \code{x} must be a numeric vector and \code{g} a grouping
#' variable of the same length.
#'
#' Each pairwise comparison is labeled \code{"B-A"}, where \code{A} precedes
#' \code{B} in the ordering of the group levels, and reports the mean rank
#' difference \eqn{\bar{R}_B - \bar{R}_A}. For one-sided alternatives,
#' \code{"greater"} tests whether \code{B} tends to have larger observations
#' than \code{A} (upper tail), and \code{"less"} tests the reverse (lower
#' tail).
#'
#' @name dunnTest
#' @aliases dunnTest dunnTest.default dunnTest.formula
#'
#' @param x a numeric vector of observations or a list of numeric vectors.
#' @param g a grouping variable corresponding to \code{x}; ignored when
#' \code{x} is a list.
#' @param method the method used to adjust the p-values for multiple
#' comparisons, one of \code{p.adjust.methods} (default is \code{"holm"}).
#' Passed directly to \code{\link{p.adjust}}.
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of \code{"two.sided"} (default), \code{"less"}
#' or \code{"greater"}. See the Details for the direction convention.
#' @param output the output format:
#'   \itemize{
#'     \item \code{"list"} pairwise comparison table.
#'     \item \code{"matrix"} lower-triangular matrix of adjusted p-values.
#'   }
#' @param alpha the significance level used to compile the groups flagged
#' as significantly different in the label attribute of the p-value matrix
#' (default is \code{0.05}).
#' @param formula a formula of the form \code{response ~ group}.
#' @param data an optional data frame containing the variables in
#' \code{formula}.
#' @param subset an optional expression specifying a subset of observations.
#' @param na.action a function indicating how missing values should be
#' handled.
#' @param \dots further arguments passed to methods.
#'
#' @return
#' An object of class \code{"rankTest"} containing:
#' \describe{
#'   \item{\code{res}}{
#'     pairwise comparison results. Depending on \code{output},
#'     either a table of mean-rank differences and adjusted p-values
#'     or a lower-triangular p-value matrix.
#'   }
#'   \item{\code{pmat}}{
#'     symmetric matrix of adjusted p-values.
#'   }
#' }
#'
#' @seealso \code{\link{kruskal.test}}, \code{\link{wilcox.test}},
#' \code{\link{p.adjust}}
#'
#' @references
#' Dunn, O. J. (1961) Multiple comparisons among means. \emph{Journal of
#' the American Statistical Association}, 56 (293), 52-64.
#'
#' Dunn, O. J. (1964) Multiple comparisons using rank sums.
#' \emph{Technometrics}, 6 (3), 241-252.
#'
#' @examples
#' ## Hollander & Wolfe (1973), p. 116
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
#' y <- c(3.8, 2.7, 4.0, 2.4)
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
#'
#' dunnTest(list(x, y, z))
#'
#' x <- c(x, y, z)
#' g <- factor(
#'   rep(1:3, c(5, 4, 5)),
#'   labels = c(
#'     "Normal subjects",
#'     "Subjects with obstructive airway disease",
#'     "Subjects with asbestosis"
#'   )
#' )
#'
#' kruskal.test(x, g)
#' dunnTest(x, g)
#'
#' ## Formula interface
#' dunnTest(Ozone ~ factor(Month), data = airquality)
#'
#' @family test.posthoc
#' @concept k-sample
#' @concept nonparametric
#'
#' @export
dunnTest <- function(x, ...)
  UseMethod("dunnTest")


# ======================================================================
# Formula method
# ======================================================================

#' @rdname dunnTest
#' @export
dunnTest.formula <- function(formula,
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

  y <- dunnTest(x = pf$x, g = pf$group, ...)

  attr(y, "data.name") <- pf$data.name

  y
}


# ======================================================================
# Default method
# ======================================================================

#' @rdname dunnTest
#' @export
dunnTest.default <- function(x,
                             g,
                             method = p.adjust.methods,
                             alternative = c("two.sided", "less",
                                             "greater"),
                             output = c("list", "matrix"),
                             alpha = 0.05,
                             ...) {

  alternative <- match.arg(alternative)
  output      <- match.arg(output)
  method      <- match.arg(method)

  dat <- resolveGroups(x, g)

  x <- dat$x
  g <- dat$groups

  N <- dat$n
  k <- dat$k
  # coerce to plain numeric vector to avoid table dimname artefacts in outer()
  n <- as.numeric(dat$group.sizes)
  names(n) <- dat$group.names
  nms <- dat$group.names

  rnk  <- rank(x)
  mrnk <- tapply(rnk, g, mean)

  tau <- table(rnk[allDuplicated(rnk)])

  tiesadj <- sum(tau^3 - tau) / (12 * (N - 1))

  # entry [B, A] is mean rank of B minus mean rank of A
  mrnkdiff <- outer(mrnk, mrnk, "-")

  z <- mrnkdiff / sqrt(
    ((N * (N + 1) / 12) - tiesadj) * outer(1 / n, 1 / n, "+")
  )

  # Comparisons are taken from the lower triangle, labeled "B-A" with
  # z = (Rbar_B - Rbar_A) / se:
  # "greater": B tends to exceed A -> upper tail of the signed z
  # "less":    A tends to exceed B -> lower tail of the signed z
  pvals <- switch(
    alternative,
    "two.sided" = 2 * pnorm(abs(z), lower.tail = FALSE),
    "greater"   = pnorm(z, lower.tail = FALSE),
    "less"      = pnorm(z)
  )

  keep <- lower.tri(pvals)

  pvals <- pvals[keep]
  pvals <- p.adjust(pvals, method = method)

  # --- p-value matrix -----------------------------------------------------

  pmat <- matrix(
    NA_real_,
    nrow = length(nms),
    ncol = length(nms),
    dimnames = list(nms, nms)
  )

  pmat[lower.tri(pmat, diag = FALSE)] <- pvals

  pmatxt <- pmat
  pmatxt[upper.tri(pmatxt)] <- t(pmatxt)[upper.tri(pmatxt)]
  diag(pmatxt) <- 1

  attr(pmatxt, "lbl") <- apply(
    pmatxt,
    1,
    function(x)
      paste(rownames(pmatxt)[x < alpha], collapse = ",")
  )

  # --- output -------------------------------------------------------------

  out <- list()

  if (output == "list") {

    dnames <- list(NULL, c("mean rank diff", "pval"))

    if (!is.null(nms)) {
      dnames[[1L]] <- outer(
        nms, nms,
        paste, sep = "-")[keep]
    }

    out$res <- array(
      c(mrnkdiff[keep], pvals),
      c(length(mrnkdiff[keep]), 2L),
      dnames
    )

  } else {
    out$res <- pmat[-1, -ncol(pmat), drop = FALSE]
  }

  out$pmat <- pmatxt

  class(out) <- "rankTest"

  attr(out, "main") <- gettextf(
    "Dunn's test of multiple comparisons using rank sums : %s",
    method
  )
  attr(out, "method")    <- method
  attr(out, "output")    <- output
  attr(out, "data.name") <- dat$data.name

  out
}
