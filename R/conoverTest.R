
#' Conover's Test for Pairwise Rank Comparisons After a Kruskal–Wallis Test
#'
#' A nonparametric post hoc test for multiple pairwise comparisons
#' following a significant Kruskal-Wallis test, based on rank data.
#'
#' `conoverTest` performs the post hoc pairwise multiple-comparison
#' procedure appropriate after rejection of the Kruskal-Wallis null
#' hypothesis. The test is based on the Conover-Iman rank-sum statistic and
#' is generally more powerful than Dunn's procedure. It is intended as a
#' post hoc procedure following a significant Kruskal-Wallis test, i.e.
#' typically for three or more groups.
#'
#' Interpretation in terms of stochastic dominance requires the additional
#' assumption that the cumulative distribution functions of the compared
#' groups do not cross.
#'
#' If `x` is a list, its elements are taken as the samples to be
#' compared and must be numeric vectors. In this case `g` is ignored.
#' Otherwise, `x` must be a numeric vector and `g` a grouping
#' variable of the same length.
#'
#' Each pairwise comparison is labeled `"B-A"`, where `A` precedes
#' `B` in the ordering of the group levels, and reports the mean rank
#' difference \eqn{\bar{R}_B - \bar{R}_A}. For one-sided alternatives,
#' `"greater"` tests whether `B` tends to have larger observations
#' than `A` (upper tail), and `"less"` tests the reverse (lower
#' tail).
#'
#' @name conoverTest
#' @aliases conoverTest conoverTest.default conoverTest.formula
#'
#' @param x a numeric vector of observations or a list of numeric vectors.
#' @param g a grouping variable corresponding to `x`; ignored when
#' `x` is a list.
#' @param method the method used to adjust the p-values for multiple
#' comparisons, one of `p.adjust.methods` (default is `"holm"`).
#' Passed directly to [p.adjust()].
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of `"two.sided"` (default), `"less"`
#' or `"greater"`. See the Details for the direction convention.
#' @param output the output format:
#'   \itemize{
#'     \item `"list"` pairwise comparison table.
#'     \item `"matrix"` lower-triangular matrix of adjusted p-values.
#'   }
#' @param alpha the significance level used to compile the groups flagged
#' as significantly different in the label attribute of the p-value matrix
#' (default is `0.05`).
#' @param formula a formula of the form `response ~ group`.
#' @param data an optional data frame containing the variables in
#' `formula`.
#' @param subset an optional expression specifying a subset of observations.
#' @param na.action a function indicating how missing values should be
#' handled.
#' @param \dots further arguments passed to methods.
#'
#' @return
#' An object of class `"rankTest"` containing:
#'   \item{`res`}{
#'     pairwise comparison results. Depending on `output`,
#'     either a table of mean-rank differences and adjusted p-values
#'     or a lower-triangular p-value matrix.
#'   }
#'   \item{`pmat`}{
#'     symmetric matrix of adjusted p-values.
#'   }
#'
#' @seealso [kruskal.test()], [wilcox.test()], [p.adjust()]
#'
#' @references
#' Conover, W. J. and Iman, R. L. (1979) On multiple-comparisons procedures.
#' *Technical Report LA-7677-MS*, Los Alamos Scientific Laboratory.
#'
#' Conover, W. J. (1999) *Practical Nonparametric Statistics*, 3rd ed.,
#' Hoboken, NJ: Wiley.
#'
#' @examples
#' ## Hollander & Wolfe (1973), p. 116
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
#' y <- c(3.8, 2.7, 4.0, 2.4)
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
#'
#' conoverTest(list(x, y, z))
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
#' conoverTest(x, g)
#'
#' ## Formula interface
#' conoverTest(Ozone ~ factor(Month), data = airquality)
#'
#' @family test.posthoc
#' @concept k-sample
#' @concept nonparametric
#'
#' @export
conoverTest <- function(x, ...)
  UseMethod("conoverTest")


# ======================================================================
# Formula method
# ======================================================================

#' @rdname conoverTest
#' @export
conoverTest.formula <- function(formula,
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

  y <- conoverTest(x = pf$x, g = pf$group, ...)

  attr(y, "data.name") <- pf$data.name

  y
}


# ======================================================================
# Default method
# ======================================================================

#' @rdname conoverTest
#' @export
conoverTest.default <- function(x,
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

  tiesadj <- 1 - sum(tau^3 - tau) / (N^3 - N)

  # entry [B, A] is mean rank of B minus mean rank of A
  mrnkdiff <- outer(mrnk, mrnk, "-")

  # Kruskal-Wallis H statistic
  H <- (12 / (N * (N + 1))) *
    sum(tapply(rnk, g, sum)^2 / n) - 3 * (N + 1)

  if (tiesadj == 1) {
    s2 <- N * (N + 1) / 12
  } else {
    s2 <- (1 / (N - 1)) * (sum(rnk^2) - N * ((N + 1)^2 / 4))
  }

  tval <- mrnkdiff / sqrt(
    s2 * ((N - 1 - H / tiesadj) / (N - k)) * outer(1 / n, 1 / n, "+")
  )

  # Comparisons are taken from the lower triangle, labeled "B-A" with
  # t = (Rbar_B - Rbar_A) / se:
  # "greater": B tends to exceed A -> upper tail of the signed t
  # "less":    A tends to exceed B -> lower tail of the signed t
  pvals <- switch(
    alternative,
    "two.sided" = 2 * pt(abs(tval), df = N - k, lower.tail = FALSE),
    "greater"   = pt(tval, df = N - k, lower.tail = FALSE),
    "less"      = pt(tval, df = N - k)
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

    dnames <- list(
      NULL,
      c("mean rank diff", "pval")
    )

    if (!is.null(nms)) {
      dnames[[1L]] <- outer(
        nms,
        nms,
        paste,
        sep = "-"
      )[keep]
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
    "Conover's test of multiple comparisons : %s",
    method
  )
  attr(out, "method")    <- method
  attr(out, "output")    <- output
  attr(out, "data.name") <- dat$data.name

  out
}
