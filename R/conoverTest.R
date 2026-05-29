
#' Conover's Test of Multiple Comparisons
#'
#' A nonparametric post hoc test for multiple pairwise comparisons
#' following a significant Kruskal-Wallis test, based on rank data.
#'
#' \code{conoverTest} performs the post hoc pairwise multiple-comparison
#' procedure appropriate after rejection of the Kruskal-Wallis null
#' hypothesis. The test is based on the Conover-Iman rank-sum statistic and
#' is generally more powerful than Dunn's procedure.
#'
#' Interpretation in terms of stochastic dominance requires the additional
#' assumption that the cumulative distribution functions of the compared
#' groups do not cross.
#'
#' If \code{x} is a list, its elements are taken as the samples to be
#' compared and must be numeric vectors. In this case \code{g} is ignored.
#'
#' Otherwise, \code{x} must be a numeric vector and \code{g} a grouping
#' variable of the same length.
#'
#' Pairwise comparisons are formed as \eqn{\bar{R}_i - \bar{R}_j} for all
#' \eqn{i < j} (lower triangle of the mean-rank difference matrix, with
#' groups ordered by their factor levels). For one-sided alternatives,
#' \code{"greater"} tests whether group \eqn{i} tends to have larger
#' observations than group \eqn{j}, and \code{"less"} tests the reverse.
#'
#' @name conoverTest
#'
#' @aliases conoverTest conoverTest.default conoverTest.formula
#'
#' @param x a numeric vector of observations or a list of numeric vectors.
#' @param g a grouping variable corresponding to \code{x}. Ignored when
#'   \code{x} is a list.
#' @param method method used to adjust p-values for multiple comparisons.
#'   Passed directly to \code{\link{p.adjust}}.
#' @param alternative a character string specifying the alternative
#'   hypothesis. Must be one of \code{"two.sided"} (default),
#'   \code{"greater"} or \code{"less"}. See Details for the direction
#'   convention.
#' @param output output format:
#'   \itemize{
#'     \item \code{"list"} pairwise comparison table
#'     \item \code{"matrix"} lower-triangular matrix of adjusted p-values
#'   }
#' @param formula a formula of the form \code{response ~ group}.
#'   Requires at least three groups (\eqn{k \ge 3}); Conover's test is
#'   intended as a post hoc procedure following a significant
#'   Kruskal-Wallis test.
#' @param data an optional data frame containing the variables in
#'   \code{formula}.
#' @param subset an optional expression specifying a subset of observations.
#' @param na.action a function indicating how missing values should be handled.
#' @param \dots further arguments passed to methods.
#'
#' @return
#' An object of class \code{c("rankTest", "htest")} containing:
#' \describe{
#'   \item{res}{
#'     Pairwise comparison results. Depending on \code{output},
#'     either a table of mean-rank differences and adjusted p-values
#'     or a lower-triangular p-value matrix.
#'   }
#'   \item{pmat}{
#'     Symmetric matrix of adjusted p-values.
#'   }
#' }
#'
#' @seealso
#' \code{\link{kruskal.test}},
#' \code{\link{wilcox.test}},
#' \code{\link{p.adjust}}
#'
#' @references
#' Conover, W. J. and Iman, R. L. (1979).
#' On multiple-comparisons procedures.
#' \emph{Technical Report LA-7677-MS},
#' Los Alamos Scientific Laboratory.
#'
#' Conover, W. J. (1999).
#' \emph{Practical Nonparametric Statistics}.
#' 3rd edition. Wiley, Hoboken, NJ.
#'
#' @examples
#'
#' ## Hollander & Wolfe (1973), p. 116
#'
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
#' conoverTest(Ozone ~ Month, data = airquality)
#'
#' @family test.posthoc
#' @concept multiple-comparisons
#' @concept nonparametric
#' @concept hypothesis-testing
#'
#'
#' @export
conoverTest <- function(x, ...)
  UseMethod("conoverTest")


# ======================================================================
# Formula method
# ======================================================================

#' @rdname conoverTest
#' @export
conoverTest.formula <- function(
    formula,
    data,
    subset,
    na.action,
    ...
) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")
  
  ## IMPORTANT!!
  ## --- capture subset / na.action HERE ---
  subset_expr <- if (!missing(subset)) substitute(subset) else NULL
  na_expr     <- if (!missing(na.action)) substitute(na.action) else NULL
  
  pf <- resolveFormula(
    formula   = formula,
    data      = data,
    subset    = subset_expr,
    na.action = na_expr,
    allowed   = "n-sample-independent"
  )
  
  y <- conoverTest(
    x = pf$x,
    g = pf$group,
    ...
  )
  
  y$data.name <- pf$data.name
  
  y
  
}


#' @rdname conoverTest
#' @export
conoverTest.default <- function(
    x,
    g,
    method = c(
      "holm", "hochberg", "hommel",
      "bonferroni", "BH", "BY", "fdr", "none"
    ),
    alternative = c("two.sided", "less", "greater"),
    output = c("list", "matrix"),
    ...
) {
  
  alternative <- match.arg(alternative)
  output      <- match.arg(output)
  method      <- match.arg(method)
  
  dat <- resolveGroups(x, g)
  
  x   <- dat$x
  g   <- dat$g
  
  N   <- dat$n
  k   <- dat$k
  # coerce to plain numeric vector to avoid table dimname artefacts in outer()
  n   <- as.numeric(dat$group.sizes)
  names(n) <- dat$group.names
  nms <- dat$group.names
  
  rnk  <- rank(x)
  mrnk <- tapply(rnk, g, mean)
  
  tau <- table(rnk[allDuplicated(rnk)])
  
  tiesadj <- 1 - sum(tau^3 - tau) / (N^3 - N)
  
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
  
  # Pairs are formed as R_i - R_j for i < j (lower triangle).
  # "greater": i tends to exceed j  -> large positive t -> upper tail
  # "less":    j tends to exceed i  -> large negative t -> lower tail of |t|
  if (alternative == "greater") {
    pvals <- pt(abs(tval), df = N - k, lower.tail = FALSE)
    
  } else if (alternative == "less") {
    pvals <- pt(abs(tval), df = N - k)
    
  } else {
    pvals <- 2 * pt(abs(tval), df = N - k, lower.tail = FALSE)
  }
  
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
      paste(rownames(pmatxt)[x < 0.05], collapse = ",")
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
  
  class(out) <- c("rankTest", "htest")
  
  attr(out, "main")      <- gettextf(
    "Conover's test of multiple comparisons : %s",
    method
  )
  attr(out, "method")    <- method
  attr(out, "output")    <- output
  attr(out, "data.name") <- dat$data.name
  
  out
  
}

