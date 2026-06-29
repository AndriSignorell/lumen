
#' Nemenyi Test of Multiple Comparisons
#'
#' A nonparametric post hoc test for multiple pairwise comparisons following
#' a significant Kruskal-Wallis test, based on differences in mean ranks.
#'
#' Performs Nemenyi's multiple comparison procedure for independent samples.
#' The test compares all pairs of groups using rank differences and controls
#' the family-wise error rate through the studentized range distribution
#' (Tukey-type approximation) or an asymptotic chi-squared approximation.
#'
#' Nemenyi's test is commonly used as a post hoc procedure after a significant
#' \code{\link{kruskal.test}} when all pairwise comparisons between groups are
#' of interest. Unlike \code{\link{dunnTest}} and
#' \code{\link{conoverTest}}, no additional p-value adjustment is applied,
#' since multiplicity control is built into the test statistic.
#'
#' If \code{x} is a list, its elements are taken as the samples to be compared,
#' and hence have to be numeric data vectors. In this case, \code{g} is
#' ignored and one can simply use \code{nemenyiTest(x)}.
#'
#' Otherwise, \code{x} must be a numeric vector and \code{g} a grouping factor
#' (or vector coercible to a factor) of the same length.
#'
#' @name nemenyiTest
#' @aliases nemenyiTest nemenyiTest.default nemenyiTest.formula
#'
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a grouping factor corresponding to \code{x}. Ignored if
#'   \code{x} is a list.
#' @param dist character string specifying the reference distribution used for
#'   the test statistic. One of \code{"tukey"} (default) or \code{"chisq"}.
#' @param output character string specifying the output format. One of
#'   \code{"list"} (default) or \code{"matrix"}.
#' @param formula a formula of the form \code{response ~ group}.
#' @param data an optional data frame containing the variables in
#'   \code{formula}.
#' @param subset an optional expression specifying a subset of observations to
#'   be used.
#' @param na.action a function specifying how missing values should be handled.
#' @param \dots further arguments passed to methods.
#'
#' @return An object of class \code{c("rankTest", "htest")} containing:
#' \describe{
#'   \item{res}{pairwise comparison results, either as a list or matrix}
#'   \item{pmat}{symmetric matrix of adjusted p-values}
#' }
#'
#' Additional information is stored in attributes:
#' \code{method}, \code{output}, \code{main}, and \code{data.name}.
#'
#' @seealso
#' \code{\link{kruskal.test}}
#'
#' @references
#' Nemenyi, P. B. (1963).
#' \emph{Distribution-Free Multiple Comparisons}.
#' PhD thesis, Princeton University.
#'
#' Hollander, M., Wolfe, D. A. and Chicken, E. (2014).
#' \emph{Nonparametric Statistical Methods}.
#' 3rd ed. Wiley.
#'
#' @examples
#' ## Hollander & Wolfe example
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
#' y <- c(3.8, 2.7, 4.0, 2.4)
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
#'
#' nemenyiTest(list(x, y, z))
#'
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)))
#'
#' nemenyiTest(x, g)
#'
#' ## Formula interface
#' nemenyiTest(Ozone ~ Month, data = airquality)
#'
#' @rdname nemenyiTest


#' @family test.posthoc  
#' @concept post-hoc  
#' @concept nonparametric
#'
#'
#' @export
nemenyiTest <- function(x, ...)
  UseMethod("nemenyiTest")


# ======================================================================
# Formula method
# ======================================================================

#' @rdname nemenyiTest
#' @export
nemenyiTest.formula <- function(formula, data, subset, na.action, ...) {
  
  if (missing(formula) || (length(formula) != 3L) ||
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
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
  
  y <- nemenyiTest(
    x = pf$x,
    g = pf$group,
    ...
  )
  
  y$data.name <- pf$data.name
  
  y
  
}


#' @rdname nemenyiTest
#' @export
nemenyiTest.default <- function(
    x,
    g,
    dist = c("tukey", "chisq"),
    output = c("list", "matrix"),
    ...
) {
  
  output <- match.arg(output)
  dist   <- match.arg(dist)
  
  gd <- resolveGroups(x, g)
  
  x <- gd$x
  g <- gd$groups
  
  N   <- gd$n
  nms <- gd$group.names
  # coerce to plain numeric vector to avoid table dimname artefacts in outer()
  n   <- as.numeric(gd$group.sizes)
  names(n) <- nms
  
  rnk  <- rank(x)
  mrnk <- tapply(rnk, g, mean)
  
  tau <- table(rnk[allDuplicated(rnk)])
  
  tiesadj <- min(1, 1 - sum(tau^3 - tau) / (N^3 - N))
  
  mrnkdiff <- outer(mrnk, mrnk, "-")
  
  if (dist == "chisq") {
    
    chi <- mrnkdiff^2 /
      ((N * (N + 1) / 12) * outer(1 / n, 1 / n, "+"))
    
    pvals <- pchisq(tiesadj * chi, df = length(n) - 1,
                    lower.tail = FALSE)
    
  } else {
    
    z <- abs(mrnkdiff) / sqrt(
      (N * (N + 1) / 12) * outer(1 / n, 1 / n, "+"))
    
    pvals <- ptukey(z * sqrt(2), nmeans = length(n),
                    df = Inf, lower.tail = FALSE)
  }
  
  keep      <- lower.tri(pvals)
  pvals_vec <- pvals[keep]
  
  # --- p-value matrix ----------------------------------------------------
  
  pmat <- matrix(
    NA_real_,
    nrow = length(nms),
    ncol = length(nms),
    dimnames = list(nms, nms)
  )
  
  pmat[lower.tri(pmat)] <- pvals_vec
  
  pmatxt <- pmat
  pmatxt[upper.tri(pmatxt)] <- t(pmatxt)[upper.tri(pmatxt)]
  diag(pmatxt) <- 1
  
  attr(pmatxt, "lbl") <- apply(
    pmatxt,
    1,
    function(x)
      paste(rownames(pmatxt)[x < 0.05], collapse = ",")
  )
  
  # --- output ------------------------------------------------------------
  
  out <- list()
  
  if (output == "list") {
    
    dnames <- list(
      NULL,
      c("mean rank diff", "pval")
    )
    
    if (!is.null(nms)) {
      dnames[[1L]] <- outer(nms, nms, paste, sep = "-")[keep]
    }
    
    out$res <- array(
      c(mrnkdiff[keep], pvals_vec),
      c(length(mrnkdiff[keep]), 2L),
      dnames
    )
    
  } else {
    out$res <- pmat[-1, -ncol(pmat), drop = FALSE]
  }
  
  out$pmat <- pmatxt
  
  class(out) <- c("rankTest", "htest")
  
  attr(out, "main")      <- gettextf(
    "Nemenyi's test of multiple comparisons for independent samples (%s)",
    dist
  )
  attr(out, "method")    <- "none"
  attr(out, "output")    <- output
  attr(out, "data.name") <- gd$data.name
  
  out
  
}

