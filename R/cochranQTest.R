
#' Cochran's Q Test for Comparing Matched Binary Responses Across Multiple Conditions
#'
#' A nonparametric test for dependent samples with dichotomous data,
#' assessing whether proportions differ across multiple conditions or
#' time points. It generalizes the McNemar test to more than two groups.
#'
#' Performs Cochran's Q test for related samples with a binary response.
#' The test is appropriate for unreplicated complete block designs (i.e.,
#' matched or paired data), where each block contains exactly one observation
#' for each group.
#'
#' Cochran's Q test is a nonparametric method for testing whether the
#' proportions of a binary outcome differ across multiple related groups,
#' while accounting for block effects. It can be regarded as an extension
#' of McNemar's test to more than two groups. The null hypothesis is that,
#' apart from block effects, the probability of a "success" (coded as 1)
#' is the same in all groups.
#'
#' If \code{y} is a matrix, groups and blocks are inferred from the columns
#' and rows, respectively. Missing values in \code{y} lead to the removal
#' of the corresponding blocks (for both methods). Missing values are not
#' allowed in \code{groups} or \code{blocks}. If all blocks give identical
#' responses, the statistic is 0 by convention and the p-value 1.
#'
#' With \code{method = "asymptotic"} the test statistic is referred to its
#' approximate chi-squared distribution. With \code{method = "approximate"}
#' a Monte Carlo permutation p-value is obtained via
#' \code{coin::symmetry_test()}, whose quadratic-form statistic coincides
#' with Cochran's Q for this design.
#'
#' Cochran's Q test is closely related to the Friedman test, but is
#' specifically designed for binary (0/1) responses.
#'
#' @name cochranQTest
#' @aliases cochranQTest cochranQTest.default cochranQTest.formula
#'
#' @param y either a numeric vector of data values, or a matrix with the
#' blocks in the rows and the groups in the columns.
#' @param method a character string specifying how the p-value is computed,
#' one of \code{"asymptotic"} (default, chi-squared approximation) or
#' \code{"approximate"} (Monte Carlo permutation via the \pkg{coin}
#' package).
#' @param nresample the number of Monte Carlo replicates used for
#' \code{method = "approximate"} (default is \code{1e4}).
#' @param groups a vector giving the group for the corresponding elements
#' of \code{y} if this is a vector; ignored if \code{y} is a matrix. If not
#' a factor object, it is coerced to one.
#' @param blocks a vector giving the block for the corresponding elements
#' of \code{y} if this is a vector; ignored if \code{y} is a matrix. If not
#' a factor object, it is coerced to one.
#' @param formula a formula of the form \code{y ~ groups | blocks}.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula. By
#' default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to
#' be used.
#' @param na.action a function which indicates what should happen when the
#' data contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' 
#' @return A list with class \code{"htest"} containing the following
#' components:
#'   \item{\code{statistic}}{the value of Cochran's chi-squared statistic.}
#'   \item{\code{parameter}}{the degrees of freedom of the approximate
#'     chi-squared distribution of the test statistic (asymptotic method
#'     only).}
#'   \item{\code{p.value}}{the p-value of the test.}
#'   \item{\code{method}}{a character string indicating the test performed and
#'     the method used to compute the p-value.}
#'   \item{\code{data.name}}{a character string giving the names of the data.}
#'
#' @references Cochran, W.G. (1950) The Comparison of Percentages in Matched
#' Samples. \emph{Biometrika}, 37 (3/4): 256-266.
#'
#' @examples
#' # example in:
#' # http://support.sas.com/documentation/cdl/en/statugfreq/63124/PDF/default/statugfreq.pdf
#' # pp. S. 1824
#'
#' # create the dataset
#' d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
#'             rep(1:8, c(6,2,2,6,16,4,4,6)), ]
#' row.names(d.frm) <- NULL
#'
#' # rearrange to long shape
#' d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[c(1:3)],
#'                   v.names="resp", direction="long")
#'
#' # after having done the hard work of data organisation,
#' # performing the test is a piece of cake....
#' cochranQTest(resp ~ time | id, data=d.long)
#'
#' # and let's perform a post hoc analysis using mcnemar's test
#' z <- split(d.long, f=d.long$time)
#' pairwise.table(function(i, j) {
#'     mcnemar.test(z[[i]]$resp, z[[j]]$resp, correct=FALSE)$p.value
#'   },
#'   level.names = names(z),
#'   p.adjust.method = "fdr"
#' )
#'
#' @family test.categorical
#' @concept k-sample
#' @concept nonparametric
#'
#' @export
cochranQTest <- function(y, ...) {
  UseMethod("cochranQTest")
}



# =========================
# DEFAULT METHOD
# =========================

#' @rdname cochranQTest
#' @export
cochranQTest.default <- function(y, groups, blocks,
                                 method    = c("asymptotic", "approximate"),
                                 nresample = 1e4,
                                 na.action = na.omit,
                                 ...) {

  method <- match.arg(method)
  DNAME  <- deparse1(substitute(y))

  if (is.matrix(y)) {

    # validate binary
    if (!all(y %in% c(0L, 1L, NA)))
      stop("'y' must contain only 0/1 values")

    y_mat <- y
    if (is.null(colnames(y_mat)))
      colnames(y_mat) <- seq_len(ncol(y_mat))

  } else {

    if (anyNA(groups) || anyNA(blocks))
      stop("NA's are not allowed in 'groups' or 'blocks'")

    if (any(diff(c(length(y), length(groups), length(blocks))) != 0L))
      stop("'y', 'groups' and 'blocks' must have the same length")

    DNAME <- paste0(DNAME, ", ", deparse1(substitute(groups)),
                    " and ", deparse1(substitute(blocks)))

    if (any(table(groups, blocks) != 1L))
      stop("not an unreplicated complete block design")

    groups <- factor(groups)
    blocks <- factor(blocks)
    o      <- order(groups, blocks)
    y      <- asBinary(y[o])
    groups <- groups[o]
    blocks <- blocks[o]

    y_mat <- matrix(y, ncol = nlevels(groups),
                    dimnames = list(NULL, levels(groups)))
  }

  # missing values in y lead to the removal of the corresponding blocks
  y_mat <- y_mat[complete.cases(y_mat), , drop = FALSE]

  if (method == "approximate")
    .cochranApproximate(y_mat, DNAME, nresample)
  else
    .cochranAsymptotic(y_mat, DNAME)
}



#' @rdname cochranQTest
#' @export
cochranQTest.formula <- function(formula,
                                 data,
                                 subset,
                                 na.action = na.pass,
                                 method    = c("asymptotic", "approximate"),
                                 nresample = 1e4,
                                 ...) {

  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")

  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "n-sample-dependent"
  )

  if (!missing(data))
    args$data <- data

  if (!missing(subset))
    args$subset <- substitute(subset)

  d <- do.call(resolveFormula, args)

  res <- cochranQTest.default(
    y         = d$response,
    groups    = d$treatment,   # resolveFormula: 'group' renamed to 'treatment' for n-sample-dependent
    blocks    = d$block,
    method    = method,
    nresample = nresample,
    ...
  )

  res$data.name <- d$data.name
  res
}




# == internal helper functions ============================================


.cochranAsymptotic <- function(x, DNAME) {

  n <- nrow(x)
  k <- ncol(x)

  if (k < 2) stop("need at least 2 groups")
  if (n < 2) stop("need at least 2 blocks")

  Sj  <- colSums(x)
  Ri  <- rowSums(x)
  Tot <- sum(x)

  denom <- k * Tot - sum(Ri^2)

  if (denom == 0) {
    # all blocks give identical responses -> Q = 0 by convention
    return(structure(list(
      statistic = c("Cochran's Q" = 0),
      parameter = c(df = k - 1L),
      p.value   = 1,
      method    = "Cochran's Q test (asymptotic)",
      data.name = DNAME
    ), class = "htest"))
  }

  Q <- (k - 1) * (k * sum(Sj^2) - Tot^2) / denom

  structure(list(
    statistic = c("Cochran's Q" = Q),
    parameter = c(df = k - 1L),
    p.value   = pchisq(Q, df = k - 1, lower.tail = FALSE),
    method    = "Cochran's Q test (asymptotic)",
    data.name = DNAME
  ), class = "htest")
}


.cochranApproximate <- function(x, DNAME, nresample) {

  if (!requireNamespace("coin", quietly = TRUE))
    stop("package 'coin' required for method = 'approximate'")

  k <- ncol(x)
  n <- nrow(x)

  df <- data.frame(
    y      = factor(as.vector(x)),
    groups = factor(rep(colnames(x), each = n)),
    blocks = factor(rep(seq_len(n), times = k))
  )

  res <- coin::symmetry_test(
    y ~ groups | blocks,
    data         = df,
    teststat     = "quad",
    distribution = coin::approximate(nresample = nresample)
  )

  structure(list(
    statistic = c("Cochran's Q" = as.numeric(coin::statistic(res))),
    p.value   = as.numeric(coin::pvalue(res)),
    method    = paste0("Cochran's Q test (approximate, B = ",
                       nresample, ")"),
    data.name = DNAME
  ), class = "htest")
}
