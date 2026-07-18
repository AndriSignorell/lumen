
#' Levene's Test for Homogeneity of Variance
#'
#' A test for homogeneity of variances across two or more groups, more
#' robust to departures from normality than Bartlett's test.
#'
#' Let \eqn{X_{ij}} be the jth observation of X for the ith group. Let
#' \eqn{Z_{ij} = |X_{ij} - X_i|}, where \eqn{X_i} is the center (by
#' default the median) of X in the ith group. Levene's test statistic is
#' \deqn{W_0 = \frac{\sum_i n_i (\bar{Z}_i - \bar{Z})^2 / (g - 1)}{\sum_i
#' \sum_j (Z_{ij} - \bar{Z}_i)^2 / \sum_i (n_i - 1)}}
#' where \eqn{n_i} is the number of observations in group i and g is the
#' number of groups. Using the median instead of the mean for \eqn{X_i}
#' yields the more robust variant known as the Brown-Forsythe test, which
#' is the default here.
#'
#' @name leveneTest
#' @aliases leveneTest leveneTest.formula leveneTest.default
#'
#' @param x a numeric vector of data values (default method), or a
#' \code{formula} (formula method).
#' @param g factor defining the groups; ignored if \code{x} is a list.
#' @param center the name of a function to compute the center of each
#' group; \code{mean} gives the original Levene's test, the default,
#' \code{median}, provides the more robust Brown-Forsythe test.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs}
#' gives the data values and \code{rhs} the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}. By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to
#' be used.
#' @param na.action a function which indicates what should happen when the
#' data contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param .centerName internal, not intended to be set by the user. Used
#' to pass the deparsed name of the \code{center} function through the
#' method dispatch chain (from \code{leveneTest.formula} to
#' \code{leveneTest.default}), since \code{substitute(center)} would
#' otherwise only resolve to the literal symbol \code{"center"} rather
#' than the original expression (e.g. \code{mean} or \code{median})
#' supplied by the caller.
#' @param \dots arguments to be passed down, e.g. \code{data} for the
#' formula method; can also be used to pass arguments to the function
#' given by \code{center} (e.g. \code{trim = 0.1} for a trimmed mean when
#' \code{center = mean}).
#'
#' @return An object of class \code{"htest"} representing the result of
#' the hypothesis test.
#'
#' @note
#' Based on \code{car::leveneTest()} by John Fox, with contributions by
#' Derek Ogle and Brian Ripley, adapted to conform to package standards.
#'
#' @references
#' Fox, J. and Weisberg, S. (2019) \emph{An R Companion to Applied
#' Regression}, 3rd ed., Thousand Oaks, CA: Sage.
#'
#' Levene, H. (1960) Robust tests for equality of variances. In: Olkin,
#' I. et al., eds.: \emph{Contributions to Probability and Statistics:
#' Essays in Honor of Harold Hotelling}, Stanford University Press,
#' pp. 278-292.
#'
#' @seealso \code{\link{fligner.test}} for a rank-based (nonparametric)
#' k-sample test for homogeneity of variances, \code{\link{bartlett.test}}
#' for a parametric alternative
#'
#' @examples
#' ## example from ansari.test:
#' ## Hollander & Wolfe (1973, p. 86f):
#' ## Serum iron determination using Hyland control sera
#' serum <- data.frame(
#'   grp = rep(c("ramsay", "jung.parekh"), each = 20),
#'   x   = c(111, 107, 100, 99, 102, 106, 109, 108, 104, 99,
#'           101, 96, 97, 102, 107, 113, 116, 113, 110, 98,
#'           107, 108, 106, 98, 105, 103, 110, 105, 104,
#'           100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99)
#' )
#'
#' leveneTest(x ~ grp, data = serum)
#' leveneTest(x ~ grp, data = serum, center = mean)
#' leveneTest(x ~ grp, data = serum, center = mean, trim = 0.1)
#'
#' leveneTest(c(rnorm(10), rnorm(10, 0, 2)),
#'            factor(rep(c("A", "B"), each = 10)))
#'
#' leveneTest(Ozone ~ factor(Month), data = airquality)
#'
#' leveneTest(count ~ spray, data = InsectSprays)
#' # Compare this to fligner.test() and bartlett.test()
#'
#' @family test.variance
#' @concept variance-test
#' @concept homogeneity
#'
#' @export
leveneTest <- function(x, ...)
  UseMethod("leveneTest")



#' @rdname leveneTest
#' @export
leveneTest.formula <- function(formula, data, subset,
                               na.action = na.pass,
                               center    = median, ...) {

  subset_expr <- if (!missing(subset)) substitute(subset) else NULL

  res <- resolveFormula(formula, data,
                        subset    = subset_expr,
                        na.action = na.action,
                        allowed   = "n-sample-independent")

  y <- leveneTest.default(x           = res$x,
                          g           = res$group,
                          center      = center,
                          .centerName = deparse1(substitute(center)),
                          ...)

  y$data.name <- res$data.name

  y
}



#' @rdname leveneTest
#' @export
leveneTest.default <- function(x, g, center = median, .centerName = NULL,
                               ...) {

  if (is.list(x)) {

    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    if (!missing(g))
      warning("'x' is a list, so ignoring argument 'g'")

    DNAME <- deparse1(substitute(x))

    x <- lapply(x, function(u) u[complete.cases(u)])
    if (!all(vapply(x, is.numeric, logical(1))))
      warning("some elements of 'x' are not numeric and will be ",
              "coerced to numeric")

    k <- length(x)
    l <- lengths(x)
    if (any(l == 0L))
      stop("all groups must contain data")

    g <- factor(rep.int(seq_len(k), l))
    x <- unlist(x)

  } else {

    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")

    DNAME <- paste(deparse1(substitute(x)), "and", deparse1(substitute(g)))

    OK <- complete.cases(x, g)
    x  <- x[OK]
    g  <- factor(g[OK])
    k  <- nlevels(g)
    if (k < 2L)
      stop("all observations are in the same group")
  }

  n <- length(x)
  if (n < 2L)
    stop("not enough observations")

  meds <- tapply(x, g, center, ...)
  resp <- abs(x - meds[g])

  ANOVA_TAB <- anova(lm(resp ~ g))
  rownames(ANOVA_TAB)[2] <- " "

  dots <- unlist(match.call(expand.dots = FALSE)$...)
  center_x <- if (!is.null(.centerName)) .centerName else
    deparse1(substitute(center))
  if (!is.null(dots))
    center_x <- paste0(
      center_x,
      gettextf("(%s)",
               paste(gettextf("%s=%s", names(dots), dots),
                     collapse = ", "))
    )

  STATISTIC <- ANOVA_TAB$`F value`[1]
  PARAMETER <- ANOVA_TAB$Df
  PVAL      <- ANOVA_TAB$`Pr(>F)`[1]

  names(STATISTIC) <- "F"
  names(PARAMETER) <- c("num df", "denom df")

  structure(
    list(
      statistic = STATISTIC,
      parameter = PARAMETER,
      p.value   = PVAL,
      method    = gettextf(
        "Levene's Test for Homogeneity of Variance (center = %s)",
        center_x),
      data.name = DNAME,
      anova_tab = ANOVA_TAB
    ),
    class = "htest"
  )
}
