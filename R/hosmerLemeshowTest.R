
#' Hosmer-Lemeshow Goodness of Fit Test
#'
#' Computes the Hosmer-Lemeshow C or H goodness-of-fit test for a logistic
#' regression model, assessing whether observed event rates match predicted
#' probabilities across grouped subsets of the data.
#'
#' The C statistic groups observations by quantiles of the fitted
#' probabilities (equal-frequency bins), while the H statistic uses
#' equal-width intervals on \eqn{[0, 1]} (e.g. 0--0.1, 0.1--0.2, ... for
#' \code{nGroups = 10}), as defined in Lemeshow and Hosmer (1982). Groups
#' that contain no observations are dropped with a warning, and the
#' degrees of freedom are based on the number of groups actually used.
#'
#' Both statistics are asymptotically compared with a chi-squared
#' distribution with \eqn{g - 2} degrees of freedom (\eqn{g} the number of
#' groups used) under the null hypothesis of adequate fit. This is the
#' convention for a model evaluated on its development data; for an
#' external validation sample, \eqn{g} degrees of freedom would be
#' appropriate. The approximation may be unreliable in small samples; a
#' warning is issued when expected cell counts fall below 5.
#'
#' @name hosmerLemeshowTest
#' @param x either a numeric vector of fitted probabilities, each in
#' \eqn{[0, 1]} and without missing values, or a fitted binomial
#' \code{\link{glm}} object, from which fitted probabilities and observed
#' outcomes are extracted.
#' @param obs a numeric vector of observed binary outcomes (0 or 1) of the
#' same length as \code{x}, without missing values; unused for the
#' \code{glm} method.
#' @param nGroups integer, the number of groups (default is \code{10}).
#' Must be at least 3.
#' @param type the type of statistic, one of \code{"C"} (default,
#' quantile-based groups) or \code{"H"} (equal-width groups on
#' \eqn{[0, 1]}).
#' @param \dots further arguments passed to methods.
#'
#' @return An object of class \code{c("HosmerLemeshowTest", "htest")},
#' which is a list with components:
#' \item{statistic}{the chi-squared test statistic.}
#' \item{parameter}{the degrees of freedom (number of groups used minus 2).}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a character string describing the test.}
#' \item{type}{the type of statistic computed (\code{"C"} or \code{"H"}).}
#' \item{nGroups}{the number of groups actually used (may be less than
#' requested if quantile breaks coincide or groups are empty).}
#' \item{observed}{matrix of observed counts for both outcome classes per
#' group.}
#' \item{expected}{matrix of expected counts for both outcome classes per
#' group.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' The \code{print} method accepts a \code{details} argument; if
#' \code{TRUE}, observed and expected counts for both outcome classes are
#' printed by group.
#'
#' @references
#' Lemeshow, S. and Hosmer, D. W. (1982) A review of goodness of fit
#' statistics for use in the development of logistic regression models.
#' \emph{American Journal of Epidemiology}, \bold{115}(1), 92--106.
#'
#' Hosmer, D. W., Lemeshow, S. and Sturdivant, R. X. (2013) \emph{Applied
#' Logistic Regression}, 3rd ed., New York: Wiley.
#'
#' @seealso \code{\link{glm}}
#'
#' @examples
#' set.seed(111)
#' x1  <- factor(sample(1:3, 50, replace = TRUE))
#' x2  <- rnorm(50)
#' obs <- sample(c(0, 1), 50, replace = TRUE)
#'
#' model <- glm(obs ~ x1 + x2, family = binomial)
#'
#' # glm method: probabilities and outcomes are extracted from the model
#' hosmerLemeshowTest(model)
#'
#' # equivalent call with explicit vectors
#' res <- hosmerLemeshowTest(fitted(model), obs, type = "C")
#' res
#'
#' print(res, details = TRUE)
#'
#' @family test.regression
#' @concept regression-diagnostics
#' @concept goodness-of-fit
#' @concept calibration
#'
#' @export
hosmerLemeshowTest <- function(x, ...)
  UseMethod("hosmerLemeshowTest")



#' @rdname hosmerLemeshowTest
#' @export
hosmerLemeshowTest.glm <- function(x, nGroups = 10, type = c("C", "H"),
                                   ...) {

  if (!(x$family$family %in% c("binomial", "quasibinomial")))
    stop("'x' must be a binomial glm")

  obs <- model.response(model.frame(x))

  if (is.matrix(obs))
    stop("matrix responses (cbind(successes, failures)) are not supported; ",
         "supply fitted probabilities and binary outcomes directly")

  if (is.factor(obs)) {
    if (nlevels(obs) != 2L)
      stop("the response must have exactly two levels")
    obs <- as.integer(obs) - 1L        # first level = failure, as in glm
  } else if (is.logical(obs)) {
    obs <- as.integer(obs)
  }

  res <- hosmerLemeshowTest.default(x = fitted(x), obs = obs,
                                    nGroups = nGroups, type = type)
  res$data.name <- deparse1(formula(x))

  res
}



#' @rdname hosmerLemeshowTest
#' @export
hosmerLemeshowTest.default <- function(x, obs, nGroups = 10,
                                       type = c("C", "H"), ...) {

  DNAME <- paste(deparse1(substitute(x)), "and", deparse1(substitute(obs)))

  ## input validation -----------------------------------------------------

  type <- match.arg(type)

  fit <- x

  if (!is.numeric(fit) || !is.numeric(obs))
    stop("'x' and 'obs' must be numeric vectors")
  if (length(fit) != length(obs))
    stop("'x' and 'obs' must have the same length")
  if (anyNA(fit) || anyNA(obs))
    stop("'x' and 'obs' must not contain missing values")
  if (any(fit < 0 | fit > 1))
    stop("'x' must contain probabilities in [0, 1]")
  if (!all(obs %in% c(0, 1)))
    stop("'obs' must be binary (0 or 1 only)")
  if (!is.numeric(nGroups) || length(nGroups) != 1L || is.na(nGroups) ||
      nGroups < 3 || nGroups != as.integer(nGroups))
    stop("'nGroups' must be a single integer >= 3")

  nGroups <- as.integer(nGroups)

  ## grouping -------------------------------------------------------------

  if (type == "C") {

    brks <- unique(quantile(fit, probs = seq(0, 1, by = 1 / nGroups)))

    if (length(brks) - 1L < 3L)
      stop("unable to construct at least 3 groups ",
           "from the fitted probabilities")
    if (length(brks) - 1L < nGroups)
      warning("found only ", length(brks) - 1L,
              " distinct groups for the Hosmer-Lemeshow C statistic")

    cutfit <- cut(fit, breaks = brks, include.lowest = TRUE)

  } else {

    # equal-width intervals on [0, 1], Lemeshow and Hosmer (1982)
    cutfit <- cut(fit, breaks = seq(0, 1, length.out = nGroups + 1L),
                  include.lowest = TRUE)
  }

  ## observed / expected --------------------------------------------------

  Obs <- xtabs(cbind("0s" = 1 - obs, "1s" = obs) ~ cutfit)
  Exp <- xtabs(cbind("0s" = 1 - fit, "1s" = fit) ~ cutfit)

  # drop groups without observations (possible for type = "H")
  empty <- rowSums(Obs) == 0
  if (any(empty)) {
    warning(sum(empty), " empty group(s) dropped; ",
            "degrees of freedom adjusted accordingly")
    Obs <- Obs[!empty, , drop = FALSE]
    Exp <- Exp[!empty, , drop = FALSE]
  }

  nGroupsUsed <- nrow(Obs)
  if (nGroupsUsed < 3L)
    stop("fewer than 3 non-empty groups; the test cannot be computed")

  if (any(Exp < 5))
    warning("some expected counts are less than 5; ",
            "the chi-squared approximation may be unreliable")

  ## test statistic -------------------------------------------------------

  STATISTIC <- sum((Obs - Exp)^2 / Exp)
  PARAMETER <- nGroupsUsed - 2L
  PVAL      <- pchisq(STATISTIC, df = PARAMETER, lower.tail = FALSE)

  structure(
    list(
      statistic = c("X-squared" = STATISTIC),
      parameter = c(df = PARAMETER),
      p.value   = PVAL,
      method    = paste("Hosmer-Lemeshow", type, "statistic"),
      type      = type,
      nGroups   = nGroupsUsed,
      observed  = Obs,
      expected  = Exp,
      data.name = DNAME
    ),
    class = c("HosmerLemeshowTest", "htest")
  )
}



#' @param digits number of significant digits to display.
#' @param details logical; if \code{TRUE}, prints observed and expected
#' counts for both outcome classes by group. Default is \code{FALSE}.
#'
#' @rdname hosmerLemeshowTest
#' @export
print.HosmerLemeshowTest <- function(x, digits = 4, details = FALSE, ...) {

  cat("\n ", x$method, "\n\n")
  cat("data: ", x$data.name, "\n")
  cat("X-squared =", format(signif(x$statistic, digits)),
      ", df =", x$parameter,
      ", p-value =", format.pval(x$p.value, digits = digits), "\n")
  cat("Number of groups:", x$nGroups, "\n")

  if (isTRUE(details)) {
    tab <- cbind(
      "Obs 0s" = x$observed[, "0s"],
      "Exp 0s" = round(x$expected[, "0s"], digits),
      "Obs 1s" = x$observed[, "1s"],
      "Exp 1s" = round(x$expected[, "1s"], digits)
    )
    rownames(tab) <- formatC(rownames(tab),
                             width = max(nchar(rownames(tab))))
    cat("\nObserved vs Expected counts by group:\n")
    print(tab)
  }

  cat("\n")
  invisible(x)
}
