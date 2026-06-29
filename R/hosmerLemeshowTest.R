
#' Hosmer-Lemeshow Goodness of Fit Test
#'
#' Computes the Hosmer-Lemeshow C or H goodness-of-fit test for a logistic
#' regression model, assessing whether observed event rates match predicted
#' probabilities across grouped subsets of the data.
#'
#' The C statistic groups observations by quantiles of the fitted probabilities
#' (equal-frequency bins), while the H statistic uses equal-width intervals.
#' Both statistics are asymptotically compared with a chi-squared distribution
#' with \code{nGroups - 2} degrees of freedom under the null hypothesis of
#' adequate fit; the approximation may be unreliable in small samples.
#'
#' @name hosmerLemeshowTest
#' @param fit numeric vector of fitted probabilities, each in \code{[0, 1]},
#'   without missing values.
#' @param obs numeric vector of observed binary outcomes (0 or 1), of the same
#'   length as \code{fit}, without missing values.
#' @param nGroups integer, number of groups. Must be \code{>= 3}.
#'   Default is \code{10}.
#' @param type character, one of \code{"C"} (quantile-based groups) or
#'   \code{"H"} (equal-width groups).
#'
#' @return An object of class \code{c("HosmerLemeshowTest", "htest")}, which
#'   is a list with components:
#'   \item{statistic}{the chi-squared test statistic, named \code{"X-squared"}.}
#'   \item{parameter}{degrees of freedom (\code{nGroups - 2}), named
#'     \code{"df"}.}
#'   \item{p.value}{p-value from the chi-squared distribution.}
#'   \item{method}{a character string describing the test.}
#'   \item{type}{the type of statistic computed (\code{"C"} or \code{"H"}).}
#'   \item{nGroups}{number of groups actually used (may be less than requested
#'     for type \code{"C"} if fewer unique quantile breaks are found).}
#'   \item{observed}{matrix of observed counts for both outcome classes per
#'     group.}
#'   \item{expected}{matrix of expected counts for both outcome classes per
#'     group.}
#'   \item{data.name}{a character string with the names of \code{fit} and
#'     \code{obs}.}
#'
#'   The \code{print} method accepts a \code{details} argument; if \code{TRUE},
#'   observed and expected counts for both outcome classes are printed by group.
#'
#' @seealso \code{\link{leCessieTest}}, \code{\link{glm}}
#'
#' @references
#' Lemeshow, S., Hosmer, D.W. (1982). A review of goodness of fit statistics
#' for use in the development of logistic regression models.
#' \emph{American Journal of Epidemiology}, \bold{115}(1), 92--106.
#'
#' @examples
#' set.seed(111)
#' x1 <- factor(sample(1:3, 50, replace = TRUE))
#' x2 <- rnorm(50)
#' obs <- sample(c(0, 1), 50, replace = TRUE)
#' fit <- glm(obs ~ x1 + x2, family = binomial)
#'
#' res <- hosmerLemeshowTest(fit = fitted(fit), obs = obs, type = "C")
#' res
#'
#' print(res, details = TRUE)
#'
#' @rdname hosmerLemeshowTest

#' @family test.regression  
#' @concept regression-diagnostics  
#' @concept goodness-of-fit  
#' @concept calibration
#'
#'
#' @export
hosmerLemeshowTest <- function(fit, obs, nGroups = 10, type = c("C", "H")) {
  
  # --- input validation -------------------------------------------------------
  
  type <- match.arg(type)
  
  if (!is.numeric(fit) || !is.numeric(obs))
    stop("'fit' and 'obs' must be numeric vectors.")
  if (length(fit) != length(obs))
    stop("'fit' and 'obs' must have the same length.")
  if (anyNA(fit))
    stop("'fit' must not contain missing values.")
  if (anyNA(obs))
    stop("'obs' must not contain missing values.")
  if (any(fit < 0 | fit > 1))
    stop("'fit' must contain probabilities in [0, 1].")
  if (!all(obs %in% c(0, 1)))
    stop("'obs' must be binary (0 or 1 only).")
  if (!is.numeric(nGroups) || length(nGroups) != 1L || nGroups < 3)
    stop("'nGroups' must be a single integer >= 3.")
  if (nGroups != as.integer(nGroups))
    stop("'nGroups' must be an integer.")
  
  nGroups <- as.integer(nGroups)
  
  # --- grouping ---------------------------------------------------------------
  
  if (type == "C") {
    brks <- unique(
      quantile(fit, probs = seq(0, 1, by = 1 / nGroups))
    )
    nGroups_actual <- length(brks) - 1L
    
    if (nGroups_actual < 3L)
      stop("Unable to construct at least 3 groups from fitted probabilities.")
    if (nGroups_actual < nGroups)
      warning(
        "Found only ", nGroups_actual,
        " distinct groups for Hosmer-Lemeshow C statistic."
      )
    
    nGroups <- nGroups_actual
    cutfit  <- cut(fit, breaks = brks, include.lowest = TRUE)
    
  } else {
    cutfit <- cut(fit, breaks = nGroups, include.lowest = TRUE)
  }
  
  # --- observed / expected ----------------------------------------------------
  
  Obs <- xtabs(cbind("0s" = 1 - obs, "1s" = obs) ~ cutfit)
  Exp <- xtabs(cbind("0s" = 1 - fit, "1s" = fit) ~ cutfit)
  
  # --- test statistic ---------------------------------------------------------
  
  if (any(Exp <= .Machine$double.eps))
    warning(
      "Some expected cell counts are near zero; ",
      "the chi-squared approximation may be unreliable."
    )
  
  chisq   <- sum((Obs - Exp)^2 / Exp, na.rm = TRUE)
  df      <- nGroups - 2L
  p.value <- pchisq(chisq, df = df, lower.tail = FALSE)
  
  data.name <- paste(
    paste(deparse(substitute(fit)), collapse = ""),
    "and",
    paste(deparse(substitute(obs)), collapse = "")
  )
  method <- paste("Hosmer-Lemeshow", type, "statistic")
  
  # --- return -----------------------------------------------------------------
  
  structure(
    list(
      statistic = c("X-squared" = chisq),
      parameter = c("df" = df),
      p.value   = p.value,
      method    = method,
      type      = type,
      nGroups   = nGroups,
      observed  = Obs,
      expected  = Exp,
      data.name = data.name
    ),
    class = c("HosmerLemeshowTest", "htest")
  )
}


#' @param x object of class \code{"HosmerLemeshowTest"}.
#' @param digits number of significant digits to display.
#' @param details logical; if \code{TRUE}, prints observed and expected counts
#'   for both outcome classes by group. Default is \code{FALSE}.
#' @param ... further arguments passed to or from methods.
#'
#' @rdname hosmerLemeshowTest
#' @export
print.HosmerLemeshowTest <- function(x, digits = 4, details = FALSE, ...) {
  cat("\n ", x$method, "\n\n")
  cat("data: ", x$data.name, "\n")
  cat(
    "X-squared =", format(signif(x$statistic, digits)),
    ", df =", x$parameter,
    ", p-value =", format.pval(x$p.value, digits = digits), "\n"
  )
  cat("Number of groups:", x$nGroups, "\n")
  
  if (isTRUE(details)) {
    tab <- cbind(
      "Obs 0s" = x$observed[, "0s"],
      "Exp 0s" = round(x$expected[, "0s"], digits),
      "Obs 1s" = x$observed[, "1s"],
      "Exp 1s" = round(x$expected[, "1s"], digits)
    )
    rownames(tab) <- formatC(rownames(tab), width = max(nchar(rownames(tab))))
    cat("\nObserved vs Expected counts by group:\n")
    print(tab)
  }
  
  cat("\n")
  invisible(x)
}

