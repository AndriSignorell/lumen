
#' le Cessie-van Houwelingen-Copas-Hosmer Global Goodness of Fit Test
#'
#' Computes the le Cessie-van Houwelingen-Copas-Hosmer unweighted sum of
#' squares test for global goodness of fit of a logistic regression model.
#'
#' The test compares the observed sum of squared residuals to its expected
#' value under the null hypothesis of correct model specification. The
#' standardised difference follows an approximate standard normal distribution.
#'
#' Unlike the Hosmer-Lemeshow tests, this test does not rely on grouping and
#' is therefore sensitive to a different class of model misspecification.
#' It requires the covariate matrix \code{X} used to fit the model.
#'
#' @name leCessieTest
#' @param fit numeric vector of fitted probabilities, each in \code{[0, 1]},
#'   without missing values.
#' @param obs numeric vector of observed binary outcomes (0 or 1), of the same
#'   length as \code{fit}, without missing values.
#' @param X numeric matrix of covariates, with \code{nrow(X) == length(fit)},
#'   without missing values. Use \code{model.matrix(fit)[, -1, drop = FALSE]}
#'   to ensure a matrix is passed even for single-predictor models.
#'
#' @return An object of class \code{c("LeCessieTest", "htest")}, which is a
#'   list with components:
#'   \item{statistic}{the standardised Z statistic, named \code{"Z"}.}
#'   \item{p.value}{two-sided p-value from the standard normal distribution.}
#'   \item{method}{a character string describing the test.}
#'   \item{sse}{observed sum of squared errors.}
#'   \item{expected}{expected sum of squared errors under H0.}
#'   \item{sd}{standard deviation of the sum of squared errors under H0.}
#'   \item{data.name}{a character string with the names of \code{fit} and
#'     \code{obs}.}
#'
#' @note
#' Adapted from code by Matthias Kohl (MKmisc) to conform to package
#' standards.
#'
#' @seealso \code{\link{hosmerLemeshowTest}}, \code{\link{glm}}
#'
#' @references
#' Hosmer, D.W., Hosmer, T., le Cessie, S., Lemeshow, S. (1997). A comparison
#' of goodness-of-fit tests for the logistic regression model.
#' \emph{Statistics in Medicine}, \bold{16}, 965--980.
#'
#' @examples
#' set.seed(111)
#' x1 <- factor(sample(1:3, 50, replace = TRUE))
#' x2 <- rnorm(50)
#' obs <- sample(c(0, 1), 50, replace = TRUE)
#' fit <- glm(obs ~ x1 + x2, family = binomial)
#'
#' leCessieTest(fit = fitted(fit), obs = obs, X = model.matrix(fit)[, -1, drop = FALSE])
#'
#' @rdname leCessieTest



#' @family test.regression  
#' @concept regression-diagnostics  
#' @concept goodness-of-fit
#'
#'
#' @export
leCessieTest <- function(fit, obs, X) {
  
  # --- input validation -------------------------------------------------------
  
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
  if (!is.matrix(X) || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
  if (nrow(X) != length(fit))
    stop("'X' must have the same number of rows as the length of 'fit'.")
  if (anyNA(X))
    stop("'X' must not contain missing values.")
  
  # --- test statistic ---------------------------------------------------------
  
  p   <- fit
  y   <- obs == 1L
  wt  <- p * (1 - p)
  sse <- sum((y - p)^2)
  ev  <- sum(wt)
  
  # weighted projection of (1 - 2p) onto X to obtain SD of SSE under H0
  d      <- 1 - 2 * p
  z      <- lm.wfit(X, d, wt, method = "qr")
  sd     <- sqrt(sum((z$residuals * sqrt(z$weights))^2))
  
  if (sd <= .Machine$double.eps)
    stop(
      "Unable to compute test statistic: ",
      "the estimated standard deviation is zero."
    )
  
  z_stat  <- (sse - ev) / sd
  p.value <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
  
  data.name <- paste(
    paste(deparse(substitute(fit)), collapse = ""),
    "and",
    paste(deparse(substitute(obs)), collapse = "")
  )
  method <- paste(
    "le Cessie-van Houwelingen-Copas-Hosmer",
    "global goodness of fit test"
  )
  
  # --- return -----------------------------------------------------------------
  
  structure(
    list(
      statistic = c("Z" = z_stat),
      p.value   = p.value,
      method    = method,
      sse       = sse,
      expected  = ev,
      sd        = sd,
      data.name = data.name
    ),
    class = c("LeCessieTest", "htest")
  )
}


#' @param x object of class \code{"LeCessieTest"}.
#' @param digits number of significant digits to display.
#' @param ... further arguments passed to or from methods.
#'

#' @rdname leCessieTest
#' @export
print.LeCessieTest <- function(x, digits = 4, ...) {
  cat("\n ", x$method, "\n\n")
  cat("data: ", x$data.name, "\n")
  cat(
    "Z =", format(signif(x$statistic, digits)),
    ", p-value =", format.pval(x$p.value, digits = digits), "\n"
  )
  cat(
    "Sum of squared errors:", format(signif(x$sse, digits)),
    "(expected:", format(signif(x$expected, digits)), ")\n\n"
  )
  invisible(x)
}



