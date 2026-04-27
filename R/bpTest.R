
#' Breusch-Pagan Test for Heteroscedasticity
#'
#' Tests the null hypothesis of homoscedasticity (constant error variance)
#' against heteroscedasticity using the Koenker variant of the Breusch-Pagan
#' test. The test statistic is \eqn{BP = n \cdot R^2} from an auxiliary
#' regression of squared residuals on fitted values, asymptotically
#' distributed as \eqn{\chi^2} with \eqn{k} degrees of freedom, where
#' \eqn{k} is the number of predictors.
#'
#' @param lm_fit a fitted \code{\link[stats]{lm}} object.
#'
#' @return An object of class \code{"htest"} with the following components:
#'   \describe{
#'     \item{\code{statistic}}{the BP test statistic.}
#'     \item{\code{parameter}}{degrees of freedom.}
#'     \item{\code{p.value}}{p-value based on the \eqn{\chi^2} distribution.}
#'     \item{\code{method}}{character string describing the test.}
#'     \item{\code{data.name}}{the formula of the fitted model.}
#'   }
#'
#' @references
#'   Breusch, T.S. and Pagan, A.R. (1979). A simple test for heteroscedasticity
#'   and random coefficient variation. \emph{Econometrica}, 47, 1287--1294.
#'
#'   Koenker, R. (1981). A note on studentizing a test for heteroscedasticity.
#'   \emph{Journal of Econometrics}, 17, 107--112.
#'
#' @seealso \code{\link[stats]{lm}}
#'
#' @family linear_models
#' @concept heteroscedasticity homoscedasticity regression diagnostics
#'
#' @examples
#' fit <- lm(Sepal.Length ~ Sepal.Width, data = iris)
#' bpTest(fit)
#'



#' @export
bpTest <- function(lm_fit) {

  if (!inherits(lm_fit, "lm"))
    stop("'lm_fit' must be a fitted lm object")

  n      <- length(residuals(lm_fit))
  resid2 <- residuals(lm_fit)^2
  aux    <- lm(resid2 ~ fitted(lm_fit))
  r2     <- summary(aux)$r.squared
  stat   <- n * r2
  df     <- length(coef(aux)) - 1L
  p_val  <- pchisq(stat, df = df, lower.tail = FALSE)

  structure(
    list(
      statistic = c(BP = stat),
      parameter = c(df = df),
      p.value   = p_val,
      method    = "Breusch-Pagan test (Koenker)",
      data.name = deparse(formula(lm_fit))
    ),
    class = "htest"
  )
}
