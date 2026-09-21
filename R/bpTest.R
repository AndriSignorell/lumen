
#' Breusch-Pagan Test  for Detecting Heteroscedasticity in Regression Models
#'
#' Tests the null hypothesis of homoscedasticity (constant error variance)
#' against heteroscedasticity using the Koenker variant of the Breusch-Pagan
#' test. The test statistic is \eqn{BP = n \cdot R^2} from an auxiliary
#' regression of squared residuals on the regressors of the model, asymptotically
#' distributed as \eqn{\chi^2} with \eqn{k} degrees of freedom, where
#' \eqn{k} is the number of predictors.
#'
#' @param fit a fitted [lm()] object.
#'
#' @return An object of class `"htest"` with the following components:
#'     \item{`statistic`}{the BP test statistic.}
#'     \item{`parameter`}{degrees of freedom.}
#'     \item{`p.value`}{p-value based on the \eqn{\chi^2} distribution.}
#'     \item{`method`}{character string describing the test.}
#'     \item{`data.name`}{the formula of the fitted model.}
#'
#' @references
#'   Breusch, T.S. and Pagan, A.R. (1979). A simple test for heteroscedasticity
#'   and random coefficient variation. *Econometrica*, 47, 1287--1294.
#'
#'   Koenker, R. (1981). A note on studentizing a test for heteroscedasticity.
#'   *Journal of Econometrics*, 17, 107--112.
#'
#' @seealso [lm()]
#'
#' @examples
#' fit <- lm(Sepal.Length ~ Sepal.Width, data = iris)
#' bpTest(fit)
#'




#' @family test.regression  
#'
#' @export
bpTest <- function(fit) {

  if (!inherits(fit, "lm") || inherits(fit, "glm"))
    stop("'fit' must be a fitted lm object")

  w <- fit$weights
  if (!is.null(w) && !isTRUE(all.equal(as.vector(w), rep(1, length(w)))))
    stop("weighted regressions are not supported", call. = FALSE)

  # fit$residuals, not residuals(): with na.exclude residuals() pads the
  # excluded rows with NA, and length() counted them into n
  e2 <- fit$residuals^2
  n  <- length(e2)

  # Koenker's studentized version regresses e^2 on the regressors of the
  # model, not on the fitted values: with more than one regressor the
  # latter is a different test with 1 df, not the documented one with k
  X   <- model.matrix(fit)
  aux <- lm.fit(X, e2)

  r2    <- 1 - sum(aux$residuals^2) / sum((e2 - mean(e2))^2)
  stat  <- n * r2
  df    <- aux$rank - 1L
  p_val <- pchisq(stat, df = df, lower.tail = FALSE)

  structure(
    list(
      statistic = c(BP = stat),
      parameter = c(df = df),
      p.value   = p_val,
      method    = "Breusch-Pagan test (Koenker)",
      data.name = deparse1(formula(fit))
    ),
    class = "htest"
  )
}
