# Durbin-Watson Test for Detecting First-Order Autocorrelation in Regression Residuals

Tests for first-order autocorrelation in the residuals of a linear
regression model, based on the ratio of successive squared residual
differences to the total residual sum of squares.

## Usage

``` r
durbinWatsonTest(x, ...)

# S3 method for class 'formula'
durbinWatsonTest(
  x,
  data = list(),
  orderBy = NULL,
  alternative = c("greater", "two.sided", "less"),
  iterations = 15,
  exact = NULL,
  tol = 0.0000000001,
  ...
)

# S3 method for class 'lm'
durbinWatsonTest(
  x,
  orderBy = NULL,
  alternative = c("greater", "two.sided", "less"),
  iterations = 15,
  exact = NULL,
  tol = 0.0000000001,
  ...
)

# S3 method for class 'numeric'
durbinWatsonTest(
  x,
  orderBy = NULL,
  alternative = c("greater", "two.sided", "less"),
  iterations = 15,
  exact = NULL,
  tol = 0.0000000001,
  ...
)

# Default S3 method
durbinWatsonTest(x, ...)
```

## Arguments

- x:

  a symbolic description of the model to be tested (a `formula`), a
  fitted `"lm"` object, or a numeric vector of residuals.

- ...:

  further arguments passed to or from other methods.

- data:

  an optional data frame containing the variables in the model. Only
  used for the `formula` method. By default the variables are taken from
  the environment which `durbinWatsonTest` is called from.

- orderBy:

  either a vector `z` or a formula with a single explanatory variable
  like `~ z`. The observations in the model are ordered by the size of
  `z`. If set to `NULL` (the default) the observations are assumed to be
  ordered (e.g., a time series).

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"greater"` (default), `"two.sided"` or `"less"`.

- iterations:

  an integer specifying the number of iterations used by the "pan"
  algorithm when computing the exact p-value.

- exact:

  logical. If `TRUE` the exact p-value is computed via the "pan"
  algorithm; if `FALSE` a normal approximation is used. The default is
  `TRUE` for sample sizes below 100 and `FALSE` otherwise.

- tol:

  numeric tolerance. Eigenvalues smaller than `tol` are treated as zero.

## Value

An object of class `"htest"` containing the following components:

- statistic:

  the Durbin-Watson test statistic.

- p.value:

  the p-value of the test.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string indicating the test performed.

- data.name:

  a character string describing the data.

## Details

The Durbin-Watson test has the null hypothesis that the autocorrelation
of the disturbances is 0. The alternative hypothesis can be specified as
greater than, not equal to, or less than 0 via the `alternative`
argument.

Under the assumption of normally distributed disturbances, the null
distribution of the Durbin-Watson statistic is the distribution of a
linear combination of chi-squared variables. The p-value is computed
using the "pan" algorithm (Farebrother, 1980, 1984), implemented via
Rcpp. For large sample sizes the algorithm might fail to compute the
p-value; in that case a warning is issued and an approximate p-value is
returned instead, computed via a normal approximation using the mean and
variance of the Durbin-Watson statistic.

Three methods are dispatched:

- `formula`:

  Fits the model from scratch using `data`.

- `lm`:

  Extracts design matrix and response from a fitted `"lm"` object. The
  `data.name` field in the result reflects the model formula, not the
  name of the object.

- `numeric`:

  Treats `x` as a vector of pre-computed residuals and fits an
  intercept-only design matrix. Note that this is *not* equivalent to
  testing residuals from a fitted model with predictors; use the `lm` or
  `formula` method when a model exists.

## Note

Based on code by Torsten Hothorn, Achim Zeileis, Clint Cummins, Giovanni
Millo and David Mitchell previously published in the lmtest package,
with an Rcpp reimplementation of the "pan" algorithm, adapted to conform
to package standards.

## References

Durbin, J. and Watson, G. S. (1950) Testing for serial correlation in
least squares regression I. *Biometrika*, 37, 409-428.

Durbin, J. and Watson, G. S. (1951) Testing for serial correlation in
least squares regression II. *Biometrika*, 38, 159-178.

Durbin, J. and Watson, G. S. (1971) Testing for serial correlation in
least squares regression III. *Biometrika*, 58, 1-19.

Farebrother, R. W. (1980) Pan's procedure for the tail probabilities of
the Durbin-Watson statistic. *Applied Statistics*, 29, 224-227.

Farebrother, R. W. (1984) AS R53: A remark on algorithms AS 106, AS 153
and AS 155. *Applied Statistics*, 33, 366-369.

Kraemer, W. and Sonnberger, H. (1986) *The Linear Regression Model under
Test*. Heidelberg: Physica.

## See also

[`lm()`](https://rdrr.io/r/stats/lm.html)

Other test.regression:
[`bpTest()`](https://andrisignorell.github.io/lumen/reference/bpTest.md),
[`breuschGodfreyTest()`](https://andrisignorell.github.io/lumen/reference/breuschGodfreyTest.md),
[`hosmerLemeshowTest()`](https://andrisignorell.github.io/lumen/reference/hosmerLemeshowTest.md),
[`leCessieTest()`](https://andrisignorell.github.io/lumen/reference/leCessieTest.md)

## Examples

``` r
## formula method
set.seed(1)
x <- rep(c(-1, 1), 50)

err1 <- rnorm(100)
durbinWatsonTest(y ~ x, data = data.frame(y = 1 + x + err1, x = x))
#> 
#>  Durbin-Watson test
#> 
#> data:  y ~ x
#> DW = 1.9762, p-value = 0.4924
#> alternative hypothesis: true autocorrelation is greater than 0
#> 

## autocorrelated errors (rho = 0.9)
err2 <- stats::filter(err1, 0.9, method = "recursive")
durbinWatsonTest(y ~ x, data = data.frame(y = 1 + x + err2, x = x))
#> 
#>  Durbin-Watson test
#> 
#> data:  y ~ x
#> DW = 0.45961, p-value = 7.862e-15
#> alternative hypothesis: true autocorrelation is greater than 0
#> 

## lm method
fit <- lm(y ~ x, data = data.frame(y = 1 + x + err1, x = x))
durbinWatsonTest(fit)
#> 
#>  Durbin-Watson test
#> 
#> data:  y ~ x
#> DW = 1.9762, p-value = 0.4924
#> alternative hypothesis: true autocorrelation is greater than 0
#> 

## numeric method (pre-computed residuals, intercept-only design assumed)
e_t <- c(-32.33, -26.603, 2.215, -16.967, -1.148, -2.512, -1.967, 11.669,
         -0.513, 27.032, -4.422, 40.032, 23.577, 33.94, -2.787, -8.606,
          0.575, 6.848, -18.971, -29.063)
durbinWatsonTest(e_t)
#> 
#>  Durbin-Watson test
#> 
#> data:  e_t
#> DW = 1.08, p-value = 0.01328
#> alternative hypothesis: true autocorrelation is greater than 0
#> 
```
