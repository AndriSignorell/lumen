# Breusch-Godfrey Test for Detecting Higher-Order Serial Correlation in Regression Residuals

A test for autocorrelation in the residuals of regression models,
generalizing the Durbin-Watson test to handle higher-order
autocorrelation and models with lagged dependent variables.

## Usage

``` r
breuschGodfreyTest(
  formula,
  data = list(),
  order = 1,
  orderBy = NULL,
  type = c("chisq", "f"),
  fill = 0
)
```

## Arguments

- formula:

  a symbolic description for the model to be tested (or a fitted `"lm"`
  object).

- data:

  an optional data frame containing the variables in the model. By
  default the variables are taken from the environment which
  `breuschGodfreyTest` is called from.

- order:

  integer, the maximal order of serial correlation to be tested.

- orderBy:

  either a vector `z` or a formula with a single explanatory variable
  like `~ z`. The observations in the model are ordered by the size of
  `z`. If set to `NULL` (the default) the observations are assumed to be
  ordered (e.g., a time series).

- type:

  the type of test statistic to be returned, either `"chisq"` (default)
  for the chi-squared test statistic or `"f"` for the F test statistic.
  Case-insensitive.

- fill:

  starting values for the lagged residuals in the auxiliary regression.
  By default `0` but can also be set to `NA`.

## Value

A list with class `"breuschGodfreyTest"` inheriting from `"htest"`
containing the following components:

- `statistic`:

  the value of the test statistic.

- `parameter`:

  the degrees of freedom.

- `p.value`:

  the p-value of the test.

- `method`:

  a character string indicating what type of test was performed.

- `data.name`:

  a character string giving the name(s) of the data.

- `coefficients`:

  coefficient estimates from the auxiliary regression.

- `vcov`:

  the corresponding covariance matrix estimate.

## Details

`breuschGodfreyTest` performs the Breusch-Godfrey test for higher-order
serial correlation.

Under \\H_0\\ the test statistic is asymptotically chi-squared with
degrees of freedom as given in `parameter`. If `type` is set to `"f"`
the function returns a finite sample version of the test statistic,
employing an \\F\\ distribution with degrees of freedom as given in
`parameter`.

By default, the starting values for the lagged residuals in the
auxiliary regression are chosen to be 0 (as in Godfrey 1978) but could
also be set to `NA` to omit them.

`breuschGodfreyTest` also returns the coefficients and estimated
covariance matrix from the auxiliary regression that includes the lagged
residuals, accessible via [`coef()`](https://rdrr.io/r/stats/coef.html)
and [`vcov()`](https://rdrr.io/r/stats/vcov.html) on the result. (Note,
however, that standard theory does not always apply to the standard
errors and t-statistics in this regression.)

## Note

Based on code by David Mitchell and Achim Zeileis previously published
as `bgtest()` in the lmtest package, adapted to conform to package
standards.

## References

Breusch, T. S. (1978) Testing for autocorrelation in dynamic linear
models. *Australian Economic Papers*, 17, 334-355.

Godfrey, L. G. (1978) Testing against general autoregressive and moving
average error models when the regressors include lagged dependent
variables. *Econometrica*, 46, 1293-1301.

## See also

Other test.regression: [`bpTest()`](bpTest.md),
[`durbinWatsonTest()`](durbinWatsonTest.md),
[`hosmerLemeshowTest()`](hosmerLemeshowTest.md),
[`leCessieTest()`](leCessieTest.md)

## Examples

``` r
## Generate a stationary and an AR(1) series
set.seed(1)
x <- rep(c(1, -1), 50)

y1 <- 1 + x + rnorm(100)

## Perform Breusch-Godfrey test for first-order serial correlation:
breuschGodfreyTest(y1 ~ x)
#> 
#>  Breusch-Godfrey test for serial correlation of order up to 1
#> 
#> data:  y1 ~ x
#> LM test = 0.0036887, df = 1, p-value = 0.9516
#> 

## or for fourth-order serial correlation
breuschGodfreyTest(y1 ~ x, order = 4)
#> 
#>  Breusch-Godfrey test for serial correlation of order up to 4
#> 
#> data:  y1 ~ x
#> LM test = 3.0822, df = 4, p-value = 0.5442
#> 

## Compare with Durbin-Watson test results:
durbinWatsonTest(y1 ~ x)
#> 
#>  Durbin-Watson test
#> 
#> data:  y1 ~ x
#> DW = 1.9762, p-value = 0.4924
#> alternative hypothesis: true autocorrelation is greater than 0
#> 

y2 <- stats::filter(y1, 0.5, method = "recursive")
breuschGodfreyTest(y2 ~ x)
#> 
#>  Breusch-Godfrey test for serial correlation of order up to 1
#> 
#> data:  y2 ~ x
#> LM test = 19.907, df = 1, p-value = 8.128e-06
#> 
```
