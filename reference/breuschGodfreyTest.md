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
  subset,
  na.action = na.omit,
  fill = 0
)
```

## Arguments

- formula:

  a symbolic description for the model to be tested (or a fitted `"lm"`
  object, in which case the model frame is taken from the fit and
  `subset` and `na.action` are ignored).

- data:

  an optional data frame containing the variables in the model. By
  default the variables are taken from the environment which
  `breuschGodfreyTest` is called from. For a fitted `"lm"` object it is
  used for `orderBy` only, as the model itself carries its own model
  frame.

- order:

  integer, the maximal order of serial correlation to be tested. Must be
  smaller than the residual degrees of freedom of the auxiliary
  regression.

- orderBy:

  either a vector `z` or a formula with a single explanatory variable
  like `~ z`. The observations in the model are ordered by the size of
  `z`; a formula with several terms is used as successive ordering keys.
  If set to `NULL` (the default) the observations are assumed to be
  ordered (e.g., a time series). `z` may be given at the length of the
  original data: rows dropped by `subset` or by `na.action` are then
  dropped from `z` as well. Missing values in `z` are ordered last.

- type:

  the type of test statistic to be returned, either `"chisq"` (default)
  for the chi-squared test statistic or `"f"` for the F test statistic.
  Case-insensitive.

- subset:

  an optional expression indicating which observations to use.

- na.action:

  a function specifying how missing values are handled. Defaults to
  [`na.omit()`](https://rdrr.io/r/stats/na.fail.html): the auxiliary
  regression is fitted by
  [`lm.fit()`](https://rdrr.io/r/stats/lmfit.html) and cannot carry
  missing values.

- fill:

  a single value used as starting value for the lagged residuals in the
  auxiliary regression. By default `0` but can also be set to `NA`, in
  which case the leading incomplete rows are dropped.

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

- `df.residual`:

  the residual degrees of freedom of the auxiliary regression, for both
  types of test statistic.

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

Unlike `bgtest()`, the residual degrees of freedom of the auxiliary
regression are reported for both types of test statistic. They are a
property of that regression and not of the statistic derived from it,
whereas `bgtest()` hands back `NULL` for the chi-squared version.
`coeftest()` therefore refers the coefficients to a \\t\\ distribution
in either case, where `bgtest()` switches to the normal one.

## References

Breusch, T. S. (1978) Testing for autocorrelation in dynamic linear
models. *Australian Economic Papers*, 17, 334-355.

Godfrey, L. G. (1978) Testing against general autoregressive and moving
average error models when the regressors include lagged dependent
variables. *Econometrica*, 46, 1293-1301.

## See also

[`durbinWatsonTest()`](durbinWatsonTest.md)

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

## finite sample F version, and dropping the leading lags instead of
## filling them with zeros
breuschGodfreyTest(y2 ~ x, order = 4, type = "f")
#> 
#>  Breusch-Godfrey test for serial correlation of order up to 4
#> 
#> data:  y2 ~ x
#> LM test = 7.2942, df1 = 4, df2 = 94, p-value = 3.682e-05
#> 
breuschGodfreyTest(y2 ~ x, order = 4, fill = NA)
#> 
#>  Breusch-Godfrey test for serial correlation of order up to 4
#> 
#> data:  y2 ~ x
#> LM test = 24.401, df = 4, p-value = 6.637e-05
#> 

## transformed terms and an explicit ordering variable
d <- data.frame(y = as.vector(y2), x = x, z = rnorm(100), tt = sample(100),
                grp = rep(c("A", "B"), each = 50))
breuschGodfreyTest(y ~ x + I(z^2), data = d, orderBy = ~ tt)
#> 
#>  Breusch-Godfrey test for serial correlation of order up to 1
#> 
#> data:  y ~ x + I(z^2)
#> LM test = 1.5791, df = 1, p-value = 0.2089
#> 

## subset and orderBy combined: tt is given at the length of d and is
## reduced to the rows the model frame kept
breuschGodfreyTest(y ~ x, data = d, subset = grp == "A", orderBy = ~ tt)
#> 
#>  Breusch-Godfrey test for serial correlation of order up to 1
#> 
#> data:  y ~ x
#> LM test = 0.83472, df = 1, p-value = 0.3609
#> 

## the test can also be applied to a fitted model
breuschGodfreyTest(lm(y1 ~ x))
#> 
#>  Breusch-Godfrey test for serial correlation of order up to 1
#> 
#> data:  lm(y1 ~ x)
#> LM test = 0.0036887, df = 1, p-value = 0.9516
#> 
```
