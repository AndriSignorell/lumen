# Augmented Dickey-Fuller Unit Root Test for Detecting Unit Roots in Time Series

Performs the augmented Dickey-Fuller test for a unit root in a time
series, testing the null hypothesis of nonstationarity against the
alternative of (trend-)stationarity. It extends the simple Dickey-Fuller
test by including lagged difference terms to account for
autocorrelation.

## Usage

``` r
adfTest(
  y,
  type = c("none", "drift", "trend"),
  lags = 1,
  selectLags = c("fixed", "aic", "bic")
)
```

## Arguments

- y:

  numeric vector or univariate time series to be tested for a unit root.

- type:

  the deterministic part of the test regression, one of `"none"`
  (default), `"drift"` or `"trend"`.

- lags:

  the number of lagged difference terms to be included. If `selectLags`
  is not `"fixed"`, this is the maximum number of lags considered in the
  lag selection.

- selectLags:

  the lag selection method, one of `"fixed"` (default, use `lags` as
  given), `"aic"` or `"bic"` (choose the lag order up to `lags` that
  minimizes the respective information criterion). Case-insensitive, so
  `"AIC"` and `"BIC"` are accepted as well.

## Value

An object of class `"htest"` containing the following components:

- `statistic`:

  the tau statistic, followed by the phi statistic(s) if `type` is
  `"drift"` or `"trend"`. The p-value refers to the tau statistic.

- `parameter`:

  the number of lagged differences included in the test regression
  (after lag selection, if requested).

- `p.value`:

  the interpolated p-value of the tau statistic.

- `critical.values`:

  a matrix with the 1%, 5% and 10% critical values of all reported
  statistics, interpolated for the effective sample size (not shown on
  screen).

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string indicating the test performed.

- `data.name`:

  a character string giving the name of the data.

## Details

If `type` is set to `"none"` neither an intercept nor a trend is
included in the test regression, `"drift"` adds an intercept and
`"trend"` adds both an intercept and a linear trend.

The reported test statistic is the t statistic of the lagged level
(`tau1`, `tau2` or `tau3`, depending on `type`). For `type = "drift"`
and `type = "trend"` the F statistics `phi1` resp. `phi2` and `phi3` of
Dickey and Fuller (1981), testing joint hypotheses on the deterministic
terms, are reported as additional statistics; their critical values are
contained in the `critical.values` component of the result.

The p-value refers to the tau statistic and is obtained by linear
interpolation in the finite sample quantiles given by Fuller (1976),
following the approach of
[`tseries::adf.test()`](https://rdrr.io/pkg/tseries/man/adf.test.html).
If the statistic falls outside the range of the table, the p-value is
reported as the respective boundary (0.01 or 0.99) and a warning is
issued. The critical values are taken from Hamilton (1994) and Dickey
and Fuller (1981).

Missing values are not allowed.

## Note

Based on code by Bernhard Pfaff previously published in the urca
package, adapted to conform to package standards.

## References

Dickey, D. A. and Fuller, W. A. (1979) Distribution of the estimators
for autoregressive time series with a unit root, *Journal of the
American Statistical Association*, **74**, 427–431.

Dickey, D. A. and Fuller, W. A. (1981) Likelihood ratio statistics for
autoregressive time series with a unit root, *Econometrica*, **49**,
1057–1072.

Fuller, W. A. (1976) *Introduction to Statistical Time Series*, New
York: Wiley.

Hamilton, J. D. (1994) *Time Series Analysis*, Princeton: Princeton
University Press.

## See also

[`kpssTest`](kpssTest.md)

Other test.timeseries: [`bartelsRankTest()`](BartelsRankTest.md),
[`kpssTest()`](kpssTest.md), [`runsTest()`](runsTest.md),
[`vonNeumannTest()`](vonNeumannTest.md)

## Examples

``` r
adfTest(AirPassengers, lags = 3, type = "trend")
#> Warning: p-value smaller than reported p-value
#> 
#>  Augmented Dickey-Fuller Test
#> 
#> data:  AirPassengers
#> tau3 = -6.9358, phi2 = 16.4042, phi3 = 24.0601, lags = 3, p-value =
#> 0.01
#> alternative hypothesis: trend-stationary
#> 

# a random walk should not be rejected
set.seed(5)
rw <- cumsum(rnorm(200))
adfTest(rw, type = "drift")
#> 
#>  Augmented Dickey-Fuller Test
#> 
#> data:  rw
#> tau2 = -2.8598, phi1 = 4.1428, lags = 1, p-value = 0.05381
#> alternative hypothesis: stationary
#> 
```
