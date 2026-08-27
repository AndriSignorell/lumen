# Kwiatkowski-Phillips-Schmidt-Shin Test for Assessing Level or Trend Stationarity

A test for stationarity in time series, complementary to unit root tests
such as the ADF test: it tests the null hypothesis of stationarity
against the alternative of a unit root.

## Usage

``` r
kpssTest(
  y,
  type = c("mu", "tau"),
  lags = c("short", "long", "nil"),
  useLag = NULL
)
```

## Arguments

- y:

  numeric vector or univariate time series to be tested for
  stationarity.

- type:

  the deterministic part of the model, one of `"mu"` (default, constant)
  or `"tau"` (constant plus linear trend).

- lags:

  the rule for the number of lags used for the error term correction,
  one of `"short"` (default), `"long"` or `"nil"`. See the Details.
  Ignored if `useLag` is given.

- useLag:

  an optional integer explicitly specifying the number of lags,
  overriding `lags`.

## Value

An object of class `"htest"` containing the following components:

- statistic:

  the value of the KPSS test statistic.

- parameter:

  the number of lags used for the error term correction.

- p.value:

  the interpolated p-value of the test.

- critical.values:

  the asymptotic critical values at the 10%, 5%, 2.5% and 1%
  significance levels (not shown on screen).

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string indicating the test performed.

- data.name:

  a character string giving the name of the data.

## Details

Performs the KPSS test (Kwiatkowski et al., 1992), where the null
hypothesis is stationarity. The test types specify as deterministic
component either a constant (`"mu"`, null hypothesis of level
stationarity) or a constant with linear trend (`"tau"`, null hypothesis
of trend stationarity).

`lags = "short"` sets the number of lags used for the long-run variance
estimation to \\4 (n/100)^{1/4}\\, whereas `lags = "long"` sets it to
\\12 (n/100)^{1/4}\\ (each truncated to an integer). With `lags = "nil"`
no error correction is made. Alternatively, an explicit number of lags
can be given via `useLag`, which then takes precedence.

The p-value is obtained by linear interpolation in the asymptotic
critical values of Kwiatkowski et al. (1992, Table 1), following the
approach of
[`tseries::kpss.test()`](https://rdrr.io/pkg/tseries/man/kpss.test.html).
If the statistic falls outside the range of the table, the p-value is
reported as the respective boundary (0.01 or 0.10) and a warning is
issued.

Missing values are silently removed.

## Note

Based on code by Bernhard Pfaff previously published in the urca
package, adapted to conform to package standards.

## References

Kwiatkowski, D., Phillips, P. C. B., Schmidt, P. and Shin, Y. (1992)
Testing the null hypothesis of stationarity against the alternative of a
unit root: How sure are we that economic time series have a unit root?
*Journal of Econometrics*, **54**, 159–178.

## See also

Other test.timeseries:
[`adfTest()`](https://andrisignorell.github.io/lumen/reference/adfTest.md),
[`bartelsRankTest()`](https://andrisignorell.github.io/lumen/reference/BartelsRankTest.md),
[`runsTest()`](https://andrisignorell.github.io/lumen/reference/runsTest.md),
[`vonNeumannTest()`](https://andrisignorell.github.io/lumen/reference/vonNeumannTest.md)

## Examples

``` r
# trend-stationary series: null hypothesis is not rejected
set.seed(1)
x <- 0.2 * seq_len(200) + rnorm(200)
kpssTest(x, type = "tau")
#> Warning: p-value greater than reported p-value
#> 
#>  KPSS test for trend stationarity
#> 
#> data:  x
#> KPSS = 0.065919, lags = 4, p-value = 0.1
#> alternative hypothesis: not trend stationary (unit root)
#> 

# random walk: null hypothesis of level stationarity is rejected
set.seed(2)
rw <- cumsum(rnorm(200))
kpssTest(rw, type = "mu")
#> 
#>  KPSS test for level stationarity
#> 
#> data:  rw
#> KPSS = 0.51701, lags = 4, p-value = 0.03784
#> alternative hypothesis: not level stationary (unit root)
#> 

kpssTest(AirPassengers, type = "tau", lags = "short")
#> Warning: p-value greater than reported p-value
#> 
#>  KPSS test for trend stationarity
#> 
#> data:  AirPassengers
#> KPSS = 0.09615, lags = 4, p-value = 0.1
#> alternative hypothesis: not trend stationary (unit root)
#> 
```
