# (Robust) Jarque-Bera Test for Assessing Normality From Skewness and Kurtosis

A goodness-of-fit test for normality based on sample skewness and
kurtosis, in its classical form or the robust modification of Gel and
Gastwirth (2008).

## Usage

``` r
jarqueBeraTest(x, robust = TRUE, method = c("chisq", "mc"), R = 1000)
```

## Arguments

- x:

  a numeric vector of data values.

- robust:

  logical, if `TRUE` (default) the robust test of Gel and
  Gastwirth (2008) is performed, otherwise the classical Jarque-Bera
  test.

- method:

  a character string specifying how the p-value is computed, one of
  `"chisq"` (default, asymptotic chi-squared approximation) or `"mc"`
  (Monte Carlo simulation).

- R:

  the number of Monte Carlo replicates used for `method = "mc"` (default
  is `1000`).

## Value

A list with class `"htest"` containing the following components:

- statistic:

  the value of the test statistic.

- parameter:

  the degrees of freedom (`method = "chisq"`), or the number of Monte
  Carlo replicates (`method = "mc"`).

- p.value:

  the p-value of the test.

- method:

  a character string indicating the test performed and how the p-value
  was computed.

- data.name:

  a character string giving the name of the data.

## Details

This function performs either the classical Jarque-Bera test or the
robust Jarque-Bera test proposed by Gel and Gastwirth (2008). The robust
version (the default) replaces the classical standard deviation by the
average absolute deviation from the median (MAAD), which makes the test
considerably less sensitive to outliers. The kurtosis component is then
scaled with the constant \\C_2 = 64\\ derived in Gel and Gastwirth
(2008), the skewness component with \\C_1 = 6\\ as in the classical
test.

With `method = "chisq"` the statistic is referred to its asymptotic
chi-squared distribution with 2 degrees of freedom. With `method = "mc"`
the p-value is estimated by Monte Carlo simulation from `R` samples of
the same size drawn from the standard normal distribution, using the
finite-sample correction \\(m + 1)/(R + 1)\\, where \\m\\ is the number
of simulated statistics at least as large as the observed one.

Missing values are silently removed.

## Note

Based on code by W. Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth and
Weiwen Miao previously published as `rjb.test()` in the lawstat package,
adapted to conform to package standards.

## References

Gel, Y. R. and Gastwirth, J. L. (2008) A robust modification of the
Jarque-Bera test of normality. *Economics Letters*, 99, 30-32.

Jarque, C. and Bera, A. (1980) Efficient tests for normality,
homoscedasticity and serial independence of regression residuals.
*Economics Letters*, 6, 255-259.

## See also

[`shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html)

Other test.normality: [`andersonDarlingTest()`](andersonDarlingTest.md),
[`cramerVonMisesTest()`](cramerVonMisesTest.md),
[`lillieTest()`](lillieTest.md), [`pearsonTest()`](pearsonTest.md),
[`shapiroFranciaTest()`](shapiroFranciaTest.md)

## Examples

``` r
set.seed(1)
x <- rnorm(100)

jarqueBeraTest(x)                   # robust version (default)
#> 
#>  Robust Jarque-Bera test (chi-squared approximation)
#> 
#> data:  x
#> X-squared = 0.086947, df = 2, p-value = 0.9575
#> 
jarqueBeraTest(x, robust = FALSE)   # classical Jarque-Bera test
#> 
#>  Jarque-Bera test (chi-squared approximation)
#> 
#> data:  x
#> X-squared = 0.087201, df = 2, p-value = 0.9573
#> 

# Monte Carlo p-value
jarqueBeraTest(x, method = "mc", R = 5000)
#> 
#>  Robust Jarque-Bera test (Monte Carlo)
#> 
#> data:  x
#> X-squared = 0.086947, R = 5000, p-value = 0.9508
#> 

# heavy-tailed alternative
jarqueBeraTest(rt(100, df = 2))
#> 
#>  Robust Jarque-Bera test (chi-squared approximation)
#> 
#> data:  rt(100, df = 2)
#> X-squared = 9299.2, df = 2, p-value < 2.2e-16
#> 
```
