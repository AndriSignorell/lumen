# Von Neumann's Successive Difference Test for Detecting Serial Dependence in Continuous Sequences

A test for randomness or autocorrelation in a sequence, based on the
mean square successive difference relative to the sample variance,
closely related to the Durbin-Watson test.

## Usage

``` r
vonNeumannTest(
  x,
  alternative = c("two.sided", "less", "greater"),
  unbiased = TRUE
)
```

## Arguments

- x:

  a numeric vector containing the observations.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"greater"` or `"less"`.

- unbiased:

  logical. If `TRUE` (default), applies the finite-sample correction
  \\n/(n-1)\\ so that VN is an unbiased estimate of the population
  value.

## Value

A list with class `"htest"` containing:

- statistic:

  the normalized z-statistic.

- parameter:

  named vector with `n`, the sample size after removal of `NA`s.

- p.value:

  the p-value of the test.

- alternative:

  a character string describing the alternative hypothesis.

- data.name:

  a character string giving the name of the data.

- vn:

  the value of the VN statistic (not printed).

## Details

The test is based on the von Neumann ratio statistic.

The VN test statistic is in the unbiased case
\$\$VN=\frac{\sum\_{i=1}^{n-1}(x_i-x\_{i+1})^2 \cdot
n}{\sum\_{i=1}^{n}\left(x_i-\bar{x}\right)^2 \cdot (n-1)}\$\$

It is known that \\(VN-\mu)/\sigma\\ is asymptotically standard normal,
where \\\mu = 2n/(n-1)\\ and \\\sigma^2 = 4 n^2 (n-2) /
\[(n+1)(n-1)^3\]\\.

The VN test statistic is in the original (biased) case
\$\$VN=\frac{\sum\_{i=1}^{n-1}(x_i-x\_{i+1})^2}{
\sum\_{i=1}^{n}\left(x_i-\bar{x}\right)^2}\$\$

The test statistic \\(VN-2)/\sigma\\ is asymptotically standard normal,
where \\\sigma^2 = 4(n-2) / \[(n+1)(n-1)\]\\.

Missing values are silently removed.

## References

von Neumann, J. (1941) Distribution of the ratio of the mean square
successive difference to the variance. *Annals of Mathematical
Statistics* **12**, 367–395.

Young, L. C. (1941) On randomness in ordered sequences. *Annals of
Mathematical Statistics* **12**, 293–300.

Bartels, R. (1982) The Rank Version of von Neumann's Ratio Test for
Randomness. *Journal of the American Statistical Association*,
**77**(377), 40–46.

## See also

Other test.timeseries: [`adfTest()`](adfTest.md),
[`bartelsRankTest()`](BartelsRankTest.md), [`kpssTest()`](kpssTest.md),
[`runsTest()`](runsTest.md)

## Examples

``` r
set.seed(2)
vonNeumannTest(runif(20))
#> 
#>  Von Neumann Successive Difference Test
#> 
#> data:  runif(20)
#> z = 0.73709, n = 20, p-value = 0.4611
#> alternative hypothesis: two.sided
#> 

# trend: small VN expected
vonNeumannTest(cumsum(rnorm(30)), alternative = "less")
#> 
#>  Von Neumann Successive Difference Test
#> 
#> data:  cumsum(rnorm(30))
#> z = -4.8849, n = 30, p-value = 5.174e-07
#> alternative hypothesis: less
#> 
```
