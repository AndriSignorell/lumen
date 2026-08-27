# Shapiro-Francia Test for Assessing Normality From Normal Scores

A goodness-of-fit test for normality based on the correlation between
the ordered sample values and the corresponding expected normal order
statistics, particularly suited for larger sample sizes than the
Shapiro-Wilk test.

## Usage

``` r
shapiroFranciaTest(x)
```

## Arguments

- x:

  a numeric vector of data values, the number of which must be between 5
  and 5000. Missing values are allowed.

## Value

A list with class “htest” containing the following components:

- statistic:

  the value of the Shapiro-Francia statistic.

- p.value :

  the p-value for the test.

- method:

  the character string “Shapiro-Francia normality test”.

- data.name:

  a character string giving the name(s) of the data.

## Details

Performs the Shapiro-Francia test for the composite hypothesis of
normality, see e.g. Thode (2002, Sec. 2.3.2).

The test statistic of the Shapiro-Francia test is simply the squared
correlation between the ordered sample values and the (approximated)
expected ordered quantiles from the standard normal distribution. The
p-value is computed from the formula given by Royston (1993).

## Note

The Shapiro-Francia test is known to perform well, see also the comments
by Royston (1993). The expected ordered quantiles from the standard
normal distribution are approximated by `qnorm(ppoints(x, a = 3/8))`,
being slightly different from the approximation
`qnorm(ppoints(x, a = 1/2))` used for the normal quantile-quantile plot
by [`qqnorm`](https://rdrr.io/r/stats/qqnorm.html) for sample sizes
greater than 10.

Based on code by Juergen Gross, adapted to conform to package standards.

## References

Royston, P. (1993): A pocket-calculator algorithm for the
Shapiro-Francia test for non-normality: an application to medicine.
Statistics in Medicine, 12, 181–184.

Thode Jr., H.C. (2002): Testing for Normality. Marcel Dekker, New York.

## See also

[`shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html) for
performing the Shapiro-Wilk test for normality

Other test.normality:
[`andersonDarlingTest()`](https://andrisignorell.github.io/lumen/reference/andersonDarlingTest.md),
[`cramerVonMisesTest()`](https://andrisignorell.github.io/lumen/reference/cramerVonMisesTest.md),
[`jarqueBeraTest()`](https://andrisignorell.github.io/lumen/reference/jarqueBeraTest.md),
[`lillieTest()`](https://andrisignorell.github.io/lumen/reference/lillieTest.md),
[`pearsonTest()`](https://andrisignorell.github.io/lumen/reference/pearsonTest.md)

## Examples

``` r

shapiroFranciaTest(rnorm(100, mean = 5, sd = 3))
#> 
#>  Shapiro-Francia normality test
#> 
#> data:  rnorm(100, mean = 5, sd = 3)
#> W = 0.98753, p-value = 0.4013
#> 
shapiroFranciaTest(runif(100, min = 2, max = 4))
#> 
#>  Shapiro-Francia normality test
#> 
#> data:  runif(100, min = 2, max = 4)
#> W = 0.96044, p-value = 0.005673
#> 
```
