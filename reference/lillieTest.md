# Lilliefors (Kolmogorov-Smirnov) Test for Assessing Normality With Estimated Parameters

A goodness-of-fit test for normality based on the Kolmogorov-Smirnov
statistic, adapted for the case where the population mean and variance
are unknown and must be estimated from the sample (Lilliefors test).

## Usage

``` r
lillieTest(x)
```

## Arguments

- x:

  a numeric vector of data values, the number of which must be at least
  5.

## Value

A list with class `"htest"` containing the following components:

- statistic:

  the value of the Lilliefors (Kolmogorov-Smirnov) statistic.

- p.value:

  the p-value of the test.

- method:

  a character string indicating the test performed.

- data.name:

  a character string giving the name of the data.

## Details

Performs the Lilliefors (Kolmogorov-Smirnov) test for the composite
hypothesis of normality, see e.g. Thode (2002, Sec. 5.1.1).

The Lilliefors (Kolmogorov-Smirnov) test is an EDF omnibus test for the
composite hypothesis of normality. The test statistic is the maximal
absolute difference between empirical and hypothetical cumulative
distribution function. It may be computed as \\D = \max\\D^{+},
D^{-}\\\\ with \$\$D^{+} = \max\_{i=1,\ldots,n}\\i/n - p\_{(i)}\\, \quad
D^{-} = \max\_{i=1,\ldots,n}\\p\_{(i)} - (i-1)/n\\,\$\$ where \\p\_{(i)}
= \Phi(\[x\_{(i)} - \bar{x}\]/s)\\. Here, \\\Phi\\ is the cumulative
distribution function of the standard normal distribution, and
\\\bar{x}\\ and \\s\\ are mean and standard deviation of the data
values. The p-value is computed from the Dallal-Wilkinson (1986)
formula, which is claimed to be only reliable when the p-value is
smaller than 0.1. If the Dallal-Wilkinson p-value turns out to be
greater than 0.1, then the p-value is computed from the distribution of
the modified statistic \\Z = D(\sqrt{n} - 0.01 + 0.85/\sqrt{n})\\, see
Stephens (1974), the actual p-value formula being obtained by a
simulation and approximation process.

Missing values are silently removed.

## Note

The Lilliefors (Kolmogorov-Smirnov) test is the most well-known EDF
omnibus test for normality. Compared to the Anderson-Darling test and
the Cramer-von Mises test it is known to perform worse. Although the
test statistic obtained from `lillieTest(x)` is the same as that
obtained from `ks.test(x, "pnorm", mean(x), sd(x))`, it is not correct
to use the p-value from the latter for the composite hypothesis of
normality (mean and variance unknown), since the distribution of the
test statistic is different when the parameters are estimated.

Based on code by Juergen Gross previously published in the nortest
package, adapted to conform to package standards.

## References

Dallal, G.E. and Wilkinson, L. (1986) An analytic approximation to the
distribution of Lilliefors' test for normality. *The American
Statistician*, 40, 294-296.

Stephens, M.A. (1974) EDF statistics for goodness of fit and some
comparisons. *Journal of the American Statistical Association*, 69,
730-737.

Thode Jr., H.C. (2002) *Testing for Normality*, New York: Marcel Dekker.

## See also

[`shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html)

Other test.normality:
[`andersonDarlingTest()`](https://andrisignorell.github.io/lumen/reference/andersonDarlingTest.md),
[`cramerVonMisesTest()`](https://andrisignorell.github.io/lumen/reference/cramerVonMisesTest.md),
[`jarqueBeraTest()`](https://andrisignorell.github.io/lumen/reference/jarqueBeraTest.md),
[`pearsonTest()`](https://andrisignorell.github.io/lumen/reference/pearsonTest.md),
[`shapiroFranciaTest()`](https://andrisignorell.github.io/lumen/reference/shapiroFranciaTest.md)

## Examples

``` r
set.seed(1)
lillieTest(rnorm(100, mean = 5, sd = 3))
#> 
#>  Lilliefors (Kolmogorov-Smirnov) normality test
#> 
#> data:  rnorm(100, mean = 5, sd = 3)
#> D = 0.047014, p-value = 0.8479
#> 
lillieTest(runif(100, min = 2, max = 4))
#> 
#>  Lilliefors (Kolmogorov-Smirnov) normality test
#> 
#> data:  runif(100, min = 2, max = 4)
#> D = 0.10079, p-value = 0.014
#> 
```
