# Cramer-Von Mises Test for Assessing Deviations From Normality

A goodness-of-fit test based on the integrated squared discrepancies
between the empirical and theoretical distribution functions. Similar to
the Anderson-Darling test, but without increased weighting of the
distribution tails.

## Usage

``` r
cramerVonMisesTest(x)
```

## Arguments

- x:

  a numeric vector of data values, the number of which must be at least
  8.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of the Cramer-von Mises statistic.

- `p.value`:

  the p-value of the test.

- `method`:

  a character string indicating the test performed.

- `data.name`:

  a character string giving the name of the data.

## Details

Performs the Cramer-von Mises test for the composite hypothesis of
normality, see e.g. Thode (2002, Sec. 5.1.3).

The Cramer-von Mises test is an EDF omnibus test for the composite
hypothesis of normality. The test statistic is \$\$W = \frac{1}{12 n} +
\sum\_{i=1}^{n} \left(p\_{(i)} - \frac{2i - 1}{2n}\right)^2,\$\$ where
\\p\_{(i)} = \Phi(\[x\_{(i)} - \bar{x}\]/s)\\. Here, \\\Phi\\ is the
cumulative distribution function of the standard normal distribution,
and \\\bar{x}\\ and \\s\\ are mean and standard deviation of the data
values. The p-value is computed from the modified statistic \\Z = W (1 +
0.5/n)\\ according to Table 4.9 in Stephens (1986).

Missing values are silently removed.

## Note

Based on code by Juergen Gross previously published in the nortest
package, adapted to conform to package standards.

## References

Stephens, M.A. (1986) Tests based on EDF statistics. In: D'Agostino,
R.B. and Stephens, M.A., eds.: *Goodness-of-Fit Techniques*. New York:
Marcel Dekker.

Thode Jr., H.C. (2002) *Testing for Normality*. New York: Marcel Dekker.

## See also

[`shapiro.test`](https://rdrr.io/r/stats/shapiro.test.html) for
performing the Shapiro-Wilk test for normality,
[`andersonDarlingTest`](https://andrisignorell.github.io/lumen/reference/andersonDarlingTest.md),
[`pharos::plotQQ()`](https://andrisignorell.github.io/pharos/reference/plotQQ.html)
for producing extended normal quantile-quantile plots

Other test.normality:
[`andersonDarlingTest()`](https://andrisignorell.github.io/lumen/reference/andersonDarlingTest.md),
[`jarqueBeraTest()`](https://andrisignorell.github.io/lumen/reference/jarqueBeraTest.md),
[`lillieTest()`](https://andrisignorell.github.io/lumen/reference/lillieTest.md),
[`pearsonTest()`](https://andrisignorell.github.io/lumen/reference/pearsonTest.md),
[`shapiroFranciaTest()`](https://andrisignorell.github.io/lumen/reference/shapiroFranciaTest.md)

## Examples

``` r
set.seed(1)
cramerVonMisesTest(rnorm(100, mean = 5, sd = 3))
#> 
#>  Cramer-von Mises normality test
#> 
#> data:  rnorm(100, mean = 5, sd = 3)
#> W = 0.026031, p-value = 0.8945
#> 
cramerVonMisesTest(runif(100, min = 2, max = 4))
#> 
#>  Cramer-von Mises normality test
#> 
#> data:  runif(100, min = 2, max = 4)
#> W = 0.27466, p-value = 0.0006342
#> 
```
