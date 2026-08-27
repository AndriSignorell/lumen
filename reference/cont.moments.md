# Mean and Variance of Continuous Distributions

Computes the mean and variance of common continuous distributions given
their parameters.

## Usage

``` r
mnorm(mean, sd)

mexp(rate)

mgamma(shape, rate)

mlnorm(meanlog, sdlog)

mbeta(shape1, shape2)

mchisq(df)

mt(df)

mf(df1, df2)

mtri(min = 0, max = 1, mode = 0.5)
```

## Arguments

- mean:

  mean of the normal distribution.

- sd:

  standard deviation of the normal distribution.

- rate:

  rate parameter (1/mean) of the exponential distribution.

- shape:

  shape parameter of the gamma distribution.

- meanlog, sdlog:

  mean and standard deviation on the log scale (log-normal
  distribution).

- shape1, shape2:

  shape parameters of the beta distribution (\\\alpha\\ and \\\beta\\).

- df:

  degrees of freedom (chi-squared and t-distribution).

- df1, df2:

  numerator and denominator degrees of freedom (F-distribution).

- min:

  lower limit of the triangular distribution.

- max:

  upper limit of the triangular distribution.

- mode:

  mode of the triangular distribution.

## Value

A named numeric vector with elements `mean` and `variance`. Returns `NA`
where moments do not exist.

## Details

|  |  |  |
|----|----|----|
| **Distribution ` `** | **Mean ` `** | **Variance** |
| Normal | \\\mu\\ | \\\sigma^2\\ |
| Exponential | \\\frac{1}{\lambda}\\ | \\\frac{1}{\lambda^2}\\ |
| Gamma | \\\frac{\alpha}{\beta}\\ | \\\frac{\alpha}{\beta^2}\\ |
| Log-normal | \\\exp(\mu\_{log} + \frac{1}{2}\sigma\_{log}^2)\\ | \\(\exp(\sigma\_{log}^2) - 1) \exp(2\mu\_{log} + \sigma\_{log}^2)\\ |
| Beta | \\\frac{\alpha}{\alpha + \beta}\\ | \\\frac{\alpha\beta} {(\alpha+\beta)^2(\alpha+\beta+1)}\\ |
| Chi-squared | \\\nu\\ | \\2\nu\\ |
| t-distribution | \\0 \quad (\nu \> 1)\\ | \\\frac{\nu}{\nu-2} \quad (\nu \> 2)\\ |
| F-distribution | \\\frac{n_2}{n_2-2} \quad (n_2 \> 2)\\ | \\\frac{2n_2^2(n_1+n_2-2)} {n_1(n_2-2)^2(n_2-4)} \quad (n_2 \> 4)\\ |
| Triangular | \\\frac{a + b + c}{3}\\ | \\\frac{a^2 + b^2 + c^2 - ab - ac - bc}{18}\\ |
|  | where \\a\\ = `min`, \\b\\ = `max` and \\c\\ = `mode` |  |

## References

Casella, G. and Berger, R. L. (2002) *Statistical Inference*. Duxbury.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) *Continuous
Univariate Distributions*, Vol. 1. Wiley.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) *Continuous
Univariate Distributions*, Vol. 2. Wiley.

## See also

[`dnorm`](https://rdrr.io/r/stats/Normal.html),
[`dexp`](https://rdrr.io/r/stats/Exponential.html),
[`dgamma`](https://rdrr.io/r/stats/GammaDist.html),
[`dlnorm`](https://rdrr.io/r/stats/Lognormal.html),
[`dbeta`](https://rdrr.io/r/stats/Beta.html),
[`dchisq`](https://rdrr.io/r/stats/Chisquare.html),
[`dt`](https://rdrr.io/r/stats/TDist.html),
[`df`](https://rdrr.io/r/stats/Fdist.html),
[distributions-overview](distributions-overview.md)

## Examples

``` r
mnorm(mean = 0, sd = 1)
#>     mean variance 
#>        0        1 
mexp(rate = 0.5)
#>     mean variance 
#>        2        4 
mgamma(shape = 2, rate = 0.5)
#>     mean variance 
#>        4        8 
mlnorm(meanlog = 0, sdlog = 1)
#>     mean variance 
#> 1.648721 4.670774 
mbeta(shape1 = 2, shape2 = 3)
#>     mean variance 
#>     0.40     0.04 
mchisq(df = 4)
#>     mean variance 
#>        4        8 
mt(df = 5)
#>     mean variance 
#> 0.000000 1.666667 
mf(df1 = 5, df2 = 10)
#>     mean variance 
#> 1.250000 1.354167 
mtri(min = 0, max = 1, mode = 0.5)
#>       mean   variance 
#> 0.50000000 0.04166667 
mtri(min = 2, max = 10, mode = 4)
#>     mean variance 
#> 5.333333 2.888889 
```
