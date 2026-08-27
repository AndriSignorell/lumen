# Mean and Variance of Extreme Value Distributions

Computes the mean and variance of common extreme value distributions
given their parameters.

## Usage

``` r
mgumbel(loc = 0, scale = 1)

mfrechet(loc = 0, scale = 1, shape = 1)

mrweibull(loc = 0, scale = 1, shape = 1)

mgev(loc = 0, scale = 1, shape = 0)

mgpd(loc = 0, scale = 1, shape = 0)

mgompertz(shape, rate = 1)
```

## Arguments

- loc:

  location parameter.

- scale:

  scale parameter.

- shape, rate:

  vector of shape and rate parameters.

## Value

A named numeric vector with elements `mean` and `variance`. Returns `NA`
where moments do not exist.

## Details

|  |  |  |
|----|----|----|
| **Distribution** | **Mean** | **Variance** |
| Gumbel | \\a + b\gamma\\ | \\\frac{\pi^2}{6}b^2\\ |
| Fréchet | \\a + b\Gamma(1 - 1/s) \quad (s \> 1)\\ | \\b^2\left\[\Gamma(1 - 2/s) - \Gamma(1 - 1/s)^2\right\] \quad (s \> 2)\\ |
| Reverse Weibull ` ` | \\a - b\Gamma(1 + 1/s)\\ | \\b^2\left\[\Gamma(1 + 2/s) - \Gamma(1 + 1/s)^2\right\]\\ |
| GEV | \\a + b\frac{\Gamma(1-s)-1}{s} \quad (s \ne 0,\\ s \< 1)\\ ` ` | \\b^2\frac{\Gamma(1-2s)-\Gamma(1-s)^2}{s^2} \quad (s \ne 0,\\ s \< 1/2)\\ |
| GPD | \\a + \frac{b}{1-s} \quad (s \< 1)\\ | \\\frac{b^2}{(1-s)^2(1-2s)} \quad (s \< 1/2)\\ |
| Gompertz | numerical integration for \\\alpha \> 0\\ | dito |
|  | \\1/\beta\\ for \\\alpha = 0\\; | \\1/\beta^2\\ for \\\alpha = 0\\; |
|  | `NA` for \\\alpha \< 0\\ | dito |

For the first five distributions, \\a\\ = `loc`, \\b\\ = `scale`, and
\\s\\ = `shape`. For the GEV with \\s = 0\\, the Gumbel moments apply.
Furthermore, \\\gamma \approx 0.5772\\ is the Euler-Mascheroni constant.
For the Gompertz distribution, \\\alpha\\ = `shape` and \\\beta\\ =
`rate`; moments for \\\alpha \> 0\\ are computed numerically by
integration.

## References

Coles, S. (2001) *An Introduction to Statistical Modeling of Extreme
Values*. Springer.

Kotz, S. and Nadarajah, S. (2000) *Extreme Value Distributions*.
Imperial College Press.

## See also

[`dgumbel`](https://andrisignorell.github.io/lumen/reference/dpqr-gumbel.md),
[`dfrechet`](https://andrisignorell.github.io/lumen/reference/dpqr-frechet.md),
[`drweibull`](https://andrisignorell.github.io/lumen/reference/dpqr-rweibull.md),
[`dgev`](https://andrisignorell.github.io/lumen/reference/dpqr-gev.md),
[`dgpd`](https://andrisignorell.github.io/lumen/reference/dpqr-gpd.md),
[distributions-overview](https://andrisignorell.github.io/lumen/reference/distributions-overview.md)

## Examples

``` r
mgumbel(loc = 0, scale = 1)
#>      mean  variance 
#> 0.5772157 1.6449341 
mfrechet(loc = 0, scale = 1, shape = 3)
#>      mean  variance 
#> 1.3541179 0.8453031 
mrweibull(loc = 0, scale = 1, shape = 2)
#>       mean   variance 
#> -0.8862269  0.2146018 
mgev(loc = 0, scale = 1, shape = 0)
#>      mean  variance 
#> 0.5772157 1.6449341 
mgev(loc = 0, scale = 1, shape = 0.3)
#>      mean  variance 
#> 0.9935178 5.9245766 
mgev(loc = 0, scale = 1, shape = -0.3)
#>      mean  variance 
#> 0.3417643 0.9784633 
mgpd(loc = 0, scale = 1, shape = 0.3)
#>     mean variance 
#> 1.428571 5.102041 
```
