# Mean and Variance of Discrete Distributions

Computes the mean and variance of common discrete distributions given
their parameters.

## Usage

``` r
mbinom(size, prob)

mpois(lambda)

mgeom(prob)

mnbinom(size, prob)

mhyper(m, n, k)

mbenford(ndigits = 1)
```

## Arguments

- size:

  number of trials (binomial, negative binomial).

- prob:

  probability of success on each trial (binomial, geometric, negative
  binomial).

- lambda:

  mean (Poisson).

- m:

  number of white balls in the urn (hypergeometric).

- n:

  number of black balls in the urn (hypergeometric).

- k:

  number of balls drawn (hypergeometric).

- ndigits:

  number of leading digits for Benford's distribution, either `1`
  (default, support {1,...,9}) or `2` (support {10,...,99}).

## Value

A named numeric vector with elements `mean` and `variance`.

## Details

|  |  |  |
|----|----|----|
| **Distribution` `** | **Mean` `** | **Variance** |
| Binomial | \\np\\ | \\np(1-p)\\ |
| Poisson | \\\lambda\\ | \\\lambda\\ |
| Geometric | \\\frac{1-p}{p}\\ | \\\frac{1-p}{p^2}\\ |
| Negative binomial | \\\frac{r(1-p)}{p}\\ | \\\frac{r(1-p)}{p^2}\\ |
| Hypergeometric | \\\frac{km}{N}\\ | \\\frac{km}{N}\frac{n}{N}\frac{N-k}{N-1}\\ |
| Benford | \\\sum_d d\log\_{10}\left(1 + \frac{1}{d}\right)\\ | \\\sum_d d^2\log\_{10}\left(1 + \frac{1}{d}\right) - \mu^2\\ |

For the binomial distribution, \\n\\ = `size`; for the negative binomial
distribution, \\r\\ = `size`; and for the hypergeometric distribution,
\\N = m + n\\. For Benford's distribution, the sum runs over \\d \in
\\1,\ldots,9\\\\ for `ndigits = 1` and \\d \in \\10,\ldots,99\\\\ for
`ndigits = 2`. As there is no closed-form solution, the moments are
computed numerically.

## References

Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011) *Statistical
Distributions*. Fourth Edition. Wiley.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) *Continuous
Univariate Distributions*, Vol. 2. Wiley.

## See also

[`Binomial`](https://rdrr.io/r/stats/Binomial.html),
[`Poisson`](https://rdrr.io/r/stats/Poisson.html),
[`Geometric`](https://rdrr.io/r/stats/Geometric.html),
[`NegBinomial`](https://rdrr.io/r/stats/NegBinomial.html),
[`Hypergeometric`](https://rdrr.io/r/stats/Hypergeometric.html),
[distributions-overview](distributions-overview.md)

## Examples

``` r
mbinom(size = 10, prob = 0.5)
#>     mean variance 
#>      5.0      2.5 
mpois(lambda = 3)
#>     mean variance 
#>        3        3 
mgeom(prob = 0.3)
#>     mean variance 
#> 2.333333 7.777778 
mnbinom(size = 5, prob = 0.3)
#>     mean variance 
#> 11.66667 38.88889 
mhyper(m = 10, n = 5, k = 4)
#>      mean  variance 
#> 2.6666667 0.6984127 
mbenford(ndigits = 1)
#>     mean variance 
#> 3.440237 6.056513 
mbenford(ndigits = 2)
#>      mean  variance 
#>  38.58976 621.83174 
```
