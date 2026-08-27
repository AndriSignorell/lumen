# Maxima of Two Gumbel Distributions

The extended Gumbel distribution models the maximum of two independent
Gumbel-distributed random variables with potentially different location
and scale parameters. It is parameterized by two pairs of location and
scale parameters, with the constraint that the first location parameter
does not exceed the second.

## Usage

``` r
dgumbelx(x, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1, log = FALSE)

pgumbelx(q, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1, lower.tail = TRUE)

qgumbelx(
  p,
  interval,
  loc1 = 0,
  scale1 = 1,
  loc2 = 0,
  scale2 = 1,
  lower.tail = TRUE,
  ...
)

rgumbelx(n, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- loc1, scale1, loc2, scale2:

  location and scale parameters of the two Gumbel distributions. The
  second location parameter must be greater than or equal to the first
  location parameter.

- log:

  logical; if `TRUE`, the log density is returned.

- lower.tail:

  logical; if `TRUE` (default), probabilities are `P[X <= x]`,
  otherwise, `P[X > x]`.

- p:

  vector of probabilities.

- interval:

  a length two vector containing the end-points of the interval to be
  searched for the quantiles, passed to the uniroot function.

- ...:

  other arguments passed to uniroot.

- n:

  number of observations.

## Value

`dgumbelx()` gives the density function, `pgumbelx()` gives the
distribution function, `qgumbelx()` gives the quantile function, and
`rgumbelx()` generates random deviates.

## Details

Density function, distribution function, quantile function and random
generation for the maxima of two Gumbel distributions, each with
different location and scale parameters.

## Note

Based on code by Alec Stephenson previously published in the evd
package, adapted to conform to package standards.

## See also

[distributions-overview](distributions-overview.md);
[`uniroot`](https://rdrr.io/r/stats/uniroot.html), which `qgumbelx()`
uses for root finding

## Examples

``` r

dgumbelx(2:4, 0, 1.1, 1, 0.5)
#> [1] 0.31056307 0.08836749 0.02808872
pgumbelx(2:4, 0, 1.1, 1, 0.5)
#> [1] 0.7425568 0.9196951 0.9715848
qgumbelx(seq(0.9, 0.6, -0.1), interval = c(0,10), 0, 1.2, 2, 0.5)
#> [1] 3.489999 2.983361 2.692015 2.478487
rgumbelx(6, 0, 1.1, 1, 0.5)
#> [1] 0.9823914 0.3906958 0.8207195 1.4664357 0.5397208 0.8659020
p <- (1:9)/10
pgumbelx(qgumbelx(p, interval = c(0,10), 0, 0.5, 1, 2), 0, 0.5, 1, 2)
#> [1] 0.1000003 0.1999999 0.2999955 0.3999998 0.5000004 0.6000001 0.7000031
#> [8] 0.8000001 0.9000000
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
```
