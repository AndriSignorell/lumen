# Maxima of Two Gumbel Distributions

The extended Gumbel distribution models the maximum of two independent
Gumbel-distributed random variables with potentially different location
and scale parameters. It is parameterized by two pairs of location and
scale parameters.

## Usage

``` r
dgumbelx(x, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1, log = FALSE)

pgumbelx(
  q,
  loc1 = 0,
  scale1 = 1,
  loc2 = 0,
  scale2 = 1,
  lower.tail = TRUE,
  log.p = FALSE
)

qgumbelx(
  p,
  loc1 = 0,
  scale1 = 1,
  loc2 = 0,
  scale2 = 1,
  lower.tail = TRUE,
  log.p = FALSE,
  interval = NULL,
  ...
)

rgumbelx(n, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- loc1, scale1, loc2, scale2:

  location and scale parameters of the two Gumbel distributions. The
  distribution is symmetric in the two margins, so their order is
  immaterial.

- log, log.p:

  logical; if `TRUE`, probabilities `p` are given as `log(p)` and the
  density is returned on the log scale.

- lower.tail:

  logical; if `TRUE` (default), probabilities are `P[X <= x]`,
  otherwise, `P[X > x]`.

- p:

  vector of probabilities.

- interval:

  a length two vector containing the end-points of the interval to be
  searched for the quantiles, passed to
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html). By default a
  bracketing interval is derived from the quantiles of the two Gumbel
  margins.

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
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html), which `qgumbelx()`
uses for root finding

## Examples

``` r

dgumbelx(2:4, 0, 1.1, 1, 0.5)
#> [1] 0.31056307 0.08836749 0.02808872
pgumbelx(2:4, 0, 1.1, 1, 0.5)
#> [1] 0.7425568 0.9196951 0.9715848
qgumbelx(seq(0.9, 0.6, -0.1), 0, 1.2, 2, 0.5)
#> [1] 3.489993 2.983368 2.692006 2.478481
rgumbelx(6, 0, 1.1, 1, 0.5)
#> [1] 0.8957238 1.3756205 2.2961160 0.8621185 2.7607791 1.4477668
p <- (1:9)/10
pgumbelx(qgumbelx(p, 0, 0.5, 1, 2), 0, 0.5, 1, 2)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
```
