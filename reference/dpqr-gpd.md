# Generalized Pareto Distribution

The Generalized Pareto Distribution (GPD) is a two-parameter family of
distributions used to model exceedances over a high threshold, commonly
applied in extreme value theory as the limiting distribution of
threshold excesses.

## Usage

``` r
dgpd(x, loc = 0, scale = 1, shape = 0, log = FALSE)

pgpd(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)

qgpd(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)

rgpd(n, loc = 0, scale = 1, shape = 0)
```

## Arguments

- x, q:

  vector of quantiles.

- loc, scale, shape:

  location, scale and shape parameters; the `shape` argument cannot be a
  vector (must have length one).

- log:

  logical; if `TRUE`, the log density is returned.

- lower.tail:

  logical; if `TRUE` (default), probabilities are `P[X <= x]`,
  otherwise, P`[X > x]`.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dgpd()` gives the density function, `pgpd()` gives the distribution
function, `qgpd()` gives the quantile function, and `rgpd()` generates
random deviates.

## Details

Density function, distribution function, quantile function and random
generation for the generalized Pareto distribution (GPD) with location,
scale and shape parameters.

The generalized Pareto distribution function (Pickands, 1975) with
parameters \\\code{loc} = a\\, \\\code{scale} = b\\ and \\\code{shape} =
s\\ is \$\$G(z) = 1 - \\1+s(z-a)/b\\^{-1/s}\$\$ for \\1+s(z-a)/b \> 0\\
and \\z \> a\\, where \\b \> 0\\. If \\s = 0\\ the distribution is
defined by continuity.

## Note

Based on code by Alec Stephenson previously published in the evd
package, adapted to conform to package standards.

## References

Pickands, J. (1975) Statistical inference using extreme order
statistics. *Annals of Statistics*, **3**, 119–131.

## See also

[distributions-overview](distributions-overview.md); `evd::fpot()` for
fitting peaks-over-threshold models

## Examples

``` r

dgpd(2:4, 1, 0.5, 0.8)
#> [1] 0.23299144 0.07919889 0.03831043
pgpd(2:4, 1, 0.5, 0.8)
#> [1] 0.6971111 0.8336823 0.8888998
qgpd(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#> [1] 5.318483 3.639936 3.012506 2.675864
rgpd(6, 1, 0.5, 0.8)
#> [1] 1.134926 1.031631 1.094098 3.241592 1.327359 1.104716
p <- (1:9)/10
pgpd(qgpd(p, 1, 2, 0.8), 1, 2, 0.8)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
```
