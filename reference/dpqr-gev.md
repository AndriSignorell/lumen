# Generalized Extreme Value Distribution

The Generalized Extreme Value (GEV) distribution unifies the three
extreme value distributions — Gumbel (Type I), Fréchet (Type II), and
Reverse Weibull (Type III) — into a single family, parameterized by
location, scale, and a shape parameter that determines which type
applies.

## Usage

``` r
dgev(x, loc = 0, scale = 1, shape = 0, log = FALSE)

pgev(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)

qgev(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)

rgev(n, loc = 0, scale = 1, shape = 0)
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

`dgev()` gives the density function, `pgev()` gives the distribution
function, `qgev()` gives the quantile function, and `rgev()` generates
random deviates.

## Details

Density function, distribution function, quantile function and random
generation for the generalized extreme value (GEV) distribution with
location, scale and shape parameters.

The GEV distribution function with parameters \\\code{loc} = a\\,
\\\code{scale} = b\\ and \\\code{shape} = s\\ is \$\$G(z) =
\exp\left\[-\\1+s(z-a)/b\\^{-1/s}\right\]\$\$ for \\1+s(z-a)/b \> 0\\,
where \\b \> 0\\. If \\s = 0\\ the distribution is defined by
continuity. If \\1+s(z-a)/b \leq 0\\, the value \\z\\ is either greater
than the upper end point (if \\s \< 0\\), or less than the lower end
point (if \\s \> 0\\).

The parametric form of the GEV encompasses that of the Gumbel, Frechet
and reverse Weibull distributions, which are obtained for \\s = 0\\, \\s
\> 0\\ and \\s \< 0\\ respectively. It was first introduced by Jenkinson
(1955).

## Note

Based on code by Alec Stephenson previously published in the evd
package, adapted to conform to package standards.

## References

Jenkinson, A. F. (1955) The frequency distribution of the annual maximum
(or minimum) of meteorological elements. *Quart. J. R. Met. Soc.*,
**81**, 158–171.

## See also

[distributions-overview](https://andrisignorell.github.io/lumen/reference/distributions-overview.md);
`evd::fgev()` for fitting the GEV to data

## Examples

``` r

dgev(2:4, 1, 0.5, 0.8)
#> [1] 0.17210639 0.06706381 0.03428205
pgev(2:4, 1, 0.5, 0.8)
#> [1] 0.7386812 0.8467772 0.8948490
qgev(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#> [1] 5.157141 3.449973 2.800811 2.444700
rgev(6, 1, 0.5, 0.8)
#> [1] 0.9494715 0.8191944 3.4902610 0.9856113 0.7009025 1.2436763
p <- (1:9)/10
pgev(qgev(p, 1, 2, 0.8), 1, 2, 0.8)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
```
