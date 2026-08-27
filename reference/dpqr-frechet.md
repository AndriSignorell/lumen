# Fréchet Distribution

The Fréchet distribution, also known as the Type II extreme value
distribution, is a continuous probability distribution for the maximum
of a sequence of independent random variables. It has a lower bound and
a heavy right tail, and is parameterized by location, scale, and shape.

## Usage

``` r
dfrechet(x, loc = 0, scale = 1, shape = 1, log = FALSE)

pfrechet(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)

qfrechet(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)

rfrechet(n, loc = 0, scale = 1, shape = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- loc, scale, shape:

  location, scale and shape parameters (can be given as vectors).

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

`dfrechet()` gives the density function, `pfrechet()` gives the
distribution function, `qfrechet()` gives the quantile function, and
`rfrechet()` generates random deviates.

## Details

Density function, distribution function, quantile function and random
generation for the Frechet distribution with location, scale and shape
parameters.

The Frechet distribution function with parameters \\\code{loc} = a\\,
\\\code{scale} = b\\ and \\\code{shape} = s\\ is \$\$G(z) =
\exp\left\\-\left(\frac{z-a}{b}\right)^{-s}\right\\\$\$ for \\z \> a\\
and zero otherwise, where \\b \> 0\\ and \\s \> 0\\.

## Note

Based on code by Alec Stephenson previously published in the evd
package, adapted to conform to package standards.

## See also

[distributions-overview](https://andrisignorell.github.io/lumen/reference/distributions-overview.md)

## Examples

``` r

dfrechet(2:4, 1, 0.5, 0.8)
#> [1] 0.25871959 0.09487423 0.05010381
pfrechet(2:4, 1, 0.5, 0.8)
#> [1] 0.5630712 0.7190122 0.7878127
qfrechet(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#> [1] 10.329571  5.260165  3.813966  3.157788
rfrechet(6, 1, 0.5, 0.8)
#> [1] 1.147215 1.265835 1.955350 1.197307 1.533234 8.566491
p <- (1:9)/10
pfrechet(qfrechet(p, 1, 2, 0.8), 1, 2, 0.8)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
```
