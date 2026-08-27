# Gumbel Distribution

The Gumbel distribution, also known as the Type I extreme value
distribution, is a continuous distribution used to model the maximum (or
minimum) of a number of samples of various distributions. It is
parameterized by location and scale and arises naturally in extreme
value theory.

## Usage

``` r
dgumbel(x, loc = 0, scale = 1, log = FALSE)

pgumbel(q, loc = 0, scale = 1, lower.tail = TRUE)

qgumbel(p, loc = 0, scale = 1, lower.tail = TRUE)

rgumbel(n, loc = 0, scale = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- loc, scale:

  location and scale parameters (can be given as vectors).

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

`dgumbel()` gives the density function, `pgumbel()` gives the
distribution function, `qgumbel()` gives the quantile function, and
`rgumbel()` generates random deviates.

## Details

Density function, distribution function, quantile function and random
generation for the Gumbel distribution with location and scale
parameters.

The Gumbel distribution function with parameters \\\code{loc} = a\\ and
\\\code{scale} = b\\ is \$\$G(z) =
\exp\left\\-\exp\left\[-\left(\frac{z-a}{b}\right)\right\]\right\\\$\$
for all real \\z\\, where \\b \> 0\\.

## Note

Based on code by Alec Stephenson previously published in the evd
package, adapted to conform to package standards.

## See also

[distributions-overview](https://andrisignorell.github.io/lumen/reference/distributions-overview.md)

## Examples

``` r

dgumbel(-1:2, -1, 0.5)
#> [1] 0.735758882 0.236409903 0.035966459 0.004945231
pgumbel(-1:2, -1, 0.5)
#> [1] 0.3678794 0.8734230 0.9818511 0.9975243
qgumbel(seq(0.9, 0.6, -0.1), 2, 0.5)
#> [1] 3.125184 2.749970 2.515465 2.335863
rgumbel(6, -1, 0.5)
#> [1] -0.33366250 -0.82888510 -1.00417239 -1.30921410  0.47954005  0.07979104
p <- (1:9)/10
pgumbel(qgumbel(p, -1, 2), -1, 2)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
```
