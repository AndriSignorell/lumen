# Reverse Weibull Distribution

The Reverse Weibull distribution, also known as the Type III extreme
value distribution, is the distribution of the negative of a
Weibull-distributed random variable. It has an upper bound and a
left-skewed density, and is parameterized by location, scale, and shape.

## Usage

``` r
drevweibull(x, loc = 0, scale = 1, shape = 1, log = FALSE)

prevweibull(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)

qrevweibull(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)

rrevweibull(n, loc = 0, scale = 1, shape = 1)

dnweibull(x, loc = 0, scale = 1, shape = 1, log = FALSE)

pnweibull(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)

qnweibull(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)

rnweibull(n, loc = 0, scale = 1, shape = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- loc, scale, shape:

  location, scale and shape parameters (can be given as vectors).

- log, log.p:

  logical; if `TRUE`, probabilities `p` are given as `log(p)` and the
  density is returned on the log scale.

- lower.tail:

  logical; if `TRUE` (default), probabilities are `P[X <= x]`,
  otherwise, P`[X > x]`.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`drevweibull()` and `dnweibull()` give the density function,
`prevweibull()` and `pnweibull()` give the distribution function,
`qrevweibull()` and `qnweibull()` give the quantile function,
`rrevweibull()` and `rnweibull()` generate random deviates.

## Details

Density function, distribution function, quantile function and random
generation for the reverse (sometimes called negative) Weibull
distribution with location, scale and shape parameters.

The reverse Weibull distribution function with parameters \\\`loc\` =
a\\, \\\`scale\` = b\\ and \\\`shape\` = s\\ is \$\$G(z) =
\exp\left\\-\left\[-\left(\frac{z-a}{b}\right)\right\]^s\right\\\$\$ for
\\z \< a\\ and one otherwise, where \\b \> 0\\ and \\s \> 0\\.

**Note:** Within extreme value theory the reverse Weibull distibution
(also known as the negative Weibull distribution) is often referred to
as the Weibull distribution. We make a distinction to avoid confusion
with the three-parameter distribution used in survival analysis, which
is related by a change of sign to the distribution given above.

## Note

Based on code by Alec Stephenson previously published in the evd
package, adapted to conform to package standards.

## See also

[distributions-overview](distributions-overview.md)

## Examples

``` r

drevweibull(-5:-3, -1, 0.5, 0.8)
#> [1] 0.005386194 0.016885315 0.058502349
prevweibull(-5:-3, -1, 0.5, 0.8)
#> [1] 0.005102464 0.015101477 0.048246445
qrevweibull(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#> [1] 1.969986 1.923317 1.862180 1.784071
rrevweibull(6, -1, 0.5, 0.8)
#> [1] -1.404639 -1.058543 -1.486602 -1.973033 -1.173827 -1.013722
p <- (1:9)/10
prevweibull(qrevweibull(p, -1, 2, 0.8), -1, 2, 0.8)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
```
