# Distributions of Maxima and Minima

Density, distribution, quantile, and random generation functions for the
maximum or minimum of a given number of independent and identically
distributed random variables from any specified distribution, derived
analytically from the underlying distribution function.

## Usage

``` r
dextreme(x, dFun, pFun, ..., distn, mlen = 1, largest = TRUE, log = FALSE)

pextreme(q, pFun, ..., distn, mlen = 1, largest = TRUE, lower.tail = TRUE)

qextreme(p, qFun, ..., distn, mlen = 1, largest = TRUE, lower.tail = TRUE)

rextreme(n, qFun, ..., distn, mlen = 1, largest = TRUE)
```

## Arguments

- x, q:

  vector of quantiles.

- dFun, pFun, qFun:

  density, distribution and quantile function of the specified
  distribution. The density function must have a `log` argument (a
  simple wrapper can always be constructed to achieve this).

- ...:

  parameters of the specified distribution.

- distn:

  a character string, optionally given as an alternative to `dFun`,
  `pFun` and `qFun` such that the density, distribution and quantile
  functions are formed upon the addition of the prefixes `d`, `p` and
  `q` respectively.

- mlen:

  the number of independent variables.

- largest:

  logical; if `TRUE` (default) use maxima, otherwise minima.

- log:

  logical; if `TRUE`, the log density is returned.

- lower.tail:

  logical; if `TRUE` (default) probabilities are `P[X <= x]`, otherwise
  P`[X > x]`.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dextreme()` gives the density function, `pextreme()` gives the
distribution function and `qextreme()` gives the quantile function of
the maximum/minimum of `mlen` independent variables from a specified
distribution. `rextreme()` generates random deviates.

## Details

Density function, distribution function, quantile function and random
generation for the maximum/minimum of a given number of independent
variables from a specified distribution.

## Note

Based on code by Alec Stephenson previously published in the evd
package, adapted to conform to package standards.

## See also

[distributions-overview](https://andrisignorell.github.io/lumen/reference/distributions-overview.md)

## Examples

``` r

dextreme(2:4, dnorm, pnorm, mean = 0.5, sd = 1.2, mlen = 5)
#> [1] 0.48689660 0.17602941 0.02346192
dextreme(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5)
#> [1] 0.48689660 0.17602941 0.02346192
dextreme(2:4, distn = "exp", mlen = 2, largest = FALSE)
#> [1] 0.0366312778 0.0049575044 0.0006709253
pextreme(2:4, distn = "exp", rate = 1.2, mlen = 2)
#> [1] 0.8267938 0.9460991 0.9836082
qextreme(seq(0.9, 0.6, -0.1), distn = "exp", rate = 1.2, mlen = 2)
#> [1] 2.474783 1.873629 1.509935 1.241553
rextreme(5, qgamma, shape = 1, mlen = 10)
#> [1] 1.089317 2.953382 1.127408 2.007502 2.624579

p <- (1:9)/10
pexp(qextreme(p, distn = "exp", rate = 1.2, mlen = 1), rate = 1.2)
#> [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9

```
