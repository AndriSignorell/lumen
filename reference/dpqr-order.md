# Distributions of Order Statistics

Density, distribution, and random generation functions for a selected
order statistic (the j-th largest or smallest value) from a sample of a
given size drawn from any specified distribution, derived analytically
using the beta distribution representation of order statistics.

## Usage

``` r
dorder(x, dFun, pFun, ..., distn, mlen = 1, j = 1, largest = TRUE, log = FALSE)

porder(q, pFun, ..., distn, mlen = 1, j = 1, largest = TRUE, lower.tail = TRUE)

rorder(n, qFun, ..., distn, mlen = 1, j = 1, largest = TRUE)
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

  a character string, optionally specified as an alternative to `dFun`,
  `pFun` and `qFun` such that the density, distribution and quantile
  functions are formed upon the addition of the prefixes `d`, `p` and
  `q` respectively.

- mlen:

  the number of independent variables.

- j:

  the order statistic, taken as the `j`th largest (default) or smallest
  of `mlen`, according to the value of `largest`.

- largest:

  logical; if `TRUE` (default) use the `j`th largest order statistic,
  otherwise use the `j`th smallest.

- log:

  logical; if `TRUE`, the log density is returned.

- lower.tail:

  logical; if `TRUE` (default) probabilities are `P[X <= x]`, otherwise
  P`[X > x]`.

- n:

  number of observations.

## Value

`dorder()` gives the density function and `porder()` gives the
distribution function of a selected order statistic from a sample of
size `mlen`, from a specified distribution. `rorder()` generates random
deviates. There is no quantile function for order statistics (`qorder()`
does not exist).

## Note

Based on code by Alec Stephenson previously published in the evd
package, adapted to conform to package standards.

## See also

[distributions-overview](distributions-overview.md)

## Examples

``` r

dorder(2:4, dnorm, pnorm, mean = 0.5, sd = 1.2, mlen = 5, j = 2)
#> [1] 0.2300687782 0.0133524232 0.0001663078
dorder(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5, j = 2)
#> [1] 0.2300687782 0.0133524232 0.0001663078
dorder(2:4, distn = "exp", mlen = 2, j = 2)
#> [1] 0.0366312778 0.0049575044 0.0006709253
porder(2:4, distn = "exp", rate = 1.2, mlen = 2, j = 2)
#> [1] 0.9917703 0.9992534 0.9999323
rorder(5, qgamma, shape = 1, mlen = 10, j = 2)
#> [1] 3.1585386 2.0492856 0.7806257 1.8950382 2.8174036
```
