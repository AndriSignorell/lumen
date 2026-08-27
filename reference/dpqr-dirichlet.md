# Dirichlet Distribution

Density, random generation, and basic utilities for the Dirichlet
distribution.

## Usage

``` r
ddirichlet(x, concentration, log = FALSE)

pdirichlet(q, concentration, R = 100000)

rdirichlet(n, concentration)

qdirichlet()
```

## Arguments

- x:

  numeric vector or matrix (rows sum to 1).

- concentration:

  numeric vector of concentration parameters (\> 0).

- log:

  logical; return log-density if TRUE.

- q:

  numeric vector of quantiles.

- R:

  number of Monte Carlo simulations used to approximate the CDF.

- n:

  number of samples.

## Value

`ddirichlet()` gives a numeric vector of densities (one per row of `x`),
`pdirichlet()` gives an approximate probability, and `rdirichlet()`
generates a matrix with `n` rows of random deviates. `qdirichlet()` only
signals an error, as no unique multivariate quantile function exists.

## Details

The Dirichlet distribution is a multivariate generalization of the Beta
distribution defined on the simplex: \$\$\sum\_{i=1}^k x_i = 1, \quad
x_i \ge 0\$\$

## See also

[distributions-overview](distributions-overview.md)

## Examples

``` r
ddirichlet(c(0.2, 0.3, 0.5), c(1,1,1))
#> [1] 2
pdirichlet(c(0.2, 0.3, 0.5), c(1,1,1))
#> [1] 0
rdirichlet(5, c(1,1,1))
#>           [,1]       [,2]       [,3]
#> [1,] 0.3098634 0.04969046 0.64044617
#> [2,] 0.3710099 0.27675250 0.35223762
#> [3,] 0.6817090 0.01999573 0.29829523
#> [4,] 0.3373393 0.01923777 0.64342291
#> [5,] 0.1445963 0.85109235 0.00431132
```
