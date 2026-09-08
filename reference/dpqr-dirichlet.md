# Dirichlet Distribution

Density, random generation, and basic utilities for the Dirichlet
distribution.

## Usage

``` r
ddirichlet(x, concentration, log = FALSE)

pdirichlet(q, concentration, R = 100000)

rdirichlet(n, concentration)

qdirichlet(p, concentration, ...)
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

  number of Monte Carlo simulations used to approximate the CDF. Must be
  at most `.Machine$integer.max`.

- n:

  number of samples.

- p:

  numeric vector of probabilities; accepted by `qdirichlet()` only to
  give a meaningful error, see below.

- ...:

  further arguments, accepted by `qdirichlet()` and ignored.

## Value

`ddirichlet()` gives a numeric vector of densities (one per row of `x`),
`pdirichlet()` gives an approximate probability, and `rdirichlet()`
generates a matrix with `n` rows of random deviates. `qdirichlet()` only
signals an error, as no unique multivariate quantile function exists; it
takes the same arguments as the others so that the message is reached
rather than an argument mismatch.

## Details

The Dirichlet distribution is a multivariate generalization of the Beta
distribution defined on the simplex: \$\$\sum\_{i=1}^k x_i = 1, \quad
x_i \ge 0\$\$

## Random number generation

`pdirichlet()` evaluates the CDF by simulation and is parallelised. Each
range of draws runs its own generator, seeded from R's stream, so that
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs the result;
the seeding does depend on how the work is split, so reproducing a value
also requires the same
[`RcppParallel::setThreadOptions()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).

## See also

[distributions-overview](distributions-overview.md)

## Examples

``` r
ddirichlet(c(0.2, 0.3, 0.5), c(1,1,1))
#> [1] 2
pdirichlet(c(0.2, 0.3, 0.5), c(1,1,1))
#> [1] 0
rdirichlet(5, c(1,1,1))
#>             [,1]       [,2]       [,3]
#> [1,] 0.007664882 0.45423051 0.53810461
#> [2,] 0.406895342 0.57088656 0.02221810
#> [3,] 0.335196674 0.62893638 0.03586694
#> [4,] 0.417697775 0.08456334 0.49773888
#> [5,] 0.005668935 0.88696625 0.10736481
```
