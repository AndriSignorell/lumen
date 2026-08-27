# Random proportions summing to 1

Generates random positive values that sum to 1 by sampling from a
Dirichlet distribution with all concentration parameters equal to 1,
i.e. Dirichlet(1, ..., 1). This corresponds to a uniform distribution
over the simplex.

## Usage

``` r
rsum1(n, digits = NULL)
```

## Arguments

- n:

  integer. Number of components (must be \>= 1).

- digits:

  optional integer. If provided, the result is rounded to the specified
  number of decimal places. A correction is applied to ensure the sum
  remains exactly 1.

## Value

Numeric vector of length `n` with non-negative entries summing to 1.

## Details

Optionally, the result can be rounded to a specified number of digits
while preserving the total sum of 1.

Values are generated via normalized Gamma draws: \$\$x_i = g_i / \sum_j
g_j, \quad g_i \sim \Gamma(1, 1)\$\$

If `digits` is specified, rounding may introduce small numerical
deviations. These are corrected by adjusting the largest component,
which minimizes relative distortion and avoids negative values.

## Examples

``` r
rsum1(5)
#> [1] 0.07209899 0.22794895 0.37010675 0.20088331 0.12896201
rsum1(5, digits = 2)
#> [1] 0.03 0.00 0.36 0.05 0.56
rsum1(1)
#> [1] 1
```
