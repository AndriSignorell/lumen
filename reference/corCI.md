# Confidence Interval for a Pearson Correlation

Computes confidence intervals for a Pearson correlation coefficient
using Fisher's \\z\\-transformation. Since the sampling distribution of
the correlation coefficient is not normally distributed, Fisher's
transformation is applied to obtain approximately normally distributed
values, from which confidence intervals can be derived.

## Usage

``` r
corCI(rho, n, conf.level = 0.95, sides = c("two.sided", "left", "right"))
```

## Arguments

- rho:

  numeric; Pearson correlation coefficient. Must be a single value in
  the interval \\\[-1, 1\]\\.

- n:

  integer; sample size used to estimate the correlation. Must be at
  least 3.

- conf.level:

  confidence level of the interval, a single number in \\(0, 1)\\.
  Defaults to `0.95`; unlike the estimators in this package it cannot be
  `NA`, since the interval is all this function computes.

- sides:

  character string specifying the sidedness of the confidence interval
  (one of `"two.sided"` (default), `"left"` or `"right"`). See
  `DescToolsX::ConfidenceIntervals`.

## Value

A named numeric vector with elements:

- `est`:

  point estimate (correlation coefficient `rho`).

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

## Details

Fisher's \\z\\-transformation is defined as \$\$z = \frac{1}{2}
\log\left(\frac{1 + r}{1 - r}\right),\$\$ which stabilizes the variance
of the correlation coefficient. The transformed values are approximately
normally distributed with standard error \\1 / \sqrt{n - 3}\\.
Confidence intervals are computed on the transformed scale and then
back-transformed.

A correlation lies in \\\[-1, 1\]\\, so the open side of a one-sided
interval is reported at that boundary rather than at \\\pm\infty\\.

The transformation has nothing to say in two situations, and both are
answered with `NA` bounds and a warning rather than with an interval
that only looks informative: at \\\|\rho\| = 1\\, where \\z\\ is
infinite and the interval would collapse onto the estimate and thereby
rule out every other value, and at \\n = 3\\, where the standard error
\\1/\sqrt{n-3}\\ is infinite and the interval would be the whole range.

## Note

Based on code by William Revelle, adapted to conform to package
standards.

## Argument name

This function used to take `alternative` with the values `"less"` and
`"greater"`. It now takes `sides`, like every other interval function in
the package, and `sides` names the side carrying the *finite* bound - so
`"left"` is the former `"greater"` and `"right"` the former `"less"`.
The values produced are unchanged.

## See also

[`fisherZ`](https://andrisignorell.github.io/lumen/reference/fisherZ.md),
[`fisherZInv`](https://andrisignorell.github.io/lumen/reference/fisherZ.md),
[`cor.test`](https://rdrr.io/r/stats/cor.test.html)

## Examples

``` r
# Confidence interval for a single correlation
corCI(0.5, n = 30)
#>       est       lci       uci 
#> 0.5000000 0.1704314 0.7289586 

# one-sided: "left" carries the finite lower bound
corCI(0.5, n = 30, sides = "left")
#>       est       lci       uci 
#> 0.5000000 0.2286399 1.0000000 

# Compare multiple correlations
r <- seq(0, 0.9, by = 0.1)
t(sapply(r, corCI, n = 30))
#>       est         lci       uci
#>  [1,] 0.0 -0.36026922 0.3602692
#>  [2,] 0.1 -0.26999636 0.4442638
#>  [3,] 0.2 -0.17271392 0.5226130
#>  [4,] 0.3 -0.06757251 0.5958674
#>  [5,] 0.4  0.04642030 0.6645084
#>  [6,] 0.5  0.17043137 0.7289586
#>  [7,] 0.6  0.30584207 0.7895902
#>  [8,] 0.7  0.45429999 0.8467329
#>  [9,] 0.8  0.61778629 0.9006795
#> [10,] 0.9  0.79870459 0.9516908
```
