# Fast Correlation Test for Testing Pairwise Correlations in a Matrix

Compute correlation coefficients, p-values, and pairwise sample sizes
for all variable pairs in a numeric matrix.

## Usage

``` r
corTest(
  x,
  method = c("pearson", "spearman", "kendall"),
  use = "pairwise.complete.obs",
  triangle = c("full", "upper", "lower"),
  maxPValue = NULL
)
```

## Arguments

- x:

  a numeric matrix or data frame.

- method:

  a character string specifying the correlation method, one of
  `"pearson"` (default), `"spearman"` or `"kendall"`. Passed to
  [`cor`](https://rdrr.io/r/stats/cor.html).

- use:

  a character string giving a method for computing correlations in the
  presence of missing values, passed to
  [`cor`](https://rdrr.io/r/stats/cor.html). Default is
  `"pairwise.complete.obs"`.

- triangle:

  a character string specifying which part of the matrices should be
  returned, one of `"full"` (default), `"upper"` or `"lower"`. The other
  triangle is set to `NA`.

- maxPValue:

  optional upper limit for reported correlations. Correlations with
  p-values larger than `maxPValue` are replaced with `NA`.

## Value

A list with three matrices:

- `cor`:

  the correlation matrix.

- `pValue`:

  the matrix of two-sided p-values (`NA` on the diagonal).

- `n`:

  the matrix of sample sizes used per pair.

## Details

Pearson, Spearman, and Kendall correlations are supported via the
`method` argument passed to [`cor`](https://rdrr.io/r/stats/cor.html).

Compared to repeatedly calling
[`cor.test`](https://rdrr.io/r/stats/cor.test.html), this implementation
is fully vectorised and substantially faster for matrices with many
variables.

For Pearson and Spearman correlations, the two-sided p-value under the
null hypothesis \\\rho = 0\\ is computed from the t statistic \$\$t = r
\sqrt{\frac{n - 2}{1 - r^2}}\$\$ with \\n - 2\\ degrees of freedom. This
is exact for the Pearson coefficient under normality and an
approximation for the Spearman coefficient (as used by `cor.test` for
larger samples). For the Kendall coefficient the normal approximation
\$\$z = \frac{3 \tau \sqrt{n (n-1)}}{\sqrt{2 (2 n + 5)}}\$\$ is used;
note that no correction for ties is applied here. Cells with fewer than
3 (Pearson/Spearman) resp. 2 (Kendall) pairwise observations get an `NA`
p-value.

If `use` requests pairwise handling of missing values (the default), the
sample sizes in `n` are the pairwise counts; otherwise all pairs share
the number of complete cases.

## See also

[`cor`](https://rdrr.io/r/stats/cor.html),
[`cor.test`](https://rdrr.io/r/stats/cor.test.html)

Other test.correlation: [`fisherZ()`](fisherZ.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200), 50, 4)

res <- corTest(X)

# only significant correlations in the upper triangle
res2 <- corTest(
  X,
  triangle = "upper",
  maxPValue = 0.05
)
```
