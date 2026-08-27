# Fisher's z-Transformation

Convert a Pearson correlation coefficient to Fisher's \\z\\ scale and
back. The transformation stabilizes the variance of the correlation
coefficient and yields approximately normally distributed values.

## Usage

``` r
fisherZ(rho)

fisherZInv(z)
```

## Arguments

- rho:

  numeric vector. Pearson correlation coefficient(s), typically in the
  interval \\\[-1, 1\]\\. Values of \\\pm 1\\ are mapped to \\\pm
  \infty\\.

- z:

  numeric vector. Fisher \\z\\-transformed values.

## Value

A numeric vector as follows:

- `fisherZ`:

  Fisher \\z\\-transformed values.

- `fisherZInv`:

  correlation coefficients.

## Details

The forward transformation is defined as \$\$ z = \tanh^{-1}(r) =
\frac{1}{2}\log\left(\frac{1 + r}{1 - r}\right), \$\$ and the inverse
transformation as \$\$ r = \tanh(z). \$\$

Fisher's \\z\\-transformation is commonly used to construct confidence
intervals and perform hypothesis tests for correlation coefficients.

## See also

[`corCI`](https://andrisignorell.github.io/lumen/reference/corCI.md),
[`cor.test`](https://rdrr.io/r/stats/cor.test.html)

Other test.correlation:
[`corTest()`](https://andrisignorell.github.io/lumen/reference/corTest.md)

## Examples

``` r
# Forward and inverse transformation
r <- seq(-0.9, 0.9, by = 0.1)
z <- fisherZ(r)
fisherZInv(z)
#>  [1] -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1  0.0  0.1  0.2  0.3  0.4  0.5
#> [16]  0.6  0.7  0.8  0.9

# Round-trip accuracy
all.equal(r, fisherZInv(fisherZ(r)))
#> [1] TRUE
```
