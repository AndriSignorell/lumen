# Null Distribution of Anderson-Darling Test Statistic

`pAD` computes the cumulative distribution function, and `qAD` computes
the quantile function, of the null distribution of the Anderson-Darling
test statistic.

## Usage

``` r
pAD(q, n = Inf, lower.tail = TRUE, fast = TRUE)

qAD(p, n = Inf, lower.tail = TRUE, fast = TRUE)
```

## Arguments

- q:

  Numeric vector of quantiles (values for which the cumulative
  probability is required).

- n:

  Integer. Sample size for the Anderson-Darling test.

- lower.tail:

  Logical. If `TRUE` (the default), probabilities are \\P(X \le q)\\,
  and otherwise they are \\P(X \> q)\\.

- fast:

  Logical value indicating whether to use a fast algorithm or a slower,
  more accurate algorithm, in the case `n=Inf`.

- p:

  Numeric vector of probabilities.

## Value

A numeric vector of the same length as `p` or `q`.

## Details

`pAD` uses the algorithms and C code described in Marsaglia and
Marsaglia (2004).

`qAD` uses [`uniroot`](https://rdrr.io/r/stats/uniroot.html) to find the
quantiles.

The argument `fast` applies only when `n=Inf` and determines whether the
asymptotic distribution is approximated using the faster algorithm
`adinf` (accurate to 4-5 places) or the slower algorithm `ADinf`
(accurate to 11 places) described in Marsaglia and Marsaglia (2004).

## Note

Original C code by G. and J. Marsaglia. R interface by Adrian Baddeley,
adapted to conform to package standards.

## References

Anderson, T.W. and Darling, D.A. (1952) Asymptotic theory of certain
'goodness-of-fit' criteria based on stochastic processes. *Annals of
Mathematical Statistics* **23**, 193–212.

Anderson, T.W. and Darling, D.A. (1954) A test of goodness of fit.
*Journal of the American Statistical Association* **49**, 765–769.

Marsaglia, G. and Marsaglia, J. (2004) Evaluating the Anderson-Darling
Distribution. *Journal of Statistical Software* **9** (2), 1–5. February
2004. <http://www.jstatsoft.org/v09/i02>

## See also

[`andersonDarlingTest()`](https://andrisignorell.github.io/lumen/reference/andersonDarlingTest.md),
[distributions-overview](https://andrisignorell.github.io/lumen/reference/distributions-overview.md)

## Examples

``` r

  pAD(1.1, n=5)
#> [1] 0.6945818
  pAD(1.1)
#> [1] 0.6911871
  pAD(1.1, fast=FALSE)
#> [1] 0.6912038

  qAD(0.5, n=5)
#> [1] 0.7640094
  qAD(0.5)
#> [1] 0.7742347
```
