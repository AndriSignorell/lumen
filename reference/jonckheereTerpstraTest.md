# Jonckheere-Terpstra Test for Detecting Ordered Differences Across Independent Groups

A nonparametric test for monotonic trends across ordered independent
groups, assessing whether observations tend to increase (or decrease)
systematically with group order.

## Usage

``` r
jonckheereTerpstraTest(x, ...)

# S3 method for class 'formula'
jonckheereTerpstraTest(formula, data, subset, na.action, ...)

# Default S3 method
jonckheereTerpstraTest(
  x,
  g,
  alternative = c("two.sided", "increasing", "decreasing"),
  method = c("auto", "exact", "permutation", "asymptotic"),
  R = NULL,
  ...
)
```

## Arguments

- x:

  a numeric vector of observations, or a list of numeric vectors.

- ...:

  further arguments passed to methods.

- formula:

  a formula of the form `response ~ group`.

- data:

  an optional data frame containing the variables in `formula`.

- subset:

  an optional expression specifying a subset of observations.

- na.action:

  a function indicating how missing values should be handled.

- g:

  a grouping variable corresponding to `x`, whose (factor) level order
  defines the hypothesised ordering; ignored when `x` is a list.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"increasing"` or `"decreasing"`.

- method:

  a character string specifying the inference method, one of `"auto"`
  (default), `"exact"`, `"permutation"` or `"asymptotic"`. `"auto"` uses
  exact inference when possible (no ties, \\n \le 100\\), otherwise the
  asymptotic approximation.

- R:

  the number of permutations, required when `method = "permutation"`.

## Value

A list with class `"htest"` containing the following components:

- statistic:

  the value of the Jonckheere-Terpstra statistic.

- parameter:

  the number of groups `k` and the total sample size `n`.

- p.value:

  the p-value of the test.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string indicating the test performed and the inference
  method used.

- data.name:

  a character string giving the names of the data.

## Details

The Jonckheere-Terpstra statistic is \$\$JT = \sum\_{k\<l} \sum\_{ij}
\left\[ I(X\_{ik} \< X\_{jl}) + \frac{1}{2} I(X\_{ik} = X\_{jl})
\right\]\$\$ where \\i, j\\ index observations from ordered groups \\k,
l\\. Large values of the statistic indicate increasing trends across
groups.

Exact p-values are computed from the exact permutation distribution
using a dynamic programming recursion implemented in C++. Exact
inference is only valid without ties and, for practical reasons, is
offered for total sample sizes \\n \le 100\\.

When ties are present or sample sizes are large, permutation p-values
can be computed by permuting group labels under the null hypothesis
(`method = "permutation"`); the number of permutations is controlled by
`R`, and the reported p-value uses the finite-sample correction \\(m +
1)/(R + 1)\\. This approach remains valid in the presence of ties.

With `method = "asymptotic"` (the fallback of `"auto"` when exact
inference does not apply), a normal approximation with the tie-corrected
variance of Hollander and Wolfe (1999, Eq. 6.19) is used.

## References

Jonckheere, A. R. (1954) A distribution-free k-sample test against
ordered alternatives. *Biometrika*, 41, 133–145.

Terpstra, T. J. (1952) The asymptotic normality and consistency of
Kendall's test against trend, when ties are present in one ranking.
*Indagationes Mathematicae*, 14, 327–333.

Hollander, M. and Wolfe, D. A. (1999) *Nonparametric Statistical
Methods*, 2nd ed., New York: Wiley.

## See also

[`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)

Other test.trend:
[`cochranArmitageTest()`](https://andrisignorell.github.io/lumen/reference/cochranArmitageTest.md),
[`coxStuartTest()`](https://andrisignorell.github.io/lumen/reference/coxStuartTest.md),
[`pageTest()`](https://andrisignorell.github.io/lumen/reference/pageTest.md)

## Examples

``` r
set.seed(1)
g <- ordered(rep(1:4, each = 10))
x <- rnorm(40) + 0.5 * as.numeric(g)

jonckheereTerpstraTest(x, g)
#> 
#>  Jonckheere-Terpstra test for ordered alternatives (exact)
#> 
#> data:  x and g
#> JT = 437, k = 4, n = 40, p-value = 0.0007562
#> alternative hypothesis: two.sided
#> 

# with ties: permutation inference
x[1:2] <- mean(x[1:2])
jonckheereTerpstraTest(x, g, method = "permutation", R = 2000)
#> 
#>  Jonckheere-Terpstra test for ordered alternatives (permutation, R =
#>  2000)
#> 
#> data:  x and g
#> JT = 438, k = 4, n = 40, p-value = 0.0004998
#> alternative hypothesis: two.sided
#> 

coffee <- list(
  c_4 = c(447, 396, 383, 410),
  c_2 = c(438, 521, 468, 391, 504, 472),
  c_0 = c(513, 543, 506, 489, 407)
)
jonckheereTerpstraTest(coffee)
#> 
#>  Jonckheere-Terpstra test for ordered alternatives (exact)
#> 
#> data:  x
#> JT = 59, k = 3, n = 15, p-value = 0.01973
#> alternative hypothesis: two.sided
#> 

# Hollander & Wolfe, Example 6.2:
# motivational effect of knowledge of performance
motiv <- list(
  no    = c(40, 35, 38, 43, 44, 41),
  rough = c(38, 40, 47, 44, 40, 42),
  acc   = c(48, 40, 45, 43, 46, 44))

jonckheereTerpstraTest(motiv, alternative = "increasing")
#> 
#>  Jonckheere-Terpstra test for ordered alternatives (asymptotic)
#> 
#> data:  x
#> JT = 79, k = 3, n = 18, p-value = 0.02071
#> alternative hypothesis: increasing
#> 
## exact one-sided p-value 0.0379 as in Hollander & Wolfe

jonckheereTerpstraTest(motiv, alternative = "increasing",
                       method = "asymptotic")
#> 
#>  Jonckheere-Terpstra test for ordered alternatives (asymptotic)
#> 
#> data:  x
#> JT = 79, k = 3, n = 18, p-value = 0.02071
#> alternative hypothesis: increasing
#> 

set.seed(42)
jonckheereTerpstraTest(motiv, alternative = "increasing",
                       method = "permutation", R = 10000)
#> 
#>  Jonckheere-Terpstra test for ordered alternatives (permutation, R =
#>  10000)
#> 
#> data:  x
#> JT = 79, k = 3, n = 18, p-value = 0.0194
#> alternative hypothesis: increasing
#> 
```
