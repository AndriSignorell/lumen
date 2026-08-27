# Siegel-Tukey Test for Comparing Dispersion Between Two Independent Groups

The Siegel-Tukey test examines the null hypothesis that the variability
(scale) of `x` and `y` is equal. Rejection indicates that the two groups
differ in spread. The test is distribution-free and does not assume
normality, but it does assume equal medians under the null hypothesis of
equal scale. If the medians differ, the test may detect that difference
rather than a difference in scale; use `adjustMedian = TRUE` to remove
median differences before testing.

Ranks are assigned to the combined sorted sample in the pattern 1, 4, 5,
8, 9, ... from the extremes inward, and 2, 3, 6, 7, ... from the
second-lowest upward. If the combined sample size is odd, the median
observation is dropped before ranking (it is taken from the larger group
when group sizes differ).

Ties receive average ranks. The p-value is computed exactly (via
[`pwilcox`](https://rdrr.io/r/stats/Wilcoxon.html)) when there are no
ties and both samples are smaller than 50 observations; otherwise a
normal approximation with tie-corrected variance is used. This behaviour
can be overridden with `exact`.

**Note:** The Siegel-Tukey test has relatively low power compared to
alternatives such as
[`ansari.test`](https://rdrr.io/r/stats/ansari.test.html) or
[`mood.test`](https://rdrr.io/r/stats/mood.test.html), and may indicate
significance due to median differences rather than scale differences
when `adjustMedian = FALSE`.

## Usage

``` r
siegelTukeyTest(x, ...)

# S3 method for class 'formula'
siegelTukeyTest(formula, data, subset, na.action = na.pass, ...)

# Default S3 method
siegelTukeyTest(
  x,
  y,
  alternative = c("two.sided", "less", "greater"),
  mu = 0,
  adjustMedian = FALSE,
  exact = NA,
  correct = TRUE,
  ...
)
```

## Arguments

- x, y:

  numeric vectors of data values.

- ...:

  further arguments passed to or from methods.

- formula:

  a formula of the form `response ~ group`, where `response` is a
  numeric vector and `group` a factor or vector with exactly two levels.

- data:

  an optional data frame (or matrix, coerced to data frame) containing
  the variables in `formula`. If not supplied, variables are taken from
  `environment(formula)`.

- subset:

  an optional vector specifying a subset of observations to use.

- na.action:

  a function indicating how to handle `NA`s in the formula interface.
  Defaults to `na.pass`; `NA`s in `x` or `y` are silently dropped in the
  default method.

- alternative:

  a character string specifying the alternative hypothesis:
  `"two.sided"` (default), `"greater"`, or `"less"`. Partial matching is
  allowed.

- mu:

  a single number specifying the location parameter under the null
  hypothesis. Default is `0`.

- adjustMedian:

  logical; if `TRUE`, the median of `x` is shifted to equal the median
  of `y` before ranking, to prevent median differences from inflating
  the test statistic. Default is `FALSE`.

- exact:

  logical; if `TRUE`, an exact p-value is computed via
  [`pwilcox`](https://rdrr.io/r/stats/Wilcoxon.html). Exact computation
  is not possible in the presence of ties; a warning is issued and the
  normal approximation is used instead. If `NA` (default), exact
  computation is used when both samples have fewer than 50 observations
  and there are no ties.

- correct:

  logical; if `TRUE` (default), a continuity correction is applied in
  the normal approximation. Ignored when `exact = TRUE` or when ties are
  present (continuity correction is not appropriate with tie-corrected
  variance).

## Value

An object of class `"htest"` with the following components:

- statistic:

  the Wilcoxon rank-sum statistic \\W\\ computed on the Siegel-Tukey
  ranks.

- p.value:

  the p-value of the test.

- null.value:

  the location parameter `mu` under the null hypothesis.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string identifying the test.

- data.name:

  a character string giving the names of the data.

- exact:

  logical indicating whether the exact distribution was used.

- ties:

  logical indicating whether ties were present in the Siegel-Tukey
  ranks.

## Details

A nonparametric test for differences in scale (variability) between two
independent groups. Ranks are assigned by alternating between the
extremes and the center of the combined sorted sample, and a
Wilcoxon-Mann-Whitney statistic is then computed on these ranks.

## Note

Originally based on a blog post by Tal Galili:  
<https://www.r-statistics.com/2010/02/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/>

## References

Siegel, S. and Tukey, J. W. (1960): A nonparametric sum of ranks
procedure for relative spread in unpaired samples. *Journal of the
American Statistical Association*, **55**(291), 429–445.

Sheskin, D. J. (2004): *Handbook of Parametric and Nonparametric
Statistical Procedures*, 3rd ed. Chapman & Hall/CRC, Boca Raton, FL.

## See also

[`ansari.test`](https://rdrr.io/r/stats/ansari.test.html),
[`mood.test`](https://rdrr.io/r/stats/mood.test.html),
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)

Other test.variance:
[`leveneTest()`](https://andrisignorell.github.io/lumen/reference/leveneTest.md),
[`mosesTest()`](https://andrisignorell.github.io/lumen/reference/mosesTest.md),
[`varCI()`](https://andrisignorell.github.io/lumen/reference/varCI.md),
[`varTest()`](https://andrisignorell.github.io/lumen/reference/varTest.md)

## Examples

``` r
# Duller, S. 183
x <- c(12, 13, 29, 30)
y <- c(15, 17, 18, 24, 25, 26)
siegelTukeyTest(x, y)
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  x and y
#> W = 45, p-value = 0.009524
#> alternative hypothesis: true mu is not equal to 0
#> 
siegelTukeyTest(x, y, alternative = "greater")
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  x and y
#> W = 45, p-value = 0.004762
#> alternative hypothesis: true mu is greater than 0
#> 

# Duller, S. 323
old <- c(870, 930, 935, 1045, 1050, 1052, 1055)
new <- c(932, 970, 980, 1001, 1009, 1030, 1032, 1040, 1046)
siegelTukeyTest(old, new, alternative = "greater")
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  old and new
#> W = 102, p-value = 0.002622
#> alternative hypothesis: true mu is greater than 0
#> 
# recommended alternatives:
mood.test(old, new, alternative = "greater")
#> 
#>  Mood two-sample test of scale
#> 
#> data:  old and new
#> Z = 2.8666, p-value = 0.002075
#> alternative hypothesis: greater
#> 
ansari.test(old, new, alternative = "greater")
#> 
#>  Ansari-Bradley test
#> 
#> data:  old and new
#> AB = 18, p-value = 0.001573
#> alternative hypothesis: true ratio of scales is greater than 1
#> 

# Bortz, S. 250 -- formula interface
x  <- c(26.3, 26.5, 26.8, 27.0, 27.0, 27.2, 27.3,
        27.3, 27.4, 27.5, 27.6, 27.8, 27.9)
id <- c(2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2) - 1
siegelTukeyTest(x ~ factor(id))
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  groups[[1L]] and groups[[2L]]
#> W = 53.5, p-value = 0.8649
#> alternative hypothesis: true mu is not equal to 0
#> 

# Sachs (2007), S. 314
A <- c(10.1, 7.3, 12.6, 2.4, 6.1, 8.5, 8.8, 9.4, 10.1, 9.8)
B <- c(15.3, 3.6, 16.5, 2.9, 3.3, 4.2, 4.9, 7.3, 11.7, 13.1)
siegelTukeyTest(A, B)
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  A and B
#> W = 75.5, p-value = 0.02825
#> alternative hypothesis: true mu is not equal to 0
#> 

# equal medians, different spread
x <- c(4, 4, 5, 5, 6, 6)
y <- c(0, 0, 1, 9, 10, 10)
siegelTukeyTest(x, y)
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  x and y
#> W = 21, p-value = 0.003601
#> alternative hypothesis: true mu is not equal to 0
#> 

# unequal group sizes
x <- c(4, 4, 5, 5, 6, 6)
y <- c(0, 0, 1, 9, 10)
siegelTukeyTest(x, y)
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  x and y
#> W = 15, p-value = 0.01141
#> alternative hypothesis: true mu is not equal to 0
#> 

# median adjustment
x <- c(177, 200, 227, 230, 232, 268, 272, 297)
y <- c(47, 105, 126, 142, 158, 172, 197, 220, 225, 230, 262, 270)
siegelTukeyTest(x, y, adjustMedian = TRUE)
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  x and y
#> W = 104, p-value = 0.09788
#> alternative hypothesis: true mu is not equal to 0
#> 

# floating-point robustness (previously affected by merge precision bug)
y2 <- c(-1, 2, 2.1, 3)
x2 <- c(-5, -9, 13, 12, 90, 100)
siegelTukeyTest(x2, y2, adjustMedian = TRUE)  
#> 
#>  Siegel-Tukey test for scale differences
#> 
#> data:  x2 and y2
#> W = 30, p-value = 0.1143
#> alternative hypothesis: true mu is not equal to 0
#> 
# p ~ 0.1143
```
