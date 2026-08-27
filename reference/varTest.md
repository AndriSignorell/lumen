# Variance Test for Testing One Variance or Comparing Two Variances

Performs a one-sample or two-sample test for variance, analogous to
[`t.test`](https://rdrr.io/r/stats/t.test.html), with support for
classical and likelihood-based lowest-density (LD) two-sided p-values.

## Usage

``` r
varTest(x, ...)

# Default S3 method
varTest(
  x,
  y = NULL,
  sigma2_0 = NULL,
  alternative = c("two.sided", "less", "greater"),
  type = c("classic", "ld"),
  ...
)

# S3 method for class 'formula'
varTest(formula, data, subset, na.action = na.pass, ...)
```

## Arguments

- x:

  a numeric vector of data values, or a formula.

- ...:

  further arguments passed to methods.

- y:

  an optional second numeric vector. If provided, a two-sample variance
  test is performed.

- sigma2_0:

  a numeric value specifying the null hypothesis variance for the
  one-sample test. Required if `y` is `NULL`.

- alternative:

  character string specifying the alternative hypothesis. Must be one of
  `"two.sided"`, `"less"`, or `"greater"`.

- type:

  character string specifying the test type:

  - `"classic"`: uses the conventional two-sided p-value \\2 \cdot
    \min(P(T \le t), P(T \ge t))\\.

  - `"ld"`: uses a likelihood-based lowest-density (LD) definition, i.e.
    the probability of observing values with density less than or equal
    to the observed density under the null distribution.

- formula:

  a formula of the form `lhs ~ rhs`, where `lhs` gives the response
  values and `rhs` the corresponding groups or explanatory variables.

- data:

  an optional matrix or data frame (or similar; see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in the formula. By default the variables are taken from
  `environment(formula)`.

- subset:

  an optional vector specifying a subset of observations to be used in
  the analysis.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

## Value

An object of class `"htest"` with components:

- statistic:

  the test statistic.

- parameter:

  degrees of freedom.

- p.value:

  the p-value of the test.

- estimate:

  sample variance(s).

- null.value:

  the null hypothesis value (one-sample only).

- alternative:

  the alternative hypothesis.

- method:

  a character string indicating the test performed.

- data.name:

  description of the data.

## Details

The null hypothesis is that the ratio of the variances of the
populations from which `x` and `y` were drawn, or in the data to which
the linear models `x` and `y` were fitted, is equal to `ratio`.

For the one-sample test, the test statistic follows a chi-squared
distribution: \$\$X^2 = (n - 1) S^2 / \sigma_0^2\$\$

For the two-sample test, the test statistic follows an F distribution:
\$\$F = S_x^2 / S_y^2\$\$

The LD method corresponds to a likelihood-ratio interpretation of
extremeness, which is particularly appropriate for asymmetric null
distributions such as chi-squared and F distributions.

The formula interface is only applicable for the 2-sample tests.

## See also

[`var.test`](https://rdrr.io/r/stats/var.test.html),
[`bartlett.test`](https://rdrr.io/r/stats/bartlett.test.html) for
testing homogeneity of variances in more than two samples from normal
distributions; [`ansari.test`](https://rdrr.io/r/stats/ansari.test.html)
and [`mood.test`](https://rdrr.io/r/stats/mood.test.html) for two rank
based (nonparametric) two-sample tests for difference in scale.

Other test.variance:
[`leveneTest()`](https://andrisignorell.github.io/lumen/reference/leveneTest.md),
[`mosesTest()`](https://andrisignorell.github.io/lumen/reference/mosesTest.md),
[`siegelTukeyTest()`](https://andrisignorell.github.io/lumen/reference/siegelTukeyTest.md),
[`varCI()`](https://andrisignorell.github.io/lumen/reference/varCI.md)

## Examples

``` r
set.seed(1)
x <- rnorm(20, sd = 3)
y <- rnorm(25, sd = 2)

# One-sample test
varTest(x, sigma2_0 = 9, type = "classic")
#> 
#>  One-sample variance test (classic)
#> 
#> data:  x
#> X-squared = 15.847, df = 19, p-value = 0.665
#> alternative hypothesis: true variance is not equal to 9
#> sample estimates:
#> variance 
#> 7.506291 
#> 
varTest(x, sigma2_0 = 9, type = "ld")
#> 
#>  One-sample variance test (ld)
#> 
#> data:  x
#> X-squared = 15.847, df = 19, p-value = 0.8411
#> alternative hypothesis: true variance is not equal to 9
#> sample estimates:
#> variance 
#> 7.506291 
#> 

# Two-sample test
varTest(x, y, type = "classic")
#> 
#>  Two-sample variance test (classic)
#> 
#> data:  x and y
#> F = 2.8526, df1 = 19, df2 = 24, p-value = 0.0165
#> alternative hypothesis: two.sided
#> sample estimates:
#>   var(x)   var(y) 
#> 7.506291 2.631395 
#> 
varTest(x, y, type = "ld")
#> 
#>  Two-sample variance test (ld)
#> 
#> data:  x and y
#> F = 2.8526, df1 = 19, df2 = 24, p-value = 0.008804
#> alternative hypothesis: two.sided
#> sample estimates:
#>   var(x)   var(y) 
#> 7.506291 2.631395 
#> 

# Formula interface
df <- data.frame(
  value = c(x, y),
  group = rep(c("A", "B"), c(length(x), length(y)))
)
varTest(value ~ group, data = df)
#> 
#>  Two-sample variance test (classic)
#> 
#> data:  value ~ group
#> F = 2.8526, df1 = 19, df2 = 24, p-value = 0.0165
#> alternative hypothesis: two.sided
#> sample estimates:
#>   var(x)   var(y) 
#> 7.506291 2.631395 
#> 
```
