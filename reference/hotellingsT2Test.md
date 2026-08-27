# Hotelling's T2 Test for Comparing Multivariate Mean Vectors

The classical test for the location of a multivariate population
(one-sample) or for the difference between the mean vectors of two
multivariate populations (two-sample). It is the multivariate
generalisation of Student's t-test.

## Usage

``` r
hotellingsT2Test(x, ...)

# S3 method for class 'formula'
hotellingsT2Test(formula, data, subset, na.action = na.pass, ...)

# Default S3 method
hotellingsT2Test(x, y = NULL, mu = NULL, test = c("f", "chi"), ...)
```

## Arguments

- x:

  a numeric matrix or data frame (observations in rows, variables in
  columns).

- ...:

  further arguments passed to or from methods.

- formula:

  a formula of the form `cbind(v1, v2, ...) ~ g` where the left-hand
  side is a numeric matrix of response variables and `g` is a factor
  with exactly two levels.

- data:

  an optional data frame (or similar, see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in `formula`. Defaults to the environment of `formula`.

- subset:

  an optional vector specifying a subset of observations.

- na.action:

  a function indicating what should happen when the data contain `NA`s.
  Defaults to `getOption("na.action")`.

- y:

  an optional numeric matrix or data frame for the two-sample test. If
  `NULL` (default) a one-sample test is performed.

- mu:

  a numeric vector of length \\p\\ giving the hypothesised mean
  (one-sample) or mean difference (two-sample). `NULL` is interpreted as
  the zero vector.

- test:

  a character string selecting the reference distribution: `"f"` (exact
  F-distribution, default) or `"chi"` (chi-squared approximation).

## Value

An object of class `"htest"` containing:

- statistic:

  the value of the T2-statistic (scaled to follow an F- or chi-squared
  distribution depending on `test`).

- parameter:

  degrees of freedom of the reference distribution.

- p.value:

  the p-value of the test.

- estimate:

  the sample mean vector (one-sample) or the difference of the sample
  mean vectors (two-sample).

- null.value:

  the hypothesised mean or mean difference.

- alternative:

  always `"two.sided"`.

- method:

  a character string describing the test variant performed.

- data.name:

  a character string giving the name(s) of the input data.

## Details

When `test = "f"` the test statistic follows an exact F-distribution
under the assumption of multivariate normality. When `test = "chi"` a
chi-squared approximation is used; it relies on large-sample asymptotic
theory and is less sensitive to departures from multivariate normality
than the F-test, but remains only asymptotically correct.

In the two-sample case both populations are assumed to share the same
covariance matrix; a pooled within-group estimate is used.

The formula interface (`cbind(v1, v2) ~ g`) is available for the
two-sample case only.

## Note

Based on code by Klaus Nordhausen, adapted to conform to package
standards.

## References

Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
Analysis* (3rd ed.). Wiley.

Nordhausen, K., Sirkia, S., Oja, H., & Tyler, D. E. (2012). *ICSNP:
Tools for Multivariate Nonparametrics*. R package version 1.0-9.
<https://cran.r-project.org/package=ICSNP>

## See also

Other test.location: [`brunnerMunzelTest()`](brunnerMunzelTest.md),
[`moodMedianTest()`](moodMedianTest.md), [`signTest()`](signTest.md),
[`tTestA()`](tTestA.md), [`vanWaerdenTest()`](vanWaerdenTest.md),
[`yuenTTest()`](yuenTTest.md), [`zTest()`](zTest.md)

## Examples

``` r
math.teach <- data.frame(
  teacher = factor(rep(1:2, c(3, 6))),
  satis   = c(1, 3, 2, 4, 6, 6, 5, 5, 4),
  know    = c(3, 7, 2, 6, 8, 8, 10, 10, 6))

hotellingsT2Test(cbind(satis, know) ~ teacher, data = math.teach)
#> 
#>  Hotelling's two-sample T-squared test
#> 
#> data:  cbind(satis, know) ~ teacher
#> T.2 = 9, df1 = 2, df2 = 6, p-value = 0.01562
#> alternative hypothesis: two.sided
#> null values:
#> location difference location difference 
#>                   0                   0 
#> sample estimates:
#> mean difference of satis  mean difference of know 
#>                       -3                       -4 
#> 

# chi-squared approximation
hotellingsT2Test(cbind(satis, know) ~ teacher, data = math.teach,
                 test = "chi")
#> 
#>  Hotelling's two-sample T-squared test
#> 
#> data:  cbind(satis, know) ~ teacher
#> T.2 = 21, df = 2, p-value = 2.754e-05
#> alternative hypothesis: two.sided
#> null values:
#> location difference location difference 
#>                   0                   0 
#> sample estimates:
#> mean difference of satis  mean difference of know 
#>                       -3                       -4 
#> 
```
