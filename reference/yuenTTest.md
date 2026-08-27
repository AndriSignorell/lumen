# Yuen T-Test for Robust Comparison of Trimmed Means

Performs Yuen's robust t-test for trimmed means. Compared with the
classical t-test, the procedure is substantially less sensitive to
outliers, heavy tails, and moderate departures from normality.

The test is based on:

- trimmed means,

- winsorized variances,

- Welch-type degrees of freedom.

For paired tests, trimming is performed on the paired differences,
following Yuen (1974).

## Usage

``` r
yuenTTest(x, ...)

# S3 method for class 'formula'
yuenTTest(formula, data, subset, na.action = na.pass, ...)

# Default S3 method
yuenTTest(
  x,
  y = NULL,
  alternative = c("two.sided", "less", "greater"),
  mu = 0,
  paired = FALSE,
  conf.level = 0.95,
  trim = 0.2,
  ...
)
```

## Arguments

- x:

  numeric vector of observations.

- ...:

  further arguments passed to methods.

- formula:

  a formula of the form `lhs ~ rhs`.

- data:

  optional data frame for the formula interface.

- subset:

  optional subset expression.

- na.action:

  NA handling function.

- y:

  optional second numeric vector.

- alternative:

  character string specifying the alternative hypothesis. One of
  `"two.sided"`, `"less"`, or `"greater"`.

- mu:

  hypothesized trimmed mean (or trimmed mean difference).

- paired:

  logical indicating whether a paired test is performed.

- conf.level:

  confidence level for the confidence interval.

- trim:

  fraction of observations trimmed from each tail. Must satisfy
  `0 <= trim < 0.5`.

## Value

An object of class `"htest"`.

## Details

Robust one-, two-, and paired-sample t-tests based on trimmed means and
winsorized variances.

## References

Wilcox, R. R. (2005). *Introduction to Robust Estimation and Hypothesis
Testing*. Academic Press.

Yuen, K. K. (1974). The two-sample trimmed t for unequal population
variances. *Biometrika*, 61, 165–170.

## See also

[`t.test()`](https://rdrr.io/r/stats/t.test.html)

Other test.location: [`brunnerMunzelTest()`](brunnerMunzelTest.md),
[`hotellingsT2Test()`](hotellingsT2Test.md),
[`moodMedianTest()`](moodMedianTest.md), [`signTest()`](signTest.md),
[`tTestA()`](tTestA.md), [`vanWaerdenTest()`](vanWaerdenTest.md),
[`zTest()`](zTest.md)

## Examples

``` r
x <- rnorm(25, 100, 5)
yuenTTest(x, mu = 99)
#> 
#>  Yuen One-Sample Trimmed Mean t-test
#> 
#> data:  x
#> t = 0.74889, df = 14.0, trim = 0.2, p-value = 0.4663
#> alternative hypothesis: true trimmed mean difference is not equal to 99
#> 95 percent confidence interval:
#>  -2.115739  4.385920
#> sample estimates:
#> trimmed mean of x 
#>          100.1351 
#> 

with(sleep,
     yuenTTest(extra[group == 1],
               extra[group == 2]))
#> 
#>  Yuen Two-Sample Trimmed Mean t-test
#> 
#> data:  extra[group == 1] and extra[group == 2]
#> t = -1.5314, df = 8.7502, trim = 0.2000, p-value = 0.161
#> alternative hypothesis: true trimmed mean difference is not equal to 0
#> 95 percent confidence interval:
#>  -1.939437  3.006103
#> sample estimates:
#> trimmed mean of x trimmed mean of y 
#>         0.5333333         2.2000000 
#> 

yuenTTest(extra ~ group, data = sleep)
#> 
#>  Yuen Two-Sample Trimmed Mean t-test
#> 
#> data:  extra ~ group
#> t = -1.5314, df = 8.7502, trim = 0.2000, p-value = 0.161
#> alternative hypothesis: true trimmed mean difference is not equal to 0
#> 95 percent confidence interval:
#>  -1.939437  3.006103
#> sample estimates:
#> trimmed mean of x trimmed mean of y 
#>         0.5333333         2.2000000 
#> 
```
