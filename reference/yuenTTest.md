# Yuen T-Test for Robust Comparison of Trimmed Means

Performs Yuen's robust t-test for trimmed means. Compared with the
classical t-test, the procedure is substantially less sensitive to
outliers, heavy tails, and moderate departures from normality.

The test is based on:

- trimmed means,

- winsorized variances,

- Welch-type degrees of freedom.

For paired tests, trimming is performed on the paired differences, i.e.
the one-sample trimmed t-test (Tukey & McLaughlin, 1963) is applied to
\\x - y\\. This tests the trimmed mean of the differences, which in
general is not the difference of the trimmed means; the latter is what
e.g. `WRS2::yuend()` compares.

## Usage

``` r
yuenTTest(x, ...)

# S3 method for class 'formula'
yuenTTest(formula, data, subset, na.action = na.pass, paired = FALSE, ...)

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

  numeric vector of observations. Non-finite values (`NA`, `NaN`, `Inf`,
  `-Inf`) are removed; in the paired case the pair is removed.

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

- paired:

  logical indicating whether a paired test is performed. Only available
  in the default method: the formula interface describes independent
  groups and does not identify pairs.

- y:

  optional second numeric vector.

- alternative:

  character string specifying the alternative hypothesis. One of
  `"two.sided"`, `"less"`, or `"greater"`.

- mu:

  hypothesized trimmed mean (or trimmed mean difference).

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

**Winsorizing.** With \\g = \lfloor \mathrm{trim} \cdot n \rfloor\\, the
\\g\\ smallest observations are set to the \\(g+1)\\-th order statistic
and the \\g\\ largest to the \\(n-g)\\-th, as in Yuen (1974) and Wilcox
(2005). The winsorized variance must be taken this way, at order
statistics rather than at interpolated quantiles: the standard error
\\\sqrt{(n-1) s_w^2 / (h(h-1))}\\ with \\h = n - 2g\\ is derived for
exactly \\g\\ replaced values in each tail, the same \\g\\ observations
that [`mean()`](https://rdrr.io/r/base/mean.html) with `trim` removes.
The results agree with `WRS2::yuen()` and `PairedData::yuen.t.test()`.

**Standard error.** In all three designs the squared standard error of a
trimmed mean is \\(n-1) s_w^2 / (h(h-1))\\, with \\s_w^2\\ the
winsorized variance and \\h\\ the number of observations left after
trimming; the degrees of freedom are \\h - 1\\ (combined by Welch's
formula in the two-sample case). With `trim = 0`, or whenever \\g = 0\\,
the three tests reduce exactly to the corresponding
[`t.test()`](https://rdrr.io/r/stats/t.test.html): one-sample, paired,
and Welch. The one-sample version in `WRS2::trimse()` uses the
asymptotically equivalent \\s_w / ((1 - 2\\\mathrm{trim}) \sqrt{n})\\;
it differs in small samples, where it inflates the standard error even
if no observation is trimmed (e.g. \\n = 4\\, `trim = 0.2`).

The confidence interval is for the estimated parameter itself (the
trimmed mean, or the difference of trimmed means), independent of `mu`.

## References

Wilcox, R. R. (2005). *Introduction to Robust Estimation and Hypothesis
Testing*. Academic Press.

Tukey, J. W., & McLaughlin, D. H. (1963). Less vulnerable confidence and
significance procedures for location based on a single sample:
trimming/winsorization 1. *Sankhya A*, 25, 331–352.

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
#> t = 0.76501, df = 14.0, trim = 0.2, p-value = 0.457
#> alternative hypothesis: true trimmed mean is not equal to 99
#> 95 percent confidence interval:
#>   96.95275 103.31743
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
#> t = -1.6168, df = 8.2647, trim = 0.2000, p-value = 0.1434
#> alternative hypothesis: true trimmed mean difference is not equal to 0
#> 95 percent confidence interval:
#>  -4.0306400  0.6973066
#> sample estimates:
#> trimmed mean of x trimmed mean of y 
#>         0.5333333         2.2000000 
#> 

yuenTTest(extra ~ group, data = sleep)
#> 
#>  Yuen Two-Sample Trimmed Mean t-test
#> 
#> data:  extra ~ group
#> t = -1.6168, df = 8.2647, trim = 0.2000, p-value = 0.1434
#> alternative hypothesis: true trimmed mean difference is not equal to 0
#> 95 percent confidence interval:
#>  -4.0306400  0.6973066
#> sample estimates:
#> trimmed mean of x trimmed mean of y 
#>         0.5333333         2.2000000 
#> 
```
