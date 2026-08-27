# Sign Test for Testing a Median or Paired Median Difference

A nonparametric test for the median of a single sample or the median
difference of paired samples, based solely on the signs of deviations
from the hypothesized median.

## Usage

``` r
signTest(
  x,
  y = NULL,
  alternative = c("two.sided", "less", "greater"),
  mu = 0,
  conf.level = 0.95
)
```

## Arguments

- x:

  numeric vector of data values. Non-finite (e.g. infinite or missing)
  values will be omitted.

- y:

  an optional numeric vector of data values: as with x non-finite values
  will be omitted.

- alternative:

  a character string, one of `"greater"`, `"less"`, or `"two.sided"`, or
  the initial letter of each, indicating the specification of the
  alternative hypothesis. For one-sample tests, `alternative` refers to
  the true median of the parent population in relation to the
  hypothesized value of the median.

- mu:

  a number specifying an optional parameter used to form the null
  hypothesis. See Details.

- conf.level:

  confidence level for the returned confidence interval, restricted to
  lie between zero and one.

## Value

A list of class `htest`, containing the following components:

- statistic:

  the S-statistic (the number of positive differences between the data
  and the hypothesized median), with names attribute “S”.

- parameter:

  the total number of valid differences.

- p.value:

  the p-value for the test.

- null.value:

  is the value of the median specified by the null hypothesis. This
  equals the input argument `mu`.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  the type of test applied.

- data.name:

  a character string giving the names of the data.

- conf.int:

  a confidence interval for the median.

- estimate:

  the sample median.

## Details

Performs one- and two-sample sign tests on vectors of data.

There is no formula interface: unlike
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)'s `paired`
argument, a one-sample or paired-samples design has no natural mapping
to `response ~ group` formula semantics (the pairing order is not
expressible that way), so none is offered here.

`signTest` computes a “Dependent-samples Sign-Test” if both `x` and `y`
are provided. If only `x` is provided, the “One-sample Sign-Test” will
be computed.

For the one-sample sign-test, the null hypothesis is that the median of
the population from which `x` is drawn is `mu`. For the two-sample
dependent case, the null hypothesis is that the median for the
differences of the populations from which `x` and `y` are drawn is `mu`.
The alternative hypothesis indicates the direction of divergence of the
population median for `x` from `mu` (i.e., `"greater"`, `"less"`,
`"two.sided"`.)

The confidence levels are exact.

## References

Gibbons, J.D. and Chakraborti, S. (1992): *Nonparametric Statistical
Inference*. Marcel Dekker Inc., New York.

Kitchens, L. J. (2003): *Basic Statistics and Data Analysis*. Duxbury.

Conover, W. J. (1980): *Practical Nonparametric Statistics, 2nd ed*.
Wiley, New York.

## See also

[`t.test()`](https://rdrr.io/r/stats/t.test.html),
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html),
[`zTest`](zTest.md),
[`binom.test`](https://rdrr.io/r/stats/binom.test.html), `SIGN.test` in
the package BSDA (reporting approximative confidence intervals).

Other test.location: [`brunnerMunzelTest()`](brunnerMunzelTest.md),
[`hotellingsT2Test()`](hotellingsT2Test.md),
[`moodMedianTest()`](moodMedianTest.md), [`tTestA()`](tTestA.md),
[`vanWaerdenTest()`](vanWaerdenTest.md), [`yuenTTest()`](yuenTTest.md),
[`zTest()`](zTest.md)

## Examples

``` r

x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

signTest(x, y)
#> 
#>  Dependent-samples Sign-Test
#> 
#> data:  x and y
#> S = 7, number of differences = 9, p-value = 0.1797
#> alternative hypothesis: true median difference is not equal to 0
#> 96.1 percent confidence interval:
#>  -0.080  0.952
#> sample estimates:
#> median of the differences 
#>                      0.49 
#> 
wilcox.test(x, y, paired = TRUE)
#> 
#>  Wilcoxon signed rank exact test
#> 
#> data:  x and y
#> V = 40, p-value = 0.03906
#> alternative hypothesis: true location shift is not equal to 0
#> 


d.light <- data.frame( 
  black = c(25.85,28.84,32.05,25.74,20.89,41.05,25.01,24.96,27.47),
  white = c(18.23,20.84,22.96,19.68,19.5,24.98,16.61,16.07,24.59),
  d     = c(7.62,8,9.09,6.06,1.39,16.07,8.4,8.89,2.88)
)

d <- d.light$d

signTest(x=d, mu = 4)
#> 
#>  One-sample Sign-Test
#> 
#> data:  d
#> S = 7, number of differences = 9, p-value = 0.1797
#> alternative hypothesis: true median is not equal to 4
#> 96.1 percent confidence interval:
#>  2.88 9.09
#> sample estimates:
#> median of the differences 
#>                         8 
#> 
wilcox.test(x=d, mu = 4, conf.int = TRUE)
#> 
#>  Wilcoxon signed rank exact test
#> 
#> data:  d
#> V = 41, p-value = 0.02734
#> alternative hypothesis: true location is not equal to 4
#> 96.1 percent confidence interval:
#>   4.505 11.845
#> sample estimates:
#> (pseudo)median 
#>           7.81 
#> 

signTest(x=d, mu = 4, alternative="less")
#> 
#>  One-sample Sign-Test
#> 
#> data:  d
#> S = 7, number of differences = 9, p-value = 0.9805
#> alternative hypothesis: true median is less than 4
#> 98 percent confidence interval:
#>  -Inf 8.89
#> sample estimates:
#> median of the differences 
#>                         8 
#> 
wilcox.test(x=d, mu = 4, conf.int = TRUE, alternative="less")
#> 
#>  Wilcoxon signed rank exact test
#> 
#> data:  d
#> V = 41, p-value = 0.9902
#> alternative hypothesis: true location is less than 4
#> 95.1 percent confidence interval:
#>  -Inf 9.09
#> sample estimates:
#> (pseudo)median 
#>           7.81 
#> 

signTest(x=d, mu = 4, alternative="greater")
#> 
#>  One-sample Sign-Test
#> 
#> data:  d
#> S = 7, number of differences = 9, p-value = 0.08984
#> alternative hypothesis: true median is greater than 4
#> 98 percent confidence interval:
#>  2.88  Inf
#> sample estimates:
#> median of the differences 
#>                         8 
#> 
wilcox.test(x=d, mu = 4, conf.int = TRUE, alternative="greater")
#> 
#>  Wilcoxon signed rank exact test
#> 
#> data:  d
#> V = 41, p-value = 0.01367
#> alternative hypothesis: true location is greater than 4
#> 95.1 percent confidence interval:
#>  5.14  Inf
#> sample estimates:
#> (pseudo)median 
#>           7.81 
#> 

with(d.light, signTest(black, white))
#> 
#>  Dependent-samples Sign-Test
#> 
#> data:  black and white
#> S = 9, number of differences = 9, p-value = 0.003906
#> alternative hypothesis: true median difference is not equal to 0
#> 96.1 percent confidence interval:
#>  2.88 9.09
#> sample estimates:
#> median of the differences 
#>                         8 
#> 
# same as:
with(d.light, signTest(black - white))
#> 
#>  One-sample Sign-Test
#> 
#> data:  black - white
#> S = 9, number of differences = 9, p-value = 0.003906
#> alternative hypothesis: true median is not equal to 0
#> 96.1 percent confidence interval:
#>  2.88 9.09
#> sample estimates:
#> median of the differences 
#>                         8 
#> 
```
