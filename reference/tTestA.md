# Student's T-Test Based on Sample Statistics for Performing T-Tests From Summary Statistics

Performs one and two sample t-tests based on user supplied summary
information instead of data as in
[`t.test()`](https://rdrr.io/r/stats/t.test.html).

## Usage

``` r
tTestA(
  mx,
  sx,
  nx,
  my = NULL,
  sy = NULL,
  ny = NULL,
  alternative = c("two.sided", "less", "greater"),
  mu = 0,
  paired = FALSE,
  var.equal = FALSE,
  conf.level = 0.95,
  ...
)
```

## Arguments

- mx:

  a single number representing the sample mean of x.

- sx:

  a single number representing the sample standard deviation of x.

- nx:

  a single number representing the sample size of x.

- my:

  an optional single number representing the sample mean of y.

- sy:

  an optional single number representing the sample standard deviation
  of y.

- ny:

  an optional single number representing the sample size of y.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"greater"` or `"less"`. You can specify
  just the initial letter.

- mu:

  a number indicating the true value of the mean (or difference in means
  if you are performing a two sample test).

- paired:

  logical; paired tests are not supported. Use a one-sample test on the
  differences instead.

- var.equal:

  a logical variable indicating whether to treat the two variances as
  being equal. If `TRUE` then the pooled variance is used to estimate
  the variance otherwise the Welch (or Satterthwaite) approximation to
  the degrees of freedom is used.

- conf.level:

  confidence level of the interval.

- ...:

  further arguments to be passed to or from methods.

## Value

A list with class `"htest"` containing the following components:

- statistic:

  the value of the t-statistic.

- parameter:

  the degrees of freedom for the t-statistic.

- p.value:

  the p-value for the test.

- conf.int:

  a confidence interval for the mean appropriate to the specified
  alternative hypothesis.

- estimate:

  the estimated mean or difference in means depending on whether it was
  a one-sample test or a two-sample test.

- null.value:

  the specified hypothesized value of the mean or mean difference
  depending on whether it was a one-sample test or a two-sample test.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string indicating what type of t-test was performed.

- data.name:

  a character string giving the name(s) of the data.

## Details

`alternative = "greater"` is the alternative that `x` has a larger mean
than `y`.

The option `paired` is not supported here, as the variance of the
differences can't be calculated on the base of the variances of the two
samples. However, for calculating the paired test we can simply supply
the mean and standard deviation of the differences and use the
one-sample test with `mu = 0`.

If `var.equal` is `TRUE` then the pooled estimate of the variance is
used. By default, if `var.equal` is `FALSE` then the variance is
estimated separately for both groups and the Welch modification to the
degrees of freedom is used.

If the input data are effectively constant (compared to the larger of
the two means) an error is generated.

## See also

[`t.test()`](https://rdrr.io/r/stats/t.test.html)

Other test.location: [`brunnerMunzelTest()`](brunnerMunzelTest.md),
[`hotellingsT2Test()`](hotellingsT2Test.md),
[`moodMedianTest()`](moodMedianTest.md), [`signTest()`](signTest.md),
[`vanWaerdenTest()`](vanWaerdenTest.md), [`yuenTTest()`](yuenTTest.md),
[`zTest()`](zTest.md)

## Examples

``` r

## Classical example: Student's sleep data
mx <- 0.75
my <- 2.33
sx <- 1.789010
sy <- 2.002249
nx <- ny <- 10
tTestA(mx=mx, my=my, sx=sx, sy=sy, nx=nx, ny=ny)
#> 
#>  Welch Two Sample t-test
#> 
#> data:  mx and my
#> t = -1.8608, df = 17.776, p-value = 0.07939
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -3.3654835  0.2054835
#> sample estimates:
#> mean of x mean of y 
#>      0.75      2.33 
#> 

# compare to
with(sleep, t.test(extra[group == 1], extra[group == 2]))
#> 
#>  Welch Two Sample t-test
#> 
#> data:  extra[group == 1] and extra[group == 2]
#> t = -1.8608, df = 17.776, p-value = 0.07939
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -3.3654832  0.2054832
#> sample estimates:
#> mean of x mean of y 
#>      0.75      2.33 
#> 

# use the one sample test for the differences instead of paired=TRUE option
x <- with(sleep, extra[group == 1])
y <- with(sleep, extra[group == 2])

tTestA(mx=mean(x-y), sx=sd(x-y), nx=length(x-y))
#> 
#>  One Sample t-test
#> 
#> data:  mean(x - y)
#> t = -4.0621, df = 9, p-value = 0.002833
#> alternative hypothesis: true mean is not equal to 0
#> 95 percent confidence interval:
#>  -2.4598858 -0.7001142
#> sample estimates:
#> mean of x 
#>     -1.58 
#> 

# compared to 
t.test(x, y, paired = TRUE)
#> 
#>  Paired t-test
#> 
#> data:  x and y
#> t = -4.0621, df = 9, p-value = 0.002833
#> alternative hypothesis: true mean difference is not equal to 0
#> 95 percent confidence interval:
#>  -2.4598858 -0.7001142
#> sample estimates:
#> mean difference 
#>           -1.58 
#> 
```
