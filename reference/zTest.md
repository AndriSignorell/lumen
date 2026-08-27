# Z-Test for Testing Means With Known Population Standard Deviations

A parametric test for the mean of a normal distribution when the
population variance is known, or for comparing two means with known
variances, based on the standard normal distribution.

## Usage

``` r
zTest(x, ...)

# S3 method for class 'formula'
zTest(formula, data, subset, na.action = na.pass, ...)

# Default S3 method
zTest(
  x,
  y = NULL,
  alternative = c("two.sided", "less", "greater"),
  paired = FALSE,
  mu = 0,
  sd_pop,
  conf.level = 0.95,
  ...
)
```

## Arguments

- x:

  numeric vector of data values. Non-finite (e.g. infinite or missing)
  values will be omitted.

- ...:

  further arguments to be passed to or from methods.

- formula:

  a formula of the form `lhs ~ rhs` where `lhs` gives the data values
  and `rhs` a factor with two levels giving the corresponding groups.

- data:

  an optional matrix or data frame (or similar: see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in the formula `formula`. By default the variables are
  taken from `environment(formula)`.

- subset:

  an optional vector specifying a subset of observations to be used.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

- y:

  an optional numeric vector of data values: as with x non-finite values
  will be omitted.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"greater"` or `"less"`. You can specify
  just the initial letter.  
  For one-sample tests, `alternative` refers to the true mean of the
  parent population in relation to the hypothesized value of the mean.

- paired:

  a logical indicating whether you want a paired z-test.

- mu:

  a number specifying the hypothesized mean of the population.

- sd_pop:

  a number specifying the known standard deviation of the population.
  For the two-sample test, this single value is assumed to be the common
  known standard deviation of both populations.

- conf.level:

  confidence level for the interval computation.

## Value

A list with class "`htest`" containing the following components:

- statistic:

  the value of the z-statistic.

- p.value:

  the p-value for the test

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

  a character string indicating what type of test was performed.

- data.name:

  a character string giving the name(s) of the data.

## Details

Compute the test of hypothesis and compute confidence interval on the
mean of a population when the standard deviation of the population is
known.

Most introductory statistical texts introduce inference by using the
z-test and z-based confidence intervals based on knowing the population
standard deviation. However statistical packages often do not include
functions to do z-tests since the t-test is usually more appropriate for
real world situations. This function is meant to be used during that
short period of learning when the student is learning about inference
using z-procedures, but has not learned the t-based procedures yet. Once
the student has learned about the t-distribution the
[`t.test()`](https://rdrr.io/r/stats/t.test.html) function should be
used instead of this one (but the syntax is very similar, so this
function should be an appropriate introductory step to learning
[`t.test()`](https://rdrr.io/r/stats/t.test.html)).

The formula interface is only applicable for the 2-sample tests.

## References

Stahel, W. (2002) *Statistische Datenanalyse, 4th ed*, vieweg

## See also

[`t.test()`](https://rdrr.io/r/stats/t.test.html),
[`print.htest`](https://rdrr.io/r/stats/print.power.htest.html)

Other test.location:
[`brunnerMunzelTest()`](https://andrisignorell.github.io/lumen/reference/brunnerMunzelTest.md),
[`hotellingsT2Test()`](https://andrisignorell.github.io/lumen/reference/hotellingsT2Test.md),
[`moodMedianTest()`](https://andrisignorell.github.io/lumen/reference/moodMedianTest.md),
[`signTest()`](https://andrisignorell.github.io/lumen/reference/signTest.md),
[`tTestA()`](https://andrisignorell.github.io/lumen/reference/tTestA.md),
[`vanWaerdenTest()`](https://andrisignorell.github.io/lumen/reference/vanWaerdenTest.md),
[`yuenTTest()`](https://andrisignorell.github.io/lumen/reference/yuenTTest.md)

## Examples

``` r

x <- rnorm(25, 100, 5)
zTest(x, mu=99, sd_pop=5)
#> 
#>  One Sample z-test
#> 
#> data:  x
#> z = -0.51593, Std. Dev. Population = 5, p-value = 0.6059
#> alternative hypothesis: true mean is not equal to 99
#> 95 percent confidence interval:
#>   96.5241 100.4440
#> sample estimates:
#> mean of x 
#>  98.48407 
#> 

# the classic interface
with(sleep, zTest(extra[group == 1], extra[group == 2], sd_pop=2))
#> 
#>  Two Sample z-test
#> 
#> data:  extra[group == 1] and extra[group == 2]
#> z = -1.7665, Std. Dev. Population = 2, p-value = 0.07731
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -3.3330451  0.1730451
#> sample estimates:
#> mean of x mean of y 
#>      0.75      2.33 
#> 

# the formula interface
zTest(extra ~ group, data = sleep, sd_pop=2)
#> 
#>  Two Sample z-test
#> 
#> data:  extra ~ group
#> z = -1.7665, Std. Dev. Population = 2, p-value = 0.07731
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -3.3330451  0.1730451
#> sample estimates:
#> mean of x mean of y 
#>      0.75      2.33 
#> 


# Stahel (2002), pp. 186, 196

d.tyres <- data.frame(A=c(44.5,55,52.5,50.2,45.3,46.1,52.1,50.5,50.6,49.2),
                      B=c(44.9,54.8,55.6,55.2,55.6,47.7,53,49.1,52.3,50.7))
with(d.tyres, zTest(A, B, sd_pop=3, paired=TRUE))
#> 
#>  Paired z-test
#> 
#> data:  A and B
#> z = -2.4139, Std. Dev. Population = 3, p-value = 0.01578
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -4.1493851 -0.4306149
#> sample estimates:
#> mean of the differences 
#>                   -2.29 
#> 


d.oxen <- data.frame(ext=c(2.7,2.7,1.1,3.0,1.9,3.0,3.8,3.8,0.3,1.9,1.9),
                     int=c(6.5,5.4,8.1,3.5,0.5,3.8,6.8,4.9,9.5,6.2,4.1))
with(d.oxen, zTest(int, ext, sd_pop=1.8, paired=FALSE))
#> 
#>  Two Sample z-test
#> 
#> data:  int and ext
#> z = 3.9324, Std. Dev. Population = 1.8, p-value = 8.411e-05
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  1.513865 4.522498
#> sample estimates:
#> mean of x mean of y 
#>  5.390909  2.372727 
#> 
```
