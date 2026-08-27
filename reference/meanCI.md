# Confidence Intervals for the Mean

Collection of several approaches to determine confidence intervals for
the mean. Both, the classical way and bootstrap intervals are
implemented for both, normal and trimmed means.

## Usage

``` r
meanCI(
  x,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("classic", "boot"),
  sd = NULL,
  trim = 0,
  na.rm = FALSE,
  ...
)
```

## Arguments

- x:

  a (non-empty) numeric vector of data values.

- conf.level:

  confidence level of the interval.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`.
  `"left"` would be analogue to a hypothesis of `"greater"` in a
  `t.test`. You can specify just the initial letter.

- method:

  A vector of character strings representing the type of intervals
  required. The value should be any subset of the values `"classic"`,
  `"boot"`. See [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html).

- sd:

  the standard deviation of x. If provided it's interpreted as sd of the
  population and the normal quantiles will be used for constructing the
  confidence intervals. If left to `NULL` (default) the sample `sd(x)`
  will be calculated and used in combination with the t-distribution.

- trim:

  the fraction (0 to 0.5) of observations to be trimmed from each end of
  `x` before the mean is computed. Values of `trim` outside that range
  are taken as the nearest endpoint.

- na.rm:

  a logical value indicating whether `NA` values should be stripped
  before the computation proceeds. Defaults to FALSE.

- ...:

  further arguments are passed to the
  [`boot`](https://rdrr.io/pkg/boot/man/boot.html) function. Supported
  arguments are `type` (`"norm"`, `"basic"`, `"stud"`, `"perc"`,
  `"bca"`), `parallel` and the number of bootstrap replicates `R`. If
  not defined those will be set to their defaults, being `"basic"` for
  `type`, option `"boot.parallel"` (and if that is not set, `"no"`) for
  `parallel` and `999` for `R`.

## Value

A named numeric vector with elements:

- `est`:

  point estimate

- `lci`:

  lower confidence interval bound

- `uci`:

  upper confidence interval bound

## Details

The confidence intervals for the trimmed means use winsorized variances
as described in the references.

The bootstrap type `"stud"` (studentized) requires a variance estimate
to be returned alongside the point estimate on every bootstrap
replicate; this is supported for both the trimmed and untrimmed mean.

## References

Wilcox, R. R., Keselman H. J. (2003) Modern robust data analysis
methods: measures of central tendency *Psychol Methods*, 8(3):254-74

Wilcox, R. R. (2005) *Introduction to robust estimation and hypothesis
testing* Elsevier Academic Press

## See also

`DescToolsX::meanX()`,
[`t.test()`](https://rdrr.io/r/stats/t.test.html), [`varCI()`](varCI.md)

Other ci.location: [`meanCIn()`](meanCIn.md),
[`meanDiffCI()`](meanDiffCI.md), [`medianCI()`](medianCI.md),
[`quantileCI()`](quantileCI.md), [`sumCI()`](sumCI.md)

## Examples

``` r

x <- mtcars$mpg

meanCI(x, na.rm=TRUE)
#>      est      lci      uci 
#> 20.09062 17.91768 22.26357 
meanCI(x, conf.level=0.99, na.rm=TRUE)
#>      est      lci      uci 
#> 20.09062 17.16706 23.01419 

meanCI(x, sides="left", na.rm=TRUE)
#>      est      lci      uci 
#> 20.09062 18.28418      Inf 
# same as:
t.test(x, alternative="greater")
#> 
#>  One Sample t-test
#> 
#> data:  x
#> t = 18.857, df = 31, p-value < 2.2e-16
#> alternative hypothesis: true mean is greater than 0
#> 95 percent confidence interval:
#>  18.28418      Inf
#> sample estimates:
#> mean of x 
#>  20.09062 
#> 

meanCI(x, sd=25, na.rm=TRUE)
#>      est      lci      uci 
#> 20.09062 11.42873 28.75252 

# the different types of bootstrap confints
meanCI(x, method="boot", type="norm", na.rm=TRUE)
#>      est      lci      uci 
#> 20.09062 18.00729 22.19484 
meanCI(x, trim=0.1, method="boot", type="norm", na.rm=TRUE)
#>      est      lci      uci 
#> 19.69615 17.55049 21.84337 
meanCI(x, trim=0.1, method="boot", type="basic", na.rm=TRUE)
#>      est      lci      uci 
#> 19.69615 17.24615 21.81154 
meanCI(x, trim=0.1, method="boot", type="stud", na.rm=TRUE)
#>      est      lci      uci 
#> 19.69615 17.56002 22.71807 
meanCI(x, trim=0.1, method="boot", type="perc", na.rm=TRUE)
#>      est      lci      uci 
#> 19.69615 17.70385 21.83462 
meanCI(x, trim=0.1, method="boot", type="bca", na.rm=TRUE)
#>      est      lci      uci 
#> 19.69615 17.87268 22.25566 

meanCI(x, trim=0.1, method="boot", type="bca", R=1999, na.rm=TRUE)
#>      est      lci      uci 
#> 19.69615 17.82699 22.22681 

# Getting the meanCI for more than 1 column
round(t(sapply(mtcars[, c("mpg", "hp")], meanCI, na.rm=TRUE)), 3)
#>         est     lci     uci
#> mpg  20.091  17.918  22.264
#> hp  146.688 121.968 171.407
```
