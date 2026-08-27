# Confidence Interval for a Variance

Computes confidence intervals for a population variance using classical
chi-square, Bonett, or bootstrap methods.

## Usage

``` r
varCI(
  x,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("classic", "bonett", "boot"),
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
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `"left"` would be analogue to a
  hypothesis of `"greater"` in a `t.test`.

- method:

  vector of character strings representing the type of intervals
  required. The value should be any subset of the values `"classic"`,
  `"bonett"`, `"norm"`, `"boot"`. Bootstrap type can be set by the ...
  arguments. See [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html).

- na.rm:

  logical. Should missing values be removed? Defaults to FALSE.

- ...:

  further arguments, can be used to provide further arguments to the
  boot function.

## Value

A named numeric vector with elements:

- `est`:

  point estimate.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

## Details

The confidence interval for the variance is very sensitive to
non-normality in the data. Bonett (2006) has proposed an interval that
is nearly exact when the data is normally distributed and provides good
performance for moderately non-normal data. See the references for the
details.

## References

Bonett (2006) Approximate Confidence Interval for Standard Deviation of
Nonnormal Distributions, *Computational Statistics and Data Analysis*,
Vol. 50, pp. 775 - 782.  
https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/sdconfli.htm
(might be outdated)

## See also

[`meanCI`](https://andrisignorell.github.io/lumen/reference/meanCI.md),
[`medianCI`](https://andrisignorell.github.io/lumen/reference/medianCI.md),
[`varTest`](https://andrisignorell.github.io/lumen/reference/varTest.md),
`DescToolsX::varX`

Other test.variance:
[`leveneTest()`](https://andrisignorell.github.io/lumen/reference/leveneTest.md),
[`mosesTest()`](https://andrisignorell.github.io/lumen/reference/mosesTest.md),
[`siegelTukeyTest()`](https://andrisignorell.github.io/lumen/reference/siegelTukeyTest.md),
[`varTest()`](https://andrisignorell.github.io/lumen/reference/varTest.md)

## Examples

``` r
x <- mtcars$mpg

varCI(x, na.rm=TRUE)
#>      est      lci      uci 
#> 36.32410 23.34653 64.20343 
varCI(x, conf.level=0.99, na.rm=TRUE)
#>      est      lci      uci 
#> 36.32410 20.47258 77.88527 

x <- c(14.816, 14.863, 14.814, 14.998, 14.965, 14.824, 14.884, 14.838, 14.916,
       15.021, 14.874, 14.856, 14.860, 14.772, 14.980, 14.919)
varCI(x, conf.level=0.9)
#>         est         lci         uci 
#> 0.005285333 0.003171734 0.010918691 

# and for the standard deviation
sqrt(varCI(x, conf.level=0.9))
#>        est        lci        uci 
#> 0.07270030 0.05631815 0.10449254 


# from Bonett's paper
# expected results:
# ------------------------------------
#  conf.lvl       sd      lci      uci
# ------------------------------------
#      90.0   0.5168   0.3592   0.9359
#      95.0   0.5168   0.3263   1.0841
#      99.0   0.5168   0.2607   1.5109

p <- c(15.83, 16.01, 16.24, 16.42, 15.33, 15.44, 16.88, 16.31)
sqrt(varCI(p, method="bonett", conf.level=0.9))
#>       est       lci       uci 
#> 0.5167965 0.3592151 0.9359420 
sqrt(varCI(p, method="bonett"))
#>       est       lci       uci 
#> 0.5167965 0.3263123 1.0840670 
sqrt(varCI(p, method="bonett", conf.level=0.99))
#>       est       lci       uci 
#> 0.5167965 0.2607127 1.5108922 

# some bootstrap intervals
varCI(x, method="boot", type="norm")
#>         est         lci         uci 
#> 0.005285333 0.002935139 0.008293970 
varCI(x, method="boot", type="perc")
#>         est         lci         uci 
#> 0.005285333 0.002081296 0.007595129 
varCI(x, method="boot", type="bca")
#>         est         lci         uci 
#> 0.005285333 0.003153558 0.009063481 
```
