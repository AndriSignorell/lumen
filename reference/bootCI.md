# Simple Bootstrap Confidence Intervals

Convenience wrapper for calculating bootstrap confidence intervals for
univariate and bivariate statistics.

## Usage

``` r
bootCI(
  x,
  y = NULL,
  FUN,
  ...,
  bci.method = c("norm", "basic", "stud", "perc", "bca"),
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  R = 999
)
```

## Arguments

- x:

  a (non-empty) numeric vector of data values.

- y:

  NULL (default) or a vector with compatible dimensions to `x`, when a
  bivariate statistic is used.

- FUN:

  the function to be used.

- ...:

  further arguments are passed to the function `FUN`.

- bci.method:

  a vector of character strings representing the type of intervals
  required. The value should be any subset of the values `"norm"`,
  `"basic"`, `"stud"`, `"perc"`, `"bca"`, as it is passed on as `method`
  to [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html).

- conf.level:

  confidence level of the interval.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `"left"` would be analogue to a
  hypothesis of `"greater"` in a `t.test`.

- R:

  number of bootstrap replicates. Usually this will be a single positive
  integer. For importance resampling, some resamples may use one set of
  weights and others use a different set of weights. In this case `R`
  would be a vector of integers where each component gives the number of
  resamples from each of the rows of weights.

## Value

A named numeric vector with three elements:

- `est`:

  the estimate calculated by `FUN`.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

## Examples

``` r

set.seed(1984)
bootCI(mtcars$mpg, FUN=mean, na.rm=TRUE, bci.method="basic")
#>      est      lci      uci 
#> 20.09062 17.99062 22.19687 
bootCI(mtcars$mpg, FUN=mean, trim=0.1, na.rm=TRUE, bci.method="basic")
#>      est      lci      uci 
#> 19.69615 17.46923 21.78846 

# bootCI(mtcars$mpg, FUN=DescToolsX::skewX, na.rm=TRUE, bci.method="basic")

# bootCI(Pizza$operator, Pizza$area, FUN=cramerV)

spearman <- function(x,y) cor(x, y, method="spearman", use="p")
bootCI(mtcars$mpg, mtcars$hp, FUN=spearman)
#>        est        lci        uci 
#> -0.8946646 -0.9981280 -0.8188144 


```
