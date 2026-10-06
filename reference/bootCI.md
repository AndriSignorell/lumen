# Simple Bootstrap Confidence Intervals

Convenience wrapper for calculating bootstrap confidence intervals for
univariate and bivariate statistics.

## Usage

``` r
bootCI(
  x,
  y = NULL,
  FUN,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  R = 999,
  ...
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

- conf.level:

  confidence level of the interval.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `"left"` would be analogue to a
  hypothesis of `"greater"` in a `t.test`.

- R:

  number of bootstrap replicates, a single positive whole number.

- ...:

  further arguments. The bootstrap options are taken out first, as in
  the other interval functions of the package: `type`, the interval type
  passed to
  [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html), one of
  `"bca"` (default), `"perc"`, `"basic"`, `"norm"` or `"stud"`, and
  `parallel` and `ncpus`, passed to
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html). Everything
  else is passed to `FUN`.

## Value

A named numeric vector with three elements:

- `est`:

  the estimate calculated by `FUN`.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

## Details

`type`, `parallel` and `ncpus` therefore cannot reach `FUN` through the
dots. A statistic that has an argument of one of these names - the
`type` of [`quantile()`](https://rdrr.io/r/stats/quantile.html), say -
is wrapped: `FUN = function(z) quantile(z, 0.9, type = 6)`.

`"stud"` needs a variance estimate for every replicate, which a general
`FUN` does not deliver;
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) then
returns no such interval and `bootCI()` stops with a message saying so.

## Examples

``` r

set.seed(1984)
bootCI(mtcars$mpg, FUN=mean, na.rm=TRUE)
#>      est      lci      uci 
#> 20.09062 18.12561 22.31055 
bootCI(mtcars$mpg, FUN=mean, trim=0.1, na.rm=TRUE, type="basic")
#>      est      lci      uci 
#> 19.69615 17.46923 21.78846 

# bootCI(mtcars$mpg, FUN=DescToolsX::skewX, na.rm=TRUE, type="basic")

# bootCI(Pizza$operator, Pizza$area, FUN=cramerV)

spearman <- function(x,y) cor(x, y, method="spearman", use="p")
bootCI(mtcars$mpg, mtcars$hp, FUN=spearman)
#>        est        lci        uci 
#> -0.8946646 -0.9592579 -0.7950177 


```
