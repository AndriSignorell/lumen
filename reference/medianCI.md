# Confidence Interval for the Median

Calculate the confidence interval for the median.

## Usage

``` r
medianCI(
  x,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("exact", "boot"),
  na.rm = FALSE,
  ...
)
```

## Arguments

- x:

  a (non-empty) numeric vector of data values.

- conf.level:

  confidence level of the interval

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `"left"` would be analogue to a
  hypothesis of `"greater"` in a `t.test`.

- method:

  defining the type of interval that should be calculated (one out of
  `"exact"`, `"boot"`). Default is `"exact"`. See Details.

- na.rm:

  logical. Should missing values be removed? Defaults to `FALSE`.

- ...:

  the dots are passed on to
  [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html). In particular,
  the type of bootstrap confidence interval can be defined via this. The
  defaults are `R=999` and `type="perc"`.

## Value

A named numeric vector with elements:

- `est`:

  point estimate

- `lci`:

  lower confidence interval bound

- `uci`:

  upper confidence interval bound

## Details

The `"exact"` method is the way SAS is said to calculate the confidence
interval. This is also implemented in
[`signTest`](https://andrisignorell.github.io/lumen/reference/signTest.md).
The boot confidence interval type is calculated by means of
[`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html) with default type
`"perc"`.  
Use [`sapply`](https://rdrr.io/r/base/lapply.html),
resp.[`apply`](https://rdrr.io/r/base/apply.html), to get the confidence
intervals from a data.frame or from a matrix.

## See also

[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html),
[`median`](https://rdrr.io/r/stats/median.html),
`DescToolsX::hodgesLehmann`

Other ci.location:
[`meanCI()`](https://andrisignorell.github.io/lumen/reference/meanCI.md),
[`meanCIn()`](https://andrisignorell.github.io/lumen/reference/meanCIn.md),
[`meanDiffCI()`](https://andrisignorell.github.io/lumen/reference/meanDiffCI.md),
[`quantileCI()`](https://andrisignorell.github.io/lumen/reference/quantileCI.md),
[`sumCI()`](https://andrisignorell.github.io/lumen/reference/sumCI.md)

## Examples

``` r

set.seed(448)
x <- c(rnorm(100), NA)

medianCI(x, na.rm=TRUE)
#>      median         lci         uci 
#> -0.01987842 -0.34886600  0.22967354 
#> attr(,"conf.level")
#> [1] 0.9647998
medianCI(x, conf.level=0.99, na.rm=TRUE)
#>      median         lci         uci 
#> -0.01987842 -0.40208643  0.32829356 
#> attr(,"conf.level")
#> [1] 0.9933629

medianCI(x, na.rm=TRUE, method="exact")
#>      median         lci         uci 
#> -0.01987842 -0.34886600  0.22967354 
#> attr(,"conf.level")
#> [1] 0.9647998
medianCI(x, na.rm=TRUE, method="boot")
#>      median         lci         uci 
#> -0.01987842 -0.33366924  0.17706651 
 
x <- x[!is.na(x)]
medianCI(x, method="boot")
#>      median         lci         uci 
#> -0.01987842 -0.33366924  0.18475037 

# ... the same as
medianCI(x, method="boot", type="bca")
#>      median         lci         uci 
#> -0.01987842 -0.31847248  0.17043184 

medianCI(x, method="boot", type="basic")
#>      median         lci         uci 
#> -0.01987842 -0.24325515  0.28676766 
medianCI(x, method="boot", type="perc")
#>      median         lci         uci 
#> -0.01987842 -0.31847248  0.22967354 
medianCI(x, method="boot", type="norm", R=499)
#>      median         lci         uci 
#> -0.01987842 -0.26235685  0.24711660 
# not supported:
medianCI(x, method="boot", type="stud")
#> Warning: bootstrap type 'stud' is not supported
#>      median         lci         uci 
#> -0.01987842          NA          NA 

medianCI(x, method="boot", sides="right")
#>      median         lci         uci 
#> -0.01987842        -Inf  0.13130127 

```
