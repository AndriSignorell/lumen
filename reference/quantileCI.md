# Confidence Interval for Any Quantile

Calculates the confidence interval for any quantile. Although
bootstrapping might be a good approach for getting senisble confidence
intervals there's sometimes need to have a nonparameteric alternative.
This function offers one.

## Usage

``` r
quantileCI(
  x,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("exact", "boot"),
  probs = seq(0, 1, 0.25),
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
  must be one of `"two.sided"` (default), `"left"` or `"right"`
  (abbreviations allowed).  
  `"left"` would be analogue to a `"greater"` hypothesis in a `t.test`.

- method:

  defining the type of interval that should be calculated (one out of
  `"exact"`, `"boot"`). Default is `"exact"`. See Details.

- probs:

  numeric vector of probabilities with values in *`[0,1]`*. Values up to
  `2e-14` outside that range are accepted and moved to the nearby
  endpoint.

- na.rm:

  logical. Should missing values be removed? Defaults to `FALSE`.

- ...:

  bootstrap arguments can be provided by the dots argument. See
  [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html) for details.

## Value

A numeric matrix with one row per element of `probs` and columns:

- `est`:

  estimated quantile.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

For the `"exact"` method, the attribute `conf.level` reports the
achieved coverage (which may differ from the requested level).

## Details

The `"exact"` method corresponds to the way the confidence interval for
the median is calculated in SAS.  
The boot confidence interval type is calculated by means of
[`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html) with default type
`"basic"`.

## See also

`DescToolsX::quantileX`,
[`quantile()`](https://rdrr.io/r/stats/quantile.html),

Other ci.location:
[`meanCI()`](https://andrisignorell.github.io/lumen/reference/meanCI.md),
[`meanCIn()`](https://andrisignorell.github.io/lumen/reference/meanCIn.md),
[`meanDiffCI()`](https://andrisignorell.github.io/lumen/reference/meanDiffCI.md),
[`medianCI()`](https://andrisignorell.github.io/lumen/reference/medianCI.md),
[`sumCI()`](https://andrisignorell.github.io/lumen/reference/sumCI.md)

## Examples

``` r

x <- mtcars$mpg
quantileCI(x, probs=0.25, na.rm=TRUE)
#>        est  lci  uci
#> 25% 15.425 13.3 17.8
#> attr(,"conf.level")
#> [1] 0.9555407

quantileCI(x, na.rm=TRUE)
#>         est  lci  uci
#> 0%   10.400   NA 10.4
#> 25%  15.425 13.3 17.8
#> 50%  19.200 15.8 21.4
#> 75%  22.800 21.0 30.4
#> 100% 33.900 30.4   NA
#> attr(,"conf.level")
#> [1] 1.0000000 0.9555407 0.9649180 0.9555407 1.0000000
quantileCI(x, conf.level=0.99, na.rm=TRUE)
#>         est  lci  uci
#> 0%   10.400   NA 10.4
#> 25%  15.425 13.3 19.2
#> 50%  19.200 15.5 22.8
#> 75%  22.800 19.2 30.4
#> 100% 33.900 30.4   NA
#> attr(,"conf.level")
#> [1] 1.0000000 0.9912891 0.9929996 0.9912891 1.0000000

# multiple probs
quantileCI(1:100, method="exact" , probs = c(0.25, 0.75, .80, 0.95))
#>       est lci uci
#> 25% 25.75  17  34
#> 75% 75.25  67  84
#> 80% 80.20  73  89
#> 95% 95.05  90  99
#> attr(,"conf.level")
#> [1] 0.9512948 0.9512948 0.9532735 0.9514464
quantileCI(1:100, method="boot" , probs = c(0.25, 0.75, .80, 0.95))
#>       est      lci   uci
#> 25% 25.75 16.38201 33.75
#> 75% 75.25 66.18216 83.25
#> 80% 80.20 72.20000 87.20
#> 95% 95.05 90.00000 98.05

```
