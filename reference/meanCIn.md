# Sample Size for a Given Width of a Confidence Interval for a Mean

Returns the required sample size to obtain a given width of a confidence
interval for the sample mean. The function uses
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) to find a numeric
solution. The t distribution is used.

## Usage

``` r
meanCIn(
  ci,
  sd,
  interval = c(2, 100000),
  conf.level = 0.95,
  norm = FALSE,
  tol = .Machine$double.eps^0.5
)
```

## Arguments

- ci:

  the left and right bound of the interval, which is presumed to be
  symmetric.

- sd:

  the standard deviation of the sample.

- interval:

  the interval for the sample size to be searched into, (default is c(2,
  100000)).

- conf.level:

  confidence level, defaults to `0.95`.

- norm:

  logical, determining if the t- or normaldistribution should be used.

- tol:

  the desired accuracy (convergence tolerance).

## Value

a numeric value

## Details

The required sample sizes for a specific width of confidence interval
for the mean depends recursively on the sample size, as the sample size
defines the degrees of freedom in the t-distribution. Although in most
practical cases it will be sufficient to use the normal distribution, we
might be interested in exact results.

## See also

[`binomCIn()`](binomCI.md)

Other ci.location: [`meanCI()`](meanCI.md),
[`meanDiffCI()`](meanDiffCI.md), [`medianCI()`](medianCI.md),
[`quantileCI()`](quantileCI.md), [`sumCI()`](sumCI.md)

## Examples

``` r

meanCIn(ci=c(25, 27), sd=5)
#> [1] 98.46626
```
