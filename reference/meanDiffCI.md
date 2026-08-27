# Confidence Interval For Difference of Means

Calculates the confidence interval for the difference of two means
either the classical way or with the bootstrap approach.

## Usage

``` r
meanDiffCI(
  x,
  y,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("classic", "boot"),
  paired = FALSE,
  var.equal = FALSE,
  na.rm = FALSE,
  ...
)
```

## Arguments

- x:

  a (non-empty) numeric vector of data values.

- y:

  a (non-empty) numeric vector of data values.

- conf.level:

  confidence level of the interval.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `"left"` would be analogue to a
  hypothesis of `"greater"` in a `t.test`.

- method:

  a vector of character strings representing the type of intervals
  required. The value should be any subset of the values `"classic"`,
  `"boot"`. Bootstrap type can be provided by the dots. See
  [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html).

- paired:

  a logical indicating whether you want confidence intervals for a
  paired design. Defaults to `FALSE`.

- var.equal:

  a logical variable indicating whether to treat the two variances as
  being equal. Default is `FALSE`. If `TRUE` then the pooled variance is
  used to estimate the variance otherwise the Welch (or Satterthwaite)
  approximation to the degrees of freedom is used. Passed on to
  [`t.test()`](https://rdrr.io/r/stats/t.test.html).

- na.rm:

  logical. Should missing values be removed? Defaults to `FALSE`.

- ...:

  further arguments, can be used to provide further arguments to the
  boot function.

## Value

A named numeric vector with elements:

- `meandiff`:

  point estimate, the difference: mean(x) - mean(y)

- `lci`:

  lower confidence interval bound

- `uci`:

  upper confidence interval bound

## Details

This function collects code from two sources. The classical confidence
interval is calculated by means of
[`t.test()`](https://rdrr.io/r/stats/t.test.html). The bootstrap
intervals are strongly based on the example in
[`boot`](https://rdrr.io/pkg/boot/man/boot.html).

The bootstrap type `"stud"` (studentized) is not supported: the
statistic functions used here return only the point estimate, not a
per-replicate variance estimate, so requesting it raises an error.

## See also

[`varCI()`](varCI.md),
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html)

Other ci.location: [`meanCI()`](meanCI.md), [`meanCIn()`](meanCIn.md),
[`medianCI()`](medianCI.md), [`quantileCI()`](quantileCI.md),
[`sumCI()`](sumCI.md)

## Examples

``` r

x <- mtcars[mtcars$am == 0, "mpg"]
y <- mtcars[mtcars$am == 1, "mpg"]

meanDiffCI(x, y, na.rm=TRUE)
#>   meandiff        lci        uci 
#>  -7.244939 -11.280194  -3.209684 
meanDiffCI(x, y, conf.level=0.99, na.rm=TRUE)
#>   meandiff        lci        uci 
#>  -7.244939 -12.769128  -1.720751 

# the different types of bootstrap confints
meanDiffCI(x, y, method="boot", type="norm", na.rm=TRUE)
#>   meandiff        lci        uci 
#>  -7.244939 -10.852192  -3.676966 
meanDiffCI(x, y, method="boot", type="basic", na.rm=TRUE)
#>   meandiff        lci        uci 
#>  -7.244939 -10.595951  -3.658704 
# type="stud" is not supported (see Details) and raises an error:
# meanDiffCI(x, y, method="boot", type="stud", na.rm=TRUE)
meanDiffCI(x, y, method="boot", type="perc", na.rm=TRUE)
#>   meandiff        lci        uci 
#>  -7.244939 -11.010121  -3.540486 
meanDiffCI(x, y, method="boot", type="bca", na.rm=TRUE)
#>   meandiff        lci        uci 
#>  -7.244939 -10.929270  -3.724819 

# for long form variables
with(mtcars, with(split(mpg, am), 
  meanDiffCI(`0`, `1`) )
)
#>   meandiff        lci        uci 
#>  -7.244939 -11.280194  -3.209684 
```
