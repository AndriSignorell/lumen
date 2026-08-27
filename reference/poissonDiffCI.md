# Confidence Interval for the Difference Between Two Poisson Rates

Estimates the difference between two independent Poisson event rates and
calculates a confidence interval.

## Usage

``` r
poissonDiffCI(
  x1,
  n1 = 1,
  x2,
  n2 = 1,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("mover", "wald")
)
```

## Arguments

- x1:

  non-negative integer event count or vector of counts for the first
  sample

- n1:

  positive exposure associated with `x1`, such as observation time,
  person-time, or population at risk; may be a vector and defaults to 1

- x2:

  non-negative integer event count or vector of counts for the second
  sample

- n2:

  positive exposure associated with `x2`; may be a vector and defaults
  to 1

- conf.level:

  numeric confidence level between 0 and 1; defaults to 0.95

- sides:

  type of confidence interval: `"two.sided"`, `"left"`, or `"right"`;
  may be abbreviated

- method:

  method used to calculate the confidence interval: `"mover"` or
  `"wald"`; may be abbreviated and defaults to `"mover"`

## Value

If the arguments identify a single result, a named numeric vector with
elements:

- `est`:

  estimated rate difference

- `lci`:

  lower confidence bound

- `uci`:

  upper confidence bound

Otherwise, a `data.frame` containing these three columns followed by the
recycled values of `x1`, `n1`, `x2`, `n2`, and `conf.level`.

## Details

The function assumes two independent counts \$\$X_i \sim
\mathrm{Poisson}(n_i\lambda_i), \quad i = 1, 2.\$\$ The parameter of
interest is the rate difference \\\Delta = \lambda_1 - \lambda_2\\,
estimated by \\\hat{\Delta} = x_1/n_1 - x_2/n_2\\.

The available confidence-interval methods are:

- `"mover"`:

  the method of variance estimates recovery (MOVER), which combines
  separate exact Garwood limits for the two rates

- `"wald"`:

  the normal-approximation interval based on the standard error
  \\\sqrt{x_1/n_1^2 + x_2/n_2^2}\\

The Wald interval has zero width when both counts are zero.

For `sides = "left"`, the function returns a lower one-sided confidence
bound and sets `uci` to `Inf`. For `sides = "right"`, it returns an
upper one-sided confidence bound and sets `lci` to `-Inf`.

The numeric arguments are recycled to a common length only when their
lengths are compatible whole multiples. Incompatible lengths produce an
error. `sides` and `method` must each identify a single choice.

## References

Zou, G. Y. and Donner, A. (2008). Construction of confidence limits
about effect measures: a general approach. *Statistics in Medicine*,
**27**(10), 1693–1702.

## See also

[`poissonCI()`](poissonCI.md), [`poissonRatioCI()`](poissonRatioCI.md),
[`binomDiffCI()`](binomDiffCI.md)

## Examples

``` r
# 15 events in 100 person-years compared with
# 6 events in 120 person-years
poissonDiffCI(15, 100, 6, 120)
#>        est        lci        uci 
#> 0.10000000 0.01155263 0.20241565 

poissonDiffCI(15, 100, 6, 120, method = "wald")
#>        est        lci        uci 
#> 0.10000000 0.01419326 0.18580674 

# A 95% lower confidence bound for the rate difference
poissonDiffCI(15, 100, 6, 120, sides = "left")
#>        est        lci        uci 
#> 0.10000000 0.02462852        Inf 

# Recycling returns one row per comparison
poissonDiffCI(x1 = c(15, 20), n1 = 100, x2 = 6, n2 = 120)
#>    est        lci       uci x1  n1 x2  n2 conf.level
#> 1 0.10 0.01155263 0.2024156 15 100  6 120       0.95
#> 2 0.15 0.05243411 0.2633907 20 100  6 120       0.95
```
