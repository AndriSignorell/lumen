# Confidence Interval for the Ratio of Two Poisson Rates

Estimates the ratio of two independent Poisson event rates and
calculates a confidence interval.

## Usage

``` r
poissonRatioCI(
  x1,
  n1 = 1,
  x2,
  n2 = 1,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("exact", "midp", "wald-log")
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

  method used to calculate the confidence interval: `"exact"`, `"midp"`,
  or `"wald-log"`; may be abbreviated and defaults to `"exact"`

## Value

If the arguments identify a single result, a named numeric vector with
elements:

- `est`:

  estimated rate ratio

- `lci`:

  lower confidence bound

- `uci`:

  upper confidence bound

Otherwise, a `data.frame` containing these three columns followed by the
recycled values of `x1`, `n1`, `x2`, `n2`, and `conf.level`.

## Details

The function assumes two independent counts \$\$X_i \sim
\mathrm{Poisson}(n_i\lambda_i), \quad i = 1, 2.\$\$ The parameter of
interest is the rate ratio \\\theta = \lambda_1/\lambda_2\\, estimated
by \\\hat{\theta} = (x_1/n_1)/(x_2/n_2)\\.

The available confidence-interval methods are:

- `"exact"`:

  the exact conditional interval obtained by conditioning on \\x_1 +
  x_2\\ and transforming a Clopper–Pearson interval for the resulting
  binomial probability; this is the construction used by
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)
  for two samples

- `"midp"`:

  the corresponding conditional mid-p interval, which is generally
  shorter but does not guarantee conservative coverage

- `"wald-log"`:

  the asymptotic Wald interval, symmetric on the logarithmic scale

For `sides = "left"`, the function returns a lower one-sided confidence
bound and sets `uci` to `Inf`. For `sides = "right"`, it returns an
upper one-sided confidence bound and sets `lci` to 0, the lower limit of
the parameter space.

The log-Wald interval cannot be calculated when either count is zero and
is then returned as `[0, Inf]`. If both counts are zero, the rate ratio
and its point estimate are undefined; the function returns `NA` for
`est` and `[0, Inf]` for the confidence interval for every method. If
only `x2` is zero, the point estimate is `Inf`.

The numeric arguments are recycled to a common length only when their
lengths are compatible whole multiples. Incompatible lengths produce an
error. `sides` and `method` must each identify a single choice.

## References

Sahai, H. and Khurshid, A. (1993). Confidence intervals for the ratio of
two Poisson means. *The Mathematical Scientist*, **18**, 43–50.

Graham, P. L., Mengersen, K. and Morton, A. P. (2003). Confidence limits
for the ratio of two rates based on likelihood scores. *Statistics in
Medicine*, **22**(12), 2071–2083.

## See also

[`poissonCI()`](https://andrisignorell.github.io/lumen/reference/poissonCI.md),
[`poissonDiffCI()`](https://andrisignorell.github.io/lumen/reference/poissonDiffCI.md),
[`binomRatioCI()`](https://andrisignorell.github.io/lumen/reference/binomRatioCI.md),
[`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)

## Examples

``` r
# 15 events in 100 person-years compared with
# 6 events in 120 person-years
poissonRatioCI(15, 100, 6, 120)
#>      est      lci      uci 
#> 3.000000 1.099947 9.437411 

# The exact interval agrees with the two-sample conditional test in stats
poisson.test(c(15, 6), c(100, 120))$conf.int
#> [1] 1.099947 9.437411
#> attr(,"conf.level")
#> [1] 0.95

poissonRatioCI(15, 100, 6, 120, method = "midp")
#>      est      lci      uci 
#> 3.000000 1.188883 8.415485 

# Zero counts are handled explicitly
poissonRatioCI(0, 100, 6, 120)
#>      est      lci      uci 
#> 0.000000 0.000000 1.019173 
```
