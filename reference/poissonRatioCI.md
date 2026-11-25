# Confidence Intervals for the Ratio of Two Poisson Rates

`poissonRatioCI()` computes confidence intervals for the ratio of two
independent Poisson event rates.

## Usage

``` r
poissonRatioCI(
  x1,
  n1 = 1,
  x2,
  n2 = 1,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("exact", "mid-p", "wald-log")
)
```

## Arguments

- x1:

  number of events in the first sample, a non-negative integer.

- n1:

  exposure associated with `x1`, such as observation time, person-time
  or population at risk; a positive number, defaults to 1.

- x2:

  number of events in the second sample, a non-negative integer.

- n2:

  exposure associated with `x2`; a positive number, defaults to 1.

- conf.level:

  confidence level, defaults to 0.95.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `sides` names the side carrying the
  finite bound: `"left"` reports the lower limit and opens the upper one
  to `Inf`, `"right"` reports the upper limit and opens the lower one to
  0, the lower limit of the parameter space. A one-sided bound puts the
  whole `1 - conf.level` into its single tail.

- method:

  character string specifying which method to use; this can be one out
  of: `"exact"` (default), `"mid-p"` and `"wald-log"`. Abbreviation of
  method is accepted. See details.

## Value

If recycling yields a single case, a named numeric vector with elements:

- `est`:

  point estimate of the rate ratio.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

If recycling yields multiple cases, a data frame with one row per case
is returned. Its first three columns are `est`, `lci`, and `uci`; the
remaining columns contain the recycled values of `x1`, `n1`, `x2`, `n2`,
and `conf.level`.

## Details

The numeric arguments are recycled to a common length only when their
lengths are compatible whole multiples; incompatible lengths produce an
error. `sides` and `method` must each identify a single choice.

The function assumes two independent counts \\X_i \sim
\mathrm{Poisson}(n_i\lambda_i),\\ i = 1, 2\\. The parameter of interest
is the rate ratio \\\theta = \lambda_1/\lambda_2\\, estimated by
\\\hat{\theta} = (x_1/n_1)/(x_2/n_2)\\.

**Exact** (default): The exact conditional interval, obtained by
conditioning on \\x_1 + x_2\\ and transforming a Clopper-Pearson
interval for the resulting binomial probability. This is the
construction used by
[`poisson.test()`](https://rdrr.io/r/stats/poisson.test.html) for two
samples.

**Mid-p**: The corresponding conditional mid-p interval, designed to
reduce the conservatism of the exact interval. It is generally shorter
but does not guarantee the nominal coverage.

**Wald-log**: The asymptotic Wald interval, symmetric on the logarithmic
scale. It cannot be calculated when either count is zero and is then
returned as \\\[0, \infty\]\\.

If both counts are zero, the rate ratio is undefined; the function
returns `NA` for `est` and \\\[0, \infty\]\\ for the confidence interval
for every method. If only \\x_2\\ is zero, the point estimate is `Inf`.

## **Which interval should be used?**

The exact interval guarantees the nominal coverage but is conservative.
The mid-p interval is a common compromise with coverage close to the
nominal level (Graham et al., 2003). The log-Wald interval requires
reasonably large counts in both samples and is not available when either
count is zero.

## References

Graham, P. L., Mengersen, K. and Morton, A. P. (2003) Confidence limits
for the ratio of two rates based on likelihood scores. *Statistics in
Medicine*, **22**(12), pp. 2071-2083.

Sahai, H. and Khurshid, A. (1993) Confidence intervals for the ratio of
two Poisson means. *The Mathematical Scientist*, **18**, pp. 43-50.

## See also

[`poissonCI()`](poissonCI.md), [`poissonDiffCI()`](poissonDiffCI.md),
[`binomRatioCI()`](binomRatioCI.md),
[`poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)

## Examples

``` r

# 15 events in 100 person-years compared with
# 6 events in 120 person-years
poissonRatioCI(15, 100, 6, 120)
#>      est      lci      uci 
#> 3.000000 1.099947 9.437411 

# the exact interval corresponds to the two-sample conditional test
poisson.test(c(15, 6), c(100, 120))$conf.int
#> [1] 1.099947 9.437411
#> attr(,"conf.level")
#> [1] 0.95

poissonRatioCI(15, 100, 6, 120, method = "mid-p")
#>      est      lci      uci 
#> 3.000000 1.188883 8.415485 

# zero counts are handled explicitly
poissonRatioCI(0, 100, 6, 120)
#>      est      lci      uci 
#> 0.000000 0.000000 1.019173 
```
