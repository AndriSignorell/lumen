# Confidence Intervals for a Poisson Rate

`poissonCI()` computes confidence intervals for a Poisson event rate
using exact and approximate methods.

## Usage

``` r
poissonCI(
  x,
  n = 1,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("exact", "score", "wald", "byar")
)
```

## Arguments

- x:

  number of events, a non-negative integer.

- n:

  exposure associated with `x`, such as observation time, person-time or
  population at risk; a positive number, defaults to 1.

- conf.level:

  confidence level, defaults to 0.95.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `sides` names the side carrying the
  finite bound: `"left"` reports the lower limit and opens the upper one
  to `Inf`, `"right"` reports the upper limit and opens the lower one
  to 0. A one-sided bound at level `conf.level` is the corresponding end
  of the two-sided interval at level `2 * conf.level - 1`. `"left"`
  corresponds to the alternative `"greater"` in a hypothesis test.

- method:

  character string specifying which method to use; this can be one out
  of: `"exact"` (default), `"score"`, `"wald"` and `"byar"`. Several
  methods may be given at once; all the methods can be asked by
  `".all"`. Abbreviation of method is accepted. See details.

## Value

If recycling yields a single case, a named numeric vector with elements:

- `est`:

  point estimate of the event rate, `x/n`.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

If recycling yields multiple cases, a data frame with one row per case
is returned. Its first three columns are `est`, `lci`, and `uci`; the
remaining columns contain the recycled argument values.

## Details

All arguments are vectorized and recycled according to standard R rules.

The function assumes \\X \sim \mathrm{Poisson}(n\lambda)\\, where \\x\\
is the observed event count, \\n\\ the exposure and \\\lambda\\ the
event rate. The point estimate is \\\hat{\lambda} = x/n\\.

**Exact** (default): The so-called exact interval, based on quantiles of
the corresponding chi-square distribution and equivalent to the Garwood
(1936) interval. The same interval is returned as `conf.int` by
[`poisson.test()`](https://rdrr.io/r/stats/poisson.test.html).

**Score**: Obtained by inverting the acceptance region of the Poisson
score test.

**Wald**: Obtained by inverting the acceptance region of the
large-sample normal (Wald) test; the interval is symmetric about
\\\hat{\lambda}\\.

**Byar**: Byar's cube-root normal approximation to the exact interval
(Rothman and Boice, 1979). The lower limit is based on \\x\\, the upper
limit on \\x + 1\\; for \\x = 0\\ the lower limit is 0.

Some methods may produce a lower limit below 0. In such cases, the bound
is truncated to remain within the valid parameter space \\\[0,
\infty)\\.

## **Which interval should be used?**

The Wald interval has poor coverage for small counts and degenerates to
a single point at \\x = 0\\. The exact interval guarantees the nominal
coverage but is conservative and often unnecessarily wide. The score
interval is a good approximate alternative when conservativeness is not
required.

## References

Garwood, F. (1936) Fiducial limits for the Poisson distribution.
*Biometrika*, **28**, pp. 437-442.

Rothman, K. J. and Boice, J. D. Jr. (1979) *Epidemiologic Analysis with
a Programmable Calculator*. NIH Publication No. 79-1649. Washington, DC:
US Government Printing Office.

## See also

[`poisson.test()`](https://rdrr.io/r/stats/poisson.test.html),
[`poissonDiffCI()`](poissonDiffCI.md),
[`poissonRatioCI()`](poissonRatioCI.md), [`binomCI()`](binomCI.md)

## Examples

``` r

# deaths from horse kicks in 280 Prussian army corps-years
count <- 0:4
corpsYears <- c(144, 91, 32, 11, 2)

x <- sum(count * corpsYears)
n <- sum(corpsYears)

poissonCI(x, n)
#>       est       lci       uci 
#> 0.7000000 0.6054271 0.8051570 

# the exact interval corresponds to the one computed by poisson.test
poisson.test(x, n)$conf.int
#> [1] 0.6054271 0.8051570
#> attr(,"conf.level")
#> [1] 0.95

# several methods at once
poissonCI(x, n, method = c("exact", "score", "wald", "byar"))
#>   est       lci       uci   x   n conf.level     sides method
#> 1 0.7 0.6054271 0.8051570 196 280       0.95 two.sided  exact
#> 2 0.7 0.6086218 0.8050977 196 280       0.95 two.sided  score
#> 3 0.7 0.6020018 0.7979982 196 280       0.95 two.sided   wald
#> 4 0.7 0.6054225 0.8051602 196 280       0.95 two.sided   byar

# a 95% lower confidence bound for the event rate
poissonCI(x, n, sides = "left")
#>       est       lci       uci 
#> 0.7000000 0.6198372       Inf 

# SMR for Welsh nickel workers
poissonCI(x = 137, n = 24.19893)
#>      est      lci      uci 
#> 5.661407 4.753125 6.692709 
```
