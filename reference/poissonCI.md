# Confidence Interval for a Poisson Rate

Estimates a Poisson event rate and calculates an exact or approximate
confidence interval.

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

  non-negative integer event count or vector of counts

- n:

  positive exposure associated with `x`, such as observation time,
  person-time, or population at risk; may be a vector and defaults to 1

- conf.level:

  numeric confidence level between 0 and 1; defaults to 0.95

- sides:

  type of confidence interval: `"two.sided"`, `"left"`, or `"right"`;
  may be abbreviated

- method:

  method used to calculate the confidence interval: `"exact"`,
  `"score"`, `"wald"`, or `"byar"`; may be abbreviated and may contain
  several methods; defaults to `"exact"`

## Value

If all arguments identify a single result, a named numeric vector with
elements:

- `est`:

  estimated event rate

- `lci`:

  lower confidence bound

- `uci`:

  upper confidence bound

Otherwise, a `data.frame` containing these three columns and the
recycled argument values that identify each result.

## Details

The function assumes \$\$X \sim \mathrm{Poisson}(n\lambda),\$\$ where
`x` is the observed event count, `n` is the exposure, and \\\lambda\\ is
the event rate. The point estimate is \\\hat{\lambda} = x/n\\.

The available confidence-interval methods are:

- `"exact"`:

  the exact Poisson interval calculated by
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html),
  equivalent to the Garwood interval

- `"score"`:

  the interval obtained by inverting the Poisson score test

- `"wald"`:

  the normal-approximation interval centred at \\\hat{\lambda}\\

- `"byar"`:

  Byar's cube-root normal approximation

The lower bound is restricted to the parameter space \\\[0, \infty)\\.
For `sides = "left"`, the function returns a lower one-sided confidence
bound and sets `uci` to `Inf`. This corresponds to the alternative
`"greater"` in a hypothesis test. For `sides = "right"`, it returns an
upper one-sided confidence bound and sets `lci` to 0.

Compatible vector lengths are recycled. Supplying several methods
therefore provides the corresponding intervals in a single `data.frame`.

## References

Garwood, F. (1936). Fiducial limits for the Poisson distribution.
*Biometrika*, **28**, 437–442.

Rothman, K. J. and Boice, J. D. Jr. (1979). *Epidemiologic Analysis with
a Programmable Calculator*. NIH Publication No. 79-1649. Washington, DC:
US Government Printing Office.

## See also

[`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)

## Examples

``` r
# Deaths from horse kicks in 280 Prussian army corps-years
count <- 0:4
corpsYears <- c(144, 91, 32, 11, 2)

x <- sum(count * corpsYears)
n <- sum(corpsYears)

poissonCI(x, n)
#>       est       lci       uci 
#> 0.7000000 0.6054271 0.8051570 
poissonCI(x, n, method = c("exact", "score", "wald", "byar"))
#>   est       lci       uci   x   n conf.level     sides method
#> 1 0.7 0.6054271 0.8051570 196 280       0.95 two.sided  exact
#> 2 0.7 0.6086218 0.8050977 196 280       0.95 two.sided  score
#> 3 0.7 0.6020018 0.7979982 196 280       0.95 two.sided   wald
#> 4 0.7 0.6070833 0.8032497 196 280       0.95 two.sided   byar

# A 95% lower confidence bound for the event rate
poissonCI(x, n, sides = "left")
#>       est       lci       uci 
#> 0.7000000 0.6198372       Inf 

# SMR for Welsh nickel workers
poissonCI(x = 137, n = 24.19893)
#>      est      lci      uci 
#> 5.661407 4.753125 6.692709 
```
