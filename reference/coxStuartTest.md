# Cox-Stuart Trend Test for Detecting Monotonic Trends in Ordered Data

Tests a sequence of observations for a monotone trend by pairing each
observation in the first half of the series with the corresponding
observation in the second half and applying a sign test to the
differences.

## Usage

``` r
coxStuartTest(x, alternative = c("two.sided", "increasing", "decreasing"))
```

## Arguments

- x:

  a numeric vector of observations in sequence order

- alternative:

  a character string specifying the alternative hypothesis, one of
  `"two.sided"` (default), `"increasing"` or `"decreasing"`

## Value

An object of class `"htest"` with components

- statistic:

  the number of positive paired differences, with names attribute `"S"`

- parameter:

  the number of untied pairs entering the test

- p.value:

  the p-value

- estimate:

  the proportion of increasing pairs

- alternative:

  a character string describing the alternative hypothesis

- method:

  the character string `"Cox-Stuart trend test"`

- data.name:

  a character string giving the name of the data

## Details

The series \\x_1, \ldots, x_N\\ is split in the middle; with an odd
\\N\\ the central observation is discarded. Writing \\m = \lfloor N/2
\rfloor\\, each of the first \\m\\ remaining observations is paired with
the one \\m\\ places later in the reduced series, that is with
\\x\_{i+m}\\ for even \\N\\ and with \\x\_{i+m+1}\\ in the original
indexing for odd \\N\\. The statistic \\S\\ counts how many of the \\m\\
paired differences are positive. Pairs that are exactly tied carry no
information about direction and are dropped.

Every observation enters at most one pair, so under independence and a
continuous common distribution the signs are independent Bernoulli
variables with probability one half. The p-value from
[`binom.test`](https://rdrr.io/r/stats/binom.test.html) is then exact at
any series length, however short. That exactness rests on the
independence assumption: with serially dependent observations the signs
need not be independent and the nominal level is no longer guaranteed.

The test is a special case of
[`signTest`](https://andrisignorell.github.io/lumen/reference/signTest.md)
and inherits both its robustness and its low power: only the sign of
each paired difference is used, and half the observations enter only as
partners. It detects a monotone drift, not curvature or oscillation: a
series that rises and then falls back can easily produce a p-value near
one.

Unlike the other members of `test.trend`, which need a grouping factor
or a contingency table, this test takes a bare series and therefore has
no formula interface. Missing values are removed before the series is
split, so the pairing refers to the observed values, not to their
original positions.

## References

Cox, D. R., Stuart, A. (1955) Some quick sign tests for trend in
location and dispersion. *Biometrika*, **42**(1/2), 80-95.

## See also

[`signTest`](https://andrisignorell.github.io/lumen/reference/signTest.md),
[`jonckheereTerpstraTest`](https://andrisignorell.github.io/lumen/reference/jonckheereTerpstraTest.md),
[`mantelTrendTest`](https://andrisignorell.github.io/lumen/reference/mantelTrendTest.md),
[`bartelsRankTest`](https://andrisignorell.github.io/lumen/reference/BartelsRankTest.md),
[`runsTest`](https://andrisignorell.github.io/lumen/reference/runsTest.md)

Other test.trend:
[`cochranArmitageTest()`](https://andrisignorell.github.io/lumen/reference/cochranArmitageTest.md),
[`jonckheereTerpstraTest()`](https://andrisignorell.github.io/lumen/reference/jonckheereTerpstraTest.md),
[`pageTest()`](https://andrisignorell.github.io/lumen/reference/pageTest.md)

## Examples

``` r
## a strictly increasing series
coxStuartTest(1:12)
#> 
#>  Cox-Stuart trend test
#> 
#> data:  1:12
#> S = 6, n = 6, p-value = 0.03125
#> alternative hypothesis: two.sided
#> sample estimates:
#> proportion of increasing pairs 
#>                              1 
#> 
## [1] S = 6, n = 6, p-value = 0.03125

coxStuartTest(1:12, alternative = "increasing")$p.value
#> [1] 0.015625
## [1] 0.015625

## no trend
set.seed(1)
coxStuartTest(rnorm(50))$p.value
#> [1] 0.1077521
```
