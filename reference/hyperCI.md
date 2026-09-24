# Confidence Intervals for the Number of Successes in a Finite Population

`hyperCI()` computes confidence intervals for the number \\M\\ of
successes in a finite population of size \\N\\, when \\x\\ successes are
observed in a sample of size \\n\\ drawn without replacement, i.e. \\X
\sim \mathrm{Hyper}(M, N - M, n)\\. It is the finite population
counterpart of [`binomCI()`](binomCI.md). The limits are counts; divide
them by `N` for the population proportion.

## Usage

``` r
hyperCI(
  x,
  n,
  N,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("wilson", "wald", "clopper-pearson", "mid-p", "blaker", "wang")
)
```

## Arguments

- x:

  number of successes in the sample, an integer between 0 and `n`.

- n:

  sample size, a positive integer not larger than `N`.

- N:

  population size, a positive integer.

- conf.level:

  confidence level, defaults to 0.95. With `NA` only the point estimate
  is returned.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `sides` names the side carrying the
  finite bound: `"left"` reports the lower limit and opens the upper one
  to `N - n + x`, `"right"` reports the upper limit and opens the lower
  one to `x`. A one-sided bound at level `conf.level` is the
  corresponding end of the two-sided interval at level
  `2 * conf.level - 1`, and therefore requires `conf.level > 0.5`. The
  exceptions are `"blaker"` and `"wang"`, which are calibrated on the
  two-sided coverage only; their one-sided bound is the Clopper-Pearson
  bound (see details).

- method:

  character string specifying which method to use; this can be one out
  of: `"wilson"` (default), `"wald"`, `"clopper-pearson"`, `"mid-p"`,
  `"blaker"` and `"wang"`. All the methods can be asked by `".all"`.
  Abbreviation of method is accepted. See details.

## Value

If recycling yields a single case, a named numeric vector with elements:

- `est`:

  point estimate of the number of successes in the population, the
  unbiased \\N x / n\\.

- `lci`:

  lower confidence interval bound, an integer.

- `uci`:

  upper confidence interval bound, an integer.

If recycling yields multiple cases, a data frame with one row per case
is returned. Its first three columns are `est`, `lci`, and `uci`; the
remaining columns contain the recycled argument values.

With `conf.level = NA` no interval is computed: the point estimate is
returned as an unnamed scalar, or as the single column `est` of the data
frame.

## Details

All arguments are vectorized and recycled according to standard R rules.

The sample already fixes \\x \le M \le N - n + x\\. All limits are kept
within this range, and a one-sided interval opens its free side to it.

**Wald**: The Wald interval for the proportion with the finite
population correction, \\\hat p \pm z \sqrt{f \hat p (1 - \hat p) / n}\\
with \\\hat p = x/n\\ and \\f = (N - n)/(N - 1)\\ (Cochran 1977),
multiplied by \\N\\.

**Wilson** (default): The Wilson score interval of
[`binomCI()`](binomCI.md) with the effective sample size \\n / f\\,
which carries the finite population correction into the score variance,
multiplied by \\N\\.

The two asymptotic intervals are rounded outwards to integers. For a
census (\\n = N\\) they return \\M = x\\.

**Clopper-Pearson**: The exact interval obtained by inverting two
one-sided hypergeometric tests (Konijn 1973): the lower limit is the
smallest \\M\\ with \\P(X \ge x \mid M) \> \alpha/2\\, the upper limit
the largest \\M\\ with \\P(X \le x \mid M) \> \alpha/2\\. It guarantees
the coverage but is conservative.

**Mid-p**: As Clopper-Pearson, with half the probability of the observed
count in each tail. Not guaranteed to reach the level, but closer to it
on average.

**Blaker**: The exact interval of Blaker (2000), obtained by inverting
the test that uses the smaller tail probability as statistic. It is
contained in the Clopper-Pearson interval and keeps the level; its
acceptance region can in rare cases have gaps, the interval spans the
smallest and largest accepted value. One-sided, the Clopper-Pearson
bound is returned: an end of the two-sided Blaker interval at level
`2 * conf.level - 1` misses the level.

**Wang**: The admissible exact interval of Wang (2015). Starting from
the Clopper-Pearson interval, the limits are shrunk pairwise (\\U_x =
N - L\_{n-x}\\) from the middle of the sample space outwards, each as
far as the coverage permits. The resulting family is monotone and
symmetric, and no limit can be moved inwards, the other intervals held
fixed, without the coverage falling below `conf.level`. It is never
wider than Clopper-Pearson. One-sided, the Clopper-Pearson bound is
already the smallest exact bound and is returned instead of an end of
the two-sided Wang interval at level `2 * conf.level - 1`, which would
miss the level. The computation proceeds from \\n/2\\ towards \\x\\; it
takes below a second for \\n = 5000\\, \\N = 10^6\\.

## References

Blaker, H. (2000) Confidence curves and improved exact confidence
intervals for discrete distributions, *Canadian Journal of Statistics*
28 (4), 783-798

Cochran, W. G. (1977) *Sampling Techniques*, 3rd ed. New York: Wiley.

Konijn, H. S. (1973) *Statistical Theory of Sample Survey Design and
Analysis*. Amsterdam: North-Holland.

Wang, W. (2015) Exact optimal confidence intervals for hypergeometric
parameters, *Journal of the American Statistical Association* 110 (512),
1491-1499,
[doi:10.1080/01621459.2014.966191](https://doi.org/10.1080/01621459.2014.966191)

Wilson, E. B. (1927) Probable inference, the law of succession, and
statistical inference, *Journal of the American Statistical Association*
22, 209-212.

## See also

[`binomCI()`](binomCI.md) for sampling with replacement or an infinite
population, [`phyper()`](https://rdrr.io/r/stats/Hypergeometric.html)

Other ci.proportion: [`binomCI()`](binomCI.md),
[`binomDiffCI()`](binomDiffCI.md), [`binomRatioCI()`](binomRatioCI.md),
[`multinomCI()`](multinomCI.md)

## Examples

``` r
# audit: 10 faulty invoices in a sample of 50 out of 2000
hyperCI(x = 10, n = 50, N = 2000, method = ".all")
#>   est lci uci  x  n    N conf.level     sides          method
#> 1 400 226 658 10 50 2000       0.95 two.sided          wilson
#> 2 400 180 620 10 50 2000       0.95 two.sided            wald
#> 3 400 203 671 10 50 2000       0.95 two.sided clopper-pearson
#> 4 400 215 651 10 50 2000       0.95 two.sided           mid-p
#> 5 400 207 654 10 50 2000       0.95 two.sided          blaker
#> 6 400 211 661 10 50 2000       0.95 two.sided            wang

# the same as proportions
hyperCI(x = 10, n = 50, N = 2000, method = "wang")[c("lci", "uci")] / 2000
#>    lci    uci 
#> 0.1055 0.3305 

# the finite population correction matters once n is a sizeable
# fraction of N, and vanishes for large N
hyperCI(x = 10, n = 50, N = c(60, 200, 1e6), method = "clopper-pearson")
#>      est    lci    uci  x  n     N conf.level     sides          method
#> 1     12     10     15 10 50 6e+01       0.95 two.sided clopper-pearson
#> 2     40     23     63 10 50 2e+02       0.95 two.sided clopper-pearson
#> 3 200000 100305 337179 10 50 1e+06       0.95 two.sided clopper-pearson
binomCI(x = 10, n = 50, method = "clopper-pearson")
#>       est       lci       uci 
#> 0.2000000 0.1003022 0.3371831 

# an upper bound for the number of faulty items after a clean sample
hyperCI(x = 0, n = 50, N = 2000, sides = "right", method = "clopper-pearson")
#> est lci uci 
#>   0   0 114 
```
