# Confidence Intervals for the Ratio of Two Binomial Proportions

Computes confidence intervals for the ratio of two independent binomial
proportions (relative risk). Several classical and modern methods are
available, which may behave quite differently for small samples or
extreme proportions.

## Usage

``` r
binomRatioCI(
  x1,
  n1,
  x2,
  n2,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("koopman", "bailey", "adj-log", "katz-log", "sinh-1", "noether"),
  tol = .Machine$double.eps^0.25
)
```

## Arguments

- x1:

  number of successes in the first group (numerator).

- n1:

  number of trials in the first group.

- x2:

  number of successes in the second group (denominator).

- n2:

  number of trials in the second group.

- conf.level:

  confidence level, default is 0.95.

- sides:

  a character string specifying the type of confidence interval:
  `"two.sided"` (default), `"left"`, or `"right"`. Partial matching is
  allowed.

- method:

  one of: `"koopman"`, `"bailey"`, `"adj-log"`, `"katz-log"`,
  `"sinh-1"`, `"noether"`.

- tol:

  desired accuracy (convergence tolerance) for the iterative
  root-finding procedure used by the Koopman interval.

## Value

If recycling yields a single case, a named numeric vector with elements:

- `est`:

  point estimate of the ratio of binomial proportions,
  `(x1/n1) / (x2/n2)`.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

If recycling yields multiple cases, a data frame with one row per case
is returned. Its first three columns are `est`, `lci`, and `uci`; the
remaining columns contain the recycled argument values.

## Details

All arguments are vectorized and recycled according to standard R rules.

The ratio of proportions is estimated by \$\$ \hat{\theta} =
\frac{\hat{p}\_1}{\hat{p}\_2} = \frac{x_1 / n_1}{x_2 / n_2}. \$\$

**Katz-log**: The classical large-sample log-transformed Wald interval
described by Katz et al. (1978). This method is simple and widely used
but may perform poorly for small sample sizes or proportions close to 0
or 1.

**Adjusted-log**: A continuity-adjusted modification of the Katz
interval using additive corrections to reduce bias and improve
performance in sparse data settings (Walter, 1975; Pettigrew et al.,
1986).

**Bailey**: A skewness-corrected interval proposed by Bailey (1987)
based on a cubic transformation. Often provides improved coverage
probabilities relative to standard Wald-type intervals.

**Koopman**: An asymptotic score interval obtained by inverting
Pearson's chi-square statistic (Koopman, 1984). Confidence limits are
determined iteratively via root-finding and the method is generally
considered among the most reliable asymptotic procedures.

**Noether**: A large-sample interval based directly on the asymptotic
variance of the estimated ratio (Noether, 1957).

**Inverse hyperbolic sine**: Based on the variance-stabilizing inverse
hyperbolic sine transformation proposed by Newcombe (2001).

Several methods require special handling in cases where \\x_1 = 0\\,
\\x_2 = 0\\, \\x_1 = n_1\\, or \\x_2 = n_2\\. Small continuity
corrections are therefore applied internally where needed.

Some methods may produce infinite limits when one observed proportion is
zero. This is expected behavior for ratio parameters.

**Which interval should be used?**  
The choice of method remains an active topic of discussion. The Koopman
interval is often recommended due to its comparatively good coverage
properties across a broad range of scenarios.

## Note

Based on code by Ken Aho, adapted to conform to package standards.

## References

Bailey BJR (1987). Confidence limits to the risk ratio. *Biometrics*,
43(1), 201-205.

Katz D, Baptista J, Azen SP, Pike MC (1978). Obtaining confidence
intervals for the risk ratio in cohort studies. *Biometrics*, 34,
469-474.

Koopman PAR (1984). Confidence intervals for the ratio of two binomial
proportions. *Biometrics*, 40, 513-517.

Newcombe RG (2001). Logit confidence intervals and the inverse sinh
transformation. *The American Statistician*, 55, 200-202.

Noether GE (1957). Sample size determination for some common
nonparametric tests. *Journal of the American Statistical Association*,
52, 645-647.

Pettigrew HM, Gart JJ, Thomas DG (1986). The bias and higher cumulants
of the logarithm of a binomial variate. *Biometrika*, 73(2), 425-435.

Walter SD (1975). The distribution of Levin's measure of attributable
risk. *Biometrika*, 62(2), 371-374.

## See also

[`binom.test`](https://rdrr.io/r/stats/binom.test.html),
[`prop.test`](https://rdrr.io/r/stats/prop.test.html)

Other ci.proportion: [`binomCI()`](binomCI.md),
[`binomDiffCI()`](binomDiffCI.md), [`multinomCI()`](multinomCI.md)

## Examples

``` r

# Example from Koopman (1984)
binomRatioCI(
  x1 = 36, n1 = 40, 
  x2 = 16, n2 = 80,   method = "katz-log")
#>      est      lci      uci 
#> 4.500000 2.868550 7.059315 

binomRatioCI(
  x1 = 36, n1 = 40, 
  x2 = 16, n2 = 80,   method = "koopman")
#>      est      lci      uci 
#> 4.500000 2.939581 7.152209 

# Compare several methods
meths <- c("katz-log", "adj-log", "bailey", "koopman", "noether", "sinh-1")
binomRatioCI(
  x1 = 25, n1 = 100, 
  x2 = 10, n2 = 100,   method = meths)
#>   est       lci      uci x1  n1 x2  n2 conf.level     sides   method
#> 1 2.5 1.2678712 4.929523 25 100 10 100       0.95 two.sided katz-log
#> 2 2.5 1.2509950 4.714614 25 100 10 100       0.95 two.sided  adj-log
#> 3 2.5 1.3076711 5.171624 25 100 10 100       0.95 two.sided   bailey
#> 4 2.5 1.2953219 4.907144 25 100 10 100       0.95 two.sided  koopman
#> 5 2.5 0.8026214 4.197379 25 100 10 100       0.95 two.sided  noether
#> 6 2.5 1.2837004 4.868737 25 100 10 100       0.95 two.sided   sinh-1
#>            tol
#> 1 0.0001220703
#> 2 0.0001220703
#> 3 0.0001220703
#> 4 0.0001220703
#> 5 0.0001220703
#> 6 0.0001220703

# Sparse data
binomRatioCI(
  x1 = 1, n1 = 20, 
  x2 = 0, n2 = 20,   method = meths )
#>   est        lci uci x1 n1 x2 n2 conf.level     sides   method          tol
#> 1 Inf 0.07103723 Inf  1 20  0 20       0.95 two.sided katz-log 0.0001220703
#> 2 Inf 0.12965188 Inf  1 20  0 20       0.95 two.sided  adj-log 0.0001220703
#> 3 Inf 0.05677804 Inf  1 20  0 20       0.95 two.sided   bailey 0.0001220703
#> 4 Inf 0.27072686 Inf  1 20  0 20       0.95 two.sided  koopman 0.0001220703
#> 5 Inf 0.00000000 Inf  1 20  0 20       0.95 two.sided  noether 0.0001220703
#> 6 Inf 0.07008511 Inf  1 20  0 20       0.95 two.sided   sinh-1 0.0001220703
  
  
```
