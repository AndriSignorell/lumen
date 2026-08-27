# Moses Test of Extreme Reactions for Comparing Extreme-Value Behaviour Between Two Groups

Performs the Moses test of extreme reactions, a distribution-free
nonparametric test for the difference in extremity of scores between two
independent groups. Scores from both groups are pooled and converted to
ranks. The test statistic is the span of the control group's ranks
(range + 1). An exact one-tailed probability is computed for the raw
span, then recomputed after trimming `h` extreme ranks from each end of
the control group.

## Usage

``` r
mosesTest(x, ...)

# S3 method for class 'formula'
mosesTest(formula, data, subset, na.action = na.pass, ...)

# Default S3 method
mosesTest(
  x,
  y,
  extreme = NULL,
  ties.method = c("first", "average", "random", "max", "min"),
  ...
)

# S3 method for class 'mosesTestResult'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- x:

  Numeric vector of observations from the control group. Non-finite
  values are removed.

- ...:

  further arguments passed to or from methods

- formula:

  A formula of the form `lhs ~ rhs`, where `lhs` gives the observations
  and `rhs` the grouping variable with exactly two levels.

- data:

  Optional data frame containing the variables in `formula`.

- subset:

  Optional vector specifying a subset of observations.

- na.action:

  Function specifying how missing values are handled.

- y:

  Numeric vector of observations from the experiment group. Non-finite
  values are removed.

- extreme:

  Non-negative integer \\h\\. Number of extreme ranks trimmed from each
  end of the control group before recomputing the span. If `NULL`,
  defaults to `max(floor(0.05 * length(x)), 1)`.

- ties.method:

  Character string passed to [`rank`](https://rdrr.io/r/base/rank.html).
  Default `"first"` preserves integer-valued ranks and exact
  combinatorial validity.

- digits:

  number of significant digits to display

## Value

An object of class `"mosesTestResult"` inheriting from `"htest"` with
components:

- `statistic`:

  Named vector containing `sRaw` and `sTrimmed`.

- `p.value`:

  Named vector containing `p_raw` and `p_trimmed`.

- `extreme`:

  Effective trimming parameter \\h\\.

- `parameter`:

  Vector containing `n_control` and `n_experiment`.

- `null.value`:

  Null hypothesis of equal extremity.

- `alternative`:

  Character string: `"greater"`.

- `method`:

  Character string describing the test.

- `data.name`:

  Character string describing the data.

## Details

A nonparametric test comparing the spread (range) of two independent
groups, assessing whether the control group shows greater variability
than the experiment group.

**Distributional derivation.** Under the null hypothesis, any assignment
of \\n_k + n_e\\ pooled ranks to the two groups is equally likely. The
exact distribution of the span \\S\\ of the control group (a subset of
size \\n_k\\ drawn without replacement from \\\\1, \ldots, N\\\\) is
given by Moses (1952) as:

\$\$P(S \le s) = \frac{1}{\binom{N}{n_k}} \sum\_{i=0}^{\min(s - n_k,\\
n_e)} \binom{i + n_k - 2}{i}\\ \binom{n_e + 1 - i}{n_e - i}\$\$

This is a cumulative distribution function (CDF): \\P(\mathrm{span} \le
s)\\.

The one-tailed p-value for the observed span \\S\_{obs}\\ is:

\$\$ p = P(S \ge S\_{obs}) = 1 - P(S \le S\_{obs} - 1) \$\$

This distinction matters because some textbook presentations write the
summation formula as \\P(S = s)\\, whereas the upper-tail probability
required for hypothesis testing is obtained from the cumulative form.

**Generalisation with trimming (\\h \> 0\\).** After removing \\h\\
extreme ranks from each end of the sorted control group ranks, the
effective control-group size becomes \\n_k^\prime = n_k - 2h\\,
yielding:

\$\$ P(S_h \le s) = \frac{1}{\binom{N}{n_k}} \sum\_{i=0}^{\min(s -
n_k^\prime,\\ n_e)} \binom{i + n_k^\prime - 2}{i}\\ \binom{n_e + 2h +
1 - i}{n_e - i} \$\$

**Tie handling.** The exact combinatorial distribution assumes a
continuous underlying distribution (no ties). By default,
`ties.method = "first"` is passed to
[`rank`](https://rdrr.io/r/base/rank.html), producing integer-valued
ranks that remain compatible with the exact formula. Tied observations
are ordered according to their occurrence in the pooled sample, matching
SPSS behaviour.

Alternative tie methods such as `"average"` may produce fractional
mid-ranks. In such cases the span is rounded upward via
[`ceiling()`](https://rdrr.io/r/base/Round.html) to retain an
integer-valued test statistic, but the resulting p-values should be
regarded as approximate rather than exact.

**Numerical stability.** The exact distribution is evaluated entirely in
log-space using [`lchoose`](https://rdrr.io/r/base/Special.html) and a
log-sum-exp transformation, avoiding overflow for large sample sizes.

## References

Moses, L.E. (1952). A two-sample test. *Psychometrika*, **17**, 239–247.
[doi:10.1007/BF02288735](https://doi.org/10.1007/BF02288735)

## See also

[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html),
[`ks.test()`](https://rdrr.io/r/stats/ks.test.html),
[`ansari.test()`](https://rdrr.io/r/stats/ansari.test.html)

Other test.variance: [`leveneTest()`](leveneTest.md),
[`siegelTukeyTest()`](siegelTukeyTest.md), [`varCI()`](varCI.md),
[`varTest()`](varTest.md)

## Examples

``` r
x <- c(0.80, 0.83, 1.89, 1.04, 1.45,
       1.38, 1.91, 1.64, 0.73, 1.46)

y <- c(1.15, 0.88, 0.90, 0.74, 1.21)

mosesTest(x, y)
#> 
#>   Moses Test of Extreme Reactions 
#> 
#> data:  x and y 
#> n_control = 10,  n_experiment = 5
#> 
#>   Without trimming:                 span =   15,  p-value = 0.4285714
#>   After trimming 1 from each end:   span =   12,  p-value = 0.4335664
#> 
#> alternative hypothesis: extreme values are more likely in the control group
#> 

set.seed(1479)

x2 <- sample(1:20, 10, replace = TRUE)
y2 <- sample(5:25, 6, replace = TRUE)

mosesTest(x2, y2, extreme = 2)
#> 
#>   Moses Test of Extreme Reactions 
#> 
#> data:  x2 and y2 
#> n_control = 10,  n_experiment = 6
#> 
#>   Without trimming:                 span =   13,  p-value = 0.9642857
#>   After trimming 2 from each end:   span =    7,  p-value = 0.9423077
#> 
#> alternative hypothesis: extreme values are more likely in the control group
#> 

df <- data.frame(
  score = c(x, y),
  group = factor(rep(
    c("control", "experiment"),
    c(length(x), length(y))
  ))
)

mosesTest(score ~ group, data = df)
#> 
#>   Moses Test of Extreme Reactions 
#> 
#> data:  groups[[1L]] and groups[[2L]] 
#> n_control = 10,  n_experiment = 5
#> 
#>   Without trimming:                 span =   15,  p-value = 0.4285714
#>   After trimming 1 from each end:   span =   12,  p-value = 0.4335664
#> 
#> alternative hypothesis: extreme values are more likely in the control group
#> 
```
