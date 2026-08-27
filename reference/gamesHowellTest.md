# Games-Howell Post-Hoc Test for Pairwise Mean Comparisons Under Unequal Variances

Performs all pairwise comparisons between group means without assuming
equal variances, using a separate Welch-type standard error and
Satterthwaite degrees of freedom for every pair and referring the result
to the studentized range distribution.

## Usage

``` r
gamesHowellTest(x, ...)

# Default S3 method
gamesHowellTest(x, g, conf.level = 0.95, ...)

# S3 method for class 'formula'
gamesHowellTest(formula, data, subset, na.action = na.pass, ...)

# S3 method for class 'aov'
gamesHowellTest(x, conf.level = 0.95, ...)
```

## Arguments

- x:

  a numeric vector of observations, an `aov` object, or a formula of the
  form `lhs ~ rhs` with a numeric `lhs` and a grouping `rhs`

- ...:

  further arguments, passed to the default method

- g:

  a vector or factor giving the group for each element of `x`

- conf.level:

  confidence level of the simultaneous intervals

- formula:

  a formula of the form `lhs ~ rhs`

- data:

  an optional data frame containing the model variables

- subset:

  an optional vector specifying a subset of observations

- na.action:

  a function indicating what should happen when the data contain `NA`s

## Value

An object of class `"PostHocTest"`: a list with one matrix, named after
the grouping variable. The matrix has columns `diff` for the observed
mean difference (second group minus first), `lci` and `uci` for the
simultaneous confidence limits, and `pval` for the simultaneous p-value.
Print and plot methods are available for class `"PostHocTest"`.

## Details

Tukey's HSD,
[`scheffeTest`](https://andrisignorell.github.io/lumen/reference/scheffeTest.md)
and the parametric methods in
[`postHoc`](https://andrisignorell.github.io/lumen/reference/postHoc.md)
all rely on a single pooled error variance. When the group variances
differ, that pooled estimate is wrong for every pair that does not
happen to match it, and the procedure loses its nominal level -
liberally when the larger variance sits in the smaller group,
conservatively in the opposite case. Games and Howell (1976) replace the
pooled term by the pairwise Welch standard error, which makes this the
post-hoc counterpart of the Welch t test and of
[`oneway.test`](https://rdrr.io/r/stats/oneway.test.html).
[`yuenTTest`](https://andrisignorell.github.io/lumen/reference/yuenTTest.md)
addresses heteroscedasticity and non-normality with trimmed means and
separate Winsorized variance estimates, but provides no all-pairs
multiple-comparison procedure, so it is not the two-sample form of this
procedure.

For groups \\i\\ and \\j\\ the statistic is \$\$q = \frac{\|\bar{x}\_j -
\bar{x}\_i\|}{\sqrt{(s_i^2/n_i + s_j^2/n_j)/2}},\$\$ referred to the
studentized range distribution with \\k\\ means and Satterthwaite
degrees of freedom. Because \\q\\ already accounts for the number of
groups, the p-values and interval widths are simultaneous over all
\\k(k-1)/2\\ comparisons; no further multiplicity adjustment is applied
or needed.

With equal group sizes and equal variances the standard error reduces
exactly to Tukey's, and the two procedures differ only in the degrees of
freedom - pairwise \\2(n-1)\\ instead of pooled \\k(n-1)\\ - which makes
Games-Howell slightly the more conservative of the two in that case.

Every group needs at least two non-missing observations, since each
contributes its own variance estimate. Two obstacles can still leave an
individual pair incomputable: both groups constant, which leaves no
scale to studentize by, and Welch degrees of freedom below two, where
the studentized range is undefined. The latter is not an exotic case -
with only two observations in each of two groups the Welch degrees of
freedom lie in \\\[1, 2\]\\ and reach two only when the two variances
are exactly equal. Affected pairs are reported as `NA` with a warning
while the remaining comparisons are kept. Missing values are removed
casewise.

## References

Games, P. A., Howell, J. F. (1976) Pairwise multiple comparison
procedures with unequal n's and/or variances: a Monte Carlo study.
*Journal of Educational Statistics*, **1**(2), 113-125.

## See also

[`TukeyHSD()`](https://rdrr.io/r/stats/TukeyHSD.html),
[yuenTTest](https://andrisignorell.github.io/lumen/reference/yuenTTest.md)

Other test.posthoc:
[`conoverTest()`](https://andrisignorell.github.io/lumen/reference/conoverTest.md),
[`dscfTest()`](https://andrisignorell.github.io/lumen/reference/dscfTest.md),
[`dunnTest()`](https://andrisignorell.github.io/lumen/reference/dunnTest.md),
[`dunnettTest()`](https://andrisignorell.github.io/lumen/reference/dunnettTest.md),
[`nemenyiTest()`](https://andrisignorell.github.io/lumen/reference/nemenyiTest.md),
[`plot.PostHocTest()`](https://andrisignorell.github.io/lumen/reference/plot.PostHocTest.md),
[`postHoc`](https://andrisignorell.github.io/lumen/reference/postHoc.md),
[`scheffeTest()`](https://andrisignorell.github.io/lumen/reference/scheffeTest.md),
[`signifDiff()`](https://andrisignorell.github.io/lumen/reference/signifDiff.md),
[`steelTest()`](https://andrisignorell.github.io/lumen/reference/steelTest.md)

## Examples

``` r
gamesHowellTest(breaks ~ tension, data = warpbreaks)
#>     95% family-wise confidence level
#> 
#> Fit: gamesHowellTest.formula(breaks ~ tension, data = warpbreaks)
#> 
#> $tension
#>           diff       lci       uci  pval signif
#> M-L -10.000000 -21.00113  1.001128 0.080    .  
#> H-L -14.722222 -25.54579 -3.898659 0.006    ** 
#> H-M  -4.722222 -11.86790  2.423460 0.251       
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

## compare with the pooled-variance procedures
gamesHowellTest(breaks ~ tension, data = warpbreaks)$tension[, "pval"]
#>         M-L         H-L         H-M 
#> 0.080082272 0.006355163 0.251298244 
## [1] 0.080082272 0.006355163 0.251298244

TukeyHSD(aov(breaks ~ tension, data = warpbreaks))
#>   Tukey multiple comparisons of means
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = breaks ~ tension, data = warpbreaks)
#> 
#> $tension
#>           diff       lwr        upr     p adj
#> M-L -10.000000 -19.55982 -0.4401756 0.0384598
#> H-L -14.722222 -24.28205 -5.1623978 0.0014315
#> H-M  -4.722222 -14.28205  4.8376022 0.4630831
#> 
```
