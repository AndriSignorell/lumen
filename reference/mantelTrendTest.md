# Mantel Linear-by-Linear Association Test for Detecting Linear Association Between Ordinal Variables

A chi-squared test for linear association between two ordinal variables
in a two-way contingency table, using row and column scores.

## Usage

``` r
mantelTrendTest(
  x,
  srow = scores(x, MARGIN = 1L, method = "table"),
  scol = scores(x, MARGIN = 2L, method = "table")
)
```

## Arguments

- x:

  a numeric matrix of counts (\\r \times c\\)

- srow:

  numeric vector of row scores; length must equal `nrow(x)`. Defaults to
  the numeric row `dimnames` of `x` if present, otherwise `1:nrow(x)`.
  See the Details.

- scol:

  numeric vector of column scores; length must equal `ncol(x)`. Defaults
  to the numeric column `dimnames` of `x` if present, otherwise
  `1:ncol(x)`. See the Details.

## Value

A list of class `"htest"` containing:

- statistic:

  the Mantel linear association chi-squared statistic

- parameter:

  degrees of freedom (always 1)

- p.value:

  the p-value

- estimate:

  the Pearson correlation coefficient \\r\\

- method:

  a character string describing the test

- data.name:

  a character string giving the name of the data

## Details

Tests for linear trend in a two-way \\r \times c\\ contingency table by
computing the Mantel linear-by-linear association statistic \$\$Q\_{MH}
= (n - 1) \cdot r^2\$\$ where \\r\\ is the Pearson correlation between
the row and column variables using the supplied scores. Under the null
hypothesis of no linear association, \\Q\\ has an asymptotic chi-squared
distribution with one degree of freedom.

This test is sometimes called the Mantel-Haenszel chi-squared test for
trend, but it is *not* the stratified Cochran-Mantel-Haenszel test for
\\2 \times 2 \times k\\ tables (see
[`mantelhaen.test`](https://rdrr.io/r/stats/mantelhaen.test.html) for
that). It is a score test for ordinal association, also known as the
linear-by-linear association test.

Both variables should be measured on an ordinal scale. If `x` has
numeric row and/or column `dimnames` (e.g. income brackets or dose
levels stored as label strings that parse as numbers), those values are
used as the default scores; otherwise the default is `1:nrow(x)` resp.
`1:ncol(x)`, i.e. the categories are assumed equally spaced. This
mirrors the scoring convention used by
[`cochranArmitageTest`](cochranArmitageTest.md). The choice of scores
affects the result: any monotone scores are permitted; non-monotone
scores (neither strictly increasing nor strictly decreasing) produce a
warning, since \\r\\ would then no longer reflect a consistent ordinal
trend.

## References

Agresti, A. (2002) *Categorical Data Analysis*, John Wiley & Sons, pp.
57, 86.

Mantel, N. (1963) Chi-square tests with one degree of freedom:
extensions of the Mantel-Haenszel procedure. *Journal of the American
Statistical Association*, 58, 690-700.

## See also

[`mantelhaen.test`](https://rdrr.io/r/stats/mantelhaen.test.html) for
the stratified Cochran-Mantel-Haenszel test,
[`chisq.test`](https://rdrr.io/r/stats/chisq.test.html) for the general
chi-squared test of independence,
[`cochranArmitageTest`](cochranArmitageTest.md) for a related trend test
with a binary response

Other test.categorical: [`barnardTest()`](barnardTest.md),
[`bhapkarTest()`](bhapkarTest.md),
[`breslowDayTest()`](breslowDayTest.md),
[`cochranQTest()`](cochranQTest.md), [`gTest()`](gTest.md),
[`lehmacherTest()`](lehmacherTest.md),
[`stuartMaxwellTest()`](stuartMaxwellTest.md),
[`woolfTest()`](woolfTest.md)

## Examples

``` r
## Agresti (2002, p. 57) Job Satisfaction
Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
              dimnames = list(
                income       = c("< 15k","15-25k","25-40k","> 40k"),
                satisfaction = c("VeryD","LittleD","ModerateS","VeryS")))

mantelTrendTest(Job)
#> 
#>  Mantel linear-by-linear association test
#> 
#> data:  Job
#> X-squared = 2.983, df = 1, p-value = 0.08414
#> sample estimates:
#>         r 
#> 0.1772001 
#> 
mantelTrendTest(Job, srow = c(7.5, 20, 32.5, 60))
#> 
#>  Mantel linear-by-linear association test
#> 
#> data:  Job
#> X-squared = 3.8075, df = 1, p-value = 0.05102
#> sample estimates:
#>         r 
#> 0.2001962 
#> 

## Automatic scores from numeric dimnames
dose <- matrix(c(10, 9, 10, 7, 0, 1, 0, 3), nrow = 4,
               dimnames = list(dose = c("0", "1", "2", "3"),
                               resp = c("no", "yes")))
mantelTrendTest(dose)  # srow taken as c(0, 1, 2, 3), not 1:4
#> 
#>  Mantel linear-by-linear association test
#> 
#> data:  dose
#> X-squared = 3.4667, df = 1, p-value = 0.06262
#> sample estimates:
#>         r 
#> 0.2981424 
#> 
```
