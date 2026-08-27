# Dunn's Test for Pairwise Rank Comparisons After a Kruskal–Wallis Test

A nonparametric post hoc test for multiple pairwise comparisons
following a significant Kruskal-Wallis test, based on rank sums with
adjustment for multiple testing.

## Usage

``` r
dunnTest(x, ...)

# S3 method for class 'formula'
dunnTest(formula, data, subset, na.action, ...)

# Default S3 method
dunnTest(
  x,
  g,
  method = p.adjust.methods,
  alternative = c("two.sided", "less", "greater"),
  output = c("list", "matrix"),
  alpha = 0.05,
  ...
)
```

## Arguments

- x:

  a numeric vector of observations or a list of numeric vectors.

- ...:

  further arguments passed to methods.

- formula:

  a formula of the form `response ~ group`.

- data:

  an optional data frame containing the variables in `formula`.

- subset:

  an optional expression specifying a subset of observations.

- na.action:

  a function indicating how missing values should be handled.

- g:

  a grouping variable corresponding to `x`; ignored when `x` is a list.

- method:

  the method used to adjust the p-values for multiple comparisons, one
  of `p.adjust.methods` (default is `"holm"`). Passed directly to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"less"` or `"greater"`. See the Details
  for the direction convention.

- output:

  the output format:

  - `"list"` pairwise comparison table.

  - `"matrix"` lower-triangular matrix of adjusted p-values.

- alpha:

  the significance level used to compile the groups flagged as
  significantly different in the label attribute of the p-value matrix
  (default is `0.05`).

## Value

An object of class `"rankTest"` containing:

- `res`:

  pairwise comparison results. Depending on `output`, either a table of
  mean-rank differences and adjusted p-values or a lower-triangular
  p-value matrix.

- `pmat`:

  symmetric matrix of adjusted p-values.

## Details

`dunnTest` performs the post hoc pairwise multiple-comparison procedure
appropriate after rejection of the Kruskal-Wallis null hypothesis. In
contrast to performing separate Wilcoxon rank-sum tests, Dunn's
procedure preserves the pooled ranking and variance estimate underlying
the Kruskal-Wallis test. It is intended as a post hoc procedure
following a significant Kruskal-Wallis test, i.e. typically for three or
more groups.

If `x` is a list, its elements are taken as the samples to be compared
and must be numeric vectors. In this case `g` is ignored. Otherwise, `x`
must be a numeric vector and `g` a grouping variable of the same length.

Each pairwise comparison is labeled `"B-A"`, where `A` precedes `B` in
the ordering of the group levels, and reports the mean rank difference
\\\bar{R}\_B - \bar{R}\_A\\. For one-sided alternatives, `"greater"`
tests whether `B` tends to have larger observations than `A` (upper
tail), and `"less"` tests the reverse (lower tail).

## References

Dunn, O. J. (1961) Multiple comparisons among means. *Journal of the
American Statistical Association*, 56 (293), 52-64.

Dunn, O. J. (1964) Multiple comparisons using rank sums.
*Technometrics*, 6 (3), 241-252.

## See also

[`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html),
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html),
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)

Other test.posthoc: [`conoverTest()`](conoverTest.md),
[`dscfTest()`](dscfTest.md), [`dunnettTest()`](dunnettTest.md),
[`gamesHowellTest()`](gamesHowellTest.md),
[`nemenyiTest()`](nemenyiTest.md),
[`plot.PostHocTest()`](plot.PostHocTest.md), [`postHoc`](postHoc.md),
[`scheffeTest()`](scheffeTest.md), [`signifDiff()`](signifDiff.md),
[`steelTest()`](steelTest.md)

## Examples

``` r
## Hollander & Wolfe (1973), p. 116
x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
y <- c(3.8, 2.7, 4.0, 2.4)
z <- c(2.8, 3.4, 3.7, 2.2, 2.0)

dunnTest(list(x, y, z))
#> 
#>  Dunn's test of multiple comparisons using rank sums : holm 
#> 
#>     mean.rank.diff   pval    
#> 2-1            1.8 1.0000    
#> 3-1           -0.6 1.0000    
#> 3-2           -2.4 1.0000    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

x <- c(x, y, z)
g <- factor(
  rep(1:3, c(5, 4, 5)),
  labels = c(
    "Normal subjects",
    "Subjects with obstructive airway disease",
    "Subjects with asbestosis"
  )
)

kruskal.test(x, g)
#> 
#>  Kruskal-Wallis rank sum test
#> 
#> data:  x and g
#> Kruskal-Wallis chi-squared = 0.77143, df = 2, p-value = 0.68
#> 
dunnTest(x, g)
#> 
#>  Dunn's test of multiple comparisons using rank sums : holm 
#> 
#>                                                                   mean.rank.diff
#> Subjects with obstructive airway disease-Normal subjects                     1.8
#> Subjects with asbestosis-Normal subjects                                    -0.6
#> Subjects with asbestosis-Subjects with obstructive airway disease           -2.4
#>                                                                     pval    
#> Subjects with obstructive airway disease-Normal subjects          1.0000    
#> Subjects with asbestosis-Normal subjects                          1.0000    
#> Subjects with asbestosis-Subjects with obstructive airway disease 1.0000    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

## Formula interface
dunnTest(Ozone ~ factor(Month), data = airquality)
#> 
#>  Dunn's test of multiple comparisons using rank sums : holm 
#> 
#>     mean.rank.diff    pval    
#> 6-5    12.02991453 1.00000    
#> 7-5    41.21153846 9.9e-05 ***
#> 8-5    38.53846154 0.00032 ***
#> 9-5    11.99734748 0.74574    
#> 7-6    29.18162393 0.14891    
#> 8-6    26.50854701 0.20743    
#> 9-6    -0.03256705 1.00000    
#> 8-7    -2.67307692 1.00000    
#> 9-7   -29.21419098 0.01036 *  
#> 9-8   -26.54111406 0.02428 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
```
