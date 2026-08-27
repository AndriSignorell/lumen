# Nemenyi Test for Nonparametric All-Pairs Rank Comparisons

A nonparametric post hoc test for multiple pairwise comparisons
following a significant Kruskal-Wallis test, based on differences in
mean ranks.

## Usage

``` r
nemenyiTest(x, ...)

# S3 method for class 'formula'
nemenyiTest(formula, data, subset, na.action, ...)

# Default S3 method
nemenyiTest(
  x,
  g,
  dist = c("tukey", "chisq"),
  output = c("list", "matrix"),
  ...
)
```

## Arguments

- x:

  a numeric vector of data values, or a list of numeric data vectors.

- ...:

  further arguments passed to methods.

- formula:

  a formula of the form `response ~ group`.

- data:

  an optional data frame containing the variables in `formula`.

- subset:

  an optional expression specifying a subset of observations to be used.

- na.action:

  a function specifying how missing values should be handled.

- g:

  a grouping factor corresponding to `x`. Ignored if `x` is a list.

- dist:

  character string specifying the reference distribution used for the
  test statistic. One of `"tukey"` (default) or `"chisq"`.

- output:

  character string specifying the output format. One of `"list"`
  (default) or `"matrix"`.

## Value

An object of class `"rankTest"` containing:

- res:

  pairwise comparison results, either as a list or matrix

- pmat:

  symmetric matrix of adjusted p-values

Additional information is stored in attributes: `method`, `output`,
`main`, and `data.name`.

## Details

Performs Nemenyi's multiple comparison procedure for independent
samples. The test compares all pairs of groups using rank differences
and controls the family-wise error rate through the studentized range
distribution (Tukey-type approximation) or an asymptotic chi-squared
approximation.

Nemenyi's test is commonly used as a post hoc procedure after a
significant
[`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html) when all
pairwise comparisons between groups are of interest. Unlike
[`dunnTest`](https://andrisignorell.github.io/lumen/reference/dunnTest.md)
and
[`conoverTest`](https://andrisignorell.github.io/lumen/reference/conoverTest.md),
no additional p-value adjustment is applied, since multiplicity control
is built into the test statistic.

If `x` is a list, its elements are taken as the samples to be compared,
and hence have to be numeric data vectors. In this case, `g` is ignored
and one can simply use `nemenyiTest(x)`.

Otherwise, `x` must be a numeric vector and `g` a grouping factor (or
vector coercible to a factor) of the same length.

## References

Nemenyi, P. B. (1963). *Distribution-Free Multiple Comparisons*. PhD
thesis, Princeton University.

Hollander, M., Wolfe, D. A. and Chicken, E. (2014). *Nonparametric
Statistical Methods*. 3rd ed. Wiley.

## See also

[`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)

Other test.posthoc:
[`conoverTest()`](https://andrisignorell.github.io/lumen/reference/conoverTest.md),
[`dscfTest()`](https://andrisignorell.github.io/lumen/reference/dscfTest.md),
[`dunnTest()`](https://andrisignorell.github.io/lumen/reference/dunnTest.md),
[`dunnettTest()`](https://andrisignorell.github.io/lumen/reference/dunnettTest.md),
[`gamesHowellTest()`](https://andrisignorell.github.io/lumen/reference/gamesHowellTest.md),
[`plot.PostHocTest()`](https://andrisignorell.github.io/lumen/reference/plot.PostHocTest.md),
[`postHoc`](https://andrisignorell.github.io/lumen/reference/postHoc.md),
[`scheffeTest()`](https://andrisignorell.github.io/lumen/reference/scheffeTest.md),
[`signifDiff()`](https://andrisignorell.github.io/lumen/reference/signifDiff.md),
[`steelTest()`](https://andrisignorell.github.io/lumen/reference/steelTest.md)

## Examples

``` r
## Hollander & Wolfe example
x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
y <- c(3.8, 2.7, 4.0, 2.4)
z <- c(2.8, 3.4, 3.7, 2.2, 2.0)

nemenyiTest(list(x, y, z))
#> 
#>  Nemenyi's test of multiple comparisons for independent samples (tukey) 
#> 
#>     mean.rank.diff   pval    
#> 2-1            1.8 0.7972    
#> 3-1           -0.6 0.9720    
#> 3-2           -2.4 0.6686    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

x <- c(x, y, z)
g <- factor(rep(1:3, c(5, 4, 5)))

nemenyiTest(x, g)
#> 
#>  Nemenyi's test of multiple comparisons for independent samples (tukey) 
#> 
#>     mean.rank.diff   pval    
#> 2-1            1.8 0.7972    
#> 3-1           -0.6 0.9720    
#> 3-2           -2.4 0.6686    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

## Formula interface
nemenyiTest(Ozone ~ factor(Month), data = airquality)
#> 
#>  Nemenyi's test of multiple comparisons for independent samples (tukey) 
#> 
#>     mean.rank.diff    pval    
#> 6-5    12.02991453 0.88737    
#> 7-5    41.21153846 9.7e-05 ***
#> 8-5    38.53846154 0.00035 ***
#> 9-5    11.99734748 0.67819    
#> 7-6    29.18162393 0.16373    
#> 8-6    26.50854701 0.24773    
#> 9-6    -0.03256705 1.00000    
#> 8-7    -2.67307692 0.99853    
#> 9-7   -29.21419098 0.01136 *  
#> 9-8   -26.54111406 0.02867 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
```
