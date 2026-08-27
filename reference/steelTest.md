# Steel Test for Nonparametric Comparisons With a Common Control

A nonparametric multiple-comparison procedure for comparing several
treatment groups with a single control group based on Wilcoxon rank-sum
statistics.

## Usage

``` r
steelTest(x, ...)

# S3 method for class 'formula'
steelTest(formula, data, subset, na.action, ...)

# Default S3 method
steelTest(
  x,
  g,
  control = NULL,
  alternative = c("two.sided", "greater", "less"),
  output = c("list", "matrix"),
  alpha = 0.05,
  ...
)
```

## Arguments

- x:

  a numeric vector of observations, or a list of numeric vectors.

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

- control:

  the level of the control group against which all treatment groups are
  compared. Defaults to the first group.

- alternative:

  character string specifying the alternative hypothesis. One of
  `"two.sided"`, `"greater"` or `"less"`.

- output:

  character string specifying the output format. One of `"list"`
  (default) or `"matrix"`.

- alpha:

  the significance level used to compile the groups flagged as
  significantly different in the label attribute of the p-value matrix
  (default is `0.05`)

## Value

An object of class `"rankTest"` containing:

- res:

  comparison results. For `output="list"` a matrix with columns `W`, `z`
  and `pval`; for `output="matrix"` a many-to-one matrix of adjusted
  p-values.

- pmat:

  matrix of adjusted p-values.

- statistic:

  observed Steel test statistic.

- p.value:

  asymptotic p-value for the global Steel test.

Additional information is stored in attributes: `method`, `alternative`,
`output`, `main`, and `data.name`.

## Details

Performs Steel's many-to-one rank test for independent samples. The
procedure is the nonparametric analogue of Dunnett's test and controls
the family-wise error rate when comparing multiple treatment groups
against a common control.

The test statistic is based on pairwise Wilcoxon rank-sum statistics
comparing each treatment group with the control group. Ties are handled
using the asymptotic variance and covariance corrections described by
Scholz (2023). Adjusted p-values are obtained from the asymptotic
multivariate normal distribution of the standardized statistics.

If `x` is a list, its elements are taken as the samples to be compared
and hence have to be numeric data vectors. In this case, `g` is ignored
and one can simply use `steelTest(x)`.

Otherwise, `x` must be a numeric vector and `g` a grouping factor (or
vector coercible to a factor) of the same length.

Steel's test is the nonparametric analogue of Dunnett's test. Unlike
Dunn's or Conover's procedures, only comparisons against the designated
control group are performed.

The asymptotic covariance structure follows Scholz (2023), including tie
corrections for the Wilcoxon rank-sum statistics. P-values are obtained
from the multivariate normal distribution corresponding to the joint
asymptotic distribution of the standardized Steel statistics.

## References

Steel, R. G. D. (1959). A multiple comparison rank sum test: Treatments
versus control. *Biometrics*, **15**, 560–572.

Scholz, F. W. (2023). Improved tie correction methods for Steel-type
rank tests. *Journal of Nonparametric Statistics*, **35**, 541–563.

## See also

[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)

Other test.posthoc:
[`conoverTest()`](https://andrisignorell.github.io/lumen/reference/conoverTest.md),
[`dscfTest()`](https://andrisignorell.github.io/lumen/reference/dscfTest.md),
[`dunnTest()`](https://andrisignorell.github.io/lumen/reference/dunnTest.md),
[`dunnettTest()`](https://andrisignorell.github.io/lumen/reference/dunnettTest.md),
[`gamesHowellTest()`](https://andrisignorell.github.io/lumen/reference/gamesHowellTest.md),
[`nemenyiTest()`](https://andrisignorell.github.io/lumen/reference/nemenyiTest.md),
[`plot.PostHocTest()`](https://andrisignorell.github.io/lumen/reference/plot.PostHocTest.md),
[`postHoc`](https://andrisignorell.github.io/lumen/reference/postHoc.md),
[`scheffeTest()`](https://andrisignorell.github.io/lumen/reference/scheffeTest.md),
[`signifDiff()`](https://andrisignorell.github.io/lumen/reference/signifDiff.md)

## Examples

``` r
## Hollander & Wolfe example
x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
y <- c(3.8, 2.7, 4.0, 2.4)
z <- c(2.8, 3.4, 3.7, 2.2, 2.0)

steelTest(list(x, y, z))
#> 
#>  Steel test for multiple comparisons with a control 
#> 
#>      W          z   pval    
#> 2-1 12  0.4898979 0.8464    
#> 3-1 12 -0.1044466 0.9924    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

x <- c(x, y, z)
g <- factor(rep(1:3, c(5, 4, 5)))

steelTest(x, g)
#> 
#>  Steel test for multiple comparisons with a control 
#> 
#>      W          z   pval    
#> 2-1 12  0.4898979 0.8464    
#> 3-1 12 -0.1044466 0.9924    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

## Specify control group
steelTest(x, g, control = "1")
#> 
#>  Steel test for multiple comparisons with a control 
#> 
#>      W          z   pval    
#> 2-1 12  0.4898979 0.8464    
#> 3-1 12 -0.1044466 0.9924    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

## Formula interface
steelTest(Ozone ~ factor(Month), data = airquality)
#> 
#>  Steel test for multiple comparisons with a control 
#> 
#>         W        z    pval    
#> 6-5 152.0 1.321795 0.50052    
#> 7-5 566.5 4.183685 0.00012 ***
#> 8-5 548.5 3.854117 0.00043 ***
#> 9-5 470.0 1.568481 0.34258    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
```
