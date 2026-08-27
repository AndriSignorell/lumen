# Dwass-Steel-Critchlow-Fligner Test for Nonparametric All-Pairs Group Comparisons

Performs the Dwass-Steel-Critchlow-Fligner (DSCF) all-pairs
multiple-comparison procedure for independent samples.

## Usage

``` r
dscfTest(x, ...)

# S3 method for class 'formula'
dscfTest(formula, data, subset, na.action, ...)

# Default S3 method
dscfTest(x, g, output = c("list", "matrix"), alpha = 0.05, ...)
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

  comparison results. For `output="list"` a matrix with columns `z` and
  `pval`; for `output="matrix"` a symmetric matrix of adjusted p-values
  with diagonal 1.

- pmat:

  symmetric matrix of adjusted p-values with diagonal 1.

Additional information is stored in attributes: `method`, `output`,
`main`, and `data.name`.

## Details

The DSCF test is a nonparametric all-pairs comparison procedure based on
pairwise Wilcoxon rank-sum statistics and provides strong control of the
family-wise error rate. It may be viewed as a nonparametric analogue of
Tukey's all-pairs procedure and is generally more powerful than the
classical Nemenyi test.

For each pair of groups, observations are ranked jointly within the two
groups being compared, yielding a pairwise Mann-Whitney statistic. Test
statistics are transformed according to the
Dwass-Steel-Critchlow-Fligner procedure and evaluated using the
studentized range distribution.

If `x` is a list, its elements are taken as the samples to be compared
and hence have to be numeric data vectors. In this case, `g` is ignored
and one can simply use `dscfTest(x)`.

Otherwise, `x` must be a numeric vector and `g` a grouping factor (or
vector coercible to a factor) of the same length.

The DSCF procedure performs all pairwise comparisons using pair-specific
rankings. For each pair of groups, a Mann-Whitney-type statistic is
computed and standardized using a tie-corrected variance estimate.
Adjusted p-values are obtained from the studentized range distribution
and control the family-wise error rate across all pairwise comparisons.

The implementation reproduces the procedure described by Dwass (1960),
Steel (1960), Critchlow and Fligner (1991), and is equivalent to the
implementation in `PMCMRplus::dscfAllPairsTest()`.

## References

Dwass, M. (1960). Some k-sample rank-order tests. In: *Contributions to
Probability and Statistics*, Stanford University Press, 198–202.

Steel, R. G. D. (1960). A rank sum test for comparing all pairs of
treatments. *Technometrics*, **2**, 197–207.

Critchlow, D. E., & Fligner, M. A. (1991). On distribution-free multiple
comparisons in the one-way analysis of variance. *Communications in
Statistics - Theory and Methods*, **20**, 127–139.

## See also

[`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)

Other test.posthoc: [`conoverTest()`](conoverTest.md),
[`dunnTest()`](dunnTest.md), [`dunnettTest()`](dunnettTest.md),
[`gamesHowellTest()`](gamesHowellTest.md),
[`nemenyiTest()`](nemenyiTest.md),
[`plot.PostHocTest()`](plot.PostHocTest.md), [`postHoc`](postHoc.md),
[`scheffeTest()`](scheffeTest.md), [`signifDiff()`](signifDiff.md),
[`steelTest()`](steelTest.md)

## Examples

``` r
## Hollander & Wolfe style example
x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
y <- c(3.8, 2.7, 4.0, 2.4)
z <- c(2.8, 3.4, 3.7, 2.2, 2.0)

dscfTest(list(x, y, z))
#> 
#>  Steel-Dwass-Critchlow-Fligner all-pairs test 
#> 
#>              z   pval    
#> 1-2 -0.6928203 0.8761    
#> 1-3 -0.1477098 0.9940    
#> 2-3 -1.3856406 0.5897    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

x <- c(x, y, z)
g <- factor(rep(1:3, c(5, 4, 5)))

dscfTest(x, g)
#> 
#>  Steel-Dwass-Critchlow-Fligner all-pairs test 
#> 
#>              z   pval    
#> 1-2 -0.6928203 0.8761    
#> 1-3 -0.1477098 0.9940    
#> 2-3 -1.3856406 0.5897    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

## Matrix output
dscfTest(x, g, output = "matrix")
#> 
#>  Steel-Dwass-Critchlow-Fligner all-pairs test 
#> 
#>   1    2    3   
#> 1 1.00 0.88 0.99
#> 2 0.88 1.00 0.59
#> 3 0.99 0.59 1.00
#> attr(,"lbl")
#> 1 2 3 
#>       
#> 

## Formula interface
dscfTest(Ozone ~ factor(Month), data = airquality)
#> 
#>  Steel-Dwass-Critchlow-Fligner all-pairs test 
#> 
#>               z    pval    
#> 5-6 -1.86973366 0.67741    
#> 5-7 -5.91550472 0.00028 ***
#> 5-8 -5.44986225 0.00110 ** 
#> 5-9 -2.21910256 0.51731    
#> 6-7 -3.49686607 0.09686 .  
#> 6-8 -3.17698764 0.16274    
#> 6-9 -0.09724737 0.99999    
#> 7-8 -0.25886212 0.99975    
#> 7-9 -4.78168939 0.00648 ** 
#> 8-9 -4.17429785 0.02624 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
```
