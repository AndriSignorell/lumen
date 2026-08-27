# Stuart-Maxwell Test for Testing Marginal Homogeneity in Paired Multicategory Data

A nonparametric test for marginal homogeneity in square contingency
tables for dependent samples, generalizing the McNemar test to more than
two categories.

## Usage

``` r
stuartMaxwellTest(x, y = NULL)
```

## Arguments

- x:

  either a 2-way \\k \times k\\ contingency table in matrix form, or a
  factor.

- y:

  a factor with the same levels as x; ignored if x is a matrix.

## Value

A list with class `"htest"` containing the following components:

- statistic:

  the value of the test statistic.

- parameter:

  the degrees of freedom.

- p.value:

  the p-value of the test.

- method:

  a character string indicating what type of test was performed.

- data.name:

  a character string giving the name of the data.

## Details

This function computes the marginal homogeneity test for a \\k \times
k\\ matrix of assignments of objects to `k` categories or two vectors
`x`, `y` of category scores for `n` data objects by two raters. The
statistic is distributed as \\\chi^2\\ with `k-1` degrees of freedom.  
It can be viewed as an extension of the McNemar test to \\k \times k\\
table.

The null is that the probabilities of being classified into cells
`[i, j]` and `[j, i]` are the same.

If `x` is a matrix, it is taken as a two-dimensional contingency table,
and hence its entries should be nonnegative integers. Otherwise, both x
and y must be vectors or factors of the same length and with the same
levels.  
Incomplete cases are removed, vectors are coerced into factors, and the
contingency table is computed from these.

If there is perfect agreement for any category k, that category must be
omitted in order to invert matrix S.

If for any category `k`, all frequencies in row `k` and column `k` are
0, except possibly for the main diagonal element (e.g., for perfect
agreement for category `k`, in such cases also the corresponding row and
column marginal frequencies would be equal), then the category is not
included in the test and should be ignored, say the Stuart-Maxwell test
is performed with respect to the remaining categories only. The degree
of freedom `df` in this case can still be considered `k - 1`, where `k`
is the number of original categories; this treats omitted categories as
if they were included but contributed 0 to the value of \\\chi^2\\ - a
reasonable view since such categories have equal row and column
marginals. (See:
<https://www.john-uebersax.com/stat/mcnemar.htm#stuart>)

## Note

Based on code from Jim Lemon, adapted to conform to package standards.

## References

Stuart, A (1955) A test for homogeneity of the marginal distributions in
a two-way classification. *Biometrika*, 42, 412-416.

Maxwell, A.E. (1970) Comparing the classification of subjects by two
independent judges. *British Journal of Psychiatry*, 116, 651-655.

Agresti, A. (2002) *Categorical Data Analysis*. John Wiley & Sons, pp 86
ff.

## See also

[`mcnemar.test()`](https://rdrr.io/r/stats/mcnemar.test.html),
[`chisq.test()`](https://rdrr.io/r/stats/chisq.test.html)

Other test.categorical: [`barnardTest()`](barnardTest.md),
[`bhapkarTest()`](bhapkarTest.md),
[`breslowDayTest()`](breslowDayTest.md),
[`cochranQTest()`](cochranQTest.md), [`gTest()`](gTest.md),
[`lehmacherTest()`](lehmacherTest.md),
[`mantelTrendTest()`](mantelTrendTest.md), [`woolfTest()`](woolfTest.md)

## Examples

``` r

# Source: https://john-uebersax.com/stat/mcnemar.htm#stuart
hyp <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow=3))
stuartMaxwellTest(hyp)
#> 
#>  Stuart-Maxwell test for marginal homogeneity
#> 
#> data:  x
#> chi-squared = 13.765, df = 2, p-value = 0.001026
#> 

# same as defined with two vectors
d.hyp <- expand.grid(c("A","B","C"), c("A","B","C"))[
                     rep(1:9, times = hyp), ]
row.names(d.hyp) <- NULL

stuartMaxwellTest(x=d.hyp[,1], y=d.hyp[,2])
#> 
#>  Stuart-Maxwell test for marginal homogeneity
#> 
#> data:  x and y
#> chi-squared = 13.765, df = 2, p-value = 0.001026
#> 


mc <- as.table(matrix(c(
         732, 1524, 1575, 1577, 1602, 837, 1554, 1437, 
         1672, 1600, 841, 1363, 1385, 1484, 1524, 791), nrow=4))

stuartMaxwellTest(mc)
#> 
#>  Stuart-Maxwell test for marginal homogeneity
#> 
#> data:  x
#> chi-squared = 0.089722, df = 3, p-value = 0.993
#> 

```
