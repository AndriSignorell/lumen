# Bhapkar Marginal Homogeneity Test for Comparing Paired Marginal Distributions Across Multiple Categories

An asymptotic chi-squared test for marginal homogeneity in square
contingency tables for dependent samples, similar to the Stuart-Maxwell
test but based on a different test statistic and generally slightly more
powerful.

## Usage

``` r
bhapkarTest(x, y = NULL)
```

## Arguments

- x:

  either a 2-way \\k \times k\\ contingency table in matrix form, or a
  factor.

- y:

  a factor with the same levels as `x`; ignored if `x` is a matrix.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of the chi-squared test statistic.

- `parameter`:

  the degrees of freedom of the approximate chi-squared distribution of
  the test statistic.

- `p.value`:

  the p-value of the test.

- `method`:

  a character string indicating the test performed.

- `data.name`:

  a character string giving the name of the data.

## Details

Bhapkar's test (Bhapkar, 1966) is used to assess marginal homogeneity in
square contingency tables. It is based on the asymptotic normality of
marginal proportions and is closely related to the generalized McNemar
test, as implemented in
[`stuartMaxwellTest`](https://andrisignorell.github.io/lumen/reference/stuartMaxwellTest.md).

The two tests differ only in the estimation of the variance-covariance
matrix of the marginal proportions and are asymptotically equivalent
(Keefe, 1982), meaning that for large sample sizes they yield the same
chi-squared statistic. For finite samples, however, the Bhapkar test is
generally more powerful and is therefore preferred in practice.

## References

Bhapkar V.P. (1966) A note on the equivalence of two test criteria for
hypotheses in categorical data. *Journal of the American Statistical
Association*, 61: 228-235.

Ireland C.T., Ku H.H., and Kullback S. (1969) Symmetry and marginal
homogeneity of an r x r contingency table. *Journal of the American
Statistical Association*, 64: 1323-1341.

Keefe T.J. (1982) On the relationship between two tests for homogeneity
of the marginal distributions in a two-way classification. *Biometrika*,
69: 683-684.

Sun X., Yang Z. (2008) Generalized McNemar's Test for Homogeneity of the
Marginal Distributions. *SAS Global Forum 2008: Statistics and Data
Analysis*, Paper 382-208.

## See also

[mcnemar.test](https://rdrr.io/r/stats/mcnemar.test.html),[chisq.test](https://rdrr.io/r/stats/chisq.test.html)

Other test.categorical:
[`barnardTest()`](https://andrisignorell.github.io/lumen/reference/barnardTest.md),
[`breslowDayTest()`](https://andrisignorell.github.io/lumen/reference/breslowDayTest.md),
[`cochranQTest()`](https://andrisignorell.github.io/lumen/reference/cochranQTest.md),
[`gTest()`](https://andrisignorell.github.io/lumen/reference/gTest.md),
[`lehmacherTest()`](https://andrisignorell.github.io/lumen/reference/lehmacherTest.md),
[`mantelTrendTest()`](https://andrisignorell.github.io/lumen/reference/mantelTrendTest.md),
[`stuartMaxwellTest()`](https://andrisignorell.github.io/lumen/reference/stuartMaxwellTest.md),
[`woolfTest()`](https://andrisignorell.github.io/lumen/reference/woolfTest.md)

## Examples

``` r
# Source: https://john-uebersax.com/stat/mcnemar.htm#bhapkar
mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow=3))

bhapkarTest(mc)
#> 
#>  Bhapkar test for marginal homogeneity
#> 
#> data:  mc
#> chi-squared = 15.423, df = 2, p-value = 0.0004476
#> 
```
