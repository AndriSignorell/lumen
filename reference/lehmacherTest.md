# Lehmacher's Test for Locating Asymmetry in Paired Multicategory Data

A nonparametric test for marginal homogeneity in square contingency
tables for dependent samples, based on a normal approximation of the
cell frequency differences.

## Usage

``` r
lehmacherTest(x, y = NULL, p.adjust.method = "hochberg")

# S3 method for class 'MHTest'
print(x, digits = 1L, ...)
```

## Arguments

- x:

  either a two-dimensional square contingency table in matrix form, or a
  factor object.

- y:

  a factor object; ignored if `x` is a matrix.

- p.adjust.method:

  the method used to adjust the per-category p-values for multiple
  comparisons, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default is
  `"hochberg"`, as recommended by Lehmacher (1980).

- digits:

  a non-null value for digits specifies the minimum number of
  significant digits to be printed. See
  [`print.default`](https://rdrr.io/r/base/print.default.html).

- ...:

  further arguments to be passed to or from other methods, ignored in
  this function.

## Value

A list with class `c("MHTest", "htest")` containing the following
components:

- statistic:

  a vector with the value of the test statistic for each category.

- parameter:

  the degrees of freedom, which is always 1.

- p.value:

  a vector with the p-values of the individual tests.

- p.value.corr:

  a vector with the adjusted p-values of the individual tests (see
  `p.adjust.method`).

- method:

  a character string indicating the test performed.

- data.name:

  a character string giving the name of the data.

## Details

Performs Lehmacher's chi-squared test for marginal homogeneity in a
square two-dimensional contingency table.

Unlike Bowker's test of symmetry, which tests whether \\P(i,j) =
P(j,i)\\ for every off-diagonal cell pair, Lehmacher's test addresses
marginal homogeneity: the null hypothesis is that, for every category
\\i\\, the row and column marginal probabilities agree, \\P(i \cdot) =
P(\cdot i)\\. One test statistic is computed per category and referred
to a chi-squared distribution with 1 degree of freedom; the resulting
p-values are adjusted for multiple comparisons.

If `x` is a matrix, it is taken as a two-dimensional contingency table,
and hence its entries should be nonnegative integers. Otherwise, both
`x` and `y` must be vectors or factors of the same length. Incomplete
cases are removed, vectors are coerced into factors, and the contingency
table is computed from these.

## References

Lehmacher, W. (1980) Simultaneous sign tests for marginal homogeneity of
square contingency tables. *Biometrical Journal*, 22 (8), 795-798.

## See also

[`mcnemar.test()`](https://rdrr.io/r/stats/mcnemar.test.html) for the
2x2 case

Other test.categorical: [`barnardTest()`](barnardTest.md),
[`bhapkarTest()`](bhapkarTest.md),
[`breslowDayTest()`](breslowDayTest.md),
[`cochranQTest()`](cochranQTest.md), [`gTest()`](gTest.md),
[`mantelTrendTest()`](mantelTrendTest.md),
[`stuartMaxwellTest()`](stuartMaxwellTest.md),
[`woolfTest()`](woolfTest.md)

## Examples

``` r
x <- matrix(c(400, 40, 20, 10,
              50, 300, 60, 20,
              10, 40, 120, 5,
              5, 90, 50, 80), nrow = 4, byrow = TRUE,
            dimnames = list(LETTERS[1:4], LETTERS[1:4]))

lehmacherTest(x)
#> 
#>  Lehmacher test for marginal homogeneity
#> 
#> data:  x
#> 
#>     Chi²   p-value     p-adj      
#> A    0.2     0.667     0.667      
#> B    5.3     0.021     0.042   *  
#> C   30.4   < 0.001   < 0.001   ***
#> D   67.2   < 0.001   < 0.001   ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
```
