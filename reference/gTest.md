# G-Test for Likelihood-Ratio Testing of Count Data

A goodness-of-fit or test of independence based on the log-likelihood
ratio (G-statistic), serving as an asymptotically equivalent alternative
to the chi-squared test.

## Usage

``` r
gTest(
  x,
  y = NULL,
  correct = c("none", "williams", "yates"),
  p = rep(1/length(x), length(x)),
  rescaleP = FALSE
)
```

## Arguments

- x:

  a numeric vector or matrix. `x` and `y` can also both be factors.

- y:

  a numeric vector; ignored if `x` is a matrix. If `x` is a factor, `y`
  should be a factor of the same length.

- correct:

  the correction to be applied, one of `"none"` (default), `"williams"`
  or `"yates"`. See the Details.

- p:

  a vector of probabilities of the same length as `x` (goodness-of-fit
  test only). An error is given if any entry of `p` is negative.

- rescaleP:

  logical; if `TRUE` then `p` is rescaled (if necessary) to sum to 1. If
  `rescaleP` is `FALSE`, and `p` does not sum to 1, an error is given.

## Value

A list with class `"htest"` containing the following components:

- statistic:

  the value of the G test statistic.

- parameter:

  the degrees of freedom of the approximate chi-squared distribution of
  the test statistic.

- p.value:

  the p-value of the test.

- method:

  a character string indicating the type of test performed, and whether
  a correction was used.

- data.name:

  a character string giving the name(s) of the data.

- observed:

  the observed counts (before any continuity correction).

- expected:

  the expected counts under the null hypothesis.

## Details

`gTest` performs log-likelihood ratio contingency table tests and
goodness-of-fit tests.

The G-test is also called "Likelihood Ratio Test" and is asymptotically
equivalent to the Pearson chi-squared test but not usually used when
analyzing 2x2 tables. It is used in logistic regression and loglinear
modeling which involves contingency tables.

If `x` is a matrix with one row or column, or if `x` is a vector and `y`
is not given, then a *goodness-of-fit test* is performed (`x` is treated
as a one-dimensional contingency table). The entries of `x` must be
non-negative integers. In this case, the hypothesis tested is whether
the population probabilities equal those in `p`, or are all equal if `p`
is not given.

If `x` is a matrix with at least two rows and columns, it is taken as a
two-dimensional contingency table: the entries of `x` must be
non-negative integers. Otherwise, `x` and `y` must be vectors or factors
of the same length; cases with missing values are removed, the objects
are coerced to factors, and the contingency table is computed from
these. Then the G-test is performed on the null hypothesis that the
joint distribution of the cell counts in a 2-dimensional contingency
table is the product of the row and column marginals.

Williams' correction (Williams, 1976) divides the statistic by a factor
\\q \> 1\\ and can be used for both test types. Yates' continuity
correction is only defined for 2x2 tables (independence) resp. two data
values (goodness-of-fit).

## Note

Based on code by Pete Hurd, adapted to conform to package standards.

## References

Agresti, A. (2007) *An Introduction to Categorical Data Analysis*, 2nd
ed., New York: John Wiley & Sons. Page 38.

Sokal, R. R. and Rohlf, F. J. (2012) *Biometry: The Principles and
Practice of Statistics in Biological Research*, 4th ed., New York: W. H.
Freeman and Co.

Williams, D. A. (1976) Improved likelihood ratio tests for complete
contingency tables. *Biometrika*, 63, 33-37.

## See also

[`chisq.test`](https://rdrr.io/r/stats/chisq.test.html)

Other test.categorical: [`barnardTest()`](barnardTest.md),
[`bhapkarTest()`](bhapkarTest.md),
[`breslowDayTest()`](breslowDayTest.md),
[`cochranQTest()`](cochranQTest.md),
[`lehmacherTest()`](lehmacherTest.md),
[`mantelTrendTest()`](mantelTrendTest.md),
[`stuartMaxwellTest()`](stuartMaxwellTest.md),
[`woolfTest()`](woolfTest.md)

## Examples

``` r
## From Agresti (2007), p. 39
M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
dimnames(M) <- list(gender = c("M", "F"),
                    party  = c("Democrat", "Independent", "Republican"))

(Xsq <- gTest(M))   # Prints test summary
#> 
#>  Log likelihood ratio (G-test) test of independence without correction
#> 
#> data:  M
#> G = 30.017, df = 2, p-value = 3.034e-07
#> 

Xsq$observed        # observed counts (same as M)
#>       party
#> gender Democrat Independent Republican
#>      M      762         327        468
#>      F      484         239        477
Xsq$expected        # expected counts under the null
#>   Democrat Independent Republican
#> M 703.6714    319.6453   533.6834
#> F 542.3286    246.3547   411.3166


## Testing for population probabilities
## Case A. Tabulated data
x <- c(A = 20, B = 15, C = 25)
gTest(x)
#> 
#>  Log likelihood ratio (G-test) goodness of fit test
#> 
#> data:  x
#> G = 2.5267, df = 2, p-value = 0.2827
#> 
gTest(as.table(x))             # the same
#> 
#>  Log likelihood ratio (G-test) goodness of fit test
#> 
#> data:  as.table(x)
#> G = 2.5267, df = 2, p-value = 0.2827
#> 
x <- c(89, 37, 30, 28, 2)
p <- c(40, 20, 20, 15, 5)
try(
gTest(x, p = p)                # gives an error
)
#> Error in gTest(x, p = p) : probabilities must sum to 1
# works
p <- c(0.40, 0.20, 0.20, 0.19, 0.01)
# Expected count in category 5
# is 1.86 < 5 ==> chi square approx.
gTest(x, p = p)                # maybe doubtful, but is ok!
#> 
#>  Log likelihood ratio (G-test) goodness of fit test
#> 
#> data:  x
#> G = 5.8414, df = 4, p-value = 0.2113
#> 

## Case B. Raw data
x <- trunc(5 * runif(100))
gTest(table(x))                # NOT 'gTest(x)'!
#> 
#>  Log likelihood ratio (G-test) goodness of fit test
#> 
#> data:  table(x)
#> G = 7.3121, df = 4, p-value = 0.1203
#> 
```
