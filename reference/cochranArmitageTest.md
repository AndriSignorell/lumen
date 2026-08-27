# Cochran-Armitage Test for Detecting Trends in Proportions Across Ordered Groups

A test for linear trend in proportions across ordered categories in \\2
\times k\\ contingency tables, typically used in epidemiological
dose-response analyses.

## Usage

``` r
cochranArmitageTest(
  x,
  alternative = c("two.sided", "increasing", "decreasing")
)
```

## Arguments

- x:

  a \\k \times 2\\ or \\2 \times k\\ frequency table or matrix with
  nonnegative counts.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"increasing"` or `"decreasing"`. You can
  specify just the initial letter. See the Details for the direction
  convention.

## Value

A list of class `"htest"`, containing the following components:

- `statistic`:

  the z-statistic of the test.

- `parameter`:

  the number of levels of the ordered variable (named `dim`).

- `p.value`:

  the p-value for the test.

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  the character string "Cochran-Armitage test for trend".

- `data.name`:

  a character string giving the names of the data.

## Details

Perform a Cochran-Armitage test for trend in binomial proportions across
the levels of a single variable. This test is appropriate only when one
variable has two levels and the other variable is ordinal. The two-level
variable represents the response, and the other represents an
explanatory variable with ordered levels. The null hypothesis is the
hypothesis of no trend, which means that the binomial proportion is the
same for all levels of the explanatory variable.

The table is oriented such that the ordered explanatory variable is in
the rows and the binary response in the columns (a \\2 \times k\\ table
is transposed accordingly). Row scores are taken from the row dimnames
if these are numeric, and are sequential integers otherwise.

The alternatives refer to the trend in the proportion of the *first*
response column: `"increasing"` tests whether this proportion increases
with the ordered levels, `"decreasing"` tests the reverse.

## Note

Based on code by Eric Lecoutre, adapted to conform to package standards.

<https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html>

Results are consistent with SAS PROC FREQ. They may differ slightly from
coin's `independence_test(..., teststat = "scalar")`, which uses a
different variance formula.

## References

Agresti, A. (2002) *Categorical Data Analysis*. John Wiley & Sons.

## See also

[`prop.trend.test`](https://rdrr.io/r/stats/prop.trend.test.html), [SAS
PROC FREQ
documentation](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.5/procstat/procstat_freq_details76.htm)

Other test.trend:
[`coxStuartTest()`](https://andrisignorell.github.io/lumen/reference/coxStuartTest.md),
[`jonckheereTerpstraTest()`](https://andrisignorell.github.io/lumen/reference/jonckheereTerpstraTest.md),
[`pageTest()`](https://andrisignorell.github.io/lumen/reference/pageTest.md)

## Examples

``` r
# http://www.lexjansen.com/pharmasug/2007/sp/sp05.pdf, pp. 4
dose <- matrix(c(10,9,10,7, 0,1,0,3), byrow=TRUE, nrow=2,
               dimnames=list(resp=0:1, dose=0:3))

cochranArmitageTest(dose)
#> 
#>  Cochran-Armitage test for trend
#> 
#> data:  dose
#> Z = -1.8856, dim = 4, p-value = 0.05935
#> alternative hypothesis: two.sided
#> 
cochranArmitageTest(dose, alternative = "increasing")
#> 
#>  Cochran-Armitage test for trend
#> 
#> data:  dose
#> Z = -1.8856, dim = 4, p-value = 0.9703
#> alternative hypothesis: increasing
#> 

# compare with coin::independence_test(..., teststat = "scalar")
# (see the Note on the variance formula)
lungtumor <- data.frame(dose = rep(c(0, 1, 2), c(40, 50, 48)),
                        tumor = c(rep(c(0, 1), c(38, 2)),
                                  rep(c(0, 1), c(43, 7)),
                                  rep(c(0, 1), c(33, 15))))
tab <- table(lungtumor$dose, lungtumor$tumor)
cochranArmitageTest(tab)
#> 
#>  Cochran-Armitage test for trend
#> 
#> data:  tab
#> Z = -3.2735, dim = 3, p-value = 0.001062
#> alternative hypothesis: two.sided
#> 

# similar to prop.trend.test (uses integer scores 1..k instead of dimnames)
prop.trend.test(tab[,1], apply(tab, 1, sum))
#> 
#>  Chi-squared Test for Trend in Proportions
#> 
#> data:  tab[, 1] out of apply(tab, 1, sum) ,
#>  using scores: 1 2 3
#> X-squared = 10.716, df = 1, p-value = 0.001062
#> 

# SAS PROC FREQ reference
# https://support.sas.com/documentation/onlinedoc/stat/142/freq.pdf, pp 2868
pain <- structure(c(26, 6, 26, 7, 23, 9, 18, 14, 9, 23),
                  dim = c(2L, 5L),
                  dimnames = list(adverse = c("No", "Yes"),
                                  dose    = c("0", "1", "2", "3", "4")),
                  class = "table")

cochranArmitageTest(pain)
#> 
#>  Cochran-Armitage test for trend
#> 
#> data:  pain
#> Z = -4.7918, dim = 5, p-value = 1.653e-06
#> alternative hypothesis: two.sided
#> 
```
