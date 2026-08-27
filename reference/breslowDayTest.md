# Breslow-Day Test for Comparing Odds Ratios Across Strata

A test for homogeneity of odds ratios across several 2x2 contingency
tables (strata), commonly used to verify whether a confounding effect is
constant across subgroups, as required by the Mantel-Haenszel method.

## Usage

``` r
breslowDayTest(x, OR = NULL, correct = FALSE)
```

## Arguments

- x:

  a \\2 \times 2 \times k\\ table.

- OR:

  the odds ratio to be tested against. If left undefined (default) the
  Mantel-Haenszel estimate will be used.

- correct:

  logical, if `TRUE` the Breslow-Day test with Tarone's adjustment is
  computed, which subtracts an adjustment factor to make the resulting
  statistic asymptotically chi-squared.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of the chi-squared test statistic.

- `parameter`:

  the degrees of freedom of the approximate chi-squared distribution of
  the test statistic (\\k-1\\).

- `p.value`:

  the p-value of the test.

- `method`:

  a character string indicating the test performed.

- `data.name`:

  a character string giving the name of the data.

- `n`:

  the total number of observations (not shown on screen).

## Details

Calculates the Breslow-Day test of homogeneity for a \\2 \times 2 \times
k\\ table, in order to investigate if all \\k\\ strata have the same OR.
If `OR` is not given, the Mantel-Haenszel estimate is used.

For the Breslow-Day test to be valid, the sample size should be
relatively large in each stratum, and at least 80% of the expected cell
counts should be greater than 5. Note that this is a stricter sample
size requirement than the requirement for the Cochran-Mantel-Haenszel
test for tables, in that each stratum sample size (not just the overall
sample size) must be relatively large. Even when the Breslow-Day test is
valid, it might not be very powerful against certain alternatives, as
discussed in Breslow and Day (1980).

The statistic is referred to a chi-squared distribution with \\k-1\\
degrees of freedom; this also applies when a prespecified `OR` is
supplied. Note that Tarone's adjustment is derived for the
Mantel-Haenszel estimate; a warning is issued if `correct = TRUE` is
combined with a user-supplied `OR`.

Alternatively, it might be better to cast the entire inference problem
into the setting of a logistic regression model. Here, the underlying
question of the Breslow-Day test can be answered by investigating
whether an interaction term with the strata variable is necessary (e.g.
using a likelihood ratio test using the `anova` function).

## Note

Based on code by Michael Hoehle, adapted to conform to package
standards.

## References

Breslow, N. E. and Day, N. E. (1980) The Analysis of Case-Control
Studies. *Statistical Methods in Cancer Research: Vol. 1*. Lyon, France,
IARC Scientific Publications.

Tarone, R.E. (1985) On heterogeneity tests based on efficient scores,
*Biometrika*, 72, pp. 91-95.

Jones, M. P., O'Gorman, T. W., Lemka, J. H., and Woolson, R. F. (1989) A
Monte Carlo Investigation of Homogeneity Tests of the Odds Ratio Under
Various Sample Size Configurations. *Biometrics*, 45, 171-181.

Breslow, N. E. (1996) Statistics in Epidemiology: The Case-Control
Study. *Journal of the American Statistical Association*, 91, 14-26.

## See also

[`mantelhaen.test`](https://rdrr.io/r/stats/mantelhaen.test.html)

Other test.categorical:
[`barnardTest()`](https://andrisignorell.github.io/lumen/reference/barnardTest.md),
[`bhapkarTest()`](https://andrisignorell.github.io/lumen/reference/bhapkarTest.md),
[`cochranQTest()`](https://andrisignorell.github.io/lumen/reference/cochranQTest.md),
[`gTest()`](https://andrisignorell.github.io/lumen/reference/gTest.md),
[`lehmacherTest()`](https://andrisignorell.github.io/lumen/reference/lehmacherTest.md),
[`mantelTrendTest()`](https://andrisignorell.github.io/lumen/reference/mantelTrendTest.md),
[`stuartMaxwellTest()`](https://andrisignorell.github.io/lumen/reference/stuartMaxwellTest.md),
[`woolfTest()`](https://andrisignorell.github.io/lumen/reference/woolfTest.md)

## Examples

``` r
migraine <- xtabs(freq ~ .,
            cbind(expand.grid(treatment=c("active", "placebo"),
                              response =c("better", "same"),
                              gender   =c("female", "male")),
                  freq=c(16, 5, 11, 20, 12, 7, 16, 19))
            )

# get rid of gender
tab <- xtabs(Freq ~ treatment + response, migraine)
tab
#>          response
#> treatment better same
#>   active      28   27
#>   placebo     12   39

# only the women
female <- migraine[,, 1]
female
#>          response
#> treatment better same
#>   active      16   11
#>   placebo      5   20

# .. and the men
male <- migraine[,, 2]
male
#>          response
#> treatment better same
#>   active      12   16
#>   placebo      7   19

breslowDayTest(migraine)
#> 
#>  Breslow-Day test for homogeneity of the odds ratios
#> 
#> data:  migraine
#> X-squared = 1.4929, df = 1, p-value = 0.2218
#> 
breslowDayTest(migraine, correct = TRUE)
#> 
#>  Breslow-Day test for homogeneity of the odds ratios (Tarone corrected)
#> 
#> data:  migraine
#> X-squared = 1.4905, df = 1, p-value = 0.2221
#> 

salary <- array(
      c(38, 12, 102, 141, 12, 9, 136, 383),
      dim=c(2, 2, 2),
      dimnames=list(exposure=c("exposed", "not"),
                    disease =c("case", "control"),
                    salary  =c("<1000", ">=1000"))
                    )

# common odds ratio = 4.028269
breslowDayTest(salary, OR = 4.02)
#> 
#>  Breslow-Day test for homogeneity of the odds ratios
#> 
#> data:  salary
#> X-squared = 0.080143, df = 1, p-value = 0.7771
#> 
```
