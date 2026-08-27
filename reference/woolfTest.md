# Woolf Test for Comparing Odds Ratios Across Stratified 2×2 Tables

A test for homogeneity of odds ratios across several 2×2 contingency
tables, similar to the Breslow-Day test but based on a different test
statistic.

## Usage

``` r
woolfTest(x)
```

## Arguments

- x:

  a \\2 \times 2 \times k\\ table, where the last dimension refers to
  the strata.

## Value

A list of class `"htest"` containing the following components:

- statistic:

  the chi-squared test statistic.

- parameter:

  degrees of freedom of the approximate chi-squared distribution of the
  test statistic.

- p.value:

  \\p\\-value for the test.

- method:

  a character string indicating the type of test performed.

- data.name:

  a character string giving the name(s) of the data.

- observed:

  the per-stratum log odds ratios (not the raw table counts).

- expected:

  the inverse-variance-weighted mean log odds ratio across strata.

## Details

Test for homogeneity on \\2 \times 2 \times k\\ tables over strata
(i.e., whether the log odds ratios are the same in all strata).

## Note

Based on code by David Meyer, Achim Zeileis, Kurt Hornik, Michael
Friendly previously published as `woolf_test()` in the vcd package,
adapted to conform to package standards.

## References

Woolf, B. 1955: On estimating the relation between blood group and
disease. *Ann. Human Genet.* (London) **19**, 251-253.

## See also

[`mantelhaen.test()`](https://rdrr.io/r/stats/mantelhaen.test.html)

Other test.categorical:
[`barnardTest()`](https://andrisignorell.github.io/lumen/reference/barnardTest.md),
[`bhapkarTest()`](https://andrisignorell.github.io/lumen/reference/bhapkarTest.md),
[`breslowDayTest()`](https://andrisignorell.github.io/lumen/reference/breslowDayTest.md),
[`cochranQTest()`](https://andrisignorell.github.io/lumen/reference/cochranQTest.md),
[`gTest()`](https://andrisignorell.github.io/lumen/reference/gTest.md),
[`lehmacherTest()`](https://andrisignorell.github.io/lumen/reference/lehmacherTest.md),
[`mantelTrendTest()`](https://andrisignorell.github.io/lumen/reference/mantelTrendTest.md),
[`stuartMaxwellTest()`](https://andrisignorell.github.io/lumen/reference/stuartMaxwellTest.md)

## Examples

``` r

migraine <- xtabs(freq ~ .,
            cbind(expand.grid(treatment=c("active","placebo"),
                               response=c("better","same"),
                               gender=c("female","male")),
                  freq=c(16,5,11,20,12,7,16,19))
            )

woolfTest(migraine)
#> 
#>  Woolf Test on Homogeneity of Odds Ratios (no 3-Way assoc.)
#> 
#> data:  migraine
#> X-squared = 1.4808, df = 1, p-value = 0.2236
#> 
```
