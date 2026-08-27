# Power Calculations for Power and Sample-Size Calculations for Chi-Square Tests

Compute power of test or determine parameters to obtain target power
(same as
[`power.anova.test`](https://rdrr.io/r/stats/power.anova.test.html)).

## Usage

``` r
powerChisqTest(
  n = NULL,
  effectSize = NULL,
  df = NULL,
  sig.level = 0.05,
  power = NULL
)
```

## Arguments

- n:

  total number of observations.

- effectSize:

  effect size.

- df:

  degrees of freedom of the chi-squared distribution, e.g.
  `(rows-1)*(cols-1)` for a test of independence. Must always be
  supplied.

- sig.level:

  significance level (Type I error probability).

- power:

  target power (1 minus Type II error probability).

## Value

Object of class "power.htest", a list of the arguments (including the
computed one) augmented with 'method' and 'note' elements.

## Details

Exactly one of the parameters `effectSize`, `n`, `power` or `sig.level`
must be passed as NULL, and this parameter is determined from the
others. Note that the last one has non-NULL default, so `NULL` must be
explicitly passed, if you want to compute it. `df` must always be
supplied; it cannot be solved for.

## Note

[`uniroot`](https://rdrr.io/r/stats/uniroot.html) is used to solve power
equation for unknowns, so you may see errors from it, notably about
inability to bracket the root when invalid arguments are given.

Based on code by Stephane Champely, and Peter Dalgaard, adapted to
conform to package standards.

## References

Cohen, J. (1988) *Statistical power analysis for the behavioral sciences
(2nd ed.)* Hillsdale, NJ: Lawrence Erlbaum.

## See also

[`power.t.test`](https://rdrr.io/r/stats/power.t.test.html)

## Examples

``` r

## Exercise 7.1 P. 249 from Cohen (1988) 
powerChisqTest(effectSize=0.289, df=(4-1), n=100, sig.level=0.05)
#> 
#>      Chi squared power calculation 
#> 
#>      effectSize = 0.289
#>               n = 100
#>              df = 3
#>       sig.level = 0.05
#>           power = 0.6750777
#> 
#> NOTE: n is the number of observations
#> 

## Exercise 7.3 p. 251
powerChisqTest(effectSize=0.346, df=(2-1)*(3-1), n=140, sig.level=0.01)
#> 
#>      Chi squared power calculation 
#> 
#>      effectSize = 0.346
#>               n = 140
#>              df = 2
#>       sig.level = 0.01
#>           power = 0.8854053
#> 
#> NOTE: n is the number of observations
#> 

## Exercise 7.8 p. 270
powerChisqTest(effectSize=0.1, df=(5-1)*(6-1), power=0.80, sig.level=0.05)
#> 
#>      Chi squared power calculation 
#> 
#>      effectSize = 0.1
#>               n = 2096.079
#>              df = 20
#>       sig.level = 0.05
#>           power = 0.8
#> 
#> NOTE: n is the number of observations
#> 
#' @family power
```
