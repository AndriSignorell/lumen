# Hosmer-Lemeshow Goodness of Fit Test for Assessing Logistic Regression Calibration

Computes the Hosmer-Lemeshow C or H goodness-of-fit test for a logistic
regression model, assessing whether observed event rates match predicted
probabilities across grouped subsets of the data.

## Usage

``` r
hosmerLemeshowTest(x, ...)

# S3 method for class 'glm'
hosmerLemeshowTest(x, nGroups = 10, type = c("C", "H"), ...)

# Default S3 method
hosmerLemeshowTest(x, obs, nGroups = 10, type = c("C", "H"), ...)

# S3 method for class 'HosmerLemeshowTest'
print(x, digits = 4, details = FALSE, ...)
```

## Arguments

- x:

  either a numeric vector of fitted probabilities, each in \\\[0, 1\]\\
  and without missing values, or a fitted binomial
  [`glm`](https://rdrr.io/r/stats/glm.html) object, from which fitted
  probabilities and observed outcomes are extracted.

- ...:

  further arguments passed to methods.

- nGroups:

  integer, the number of groups (default is `10`). Must be at least 3.

- type:

  the type of statistic, one of `"C"` (default, quantile-based groups)
  or `"H"` (equal-width groups on \\\[0, 1\]\\).

- obs:

  a numeric vector of observed binary outcomes (0 or 1) of the same
  length as `x`, without missing values; unused for the `glm` method.

- digits:

  number of significant digits to display.

- details:

  logical; if `TRUE`, prints observed and expected counts for both
  outcome classes by group. Default is `FALSE`.

## Value

An object of class `c("HosmerLemeshowTest", "htest")`, which is a list
with components:

- statistic:

  the chi-squared test statistic.

- parameter:

  the degrees of freedom (number of groups used minus 2).

- p.value:

  the p-value of the test.

- method:

  a character string describing the test.

- type:

  the type of statistic computed (`"C"` or `"H"`).

- nGroups:

  the number of groups actually used (may be less than requested if
  quantile breaks coincide or groups are empty).

- observed:

  matrix of observed counts for both outcome classes per group.

- expected:

  matrix of expected counts for both outcome classes per group.

- data.name:

  a character string giving the names of the data.

The `print` method accepts a `details` argument; if `TRUE`, observed and
expected counts for both outcome classes are printed by group.

## Details

The C statistic groups observations by quantiles of the fitted
probabilities (equal-frequency bins), while the H statistic uses
equal-width intervals on \\\[0, 1\]\\ (e.g. 0–0.1, 0.1–0.2, ... for
`nGroups = 10`), as defined in Lemeshow and Hosmer (1982). Groups that
contain no observations are dropped with a warning, and the degrees of
freedom are based on the number of groups actually used.

Both statistics are asymptotically compared with a chi-squared
distribution with \\g - 2\\ degrees of freedom (\\g\\ the number of
groups used) under the null hypothesis of adequate fit. This is the
convention for a model evaluated on its development data; for an
external validation sample, \\g\\ degrees of freedom would be
appropriate. The approximation may be unreliable in small samples; a
warning is issued when expected cell counts fall below 5.

## References

Lemeshow, S. and Hosmer, D. W. (1982) A review of goodness of fit
statistics for use in the development of logistic regression models.
*American Journal of Epidemiology*, **115**(1), 92–106.

Hosmer, D. W., Lemeshow, S. and Sturdivant, R. X. (2013) *Applied
Logistic Regression*, 3rd ed., New York: Wiley.

## See also

[`glm`](https://rdrr.io/r/stats/glm.html)

Other test.regression:
[`bpTest()`](https://andrisignorell.github.io/lumen/reference/bpTest.md),
[`breuschGodfreyTest()`](https://andrisignorell.github.io/lumen/reference/breuschGodfreyTest.md),
[`durbinWatsonTest()`](https://andrisignorell.github.io/lumen/reference/durbinWatsonTest.md),
[`leCessieTest()`](https://andrisignorell.github.io/lumen/reference/leCessieTest.md)

## Examples

``` r
set.seed(111)
x1  <- factor(sample(1:3, 50, replace = TRUE))
x2  <- rnorm(50)
obs <- sample(c(0, 1), 50, replace = TRUE)

model <- glm(obs ~ x1 + x2, family = binomial)

# glm method: probabilities and outcomes are extracted from the model
hosmerLemeshowTest(model)
#> Warning: some expected counts are less than 5; the chi-squared approximation may be unreliable
#> 
#>   Hosmer-Lemeshow C statistic 
#> 
#> data:  obs ~ x1 + x2 
#> X-squared = 4.85 , df = 8 , p-value = 0.7735 
#> Number of groups: 10 
#> 

# equivalent call with explicit vectors
res <- hosmerLemeshowTest(fitted(model), obs, type = "C")
#> Warning: some expected counts are less than 5; the chi-squared approximation may be unreliable
res
#> 
#>   Hosmer-Lemeshow C statistic 
#> 
#> data:  fitted(model) and obs 
#> X-squared = 4.85 , df = 8 , p-value = 0.7735 
#> Number of groups: 10 
#> 

print(res, details = TRUE)
#> 
#>   Hosmer-Lemeshow C statistic 
#> 
#> data:  fitted(model) and obs 
#> X-squared = 4.85 , df = 8 , p-value = 0.7735 
#> Number of groups: 10 
#> 
#> Observed vs Expected counts by group:
#>               Obs 0s Exp 0s Obs 1s Exp 1s
#> [0.262,0.325]      3 3.5250      2 1.4750
#> (0.325,0.374]      4 3.2600      1 1.7400
#> (0.374,0.514]      2 2.7142      3 2.2858
#> (0.514,0.526]      3 2.4012      2 2.5988
#> (0.526,0.539]      4 2.3274      1 2.6726
#> (0.539,0.551]      2 2.2662      3 2.7338
#>  (0.551,0.56]      2 2.2145      3 2.7855
#>  (0.56,0.567]      2 2.1814      3 2.8186
#> (0.567,0.584]      1 2.1236      4 2.8764
#> (0.584,0.616]      2 1.9866      3 3.0134
#> 
```
