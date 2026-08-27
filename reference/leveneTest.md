# Levene's Test for Comparing Variances Across Groups

A test for homogeneity of variances across two or more groups, more
robust to departures from normality than Bartlett's test.

## Usage

``` r
leveneTest(x, ...)

# S3 method for class 'formula'
leveneTest(formula, data, subset, na.action = na.pass, center = median, ...)

# Default S3 method
leveneTest(x, g, center = median, .centerName = NULL, ...)
```

## Arguments

- x:

  a numeric vector of data values (default method), or a `formula`
  (formula method).

- ...:

  arguments to be passed down, e.g. `data` for the formula method; can
  also be used to pass arguments to the function given by `center` (e.g.
  `trim = 0.1` for a trimmed mean when `center = mean`).

- formula:

  a formula of the form `lhs ~ rhs` where `lhs` gives the data values
  and `rhs` the corresponding groups.

- data:

  an optional matrix or data frame (or similar: see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in the formula `formula`. By default the variables are
  taken from `environment(formula)`.

- subset:

  an optional vector specifying a subset of observations to be used.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

- center:

  the name of a function to compute the center of each group; `mean`
  gives the original Levene's test, the default, `median`, provides the
  more robust Brown-Forsythe test.

- g:

  factor defining the groups; ignored if `x` is a list.

- .centerName:

  internal, not intended to be set by the user. Used to pass the
  deparsed name of the `center` function through the method dispatch
  chain (from `leveneTest.formula` to `leveneTest.default`), since
  `substitute(center)` would otherwise only resolve to the literal
  symbol `"center"` rather than the original expression (e.g. `mean` or
  `median`) supplied by the caller.

## Value

An object of class `"htest"` representing the result of the hypothesis
test.

## Details

Let \\X\_{ij}\\ be the jth observation of X for the ith group. Let
\\Z\_{ij} = \|X\_{ij} - X_i\|\\, where \\X_i\\ is the center (by default
the median) of X in the ith group. Levene's test statistic is \$\$W_0 =
\frac{\sum_i n_i (\bar{Z}\_i - \bar{Z})^2 / (g - 1)}{\sum_i \sum_j
(Z\_{ij} - \bar{Z}\_i)^2 / \sum_i (n_i - 1)}\$\$ where \\n_i\\ is the
number of observations in group i and g is the number of groups. Using
the median instead of the mean for \\X_i\\ yields the more robust
variant known as the Brown-Forsythe test, which is the default here.

## Note

Based on `car::leveneTest()` by John Fox, with contributions by Derek
Ogle and Brian Ripley, adapted to conform to package standards.

## References

Fox, J. and Weisberg, S. (2019) *An R Companion to Applied Regression*,
3rd ed., Thousand Oaks, CA: Sage.

Levene, H. (1960) Robust tests for equality of variances. In: Olkin, I.
et al., eds.: *Contributions to Probability and Statistics: Essays in
Honor of Harold Hotelling*, Stanford University Press, pp. 278-292.

## See also

[`fligner.test`](https://rdrr.io/r/stats/fligner.test.html) for a
rank-based (nonparametric) k-sample test for homogeneity of variances,
[`bartlett.test`](https://rdrr.io/r/stats/bartlett.test.html) for a
parametric alternative

Other test.variance: [`mosesTest()`](mosesTest.md),
[`siegelTukeyTest()`](siegelTukeyTest.md), [`varCI()`](varCI.md),
[`varTest()`](varTest.md)

## Examples

``` r
## example from ansari.test:
## Hollander & Wolfe (1973, p. 86f):
## Serum iron determination using Hyland control sera
serum <- data.frame(
  grp = rep(c("ramsay", "jung.parekh"), each = 20),
  x   = c(111, 107, 100, 99, 102, 106, 109, 108, 104, 99,
          101, 96, 97, 102, 107, 113, 116, 113, 110, 98,
          107, 108, 106, 98, 105, 103, 110, 105, 104,
          100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99)
)

leveneTest(x ~ grp, data = serum)
#> 
#>  Levene's Test for Homogeneity of Variance (center = median)
#> 
#> data:  x ~ grp
#> F = 1.7865, num df = 1, denom df = 38, p-value = 0.1893
#> 
leveneTest(x ~ grp, data = serum, center = mean)
#> 
#>  Levene's Test for Homogeneity of Variance (center = mean)
#> 
#> data:  x ~ grp
#> F = 1.7879, num df = 1, denom df = 38, p-value = 0.1891
#> 
leveneTest(x ~ grp, data = serum, center = mean, trim = 0.1)
#> 
#>  Levene's Test for Homogeneity of Variance (center = mean(trim=0.1))
#> 
#> data:  x ~ grp
#> F = 1.7854, num df = 1, denom df = 38, p-value = 0.1894
#> 

leveneTest(c(rnorm(10), rnorm(10, 0, 2)),
           factor(rep(c("A", "B"), each = 10)))
#> 
#>  Levene's Test for Homogeneity of Variance (center = median)
#> 
#> data:  c(rnorm(10), rnorm(10, 0, 2)) and factor(rep(c("A", "B"), each = 10))
#> F = 1.0614, num df = 1, denom df = 18, p-value = 0.3165
#> 

leveneTest(Ozone ~ factor(Month), data = airquality)
#> 
#>  Levene's Test for Homogeneity of Variance (center = median)
#> 
#> data:  Ozone ~ factor(Month)
#> F = 3.9558, num df = 4, denom df = 111, p-value = 0.004863
#> 

leveneTest(count ~ spray, data = InsectSprays)
#> 
#>  Levene's Test for Homogeneity of Variance (center = median)
#> 
#> data:  count ~ spray
#> F = 3.8214, num df = 5, denom df = 66, p-value = 0.004223
#> 
# Compare this to fligner.test() and bartlett.test()
```
