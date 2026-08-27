# Pearson Chi-Square Test for Assessing Normality From Grouped Frequencies

Performs the Pearson chi-square test for the composite hypothesis of
normality, see e.g. Thode (2002, Sec. 5.2).

## Usage

``` r
pearsonTest(x, nClasses = ceiling(2 * (n^(2/5))), adjust = TRUE)
```

## Arguments

- x:

  a numeric vector of data values. Missing values are allowed.

- nClasses:

  the number of classes. The default is due to Moore (1986).

- adjust:

  logical; if `TRUE` (default), the p-value is computed from a
  chi-square distribution with `nClasses`-3 degrees of freedom,
  otherwise from a chi-square distribution with `nClasses`-1 degrees of
  freedom.

## Value

A list with class “htest” containing the following components:

- statistic:

  the value of the Pearson chi-square statistic.

- p.value :

  the p-value for the test.

- method:

  the character string “Pearson chi-square normality test”.

- data.name:

  a character string giving the name(s) of the data.

- nClasses:

  the number of classes used for the test.

- df:

  the degress of freedom of the chi-square distribution used to compute
  the p-value.

## Details

The Pearson test statistic is \\P=\sum (C\_{i} - E\_{i})^{2}/E\_{i}\\,
where \\C\_{i}\\ is the number of counted and \\E\_{i}\\ is the number
of expected observations (under the hypothesis) in class \\i\\. The
classes are build is such a way that they are equiprobable under the
hypothesis of normality. The p-value is computed from a chi-square
distribution with `nClasses`-3 degrees of freedom if `adjust` is `TRUE`
and from a chi-square distribution with `nClasses`-1 degrees of freedom
otherwise. In both cases this is not (!) the correct p-value, lying
somewhere between the two, see also Moore (1986).

## Note

The Pearson chi-square test is usually not recommended for testing the
composite hypothesis of normality due to its inferior power properties
compared to other tests. It is common practice to compute the p-value
from the chi-square distribution with `nClasses` - 3 degrees of freedom,
in order to adjust for the additional estimation of two parameters. (For
the simple hypothesis of normality (mean and variance known) the test
statistic is asymptotically chi-square distributed with `nClasses` - 1
degrees of freedom.) This is, however, not correct as long as the
parameters are estimated by `mean(x)` and `var(x)` (or `sd(x)`), as it
is usually done, see Moore (1986) for details. Since the true p-value is
somewhere between the two, it is suggested to run `pearsonTest` twice,
with `adjust = TRUE` (default) and with `adjust = FALSE`. It is also
suggested to slightly change the default number of classes, in order to
see the effect on the p-value. Eventually, it is suggested not to rely
upon the result of the test.

The function call `pearsonTest(x)` essentially produces the same result
as the S-PLUS function call
`chisq.gof((x-mean(x))/sqrt(var(x)), n.param.est=2)`.

Based on code by Juergen Gross, adapted to conform to package standards.

## References

Moore, D.S. (1986): Tests of the chi-squared type. In: D'Agostino, R.B.
and Stephens, M.A., eds.: Goodness-of-Fit Techniques. Marcel Dekker, New
York.

Thode Jr., H.C. (2002): Testing for Normality. Marcel Dekker, New York.

## See also

[`shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html) for
performing the Shapiro-Wilk test for normality

Other test.normality: [`andersonDarlingTest()`](andersonDarlingTest.md),
[`cramerVonMisesTest()`](cramerVonMisesTest.md),
[`jarqueBeraTest()`](jarqueBeraTest.md),
[`lillieTest()`](lillieTest.md),
[`shapiroFranciaTest()`](shapiroFranciaTest.md)

## Examples

``` r

pearsonTest(rnorm(100, mean = 5, sd = 3))
#> 
#>  Pearson chi-square normality test
#> 
#> data:  rnorm(100, mean = 5, sd = 3)
#> P = 6.6, p-value = 0.7626
#> 
pearsonTest(runif(100, min = 2, max = 4))
#> 
#>  Pearson chi-square normality test
#> 
#> data:  runif(100, min = 2, max = 4)
#> P = 12.06, p-value = 0.2811
#> 
```
