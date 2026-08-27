# Mood's Median Test for Comparing Medians Across Independent Groups

Tests whether two or more independent samples come from populations with
the same median. All observations are dichotomised at the pooled median
and the resulting 2 x k table is tested for independence.

## Usage

``` r
moodMedianTest(x, ...)

# Default S3 method
moodMedianTest(
  x,
  g,
  ties = c("below", "above", "drop"),
  method = c("chisq", "exact"),
  correct = TRUE,
  ...
)

# S3 method for class 'formula'
moodMedianTest(formula, data, subset, na.action = na.pass, ...)
```

## Arguments

- x:

  a numeric vector of observations, or a formula of the form `lhs ~ rhs`
  with a numeric `lhs` and a grouping `rhs`

- ...:

  further arguments, passed to the default method

- g:

  a vector or factor giving the group for each element of `x`

- ties:

  how to treat observations exactly equal to the pooled median, one of
  `"below"` (default), `"above"` or `"drop"`

- method:

  the test applied to the 2 x k table, either `"chisq"` (default) or
  `"exact"`

- correct:

  logical, whether to apply the continuity correction; ignored unless
  the table is 2 x 2

- formula:

  a formula of the form `lhs ~ rhs`

- data:

  an optional data frame containing the model variables

- subset:

  an optional vector specifying a subset of observations

- na.action:

  a function indicating what should happen when the data contain `NA`s

## Value

An object of class `"htest"` with components

- statistic:

  Pearson's chi-squared statistic, `NULL` for `method = "exact"`

- parameter:

  the degrees of freedom, `NULL` for `method = "exact"`

- p.value:

  the p-value

- estimate:

  the group medians

- observed:

  the 2 x k table of counts above and below the pooled median

- expected:

  the expected counts under independence

- grand.median:

  the pooled median at which the data were dichotomised

- method:

  a character string describing the test

- data.name:

  a character string giving the names of the data

## Details

Not to be confused with
[`mood.test`](https://rdrr.io/r/stats/mood.test.html), which is Mood's
two-sample test for a difference in *scale*. The median test described
here has no base R implementation.

The procedure is the k-sample counterpart of [`signTest`](signTest.md)
and shares its robustness and its modest power: only the position of
each observation relative to the pooled median enters the statistic, so
all information about distance is discarded. Its asymptotic relative
efficiency against the F test under normality is \\2/\pi\\.

That low efficiency is the price of asking a narrow question, and the
alternatives ask different ones rather than the same one better:
[`brunnerMunzelTest`](brunnerMunzelTest.md) tests the relative effect
\\P(X \< Y) + \frac{1}{2}P(X = Y) = \frac{1}{2}\\ and
[`vanWaerdenTest`](vanWaerdenTest.md) tests equality of the
distributions against normal-score location alternatives. Neither is a
test of equal medians, so they are not drop-in replacements. Use the
median test when the median is genuinely the quantity of interest, or
when only the side of a threshold is trustworthy, as with coarsely
recorded or thresholded data.

**Observations equal to the pooled median.** With an even total sample
size the pooled median usually falls between two observations and the
question does not arise. Otherwise `ties` decides: `"below"` (default)
counts them with the lower group, following Conover, `"above"` with the
upper group, and `"drop"` removes them, which retains a symmetric
above-versus-below classification at the cost of a smaller effective
sample size. Dropping them does not systematically enlarge the p-value;
it can move it either way. `"below"` and `"above"` correspond to
`mid.score` values of `"0"` and `"1"` in
[`coin::median_test`](https://rdrr.io/pkg/coin/man/LocationTests.html);
`"drop"` has no counterpart there, since `mid.score = "0.5"` scores the
median observations rather than removing them.

`method = "exact"` replaces the chi-squared approximation by Fisher's
exact test on the same table and is advisable when expected counts are
small; no test statistic is reported in that case. The continuity
correction applies only to a 2 x 2 table, as in
[`chisq.test`](https://rdrr.io/r/stats/chisq.test.html).

Missing values are removed casewise.

## References

Mood, A. M. (1950) *Introduction to the Theory of Statistics*.
McGraw-Hill, New York, pp. 394-399.

Conover, W. J. (1999) *Practical Nonparametric Statistics*, 3rd edition.
Wiley, New York, pp. 218-223.

## See also

[`mood.test()`](https://rdrr.io/r/stats/mood.test.html)

Other test.location: [`brunnerMunzelTest()`](brunnerMunzelTest.md),
[`hotellingsT2Test()`](hotellingsT2Test.md),
[`signTest()`](signTest.md), [`tTestA()`](tTestA.md),
[`vanWaerdenTest()`](vanWaerdenTest.md), [`yuenTTest()`](yuenTTest.md),
[`zTest()`](zTest.md)

## Examples

``` r
moodMedianTest(breaks ~ tension, data = warpbreaks)
#> 
#>  Mood's median test (3 groups, chi-squared approximation)
#> 
#> data:  breaks ~ tension
#> X-squared = 7.2993, df = 2, p-value = 0.026
#> sample estimates:
#>    L    M    H 
#> 29.5 27.0 20.5 
#> 

## the pooled median at which the data were split, and the resulting table
r <- moodMedianTest(breaks ~ tension, data = warpbreaks)
r$grand.median
#> [1] 26
## [1] 26
r$observed
#>        
#>          L  M  H
#>   above 12  9  4
#>   below  6  9 14

## small tables: use the exact test on the same 2 x k table
moodMedianTest(breaks ~ tension, data = warpbreaks, method = "exact")$p.value
#> [1] 0.03153902
## [1] 0.03153902
```
