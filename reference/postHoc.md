# Post Hoc Tests for Multiple Comparisons Following Anova

Provides a unified interface for several parametric post hoc tests
following a significant ANOVA. The function computes pairwise
comparisons of group means with different methods for controlling the
family-wise error rate.

## Usage

``` r
postHocTest(x, ...)

# S3 method for class 'aov'
postHocTest(
  x,
  which = NULL,
  method = c("hsd", "bonferroni", "lsd", "scheffe", "newmankeuls", "duncan"),
  conf.level = 0.95,
  ordered = FALSE,
  ...
)

# S3 method for class 'matrix'
postHocTest(
  x,
  method = c("none", "fdr", "BH", "BY", "bonferroni", "holm", "hochberg", "hommel"),
  conf.level = 0.95,
  ...
)

# S3 method for class 'table'
postHocTest(
  x,
  method = c("none", "fdr", "BH", "BY", "bonferroni", "holm", "hochberg", "hommel"),
  conf.level = 0.95,
  ...
)

# S3 method for class 'PostHocTest'
print(x, digits = getOption("digits", 3), ...)
```

## Arguments

- x:

  an object of class `aov`.

- ...:

  further arguments, not used so far.

- which:

  a character vector listing terms in the fitted model for which the
  intervals should be calculated. Defaults to all the terms.

- method:

  one of `"hsd"`, `"bonferroni"`, `"lsd"`, `"scheffe"`, `"newmankeuls"`,
  `"duncan"`, defining the method for the pairwise comparisons (may be
  abbreviated).  
  For the post hoc test of tables the methods of
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) can be supplied.
  See the detail there.

- conf.level:

  a numeric value between zero and one giving the family-wise confidence
  level to use. If this is set to NA, just a matrix with the p-values
  will be returned.

- ordered:

  a logical value indicating if the levels of the factor should be
  ordered according to increasing average in the sample before taking
  differences. If ordered is `TRUE` then the calculated differences in
  the means will all be positive. The significant differences will be
  those for which the lower end point is positive.  
  This argument will be ignored if method is not either `hsd` or
  `newmankeuls`.

- digits:

  controls the number of fixed digits to print.

## Value

An object of class `"PostHocTest"`: a list of data frames containing the
mean difference, lower and upper confidence interval bounds, and p-value
when `conf.level` is not `NA`; otherwise, a list of p-value matrices.

## Details

Post hoc tests differ in how strongly they control type I error.
Conservative methods reduce false positives but have lower power, while
more liberal methods increase power at the cost of a higher false
positive rate.

**LSD (Fisher)**: ` `No adjustment for multiple testing; highest power
but inflated type I error. Mainly suitable for a small number of groups.

**Bonferroni** (Dunn's (Bonferroni) t-test): ` `Adjusts p-values by the
number of comparisons; simple and robust, but often overly conservative.

**Tukey HSD**: ` `Controls the family-wise error rate for all pairwise
comparisons; widely used and generally recommended for balanced designs.

**Newman-Keuls**: ` `Stepwise procedure with more power than Tukey, but
weaker error control; may inflate type I error.

**Duncan**: ` `Similar to Newman-Keuls but more liberal; provides higher
power at the cost of increased false positives.

**Scheffé**: ` `Very conservative; suitable for both pairwise and
complex (contrast-based) comparisons.

**Guidance**  
Tukey HSD is typically a good default for pairwise comparisons.
Bonferroni is useful when strict error control is required. Scheffé is
appropriate for more complex contrasts.

**Tables**  
For tables pairwise chi-square tests can be performed, either without
correction or with correction for multiple testing following the logic
in [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

## See also

[`TukeyHSD`](https://rdrr.io/r/stats/TukeyHSD.html),
[`aov`](https://rdrr.io/r/stats/aov.html),
[`pairwise.t.test`](https://rdrr.io/r/stats/pairwise.t.test.html),
[`scheffeTest`](scheffeTest.md)

Other test.posthoc: [`conoverTest()`](conoverTest.md),
[`dscfTest()`](dscfTest.md), [`dunnTest()`](dunnTest.md),
[`dunnettTest()`](dunnettTest.md),
[`gamesHowellTest()`](gamesHowellTest.md),
[`nemenyiTest()`](nemenyiTest.md),
[`plot.PostHocTest()`](plot.PostHocTest.md),
[`scheffeTest()`](scheffeTest.md), [`signifDiff()`](signifDiff.md),
[`steelTest()`](steelTest.md)

## Examples

``` r

postHocTest(aov(breaks ~ tension, data = warpbreaks), method = "lsd")
#> 
#>   Posthoc multiple comparisons of means : Fisher LSD 
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = breaks ~ tension, data = warpbreaks)
#> 
#> $tension
#>           diff       lci       uci    pval signif
#> M-L -10.000000 -17.95042 -2.049581   0.015    *  
#> H-L -14.722222 -22.67264 -6.771803 < 0.001    ***
#> H-M  -4.722222 -12.67264  3.228197   0.239       
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
postHocTest(aov(breaks ~ tension, data = warpbreaks), method = "hsd")
#> 
#>   Posthoc multiple comparisons of means : Tukey HSD 
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = breaks ~ tension, data = warpbreaks)
#> 
#> $tension
#>           diff       lci        uci  pval signif
#> M-L -10.000000 -19.55982 -0.4401756 0.038    *  
#> H-L -14.722222 -24.28205 -5.1623978 0.001    ** 
#> H-M  -4.722222 -14.28205  4.8376022 0.463       
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
postHocTest(aov(breaks ~ tension, data = warpbreaks), method = "scheffe")
#> 
#>   Posthoc multiple comparisons of means : Scheffé 
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = breaks ~ tension, data = warpbreaks)
#> 
#> $tension
#>           diff       lci         uci  pval signif
#> M-L -10.000000 -19.98534 -0.01465926 0.050    *  
#> H-L -14.722222 -24.70756 -4.73688148 0.002    ** 
#> H-M  -4.722222 -14.70756  5.26311852 0.496       
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

r.aov <- aov(breaks ~ tension, data = warpbreaks)

# compare p-values:
round(cbind(
    lsd= postHocTest(r.aov, method="lsd")$tension[,"pval"]
  , bonf=postHocTest(r.aov, method="bonf")$tension[,"pval"]
), 4)
#>        lsd   bonf
#> M-L 0.0147 0.0442
#> H-L 0.0005 0.0015
#> H-M 0.2386 0.7158

# only p-values by setting conf.level to NA
postHocTest(aov(breaks ~ tension, data = warpbreaks), method = "hsd",
            conf.level=NA)
#> 
#>   Posthoc multiple comparisons of means : Tukey HSD 
#> 
#> Fit: aov(formula = breaks ~ tension, data = warpbreaks)
#> 
#> $tension
#>   L     M    
#> M 0.038     -
#> H 0.001 0.463
#> 
#> 
```
