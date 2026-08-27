# Scheffé Test for Testing Arbitrary Contrasts Following Anova

A parametric post hoc test for all possible pairwise and complex
comparisons following a significant ANOVA, controlling the familywise
error rate conservatively for any number of contrasts.

## Usage

``` r
scheffeTest(x, ...)

# Default S3 method
scheffeTest(
  x,
  g = NULL,
  which = NULL,
  contrasts = NULL,
  conf.level = 0.95,
  ...
)

# S3 method for class 'formula'
scheffeTest(formula, data, subset, na.action, ...)

# S3 method for class 'aov'
scheffeTest(x, which = NULL, contrasts = NULL, conf.level = 0.95, ...)
```

## Arguments

- x:

  either a fitted model object, usually an
  [`aov`](https://rdrr.io/r/stats/aov.html) fit, when g is left to
  `NULL` or a response variable to be evalutated by g (which mustn't be
  `NULL` then).

- ...:

  further arguments, currently not used.

- g:

  the grouping variable.

- which:

  character vector listing terms in the fitted model for which the
  intervals should be calculated. Defaults to all the terms.

- contrasts:

  a \\r \times c\\ matrix containing the contrasts to be computed, while
  `r` is the number of factor levels and `c` the number of contrasts.
  Each column must contain a full contrast ("sum") adding up to 0. Note
  that the argument `which` must be defined, when non default contrasts
  are used. Default value of `contrasts` is `NULL`. In this case all
  pairwise contrasts will be reported.

- conf.level:

  numeric value between zero and one giving the confidence level to use.
  If this is set to NA, just a matrix with the p-values will be
  returned.

- formula:

  a formula of the form `lhs ~ rhs`, where `lhs` gives the response
  values and `rhs` the corresponding groups or explanatory variables.

- data:

  an optional matrix or data frame (or similar; see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in the formula. By default the variables are taken from
  `environment(formula)`.

- subset:

  an optional vector specifying a subset of observations to be used in
  the analysis.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

## Value

A list of classes `c("PostHocTest")`, with one component for each term
requested in `which`. Each component is a matrix with columns `diff`
giving the difference in the observed means, `lci` giving the lower end
point of the interval, `uci` giving the upper end point and `pval`
giving the p-value after adjustment for the multiple comparisons.

There are print and plot methods for class `"PostHocTest"`. The plot
method does not accept `xlab`, `ylab` or `main` arguments and creates
its own values for each plot.

## Details

Scheffé's method applies to the set of estimates of all possible
contrasts among the factor level means, not just the pairwise
differences considered by Tukey's method.

## References

Robert O. Kuehl, Steel R. (2000) *Design of experiments*. Duxbury

Steel R.G.D., Torrie J.H., Dickey, D.A. (1997) *Principles and
Procedures of Statistics, A Biometrical Approach*. McGraw-Hill

## See also

[`pairwise.t.test`](https://rdrr.io/r/stats/pairwise.t.test.html),
[`TukeyHSD`](https://rdrr.io/r/stats/TukeyHSD.html)

Other test.posthoc:
[`conoverTest()`](https://andrisignorell.github.io/lumen/reference/conoverTest.md),
[`dscfTest()`](https://andrisignorell.github.io/lumen/reference/dscfTest.md),
[`dunnTest()`](https://andrisignorell.github.io/lumen/reference/dunnTest.md),
[`dunnettTest()`](https://andrisignorell.github.io/lumen/reference/dunnettTest.md),
[`gamesHowellTest()`](https://andrisignorell.github.io/lumen/reference/gamesHowellTest.md),
[`nemenyiTest()`](https://andrisignorell.github.io/lumen/reference/nemenyiTest.md),
[`plot.PostHocTest()`](https://andrisignorell.github.io/lumen/reference/plot.PostHocTest.md),
[`postHoc`](https://andrisignorell.github.io/lumen/reference/postHoc.md),
[`signifDiff()`](https://andrisignorell.github.io/lumen/reference/signifDiff.md),
[`steelTest()`](https://andrisignorell.github.io/lumen/reference/steelTest.md)

## Examples

``` r

fm1 <- aov(breaks ~ wool + tension, data = warpbreaks)
scheffeTest(x=fm1)
#> 
#>   Posthoc multiple comparisons of means: Scheffé Test 
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = breaks ~ wool + tension, data = warpbreaks)
#> 
#> $wool
#>          diff       lci       uci  pval signif
#> B-A -5.777778 -12.12841 0.5728505 0.074    .  
#> 
#> $tension
#>           diff       lci        uci  pval signif
#> M-L -10.000000 -19.76977 -0.2302286 0.044    *  
#> H-L -14.722222 -24.49199 -4.9524508 0.002    ** 
#> H-M  -4.722222 -14.49199  5.0475492 0.481       
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
scheffeTest(x=fm1, which="tension")
#> 
#>   Posthoc multiple comparisons of means: Scheffé Test 
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = breaks ~ wool + tension, data = warpbreaks)
#> 
#> $tension
#>           diff       lci        uci  pval signif
#> M-L -10.000000 -19.76977 -0.2302286 0.044    *  
#> H-L -14.722222 -24.49199 -4.9524508 0.002    ** 
#> H-M  -4.722222 -14.49199  5.0475492 0.481       
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

# some special contrasts
y <- c(7,33,26,27,21,6,14,19,6,11,11,18,14,18,19,14,9,12,6,
       24,7,10,1,10,42,25,8,28,30,22,17,32,28,6,1,15,9,15,
       2,37,13,18,23,1,3,4,6,2)
group <- factor(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,
       3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6))

r.aov <- aov(y ~ group)
scheffeTest(r.aov, contrasts=matrix( c(1,-0.5,-0.5,0,0,0,
                                       0,0,0,1,-0.5,-0.5), ncol=2) )
#> 
#>   Posthoc multiple comparisons of means: Scheffé Test 
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = y ~ group)
#> 
#> $group
#>          diff       lci      uci  pval signif
#> 1-2,3  7.2500 -6.417446 20.91745 0.637       
#> 4-5,6 14.0625  0.395054 27.72995 0.040    *  
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

# just p-values:
scheffeTest(r.aov, conf.level=NA)
#> 
#>   Posthoc multiple comparisons of means: Scheffé Test 
#> 
#> Fit: aov(formula = y ~ group)
#> 
#> $group
#>   1     2     3     4     5    
#> 2 0.927     -     -     -     -
#> 3 0.531 0.977     -     -     -
#> 4 0.848 0.273 0.054     -     -
#> 5 0.940 1.000 0.970 0.296     -
#> 6 0.400 0.934 1.000 0.031 0.920
#> 
#> 
```
