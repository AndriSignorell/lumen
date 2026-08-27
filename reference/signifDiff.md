# Significantly Different Levels per Group

Extracts, for every factor level, the set of levels it differs
significantly from, based on the pairwise comparisons of a post hoc
test.

## Usage

``` r
signifDiff(x, ...)

# S3 method for class 'PostHocTest'
signifDiff(
  x,
  alpha = NULL,
  direction = TRUE,
  labels = "numbers",
  minlength = 3L,
  sep = ", ",
  ...
)

# S3 method for class 'pairwise.htest'
signifDiff(
  x,
  alpha = 0.05,
  direction = FALSE,
  labels = "numbers",
  minlength = 3L,
  sep = ", ",
  ...
)

# S3 method for class 'signifDiff'
print(x, legend = TRUE, ...)
```

## Arguments

- x:

  an object of class `PostHocTest` as returned by
  [`postHocTest`](postHoc.md), or of class `pairwise.htest` as returned
  by [`pairwise.t.test`](https://rdrr.io/r/stats/pairwise.t.test.html)
  and friends

- ...:

  further arguments, not used so far

- alpha:

  the significance level; defaults to `1 - conf.level` of the object, or
  to 0.05 where no confidence level is stored

- direction:

  logical; if `TRUE`, a sign is appended to every label, giving the
  direction of the row level's mean relative to the listed level.
  Requires the differences to be present in the object, which is not the
  case for the p-value branch and for `pairwise.htest` objects.

- labels:

  either a keyword defining how the levels are abbreviated, one of
  `"numbers"` (the default), `"letters"`, `"LETTERS"`, `"abbreviate"` or
  `"names"`, or a character vector of labels with one entry per level

- minlength:

  the minimum length of the abbreviations, used for
  `labels = "abbreviate"` only

- sep:

  the separator between the labels

- legend:

  logical; if `TRUE`, the meaning of the signs is printed as a footer,
  separated by a rule

## Value

a list with one data frame per term, each holding the columns `label`
(the label of the level itself) and `diff` (the labels of the levels it
differs significantly from), with the levels as row names. The class is
`"signifDiff"`, the significance level is kept in the attribute `alpha`.

## Details

The result answers the question a reader of a post hoc table usually
has: not which of the pairs are significant, but which levels a given
level is distinguishable from, and in which direction.

The labels are numbers by default. Letters are available, but carry the
opposite meaning of a compact letter display, where levels sharing a
letter are those that could *not* be distinguished; here a label marks a
significant difference.

## See also

[`postHocTest`](postHoc.md),
[`pairwise.t.test`](https://rdrr.io/r/stats/pairwise.t.test.html)

Other test.posthoc: [`conoverTest()`](conoverTest.md),
[`dscfTest()`](dscfTest.md), [`dunnTest()`](dunnTest.md),
[`dunnettTest()`](dunnettTest.md),
[`gamesHowellTest()`](gamesHowellTest.md),
[`nemenyiTest()`](nemenyiTest.md),
[`plot.PostHocTest()`](plot.PostHocTest.md), [`postHoc`](postHoc.md),
[`scheffeTest()`](scheffeTest.md), [`steelTest()`](steelTest.md)

## Examples

``` r
r.aov <- aov(breaks ~ tension, data = warpbreaks)
res <- postHocTest(r.aov, method = "hsd")

signifDiff(res)
#> 
#>   Posthoc multiple comparisons of means : Tukey HSD 
#>     levels a level differs from, at alpha = 0.05
#> 
#> $tension
#>   label diff  
#> L 1     2+, 3+
#> M 2     1-    
#> H 3     1-    
#> 
#> ---
#> Sign codes:  '+' mean of the row level is higher than the listed one, '-' lower
#> 

# stricter level, without recomputing the test
signifDiff(res, alpha = 0.01)
#> 
#>   Posthoc multiple comparisons of means : Tukey HSD 
#>     levels a level differs from, at alpha = 0.01
#> 
#> $tension
#>   label diff
#> L 1     3+  
#> M 2         
#> H 3     1-  
#> 
#> ---
#> Sign codes:  '+' mean of the row level is higher than the listed one, '-' lower
#> 

# readable labels instead of numbers
signifDiff(res, labels = "abbreviate", minlength = 4)
#> 
#>   Posthoc multiple comparisons of means : Tukey HSD 
#>     levels a level differs from, at alpha = 0.05
#> 
#> $tension
#>   label diff  
#> L L     M+, H+
#> M M     L-    
#> H H     L-    
#> 
#> ---
#> Sign codes:  '+' mean of the row level is higher than the listed one, '-' lower
#> 

# the p-value branch carries no differences, hence no signs
signifDiff(postHocTest(r.aov, method = "hsd", conf.level = NA))
#> 
#>   Posthoc multiple comparisons of means : Tukey HSD 
#>     levels a level differs from, at alpha = 0.05
#> 
#> $tension
#>   label diff
#> L 1     2, 3
#> M 2     1   
#> H 3     1   
#> 

signifDiff(pairwise.t.test(warpbreaks$breaks, warpbreaks$tension))
#> 
#>   Pairwise comparisons : t tests with pooled SD 
#>     levels a level differs from, at alpha = 0.05
#> 
#> $warpbreaks$breaks and warpbreaks$tension
#>   label diff
#> L 1     2, 3
#> M 2     1   
#> H 3     1   
#> 
```
