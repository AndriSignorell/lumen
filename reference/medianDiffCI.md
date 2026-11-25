# Confidence Interval for the Difference of Two Medians

Calculates the confidence interval for the difference of two medians,
`median(x) - median(y)`, either with the distribution-free method of
Bonett and Price or with the bootstrap.

## Usage

``` r
medianDiffCI(
  x,
  y,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("bonett-price", "boot"),
  paired = FALSE,
  na.rm = FALSE,
  ...
)
```

## Arguments

- x:

  a (non-empty) numeric vector of data values.

- y:

  a (non-empty) numeric vector of data values.

- conf.level:

  confidence level of the interval.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `"left"` would be analogue to a
  hypothesis of `"greater"` in a `t.test`.

- method:

  a character string, either `"bonett-price"` (default) or `"boot"`. The
  bootstrap type can be provided by the dots, see
  [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html).

- paired:

  a logical indicating whether `x` and `y` are paired. Defaults to
  `FALSE`. Requires `method = "boot"`.

- na.rm:

  logical. Should missing values be removed? Defaults to `FALSE`. For
  paired data, incomplete pairs are removed.

- ...:

  further arguments passed to the bootstrap, e.g. `R` and `type`.

## Value

A named numeric vector with elements:

- `mediandiff`:

  point estimate, the difference `median(x) - median(y)`

- `lci`:

  lower confidence interval bound

- `uci`:

  upper confidence interval bound

## Details

`method = "bonett-price"` (default) estimates the standard error of each
sample median from a pair of order statistics, \$\$\widehat{se}\_j =
\frac{Y\_{(n_j - c_j + 1)} - Y\_{(c_j)}}{2 z_j}, \quad c_j =
\mathrm{round}\left(\frac{n_j}{2} - \sqrt{n_j}\right),\$\$ where \\z_j\\
is the normal quantile matching the exact binomial coverage of
\\\[Y\_{(c_j)}, Y\_{(n_j - c_j + 1)}\]\\, and combines both into a
normal interval. No distributional assumption is made, and the variances
of the two groups may differ. It is available for independent samples
only.

`method = "boot"` resamples within each group (stratified), or whole
pairs if `paired = TRUE`. The bootstrap distribution of a median is
discrete and lumpy in small samples; `type = "perc"` is the robust
choice there. The bootstrap type `"stud"` is not supported, as the
statistic returns no per-replicate variance estimate.

Note that for paired data the difference of the medians is not the
median of the differences. For the latter, use `medianCI(x - y)`.
Neither is it the Hodges-Lehmann shift estimate reported by
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html), which
coincides with the median difference only under a pure location shift.

## References

Bonett, D. G. and Price, R. M. (2002). Statistical inference for a
linear function of medians: confidence intervals, hypothesis testing,
and sample size requirements. *Psychological Methods* **7**(3), 370–383.
[doi:10.1037//1082-989x.7.3.370](https://doi.org/10.1037//1082-989x.7.3.370)

Price, R. M. and Bonett, D. G. (2002). Distribution-free confidence
intervals for difference and ratio of medians. *Journal of Statistical
Computation and Simulation* **72**(2), 119–124.
[doi:10.1080/00949650212140](https://doi.org/10.1080/00949650212140)

## See also

[`meanDiffCI()`](meanDiffCI.md), [`medianCI()`](medianCI.md),
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html)

Other ci.location: [`meanCI()`](meanCI.md), [`meanCIn()`](meanCIn.md),
[`meanDiffCI()`](meanDiffCI.md), [`medianCI()`](medianCI.md),
[`quantileCI()`](quantileCI.md), [`sumCI()`](sumCI.md)

## Examples

``` r
x <- mtcars[mtcars$am == 0, "mpg"]
y <- mtcars[mtcars$am == 1, "mpg"]

medianDiffCI(x, y)
#>  mediandiff         lci         uci 
#>  -5.5000000 -10.4651329  -0.5348671 
medianDiffCI(x, y, conf.level = 0.99)
#> mediandiff        lci        uci 
#>  -5.500000 -12.025291   1.025291 
medianDiffCI(x, y, sides = "left")
#> mediandiff        lci        uci 
#>  -5.500000  -9.666871        Inf 

# bootstrap
set.seed(1)
medianDiffCI(x, y, method = "boot", type = "perc")
#> mediandiff        lci        uci 
#>       -5.5      -14.0       -2.2 
medianDiffCI(x, y, method = "boot", type = "bca")
#> mediandiff        lci        uci 
#>       -5.5      -12.6       -1.8 

# paired: difference of medians vs. median of differences
set.seed(2)
drug2 <- sleep$extra[sleep$group == 2]
drug1 <- sleep$extra[sleep$group == 1]
medianDiffCI(drug2, drug1, method = "boot", paired = TRUE, type = "perc")
#> mediandiff        lci        uci 
#>       1.40       0.75       3.70 
medianCI(drug2 - drug1)
#> median    lci    uci 
#>    1.3    0.8    2.4 
#> attr(,"conf.level")
#> [1] 0.9785156
```
