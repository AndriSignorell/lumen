# Confidence Intervals for Median Absolute Deviations

Confidence intervals for the median absolute deviation (MAD) of a single
sample (`madCI`), for the difference of two MADs (`madDiffCI`), and for
the squared ratio of two MADs (`madRatioCI`). Two methods are available
throughout: an asymptotic interval based on the generalized lambda
distribution (GLD), and a parallel bootstrap interval.

## Usage

``` r
madCI(
  x,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("classic", "boot"),
  gldMethod = "TM",
  na.rm = FALSE,
  ...
)

madDiffCI(
  x,
  y,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("classic", "boot"),
  gldMethod = "TM",
  na.rm = FALSE,
  ...
)

madRatioCI(
  x,
  y,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("classic", "boot"),
  gldMethod = "TM",
  na.rm = FALSE,
  ...
)
```

## Arguments

- x:

  A non-empty numeric vector (first or only sample).

- conf.level:

  Confidence level of the interval. A single numeric value in \\(0,
  1)\\. Default `0.95`.

- sides:

  A character string specifying the side of the interval: `"two.sided"`
  (default), `"left"`, or `"right"`. Partial matching is supported.
  `"left"` sets `uci = Inf`; `"right"` sets `lci = -Inf`.

- method:

  A character string selecting the CI method: `"classic"` (asymptotic
  GLD-based, default) or `"boot"` (parallel bootstrap).

- gldMethod:

  A character string passed to `.asv.mad()` selecting the GLD estimation
  method. One of `"ML"`, `"MPS"`, `"TM"` (default), `"SM"`, `"TL"`,
  `"Lmom"`, `"DLA"`, or `"Mom"`. See
  [`fit.fkml()`](https://rdrr.io/pkg/gld/man/fit.fkml.html). Used only
  when `method = "classic"`.

- na.rm:

  Logical. Should missing values be removed before computation? Default
  `FALSE`.

- ...:

  Further arguments passed to the bootstrap engine when
  `method = "boot"`: `R`, `type`, `parallel`, `ncpus`. See Details.

- y:

  A non-empty numeric vector (second sample). Required for `madDiffCI`
  and `madRatioCI`.

## Value

A named numeric vector with three elements:

- `est`: point estimate:  
  ` ` \\\mathrm{mad}(x)\\ for `madCI`  
  ` ` \\\mathrm{mad}(x) - \mathrm{mad}(y)\\ for `madDiffCI`  
  ` ` \\(\mathrm{mad}(x)/\mathrm{mad}(y))^2\\ for `madRatioCI`

- `lci`: lower confidence bound.

- `uci`: upper confidence bound.

## Details

**Classic method** (`method = "classic"`)

All three functions follow Arachchige & Prendergast (2019) and base the
interval on the asymptotic variance of the MAD, approximated by fitting
a GLD to the data via `.asv.mad()`. The GLD estimation method is
selected with `gldMethod`.

For `madDiffCI` the asymptotic variances of the two samples are combined
as \\\widehat{\mathrm{ASV}}(x)/n_x + \widehat{\mathrm{ASV}}(y)/n_y\\.

For `madRatioCI` the interval is constructed on the log scale via the
delta method and back-transformed to guarantee positivity: \$\$
\exp\\\Bigl(\log\hat\theta \\\pm\\ z\_{\alpha/2}
\\\sqrt{\widehat{\mathrm{Var}}(\hat\theta)}\\/\\\hat\theta\Bigr), \quad
\hat\theta = \bigl(\mathrm{MAD}(x)/\mathrm{MAD}(y)\bigr)^2. \$\$

The classic method is fast and accurate for large samples but may
undercover for small or heavy-tailed distributions.

**Bootstrap method** (`method = "boot"`)

Data are resampled \\R\\ times using a parallel Rcpp worker and a
percentile or BCa interval is returned. For the two-sample functions the
two samples are resampled independently. Bootstrap arguments are passed
through `...` and extracted via `.extractBootArgs()`:

- `R`:

  Number of bootstrap replicates (default `999`).

- `type`:

  CI type: `"perc"` or `"bca"` (default).

- `parallel`:

  Parallelisation: `"no"`, `"multicore"`, or `"snow"` (default `"no"`).

- `ncpus`:

  Number of CPUs for parallel bootstrap (default
  `getOption("boot.ncpus", 1L)`).

## Note

Based on code by Arachchige Chandima N. P. G. and Prendergast Luke A.,
adapted to conform to package standards.

## References

Arachchige, C. N. P. G., & Prendergast, L. A. (2019). Confidence
intervals for median absolute deviations. *arXiv:1910.00229*
`[math.ST]`.

## See also

[`mad`](https://rdrr.io/r/stats/mad.html), `DescToolsX::madX`

## Examples

``` r
set.seed(1)
x <- rlnorm(100)
y <- rlnorm(200, meanlog = 1.2)

# single sample
madCI(x)
#>       est       lci       uci 
#> 0.9264988 0.7075985 1.1453991 
madCI(x, sides = "left")
#>       est       lci       uci 
#> 0.9264988 0.7427919       Inf 
madCI(x, method = "boot", R = 499, type = "bca")
#>       est       lci       uci 
#> 0.9264988 0.6922272 1.2100449 

# two-sample difference
madDiffCI(x, y)
#>       est       lci       uci 
#> -1.819975 -2.372101 -1.267849 
madDiffCI(x, y, method = "boot", R = 499, type = "perc")
#>       est       lci       uci 
#> -1.819975 -2.328010 -1.104668 

# two-sample squared ratio
madRatioCI(x, y)
#>        est        lci        uci 
#> 0.11379909 0.06247871 0.20727434 
madRatioCI(x, y, method = "boot", R = 499, type = "bca")
#>        est        lci        uci 
#> 0.11379909 0.05335922 0.25521800 
```
