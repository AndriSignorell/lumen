# lumen

`lumen` provides hypothesis tests, confidence intervals, and selected
statistical distributions used across the DescToolsX ecosystem. It is
designed as a focused statistical companion package: methodologically
transparent, API-consistent, and suitable for use in applied statistical
workflows.

Every function ships with a documented reference — a textbook formula,
an original paper, or a well-established R implementation it was checked
against — rather than a bare implementation.

The package is currently under active development.

## Installation

You can install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("AndriSignorell/lumen")
```

The package requires R \>= 4.2.0.

## Scope

`lumen` collects statistical procedures that are commonly useful in
exploratory analysis, inference, and methodological comparison:

- hypothesis tests for goodness-of-fit, normality, stationarity,
  randomness, contingency tables, marginal homogeneity, and
  nonparametric group comparisons
- post-hoc procedures and multiple-comparison methods, both parametric
  (Dunnett, Scheffé, Tukey via
  [`postHoc()`](https://andrisignorell.github.io/lumen/reference/postHoc.md))
  and nonparametric (Conover, Dunn, Nemenyi, Steel, DSCF)
- confidence intervals for proportions, differences and ratios of
  proportions, means, medians, variances, sums, correlations, and
  quantiles, each with classical, robust, and bootstrap variants where
  applicable
- selected probability distributions, including extreme-value (GEV, GPD,
  Gumbel, Fréchet, reversed Weibull), Dirichlet, Gompertz, triangular,
  Benford, and order-statistic distributions, plus their mean/variance
  companions
- generic bootstrap confidence interval and power-analysis helpers

## Examples

### Binomial confidence intervals

``` r

library(lumen)

binomCI(x = 37, n = 43, method = "wilson")
binomCI(x = 42, n = 43, method = "clopper-pearson")
```

### Goodness-of-fit testing

``` r

x <- rnorm(50)
andersonDarlingTest(x, null = "pnorm", mean = mean(x), sd = sd(x), estimated = TRUE)
```

### Nonparametric tests

``` r

x <- c(1.1, 1.4, 1.6, 2.0, 2.2)
y <- c(1.0, 1.2, 1.3, 1.7, 1.9)

siegelTukeyTest(x, y)
```

### Contingency-table tests

``` r

tab <- matrix(c(8, 14, 1, 3), nrow = 2)
barnardTest(tab)
```

### Post hoc comparisons after ANOVA

``` r

fit <- aov(breaks ~ tension, data = warpbreaks)
postHoc(fit, method = "hsd")
```

### Bootstrap confidence intervals

``` r

set.seed(1984)
bootCI(mtcars$mpg, FUN = mean, na.rm = TRUE, bci.method = "basic")
```

## Design principles

`lumen` follows the broader DescToolsX design philosophy:

- predictable lowerCamelCase function names, with
  `.default`/`.formula`/`.glm` S3 methods offered wherever a formula or
  model-object interface makes sense
- explicit argument validation with informative error messages
- transparent method choices — where a test or interval admits more than
  one established computation (e.g. exact vs. asymptotic, classical
  vs. bootstrap), the alternatives are exposed as an explicit
  `method`/`type` argument rather than hidden behind a single default
- clean separation between user-facing interfaces and computational
  engines, with performance-critical algorithms (exact permutation and
  rank distributions, Wellek’s exact Page test tables, run-length
  distributions) implemented in C++ via Rcpp/RcppParallel/RcppArmadillo
- compatibility with familiar base R idioms: most hypothesis tests
  return `htest`-compatible objects usable with
  [`print()`](https://rdrr.io/r/base/print.html), `$statistic`,
  `$p.value`, and friends

## Function reference

The sections below group the most commonly used functions by purpose.
The full alphabetical reference, including internal helpers and every
distribution variant, is available via `?lumen` and the [pkgdown
site](https://andrisignorell.github.io/lumen/reference/index.html).

### Goodness-of-fit and normality

[`andersonDarlingTest()`](https://andrisignorell.github.io/lumen/reference/andersonDarlingTest.md)
·
[`cramerVonMisesTest()`](https://andrisignorell.github.io/lumen/reference/cramerVonMisesTest.md)
·
[`jarqueBeraTest()`](https://andrisignorell.github.io/lumen/reference/jarqueBeraTest.md)
·
[`lillieTest()`](https://andrisignorell.github.io/lumen/reference/lillieTest.md)
·
[`pearsonTest()`](https://andrisignorell.github.io/lumen/reference/pearsonTest.md)
·
[`shapiroFranciaTest()`](https://andrisignorell.github.io/lumen/reference/shapiroFranciaTest.md)

### Time series: stationarity, autocorrelation, randomness

[`adfTest()`](https://andrisignorell.github.io/lumen/reference/adfTest.md)
·
[`kpssTest()`](https://andrisignorell.github.io/lumen/reference/kpssTest.md)
·
[`durbinWatsonTest()`](https://andrisignorell.github.io/lumen/reference/durbinWatsonTest.md)
·
[`breuschGodfreyTest()`](https://andrisignorell.github.io/lumen/reference/breuschGodfreyTest.md)
·
[`bpTest()`](https://andrisignorell.github.io/lumen/reference/bpTest.md)
·
[`bartelsRankTest()`](https://andrisignorell.github.io/lumen/reference/BartelsRankTest.md)
·
[`runsTest()`](https://andrisignorell.github.io/lumen/reference/runsTest.md)
·
[`vonNeumannTest()`](https://andrisignorell.github.io/lumen/reference/vonNeumannTest.md)

### Categorical and contingency-table tests

[`barnardTest()`](https://andrisignorell.github.io/lumen/reference/barnardTest.md)
·
[`bhapkarTest()`](https://andrisignorell.github.io/lumen/reference/bhapkarTest.md)
·
[`breslowDayTest()`](https://andrisignorell.github.io/lumen/reference/breslowDayTest.md)
·
[`stuartMaxwellTest()`](https://andrisignorell.github.io/lumen/reference/stuartMaxwellTest.md)
·
[`woolfTest()`](https://andrisignorell.github.io/lumen/reference/woolfTest.md)
·
[`cochranArmitageTest()`](https://andrisignorell.github.io/lumen/reference/cochranArmitageTest.md)
·
[`cochranQTest()`](https://andrisignorell.github.io/lumen/reference/cochranQTest.md)
· [`gTest()`](https://andrisignorell.github.io/lumen/reference/gTest.md)
·
[`lehmacherTest()`](https://andrisignorell.github.io/lumen/reference/lehmacherTest.md)
·
[`mantelTrendTest()`](https://andrisignorell.github.io/lumen/reference/mantelTrendTest.md)
·
[`hosmerLemeshowTest()`](https://andrisignorell.github.io/lumen/reference/hosmerLemeshowTest.md)
·
[`leCessieTest()`](https://andrisignorell.github.io/lumen/reference/leCessieTest.md)

### Location and scale tests

[`hotellingsT2Test()`](https://andrisignorell.github.io/lumen/reference/hotellingsT2Test.md)
·
[`signTest()`](https://andrisignorell.github.io/lumen/reference/signTest.md)
·
[`tTestA()`](https://andrisignorell.github.io/lumen/reference/tTestA.md)
·
[`yuenTTest()`](https://andrisignorell.github.io/lumen/reference/yuenTTest.md)
· [`zTest()`](https://andrisignorell.github.io/lumen/reference/zTest.md)
·
[`varTest()`](https://andrisignorell.github.io/lumen/reference/varTest.md)
·
[`leveneTest()`](https://andrisignorell.github.io/lumen/reference/leveneTest.md)
·
[`siegelTukeyTest()`](https://andrisignorell.github.io/lumen/reference/siegelTukeyTest.md)

### Post hoc and multiple comparisons

[`postHoc()`](https://andrisignorell.github.io/lumen/reference/postHoc.md)
(Tukey HSD, Bonferroni, LSD, Scheffé, Newman-Keuls, Duncan) ·
[`scheffeTest()`](https://andrisignorell.github.io/lumen/reference/scheffeTest.md)
·
[`conoverTest()`](https://andrisignorell.github.io/lumen/reference/conoverTest.md)
·
[`dunnTest()`](https://andrisignorell.github.io/lumen/reference/dunnTest.md)
·
[`dunnettTest()`](https://andrisignorell.github.io/lumen/reference/dunnettTest.md)
·
[`nemenyiTest()`](https://andrisignorell.github.io/lumen/reference/nemenyiTest.md)
·
[`steelTest()`](https://andrisignorell.github.io/lumen/reference/steelTest.md)
·
[`dscfTest()`](https://andrisignorell.github.io/lumen/reference/dscfTest.md)
·
[`vanWaerdenTest()`](https://andrisignorell.github.io/lumen/reference/vanWaerdenTest.md)
·
[`jonckheereTerpstraTest()`](https://andrisignorell.github.io/lumen/reference/jonckheereTerpstraTest.md)
·
[`pageTest()`](https://andrisignorell.github.io/lumen/reference/pageTest.md)

### Power analysis

[`powerChisqTest()`](https://andrisignorell.github.io/lumen/reference/powerChisqTest.md)

### Confidence intervals

[`binomCI()`](https://andrisignorell.github.io/lumen/reference/binomCI.md)
·
[`binomDiffCI()`](https://andrisignorell.github.io/lumen/reference/binomDiffCI.md)
·
[`binomRatioCI()`](https://andrisignorell.github.io/lumen/reference/binomRatioCI.md)
·
[`poissonCI()`](https://andrisignorell.github.io/lumen/reference/poissonCI.md)
·
[`multinomCI()`](https://andrisignorell.github.io/lumen/reference/multinomCI.md)
·
[`meanCI()`](https://andrisignorell.github.io/lumen/reference/meanCI.md)
·
[`meanCIn()`](https://andrisignorell.github.io/lumen/reference/meanCIn.md)
·
[`meanDiffCI()`](https://andrisignorell.github.io/lumen/reference/meanDiffCI.md)
·
[`medianCI()`](https://andrisignorell.github.io/lumen/reference/medianCI.md)
·
[`quantileCI()`](https://andrisignorell.github.io/lumen/reference/quantileCI.md)
· [`varCI()`](https://andrisignorell.github.io/lumen/reference/varCI.md)
· [`sumCI()`](https://andrisignorell.github.io/lumen/reference/sumCI.md)
· [`corCI()`](https://andrisignorell.github.io/lumen/reference/corCI.md)
·
[`fisherZ()`](https://andrisignorell.github.io/lumen/reference/fisherZ.md)
·
[`bootCI()`](https://andrisignorell.github.io/lumen/reference/bootCI.md)

### Distributions and moments

[`dgev()`](https://andrisignorell.github.io/lumen/reference/dpqr-gev.md)/[`pgev()`](https://andrisignorell.github.io/lumen/reference/dpqr-gev.md)/[`qgev()`](https://andrisignorell.github.io/lumen/reference/dpqr-gev.md)/[`rgev()`](https://andrisignorell.github.io/lumen/reference/dpqr-gev.md)
— generalized extreme value ·
[`dgpd()`](https://andrisignorell.github.io/lumen/reference/dpqr-gpd.md)/[`pgpd()`](https://andrisignorell.github.io/lumen/reference/dpqr-gpd.md)/[`qgpd()`](https://andrisignorell.github.io/lumen/reference/dpqr-gpd.md)/[`rgpd()`](https://andrisignorell.github.io/lumen/reference/dpqr-gpd.md)
— generalized Pareto ·
[`dgumbel()`](https://andrisignorell.github.io/lumen/reference/dpqr-gumbel.md)/[`pgumbel()`](https://andrisignorell.github.io/lumen/reference/dpqr-gumbel.md)/[`qgumbel()`](https://andrisignorell.github.io/lumen/reference/dpqr-gumbel.md)/[`rgumbel()`](https://andrisignorell.github.io/lumen/reference/dpqr-gumbel.md)
— Gumbel ·
[`dgumbelx()`](https://andrisignorell.github.io/lumen/reference/dpqr-gumbelx.md)
— reversed Gumbel extremes ·
[`dfrechet()`](https://andrisignorell.github.io/lumen/reference/dpqr-frechet.md)
— Fréchet ·
[`drweibull()`](https://andrisignorell.github.io/lumen/reference/dpqr-rweibull.md)
— reversed Weibull ·
[`dgompertz()`](https://andrisignorell.github.io/lumen/reference/dpqr-gompertz.md)
— Gompertz ·
[`ddirichlet()`](https://andrisignorell.github.io/lumen/reference/dpqr-dirichlet.md)/[`rdirichlet()`](https://andrisignorell.github.io/lumen/reference/dpqr-dirichlet.md)
— Dirichlet ·
[`dtri()`](https://andrisignorell.github.io/lumen/reference/dpqr-tri.md)
— triangular ·
[`dbenford()`](https://andrisignorell.github.io/lumen/reference/dpqr-benford.md)
— Benford’s law ·
[`dorder()`](https://andrisignorell.github.io/lumen/reference/dpqr-order.md)
— order statistics · `mcontinuous`/`mdiscrete`/`mextreme` —
mean/variance for common continuous, discrete, and extreme-value
distributions ·
[`pAD()`](https://andrisignorell.github.io/lumen/reference/pAD.md) —
null distribution of the Anderson-Darling statistic

See `?distributions-overview` for the full family index and conversion
table.

### Utilities

[`scores()`](https://andrisignorell.github.io/lumen/reference/scores.md)
— score transformations (table, rank, ridit, modified ridit) for ordinal
contingency tables ·
[`sphericityEps()`](https://andrisignorell.github.io/lumen/reference/sphericityEps.md)
— Greenhouse-Geisser and Huynh-Feldt sphericity corrections

## Dependencies

`lumen` imports `pharos`, `bedrock`, `mvtnorm`, `stats`, `boot`, and
`gld`. C++ support is provided through `Rcpp`, `RcppParallel`, and
`RcppArmadillo`. `Exact`, `goftest`, `coin`, `flexsurv`, `randtests`,
`lmtest`, `nortest`, `multcomp`, and `tseries` are used in `Suggests`,
for optional methods and test-suite cross-validation.

## Documentation

The development documentation is available at:

<https://andrisignorell.github.io/lumen/>

The source repository is available at:

<https://github.com/AndriSignorell/lumen/>

Issues and feature requests can be submitted at:

<https://github.com/AndriSignorell/lumen/issues>

## License

`lumen` is released under GPL (\>= 2).

## Status

This package is experimental and versioned as `0.0.0.924`. Interfaces
may still change before a stable release.
