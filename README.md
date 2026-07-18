# lumen

`lumen` provides hypothesis tests, confidence intervals, and selected statistical
distributions used across the DescToolsX ecosystem. It is designed as a focused
statistical companion package: methodologically transparent, API-consistent, and
suitable for use in applied statistical workflows.

Every function ships with a documented reference — a textbook formula, an original
paper, or a well-established R implementation it was checked against — rather than
a bare implementation.

The package is currently under active development.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("AndriSignorell/lumen")
```

The package requires R >= 4.2.0.

## Scope

`lumen` collects statistical procedures that are commonly useful in exploratory
analysis, inference, and methodological comparison:

- hypothesis tests for goodness-of-fit, normality, stationarity, randomness,
  contingency tables, marginal homogeneity, and nonparametric group comparisons
- post-hoc procedures and multiple-comparison methods, both parametric
  (Dunnett, Scheffé, Tukey via `postHoc()`) and nonparametric (Conover, Dunn,
  Nemenyi, Steel, DSCF)
- confidence intervals for proportions, differences and ratios of proportions,
  means, medians, variances, sums, correlations, and quantiles, each with
  classical, robust, and bootstrap variants where applicable
- selected probability distributions, including extreme-value (GEV, GPD,
  Gumbel, Fréchet, reversed Weibull), Dirichlet, Gompertz, triangular,
  Benford, and order-statistic distributions, plus their mean/variance
  companions
- generic bootstrap confidence interval and power-analysis helpers

## Examples

### Binomial confidence intervals

```r
library(lumen)

binomCI(x = 37, n = 43, method = "wilson")
binomCI(x = 42, n = 43, method = "clopper-pearson")
```

### Goodness-of-fit testing

```r
x <- rnorm(50)
andersonDarlingTest(x, null = "pnorm", mean = mean(x), sd = sd(x), estimated = TRUE)
```

### Nonparametric tests

```r
x <- c(1.1, 1.4, 1.6, 2.0, 2.2)
y <- c(1.0, 1.2, 1.3, 1.7, 1.9)

siegelTukeyTest(x, y)
```

### Contingency-table tests

```r
tab <- matrix(c(8, 14, 1, 3), nrow = 2)
barnardTest(tab)
```

### Post hoc comparisons after ANOVA

```r
fit <- aov(breaks ~ tension, data = warpbreaks)
postHoc(fit, method = "hsd")
```

### Bootstrap confidence intervals

```r
set.seed(1984)
bootCI(mtcars$mpg, FUN = mean, na.rm = TRUE, bci.method = "basic")
```

## Design principles

`lumen` follows the broader DescToolsX design philosophy:

- predictable lowerCamelCase function names, with `.default`/`.formula`/`.glm`
  S3 methods offered wherever a formula or model-object interface makes sense
- explicit argument validation with informative error messages
- transparent method choices — where a test or interval admits more than one
  established computation (e.g. exact vs. asymptotic, classical vs. bootstrap),
  the alternatives are exposed as an explicit `method`/`type` argument rather
  than hidden behind a single default
- clean separation between user-facing interfaces and computational engines,
  with performance-critical algorithms (exact permutation and rank
  distributions, Wellek's exact Page test tables, run-length distributions)
  implemented in C++ via Rcpp/RcppParallel/RcppArmadillo
- compatibility with familiar base R idioms: most hypothesis tests return
  `htest`-compatible objects usable with `print()`, `$statistic`, `$p.value`,
  and friends

## Function reference

The sections below group the most commonly used functions by purpose. The full
alphabetical reference, including internal helpers and every distribution
variant, is available via `?lumen` and the
[pkgdown site](https://andrisignorell.github.io/lumen/reference/index.html).

### Goodness-of-fit and normality

`andersonDarlingTest()` · `cramerVonMisesTest()` · `jarqueBeraTest()` ·
`lillieTest()` · `pearsonTest()` · `shapiroFranciaTest()`

### Time series: stationarity, autocorrelation, randomness

`adfTest()` · `kpssTest()` · `durbinWatsonTest()` · `breuschGodfreyTest()` ·
`bpTest()` · `bartelsRankTest()` · `runsTest()` · `vonNeumannTest()`

### Categorical and contingency-table tests

`barnardTest()` · `bhapkarTest()` · `breslowDayTest()` · `stuartMaxwellTest()` ·
`woolfTest()` · `cochranArmitageTest()` · `cochranQTest()` · `gTest()` ·
`lehmacherTest()` · `mantelTrendTest()` · `hosmerLemeshowTest()` ·
`leCessieTest()`

### Location and scale tests

`hotellingsT2Test()` · `signTest()` · `tTestA()` · `yuenTTest()` · `zTest()` ·
`varTest()` · `leveneTest()` · `siegelTukeyTest()`

### Post hoc and multiple comparisons

`postHoc()` (Tukey HSD, Bonferroni, LSD, Scheffé, Newman-Keuls, Duncan) ·
`scheffeTest()` · `conoverTest()` · `dunnTest()` · `dunnettTest()` ·
`nemenyiTest()` · `steelTest()` · `sdcfTest()` · `vanWaerdenTest()` ·
`jonckheereTerpstraTest()` · `pageTest()`

### Power analysis

`powerChisqTest()`

### Confidence intervals

`binomCI()` · `binomDiffCI()` · `binomRatioCI()` · `poissonCI()` ·
`multinomCI()` · `meanCI()` · `meanCIn()` · `meanDiffCI()` · `medianCI()` ·
`quantileCI()` · `varCI()` · `sumCI()` · `corCI()` · `fisherZ()` ·
`bootCI()`

### Distributions and moments

`dgev()`/`pgev()`/`qgev()`/`rgev()` — generalized extreme value ·
`dgpd()`/`pgpd()`/`qgpd()`/`rgpd()` — generalized Pareto ·
`dgumbel()`/`pgumbel()`/`qgumbel()`/`rgumbel()` — Gumbel ·
`dgumbelx()` — reversed Gumbel extremes ·
`dfrechet()` — Fréchet · `drweibull()` — reversed Weibull ·
`dgompertz()` — Gompertz · `ddirichlet()`/`rdirichlet()` — Dirichlet ·
`dtri()` — triangular · `dbenford()` — Benford's law ·
`dorder()` — order statistics ·
`mcontinuous`/`mdiscrete`/`mextreme` — mean/variance for common continuous,
discrete, and extreme-value distributions ·
`pAD()` — null distribution of the Anderson-Darling statistic

See `?distributions-overview` for the full family index and conversion table.

### Utilities

`scores()` — score transformations (table, rank, ridit, modified ridit) for
ordinal contingency tables · `sphericityEps()` — Greenhouse-Geisser and
Huynh-Feldt sphericity corrections

## Dependencies

`lumen` imports `pharos`, `bedrock`, `mvtnorm`, `stats`, `boot`, and `gld`.
C++ support is provided through `Rcpp`, `RcppParallel`, and `RcppArmadillo`.
`Exact`, `goftest`, `coin`, `flexsurv`, `randtests`, `lmtest`, `nortest`,
`multcomp`, and `tseries` are used in `Suggests`, for optional methods and
test-suite cross-validation.

## Documentation

The development documentation is available at:

<https://andrisignorell.github.io/lumen/>

The source repository is available at:

<https://github.com/AndriSignorell/lumen/>

Issues and feature requests can be submitted at:

<https://github.com/AndriSignorell/lumen/issues>

## License

`lumen` is released under GPL (>= 2).

## Status

This package is experimental and versioned as `0.0.0.924`. Interfaces may
still change before a stable release.
