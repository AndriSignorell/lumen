# lumen

`lumen` provides hypothesis tests, confidence intervals, and selected statistical distribution utilities used in the DescToolsX ecosystem. It is designed as a focused statistical companion package: methodologically transparent, API-consistent, and suitable for use in applied statistical workflows.

The package is currently under active development.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("AndriSignorell/lumen")
```

The package requires R >= 4.2.0.

## Scope

`lumen` collects statistical procedures that are commonly useful in exploratory analysis, inference, and methodological comparison. The current package includes:

- hypothesis tests for goodness-of-fit, normality, stationarity, randomness, contingency tables, marginal homogeneity, and nonparametric group comparisons
- confidence intervals for proportions, differences and ratios of proportions, means, medians, variances, correlations, regression coefficients, and related quantities
- post-hoc procedures and multiple-comparison helpers
- selected distribution functions, including extreme-value, Dirichlet, Gompertz, triangular, Benford, and order-statistic distributions
- bootstrap confidence interval helpers

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

### Bootstrap confidence intervals

```r
set.seed(1984)
bootCI(mtcars$mpg, FUN = mean, na.rm = TRUE, bci.method = "basic")
```

## Design principles

`lumen` follows the broader DescToolsX design philosophy:

- predictable lowerCamelCase function names
- explicit argument validation
- transparent method choices
- clean separation between user-facing interfaces and computational engines
- compatibility with familiar base R idioms where appropriate

Most hypothesis tests return objects compatible with the standard `htest` interface. Confidence interval functions generally return compact vectors or matrices with estimates and interval bounds.

## Selected functions

### Hypothesis tests

- `andersonDarlingTest()` — Anderson-Darling goodness-of-fit test
- `bartelsRankTest()` — Bartels rank test for randomness
- `barnardTest()` — Barnard's unconditional test for 2 x 2 tables
- `bhapkarTest()` — Bhapkar marginal homogeneity test
- `breslowDayTest()` — Breslow-Day test for homogeneity of odds ratios
- `cochranArmitageTest()` — Cochran-Armitage trend test
- `cochranQTest()` — Cochran's Q test
- `cramerVonMisesTest()` — Cramer-von Mises goodness-of-fit test
- `durbinWatsonTest()` — Durbin-Watson test
- `jarqueBeraTest()` — Jarque-Bera normality test
- `kpssTest()` — KPSS stationarity test
- `leveneTest()` — Levene test for equality of variances
- `lillieTest()` — Lilliefors normality test
- `siegelTukeyTest()` — Siegel-Tukey test for scale differences
- `stuartMaxwellTest()` — Stuart-Maxwell marginal homogeneity test
- `woolfTest()` — Woolf test for homogeneity of odds ratios

### Confidence intervals

- `binomCI()` — confidence intervals for binomial proportions
- `binomDiffCI()` — confidence intervals for differences of binomial proportions
- `binomRatioCI()` — confidence intervals for ratios of binomial proportions
- `bootCI()` — bootstrap confidence intervals
- `corCI()` — confidence intervals for correlations
- `meanCI()` — confidence intervals for means
- `medianCI()` — confidence intervals for medians
- `multinomCI()` — confidence intervals for multinomial proportions
- `poissonCI()` — confidence intervals for Poisson rates
- `quantileCI()` — confidence intervals for quantiles
- `varCI()` — confidence intervals for variances

### Distributions and utilities

- `dgev()`, `pgev()`, `qgev()`, `rgev()` — generalized extreme value distribution
- `dgpd()`, `pgpd()`, `qgpd()`, `rgpd()` — generalized Pareto distribution
- `dgumbel()`, `pgumbel()`, `qgumbel()`, `rgumbel()` — Gumbel distribution
- `ddirichlet()`, `pdirichlet()`, `qdirichlet()`, `rdirichlet()` — Dirichlet distribution
- `dtri()` — triangular distribution
- `dbenford()` — Benford distribution
- `scores()` — score generation helper

## Dependencies

`lumen` imports several packages used for statistical computation and infrastructure, including `boot`, `aurora`, `bedrock`, `Exact`, `mvtnorm`, `stats`, `withr`, and `gld`. C++ support is provided through `Rcpp`, `RcppParallel`, and `RcppArmadillo`.

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

This package is experimental and versioned as `0.0.0.907`. Interfaces may still change before a stable release.

