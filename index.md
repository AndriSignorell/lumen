# 📦 lumen

**Title:** Statistical Tests, Confidence Intervals, and Distributions  
**License:** GPL (≥ 2)

## 🧩 Overview

`lumen` is the inferential layer of the **DescToolsX ecosystem**. It
provides classical, robust, nonparametric and specialised hypothesis
tests, confidence intervals for the common estimands, and a set of
probability distributions that base R does not carry.

Tests return ordinary `"htest"` objects and intervals follow one
argument convention throughout, so procedures from very different
literatures can be combined without translating between interfaces.

The package is self-contained and does not require the higher-level
packages of the suite.

📖 **Documentation:** <https://andrisignorell.github.io/lumen/>

## ⚙️ Installation

``` r

install.packages("lumen")
```

Or the development version from GitHub:

``` r

remotes::install_github("AndriSignorell/lumen")
```

## 📚 Core Features

### 🔹 Location Tests

- [`tTestA()`](reference/tTestA.md), [`zTest()`](reference/zTest.md),
  [`yuenTTest()`](reference/yuenTTest.md),
  [`signTest()`](reference/signTest.md)
- [`brunnerMunzelTest()`](reference/brunnerMunzelTest.md),
  [`siegelTukeyTest()`](reference/siegelTukeyTest.md),
  [`mosesTest()`](reference/mosesTest.md)
- [`moodMedianTest()`](reference/moodMedianTest.md),
  [`vanWaerdenTest()`](reference/vanWaerdenTest.md),
  [`dscfTest()`](reference/dscfTest.md)

### 🔹 Categorical Data

- [`barnardTest()`](reference/barnardTest.md),
  [`gTest()`](reference/gTest.md),
  [`cochranArmitageTest()`](reference/cochranArmitageTest.md),
  [`cochranQTest()`](reference/cochranQTest.md)
- [`breslowDayTest()`](reference/breslowDayTest.md),
  [`woolfTest()`](reference/woolfTest.md),
  [`mantelTrendTest()`](reference/mantelTrendTest.md)
- [`bhapkarTest()`](reference/bhapkarTest.md),
  [`stuartMaxwellTest()`](reference/stuartMaxwellTest.md),
  [`lehmacherTest()`](reference/lehmacherTest.md)

### 🔹 Goodness of Fit and Normality

- [`andersonDarlingTest()`](reference/andersonDarlingTest.md),
  [`cramerVonMisesTest()`](reference/cramerVonMisesTest.md),
  [`lillieTest()`](reference/lillieTest.md)
- [`shapiroFranciaTest()`](reference/shapiroFranciaTest.md),
  [`jarqueBeraTest()`](reference/jarqueBeraTest.md),
  [`pearsonTest()`](reference/pearsonTest.md)
- [`hosmerLemeshowTest()`](reference/hosmerLemeshowTest.md),
  [`leCessieTest()`](reference/leCessieTest.md) — for logistic models

### 🔹 Multiple Comparisons and Post-hoc

- [`postHoc()`](reference/postHoc.md),
  [`plot.PostHocTest()`](reference/plot.PostHocTest.md)
- [`dunnTest()`](reference/dunnTest.md),
  [`dunnettTest()`](reference/dunnettTest.md),
  [`conoverTest()`](reference/conoverTest.md),
  [`nemenyiTest()`](reference/nemenyiTest.md)
- [`gamesHowellTest()`](reference/gamesHowellTest.md),
  [`scheffeTest()`](reference/scheffeTest.md),
  [`steelTest()`](reference/steelTest.md)

### 🔹 Trend, Randomness and Time Series

- [`adfTest()`](reference/adfTest.md),
  [`kpssTest()`](reference/kpssTest.md),
  [`runsTest()`](reference/runsTest.md),
  [`vonNeumannTest()`](reference/vonNeumannTest.md)
- [`bartelsRankTest()`](reference/BartelsRankTest.md),
  [`coxStuartTest()`](reference/coxStuartTest.md),
  [`jonckheereTerpstraTest()`](reference/jonckheereTerpstraTest.md)
- [`pageTest()`](reference/pageTest.md),
  [`durbinWatsonTest()`](reference/durbinWatsonTest.md),
  [`breuschGodfreyTest()`](reference/breuschGodfreyTest.md),
  [`bpTest()`](reference/bpTest.md)

### 🔹 Confidence Intervals

- Proportions: [`binomCI()`](reference/binomCI.md),
  [`binomDiffCI()`](reference/binomDiffCI.md),
  [`binomRatioCI()`](reference/binomRatioCI.md),
  [`multinomCI()`](reference/multinomCI.md)
- Counts and rates: [`poissonCI()`](reference/poissonCI.md),
  [`poissonDiffCI()`](reference/poissonDiffCI.md),
  [`poissonRatioCI()`](reference/poissonRatioCI.md)
- Location and scale: [`meanCI()`](reference/meanCI.md),
  [`meanDiffCI()`](reference/meanDiffCI.md),
  [`medianCI()`](reference/medianCI.md),
  [`quantileCI()`](reference/quantileCI.md),
  [`varCI()`](reference/varCI.md), [`sumCI()`](reference/sumCI.md)
- Correlation: [`corCI()`](reference/corCI.md),
  [`fisherZ()`](reference/fisherZ.md)
- General: [`bootCI()`](reference/bootCI.md) for arbitrary statistics

### 🔹 Distributions

`dpqr` families for the Benford, Dirichlet, extreme value, Fréchet, GEV,
GPD, Gompertz, Gumbel, reversed Gumbel, order-statistic, reverse Weibull
and triangular distributions, plus
[`cont.moments()`](reference/cont.moments.md),
[`disc.moments()`](reference/disc.moments.md) and
`extreme-value-moments()`.

### 🔹 Power and Sample Size

- [`meanCIn()`](reference/meanCIn.md),
  [`powerChisqTest()`](reference/powerChisqTest.md),
  [`sphericityEps()`](reference/sphericityEps.md),
  [`scores()`](reference/scores.md)

## 🚀 Design Principles

- **Consistent** — lowerCamelCase API, uniform argument names and
  ordering across the whole DescToolsX suite
- **Standard output** — tests return `"htest"` objects and print like
  the ones in base R
- **Fast** — bootstrap and permutation routines implemented in C++ via
  Rcpp, RcppArmadillo and RcppParallel
- **Documented** — references to the original literature on every
  procedure

## 🧪 Example

``` r

library(lumen)

# a confidence interval for a proportion, several methods
binomCI(37, 100, method = "wilson")

# nonparametric alternative to the t test
brunnerMunzelTest(extra ~ group, data = sleep)

# post-hoc comparisons after an ANOVA
fit <- aov(breaks ~ tension, data = warpbreaks)
postHoc(fit, method = "hsd")

# goodness of fit
andersonDarlingTest(rnorm(100), null = "pnorm")
```

## 🧱 The Suite

`lumen` builds on `bedrock` (base utilities) and `pharos` (graphics).
`DescToolsX` (descriptive statistics), `alloy` (modelling), `pons`
(MS-Office) and `swissValet` (RStudio addins) complete the family.

## 🙏 Acknowledgements

Parts of the code and documentation were reviewed with the help of large
language models (OpenAI Codex, Anthropic Claude). Every suggestion was
assessed, edited and verified by the maintainer, who remains solely
responsible for the content of this package.

## 📜 License

GPL (≥ 2)
