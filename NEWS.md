# lumen 0.1.2

First release on CRAN.

lumen is the inferential layer of the DescToolsX package suite. It brings
together classical, robust, nonparametric and specialised hypothesis tests,
confidence intervals for the common estimands, and probability
distributions not available in base R. The package can be used on its own;
it does not require the higher-level packages of the suite.

## One interface throughout

- All tests return ordinary `"htest"` objects, so they print, extract and
  combine like the procedures in base R.
- Confidence intervals follow one argument convention (`conf.level`,
  `sides`, `method`) across proportions, counts, location, scale and
  correlation, and return the estimate with its bounds `lci` and `uci`.
- Tests with group comparisons offer both a default and a formula
  interface.

## Highlights

- **Location and scale:** Brunner-Munzel with exact and Monte Carlo
  permutation inference (`brunnerMunzelTest()`), Mood's median test,
  Yuen's trimmed t-test, Moses, Siegel-Tukey, Levene and variance tests.
  Confidence intervals for the difference of two medians
  (`medianDiffCI()`) by the distribution-free method of Bonett and Price,
  which does not assume equally shaped distributions, or by bootstrap.
- **Trend and order:** Jonckheere-Terpstra with exact inference in the
  presence of ties (`jonckheereTerpstraTest()`), Cox-Stuart, Page,
  Cochran-Armitage and Mantel trend tests.
- **Categorical data:** Barnard's unconditional exact test with its own
  C++ engine (`barnardTest()`), G-test, Breslow-Day, Woolf,
  Stuart-Maxwell, Bhapkar, Cochran's Q and Lehmacher tests.
- **Proportions and counts:** `binomCI()` with 18 methods, including
  Blaker and Wang's admissible exact interval, plus intervals for
  differences, ratios and paired designs, hypergeometric (`hyperCI()`),
  multinomial and Poisson rates.
- **Post-hoc comparisons:** `postHocTest()` collects the procedures behind
  one interface, with a plot method; also Games-Howell, Dunn, Conover,
  Nemenyi, Dunnett, Scheffe, Steel and Dwass-Steel-Critchlow-Fligner
  tests, and `signifDiff()` for a compact per-level summary.
- **Goodness of fit and time series:** Anderson-Darling, Cramer-von Mises,
  Lilliefors, Shapiro-Francia, Jarque-Bera; ADF, KPSS, Durbin-Watson,
  Breusch-Godfrey, Breusch-Pagan, runs and von Neumann tests.
- **Bootstrap:** `bootCI()` for arbitrary statistics; confidence intervals
  for the MAD and its two-sample difference and ratio run in parallel C++
  (Rcpp, RcppArmadillo, RcppParallel).
- **Distributions:** Benford, Dirichlet, Frechet, GEV, GPD, Gompertz,
  Gumbel, reverse Weibull, triangular, order statistics and more, with
  `m*()` functions for the moments of common distributions.

## Acknowledgements

lumen builds on code and ideas from many authors, listed as contributors
in the DESCRIPTION file.

Parts of the code and documentation were reviewed with the help of large
language models (OpenAI Codex, Anthropic Claude). Every suggestion was
assessed, edited and verified by the maintainer, who remains solely
responsible for the content of this package.
