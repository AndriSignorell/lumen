
#' lumen: Tests and Distributions for DescToolsX
#'
#' @description
#' The **lumen** package provides a comprehensive collection of
#' statistical hypothesis tests and probability distributions designed to
#' complement the functionality of \pkg{DescToolsX}.
#'
#' The package focuses on a consistent user interface, extended methodological
#' coverage, and seamless integration of inferential procedures and distribution
#' functions within a unified framework.
#'
#' @section Hypothesis Tests:
#'
#' The package implements a wide range of statistical tests:
#'
#' **Goodness-of-Fit Tests**
#' \itemize{
#'   \item `andersonDarlingTest`, `cramerVonMisesTest`
#'   \item `lillieTest`, `jarqueBeraTest`, `shapiroFranciaTest`
#'   \item `pearsonTest`
#' }
#'
#' **Nonparametric Tests**
#' \itemize{
#'   \item `signTest`, `jonckheereTerpstraTest`, `pageTest`
#'   \item `mosesTest`, `siegelTukeyTest`, `vanWaerdenTest`
#' }
#'
#' **Post-hoc Procedures**
#' \itemize{
#'   \item `dunnTest`, `conoverTest`, `nemenyiTest`
#'   \item `scheffeTest`, `postHoc`
#' }
#'
#' **Parametric Tests**
#' \itemize{
#'   \item `tTestA`, `yuenTTest`, `zTest`, `varTest`
#'   \item `hotellingsT2Test`, `leveneTest`
#' }
#'
#' **Contingency Table Tests**
#' \itemize{
#'   \item `barnardTest`, `bhapkarTest`, `breslowDayTest`
#'   \item `cochranArmitageTest`, `cochranQTest`
#'   \item `mantelTrendTest`, `woolfTest`, `stuartMaxwellTest`
#' }
#'
#' **Time Series Tests**
#' \itemize{
#'   \item `adfTest`, `kpssTest`
#'   \item `durbinWatsonTest`, `breuschGodfreyTest`
#' }
#'
#' **Randomness and Independence Tests**
#' \itemize{
#'   \item `runsTest`, `BartelsRankTest`, `vonNeumannTest`
#' }
#'
#' @section Probability Distributions:
#'
#' The package provides density (`d*`), distribution (`p*`),
#' quantile (`q*`) and random generation (`r*`) functions for
#' a range of distributions, following the standard R conventions.
#'
#' **Extreme Value Distributions**
#' \itemize{
#'   \item Generalized Extreme Value: `dgev`, `pgev`, `qgev`, `rgev`
#'   \item Generalized Pareto: `dgpd`, `pgpd`, `qgpd`, `rgpd`
#'   \item Gumbel and extended Gumbel: `dgumbel`, `dgumbelx`, ...
#'   \item Frechet and reverse Weibull: `dfrechet`, `drweibull`, ...
#'   \item Maxima/minima distributions: `dextreme`, `pextreme`, ...
#' }
#'
#' **Special Distributions**
#' \itemize{
#'   \item Benford distribution: `dbenford`, `pbenford`, `qbenford`, `rbenford`
#'   \item Order and triangular distributions: `dorder`, `dtri`, ...
#' }
#'
#' @section Utilities:
#'
#' \itemize{
#'   \item `scores` – Score generation for ordinal contingency tables
#'   \item `corTest` – Fast correlation testing for matrices
#' }
#'
#' @section Design Principles:
#'
#' \itemize{
#'   \item Unified interface across hypothesis tests and distributions
#'   \item Standard `d/p/q/r` naming for distributions
#'   \item Clear classification via `@family` and `@concept`
#'   \item Separation of statistical procedures and data transformation utilities
#' }
#'
#' @keywords internal
#_PACKAGE"
#' @name lumen-package
NULL
