#' Distribution Functions in lumen
#'
#' @description
#' lumen provides density (`d`), CDF (`p`), quantile (`q`),
#' and random-generation (`r`) functions for several distributions
#' not covered by base R, plus moment (`m`) functions -- mean and
#' variance -- for those and for several distributions base R already
#' provides d-p-q-r for. A `-` marks a combination with no function.
#'
#' @section Extreme value distributions:
#'
#' | **Distribution** | **d-p-q-r** | **Moments** |
#' |---|---|---|
#' | Generic extreme | [dpqr-extreme] | `-` |
#' | Gen. Extreme Value | [dpqr-gev] | [mgev()] |
#' | Gumbel | [dpqr-gumbel] | [mgumbel()] |
#' | Gumbel, maximum of two \verb{  } | [dpqr-gumbelx] | `-` |
#' | Fréchet | [dpqr-frechet] | [mfrechet()] |
#' | Reverse Weibull | [dpqr-revweibull] | [mrevweibull()] |
#' | Reverse Gumbel | [dpqr-revgumbel] \verb{  } | [mrevgumbel()] |
#' | Order statistics | [dpqr-order] | `-` |
#'
#' Order statistics have no quantile function (`qorder()` does not
#' exist). [qrevgumbelExp()] is a specialized quantile for the
#' exponential parametrization of the reverse Gumbel distribution, not a
#' general-purpose `q` slot -- use [qrevgumbel()] for that.
#'
#' @section Other distributions:
#'
#' | **Distribution** | **d-p-q-r** | **Moments** |
#' |---|---|---|
#' | Benford | [dpqr-benford] | [mbenford()] |
#' | Dirichlet | [dpqr-dirichlet] | `-` |
#' | Gompertz | [dpqr-gompertz] \verb{  } | [mgompertz()] |
#' | Gen. Pareto \verb{  } | [dpqr-gpd] | [mgpd()] |
#' | Triangular | [dpqr-tri] | [mtri()] |
#'
#' @section Moments for base R distributions:
#' For these, `d`-`p`-`q`-`r` already exist in
#' \pkg{stats} -- lumen only adds the moments function.
#'
#' | **Distribution** | **d-p-q-r** | **Moments** |
#' |---|---|---|
#' | Beta | [Beta] | [mbeta()] |
#' | Binomial | [Binomial] | [mbinom()] |
#' | Chi-squared | [Chisquare] | [mchisq()] |
#' | Exponential | [Exponential] | [mexp()] |
#' | F | [FDist] | [mf()] |
#' | Gamma | [GammaDist] | [mgamma()] |
#' | Geometric | [Geometric] | [mgeom()] |
#' | Hypergeometric | [Hypergeometric] \verb{  } | [mhyper()] |
#' | Log-normal | [Lognormal] | [mlnorm()] |
#' | Negative binomial \verb{  }| [NegBinomial] | [mnbinom()] |
#' | Normal | [Normal] | [mnorm()] |
#' | Poisson | [Poisson] | [mpois()] |
#' | Student's t | [TDist] | [mt()] |
#'
#' @section Standalone:
#' [rsum1()] generates a Dirichlet-distributed sample (which sums to 1
#' by construction), related to but not part of the
#' [ddirichlet()]/[pdirichlet()]/[qdirichlet()]/[rdirichlet()]
#' suite above.
#'
#' [pAD()]/[qAD()] are the CDF and quantile function of the Anderson-Darling
#' test statistic's null distribution -- supporting functions for
#' [andersonDarlingTest()], not a distribution family of their
#' own (no `d`/`r` counterparts).
#'
#' @name distributions-overview
#' @seealso [Distributions]
NULL
