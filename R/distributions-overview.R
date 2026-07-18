#' Distribution Functions in lumen
#'
#' @description
#' lumen provides density (\code{d}), CDF (\code{p}), quantile (\code{q}),
#' and random-generation (\code{r}) functions for several distributions
#' not covered by base R, plus moment (\code{m}) functions -- mean and
#' variance -- for those and for several distributions base R already
#' provides d-p-q-r for. A \code{-} marks a combination with no function.
#'
#' @section Extreme value distributions:
#'
#' | \strong{Distribution} | \strong{d-p-q-r} | \strong{Moments} |
#' |---|---|---|
#' | Generic extreme | [dpqr-extreme] | \code{-} |
#' | Gen. Extreme Value | [dpqr-gev] | [mgev()] |
#' | Gumbel | [dpqr-gumbel] | [mgumbel()] |
#' | Gumbel, maximum of two \verb{  } | [dpqr-gumbelx] | \code{-} |
#' | Fréchet | [dpqr-frechet] | [mfrechet()] |
#' | Reversed Weibull | [dpqr-rweibull] | [mrweibull()] |
#' | Reverse Gumbel | [dpqr-RevGumbel] \verb{  } | \code{-} |
#' | Order statistics | [dpqr-order] | \code{-} |
#'
#' Order statistics have no quantile function (\code{qorder()} does not
#' exist). \code{\link{qRevGumbelExp}} is a specialized quantile for the
#' exponential parametrization of the reverse Gumbel distribution, not a
#' general-purpose \code{q} slot -- use \code{\link{qRevGumbel}} for that.
#'
#' @section Other distributions:
#'
#' | \strong{Distribution} | \strong{d-p-q-r} | \strong{Moments} |
#' |---|---|---|
#' | Benford | [dpqr-benford] | [mbenford()] |
#' | Dirichlet | [dpqr-dirichlet] | \code{-} |
#' | Gompertz | [dpqr-gompertz] \verb{  } | [mgompertz()] |
#' | Gen. Pareto \verb{  } | [dpqr-gpd] | [mgpd()] |
#' | Triangular | [dpqr-tri] | [mtri()] |
#'
#' @section Moments for base R distributions:
#' For these, \code{d}-\code{p}-\code{q}-\code{r} already exist in
#' \pkg{stats} -- lumen only adds the moments function.
#'
#' | \strong{Distribution} | \strong{d-p-q-r} | \strong{Moments} |
#' |---|---|---|
#' | Beta | [stats::Beta] | [mbeta()] |
#' | Binomial | [stats::Binomial] | [mbinom()] |
#' | Chi-squared | [stats::Chisquare] | [mchisq()] |
#' | Exponential | [stats::Exponential] | [mexp()] |
#' | F | [stats::FDist] | [mf()] |
#' | Gamma | [stats::GammaDist] | [mgamma()] |
#' | Geometric | [stats::Geometric] | [mgeom()] |
#' | Hypergeometric | [stats::Hypergeometric] \verb{  } | [mhyper()] |
#' | Log-normal | [stats::Lognormal] | [mlnorm()] |
#' | Negative binomial \verb{  }| [stats::NegBinomial] | [mnbinom()] |
#' | Normal | [stats::Normal] | [mnorm()] |
#' | Poisson | [stats::Poisson] | [mpois()] |
#' | Student's t | [stats::TDist] | [mt()] |
#'
#' @section Standalone:
#' [rsum1()] generates a Dirichlet-distributed sample (which sums to 1
#' by construction), related to but not part of the
#' \code{\link{ddirichlet}}/\code{\link{pdirichlet}}/\code{\link{qdirichlet}}/\code{\link{rdirichlet}}
#' suite above.
#'
#' [pAD()]/[qAD()] are the CDF and quantile function of the Anderson-Darling
#' test statistic's null distribution -- supporting functions for
#' \code{\link{andersonDarlingTest}}, not a distribution family of their
#' own (no \code{d}/\code{r} counterparts).
#'
#' @name distributions-overview
#' @seealso [stats::Distributions]
NULL
