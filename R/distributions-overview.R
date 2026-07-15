#' Distribution Functions in lumen
#'
#' @description
#' lumen provides density (\code{d}), CDF (\code{p}), quantile (\code{q}),
#' and random-generation (\code{r}) functions for several distributions
#' not covered by base R, plus moment (\code{m}) functions -- mean and
#' variance -- for those and for several distributions base R already
#' provides d/p/q/r for. A \code{.} marks a combination with no function.
#'
#' @section Extreme value distributions:
#'
#' | Distribution | d | p | q | r | m |
#' |---|---|---|---|---|---|
#' | Generic extreme | [dextreme()] | [pextreme()] | [qextreme()] | [rextreme()] | . |
#' | Gen. Extreme Value | [dgev()] | [pgev()] | [qgev()] | [rgev()] | [mgev()] |
#' | Gumbel | [dgumbel()] | [pgumbel()] | [qgumbel()] | [rgumbel()] | [mgumbel()] |
#' | Gumbel, maximum of two | [dgumbelx()] | [pgumbelx()] | [qgumbelx()] | [rgumbelx()] | . |
#' | Fréchet | [dfrechet()] | [pfrechet()] | [qfrechet()] | [rfrechet()] | [mfrechet()] |
#' | Reversed Weibull | [drweibull()] | [prweibull()] | [qrweibull()] | [rrweibull()] | [mrweibull()] |
#' | Reverse Gumbel | [dRevGumbel()] | [pRevGumbel()] | [qRevGumbel()] | [rRevGumbel()] | . |
#' | Order statistics | [dorder()] | [porder()] | . | [rorder()] | . |
#'
#' Order statistics have no quantile function (\code{qorder()} does not
#' exist). \code{\link{qRevGumbelExp}} is a specialized quantile for the
#' exponential parametrization of the reverse Gumbel distribution, not a
#' general-purpose \code{q} slot -- use \code{\link{qRevGumbel}} for that.
#'
#' @section Other distributions:
#'
#' | Distribution | d | p | q | r | m |
#' |---|---|---|---|---|---|
#' | Benford | [dpqr-benford] | [mbenford()] |
#' | Dirichlet | [dpqr-dirichlet] | . |
#' | Gompertz | [dpqr-gompertz] | [mgompertz()] |
#' | Gen. Pareto | [dpqr-gpd] | [mgpd()] |
#' | Triangular | [dpqr-tri] | [mtri()] |
#'
#' @section Moments for base R distributions:
#' For these, \code{d}/\code{p}/\code{q}/\code{r} already exist in
#' \pkg{stats} -- lumen only adds the moments function.
#'
#' | Distribution | stats:: | Moments |
#' |---|---|---|
#' | Beta | [stats::Beta] | [mbeta()] |
#' | Binomial | \link[stats:Binomial]{Binomial} | [mbinom()] |
#' 
#' 
#' | Chi-squared | \code{dchisq}/\code{pchisq}/\code{qchisq}/\code{rchisq} | [mchisq()] |
#' | Exponential | \code{dexp}/\code{pexp}/\code{qexp}/\code{rexp} | [mexp()] |
#' | F | \code{df}/\code{pf}/\code{qf}/\code{rf} | [mf()] |
#' | Gamma | \code{dgamma}/\code{pgamma}/\code{qgamma}/\code{rgamma} | [mgamma()] |
#' | Geometric | \code{dgeom}/\code{pgeom}/\code{qgeom}/\code{rgeom} | [mgeom()] |
#' | Hypergeometric | \code{dhyper}/\code{phyper}/\code{qhyper}/\code{rhyper} | [mhyper()] |
#' | Log-normal | \code{dlnorm}/\code{plnorm}/\code{qlnorm}/\code{rlnorm} | [mlnorm()] |
#' | Negative binomial | \code{dnbinom}/\code{pnbinom}/\code{qnbinom}/\code{rnbinom} | [mnbinom()] |
#' | Normal | \code{dnorm}/\code{pnorm}/\code{qnorm}/\code{rnorm} | [mnorm()] |
#' | Poisson | \code{dpois}/\code{ppois}/\code{qpois}/\code{rpois} | [mpois()] |
#' | Student's t | \code{dt}/\code{pt}/\code{qt}/\code{rt} | [mt()] |
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
