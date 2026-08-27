# Distribution Functions in lumen

lumen provides density (`d`), CDF (`p`), quantile (`q`), and
random-generation (`r`) functions for several distributions not covered
by base R, plus moment (`m`) functions – mean and variance – for those
and for several distributions base R already provides d-p-q-r for. A `-`
marks a combination with no function.

## Extreme value distributions

|  |  |  |
|----|----|----|
| **Distribution** | **d-p-q-r** | **Moments** |
| Generic extreme | [dpqr-extreme](dpqr-extreme.md) | `-` |
| Gen. Extreme Value | [dpqr-gev](dpqr-gev.md) | [`mgev()`](extreme-value-moments.md) |
| Gumbel | [dpqr-gumbel](dpqr-gumbel.md) | [`mgumbel()`](extreme-value-moments.md) |
| Gumbel, maximum of two ` ` | [dpqr-gumbelx](dpqr-gumbelx.md) | `-` |
| Fréchet | [dpqr-frechet](dpqr-frechet.md) | [`mfrechet()`](extreme-value-moments.md) |
| Reversed Weibull | [dpqr-rweibull](dpqr-rweibull.md) | [`mrweibull()`](extreme-value-moments.md) |
| Reverse Gumbel | [dpqr-RevGumbel](dpqr-RevGumbel.md) ` ` | `-` |
| Order statistics | [dpqr-order](dpqr-order.md) | `-` |

Order statistics have no quantile function (`qorder()` does not exist).
[`qRevGumbelExp`](dpqr-RevGumbel.md) is a specialized quantile for the
exponential parametrization of the reverse Gumbel distribution, not a
general-purpose `q` slot – use [`qRevGumbel`](dpqr-RevGumbel.md) for
that.

## Other distributions

|  |  |  |
|----|----|----|
| **Distribution** | **d-p-q-r** | **Moments** |
| Benford | [dpqr-benford](dpqr-benford.md) | [`mbenford()`](disc.moments.md) |
| Dirichlet | [dpqr-dirichlet](dpqr-dirichlet.md) | `-` |
| Gompertz | [dpqr-gompertz](dpqr-gompertz.md) ` ` | [`mgompertz()`](extreme-value-moments.md) |
| Gen. Pareto ` ` | [dpqr-gpd](dpqr-gpd.md) | [`mgpd()`](extreme-value-moments.md) |
| Triangular | [dpqr-tri](dpqr-tri.md) | [`mtri()`](cont.moments.md) |

## Moments for base R distributions

For these, `d`-`p`-`q`-`r` already exist in stats – lumen only adds the
moments function.

|  |  |  |
|----|----|----|
| **Distribution** | **d-p-q-r** | **Moments** |
| Beta | [stats::Beta](https://rdrr.io/r/stats/Beta.html) | [`mbeta()`](cont.moments.md) |
| Binomial | [stats::Binomial](https://rdrr.io/r/stats/Binomial.html) | [`mbinom()`](disc.moments.md) |
| Chi-squared | [stats::Chisquare](https://rdrr.io/r/stats/Chisquare.html) | [`mchisq()`](cont.moments.md) |
| Exponential | [stats::Exponential](https://rdrr.io/r/stats/Exponential.html) | [`mexp()`](cont.moments.md) |
| F | [stats::FDist](https://rdrr.io/r/stats/Fdist.html) | [`mf()`](cont.moments.md) |
| Gamma | [stats::GammaDist](https://rdrr.io/r/stats/GammaDist.html) | [`mgamma()`](cont.moments.md) |
| Geometric | [stats::Geometric](https://rdrr.io/r/stats/Geometric.html) | [`mgeom()`](disc.moments.md) |
| Hypergeometric | [stats::Hypergeometric](https://rdrr.io/r/stats/Hypergeometric.html) ` ` | [`mhyper()`](disc.moments.md) |
| Log-normal | [stats::Lognormal](https://rdrr.io/r/stats/Lognormal.html) | [`mlnorm()`](cont.moments.md) |
| Negative binomial ` ` | [stats::NegBinomial](https://rdrr.io/r/stats/NegBinomial.html) | [`mnbinom()`](disc.moments.md) |
| Normal | [stats::Normal](https://rdrr.io/r/stats/Normal.html) | [`mnorm()`](cont.moments.md) |
| Poisson | [stats::Poisson](https://rdrr.io/r/stats/Poisson.html) | [`mpois()`](disc.moments.md) |
| Student's t | [stats::TDist](https://rdrr.io/r/stats/TDist.html) | [`mt()`](cont.moments.md) |

## Standalone

[`rsum1()`](rsum1.md) generates a Dirichlet-distributed sample (which
sums to 1 by construction), related to but not part of the
[`ddirichlet`](dpqr-dirichlet.md)/[`pdirichlet`](dpqr-dirichlet.md)/[`qdirichlet`](dpqr-dirichlet.md)/[`rdirichlet`](dpqr-dirichlet.md)
suite above.

[`pAD()`](pAD.md)/[`qAD()`](pAD.md) are the CDF and quantile function of
the Anderson-Darling test statistic's null distribution – supporting
functions for [`andersonDarlingTest`](andersonDarlingTest.md), not a
distribution family of their own (no `d`/`r` counterparts).

## See also

[stats::Distributions](https://rdrr.io/r/stats/Distributions.html)
