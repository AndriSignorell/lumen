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
| Generic extreme | [dpqr-extreme](https://andrisignorell.github.io/lumen/reference/dpqr-extreme.md) | `-` |
| Gen. Extreme Value | [dpqr-gev](https://andrisignorell.github.io/lumen/reference/dpqr-gev.md) | [`mgev()`](https://andrisignorell.github.io/lumen/reference/extreme-value-moments.md) |
| Gumbel | [dpqr-gumbel](https://andrisignorell.github.io/lumen/reference/dpqr-gumbel.md) | [`mgumbel()`](https://andrisignorell.github.io/lumen/reference/extreme-value-moments.md) |
| Gumbel, maximum of two ` ` | [dpqr-gumbelx](https://andrisignorell.github.io/lumen/reference/dpqr-gumbelx.md) | `-` |
| Fréchet | [dpqr-frechet](https://andrisignorell.github.io/lumen/reference/dpqr-frechet.md) | [`mfrechet()`](https://andrisignorell.github.io/lumen/reference/extreme-value-moments.md) |
| Reversed Weibull | [dpqr-rweibull](https://andrisignorell.github.io/lumen/reference/dpqr-rweibull.md) | [`mrweibull()`](https://andrisignorell.github.io/lumen/reference/extreme-value-moments.md) |
| Reverse Gumbel | [dpqr-RevGumbel](https://andrisignorell.github.io/lumen/reference/dpqr-RevGumbel.md) ` ` | `-` |
| Order statistics | [dpqr-order](https://andrisignorell.github.io/lumen/reference/dpqr-order.md) | `-` |

Order statistics have no quantile function (`qorder()` does not exist).
[`qRevGumbelExp`](https://andrisignorell.github.io/lumen/reference/dpqr-RevGumbel.md)
is a specialized quantile for the exponential parametrization of the
reverse Gumbel distribution, not a general-purpose `q` slot – use
[`qRevGumbel`](https://andrisignorell.github.io/lumen/reference/dpqr-RevGumbel.md)
for that.

## Other distributions

|  |  |  |
|----|----|----|
| **Distribution** | **d-p-q-r** | **Moments** |
| Benford | [dpqr-benford](https://andrisignorell.github.io/lumen/reference/dpqr-benford.md) | [`mbenford()`](https://andrisignorell.github.io/lumen/reference/disc.moments.md) |
| Dirichlet | [dpqr-dirichlet](https://andrisignorell.github.io/lumen/reference/dpqr-dirichlet.md) | `-` |
| Gompertz | [dpqr-gompertz](https://andrisignorell.github.io/lumen/reference/dpqr-gompertz.md) ` ` | [`mgompertz()`](https://andrisignorell.github.io/lumen/reference/extreme-value-moments.md) |
| Gen. Pareto ` ` | [dpqr-gpd](https://andrisignorell.github.io/lumen/reference/dpqr-gpd.md) | [`mgpd()`](https://andrisignorell.github.io/lumen/reference/extreme-value-moments.md) |
| Triangular | [dpqr-tri](https://andrisignorell.github.io/lumen/reference/dpqr-tri.md) | [`mtri()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |

## Moments for base R distributions

For these, `d`-`p`-`q`-`r` already exist in stats – lumen only adds the
moments function.

|  |  |  |
|----|----|----|
| **Distribution** | **d-p-q-r** | **Moments** |
| Beta | [stats::Beta](https://rdrr.io/r/stats/Beta.html) | [`mbeta()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |
| Binomial | [stats::Binomial](https://rdrr.io/r/stats/Binomial.html) | [`mbinom()`](https://andrisignorell.github.io/lumen/reference/disc.moments.md) |
| Chi-squared | [stats::Chisquare](https://rdrr.io/r/stats/Chisquare.html) | [`mchisq()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |
| Exponential | [stats::Exponential](https://rdrr.io/r/stats/Exponential.html) | [`mexp()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |
| F | [stats::FDist](https://rdrr.io/r/stats/Fdist.html) | [`mf()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |
| Gamma | [stats::GammaDist](https://rdrr.io/r/stats/GammaDist.html) | [`mgamma()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |
| Geometric | [stats::Geometric](https://rdrr.io/r/stats/Geometric.html) | [`mgeom()`](https://andrisignorell.github.io/lumen/reference/disc.moments.md) |
| Hypergeometric | [stats::Hypergeometric](https://rdrr.io/r/stats/Hypergeometric.html) ` ` | [`mhyper()`](https://andrisignorell.github.io/lumen/reference/disc.moments.md) |
| Log-normal | [stats::Lognormal](https://rdrr.io/r/stats/Lognormal.html) | [`mlnorm()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |
| Negative binomial ` ` | [stats::NegBinomial](https://rdrr.io/r/stats/NegBinomial.html) | [`mnbinom()`](https://andrisignorell.github.io/lumen/reference/disc.moments.md) |
| Normal | [stats::Normal](https://rdrr.io/r/stats/Normal.html) | [`mnorm()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |
| Poisson | [stats::Poisson](https://rdrr.io/r/stats/Poisson.html) | [`mpois()`](https://andrisignorell.github.io/lumen/reference/disc.moments.md) |
| Student's t | [stats::TDist](https://rdrr.io/r/stats/TDist.html) | [`mt()`](https://andrisignorell.github.io/lumen/reference/cont.moments.md) |

## Standalone

[`rsum1()`](https://andrisignorell.github.io/lumen/reference/rsum1.md)
generates a Dirichlet-distributed sample (which sums to 1 by
construction), related to but not part of the
[`ddirichlet`](https://andrisignorell.github.io/lumen/reference/dpqr-dirichlet.md)/[`pdirichlet`](https://andrisignorell.github.io/lumen/reference/dpqr-dirichlet.md)/[`qdirichlet`](https://andrisignorell.github.io/lumen/reference/dpqr-dirichlet.md)/[`rdirichlet`](https://andrisignorell.github.io/lumen/reference/dpqr-dirichlet.md)
suite above.

[`pAD()`](https://andrisignorell.github.io/lumen/reference/pAD.md)/[`qAD()`](https://andrisignorell.github.io/lumen/reference/pAD.md)
are the CDF and quantile function of the Anderson-Darling test
statistic's null distribution – supporting functions for
[`andersonDarlingTest`](https://andrisignorell.github.io/lumen/reference/andersonDarlingTest.md),
not a distribution family of their own (no `d`/`r` counterparts).

## See also

[stats::Distributions](https://rdrr.io/r/stats/Distributions.html)
