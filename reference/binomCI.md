# Confidence Intervals for Binomial Proportions

`binomCI()` computes confidence intervals for binomial proportions using
a wide range of commonly proposed methods.

`binomCIn()` computes the required sample size to obtain a binomial
confidence interval of a specified width, as calculated by `binomCI()`.
The function uses [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) to
numerically solve for the corresponding sample size.

## Usage

``` r
binomCI(
  x,
  n,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("wilson", "wilson-cc", "wilson-mod", "wald", "wald-cc", "jeffreys",
    "jeffreys-mod", "clopper-pearson", "agresti-coull", "pratt", "arcsine", "logit",
    "witting", "mid-p", "blaker", "likelihood", "khouadji"),
  stdEst = TRUE
)

binomCIn(
  p = 0.5,
  width,
  interval = c(1, 100000),
  conf.level = 0.95,
  sides = "two.sided",
  method = "wilson"
)
```

## Arguments

- x:

  number of successes.

- n:

  number of trials.

- conf.level:

  confidence level, defaults to 0.95.

- sides:

  a character string specifying the side of the confidence interval,
  must be one of `"two.sided"` (default), `"left"` or `"right"`. You can
  specify just the initial letter. `"left"` would be analogue to a
  hypothesis of `"greater"` in a `t.test`.

- method:

  character string specifying which method to use; this can be one out
  of: `"wald"`, `"wald-cc"`,`"wilson"` (default), `"wilson-cc"`,
  `"agresti-coull"`, `"jeffreys"`, `"wilson-mod"`, `"jeffreys-mod"`,
  `"clopper-pearson"`, `"arcsine"`, `"logit"`, `"witting"`, `"pratt"`,
  `"mid-p"`, `"likelihood"` and `"blaker"`. All the methods can be asked
  by `".all"`. Abbreviation of method is accepted. See details.

- stdEst:

  logical, specifying if the standard point estimator for the proportion
  value `x/n` should be returned (`TRUE`, default) or the
  method-specific internally used alternative point estimate (`FALSE`).

- p:

  probability for success, defaults to `0.5` as worst case.

- width:

  the width of the confidence interval.

- interval:

  a vector containing the end-points of the interval to be searched for
  the root. The defaults are set to `c(1, 100000)`.

## Value

If recycling yields a single case, a named numeric vector with elements:

- `est`:

  point estimate of the binomial proportion; `x/n` if `stdEst = TRUE`,
  otherwise the method-specific estimate.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

If recycling yields multiple cases, a data frame with one row per case
is returned. Its first three columns are `est`, `lci`, and `uci`; the
remaining columns contain the recycled argument values.

`binomCIn` returns a single numeric value giving the required sample
size.

## Details

All arguments are vectorized and recycled according to standard R rules.

**Wald**: Obtained by inverting the acceptance region of the
large-sample normal (Wald) test.

**Wald with continuity correction**: A continuity-corrected version of
the Wald interval, obtained by adding 1/(2n) to the standard Wald
limits.

**Wilson** (default): Introduced by Wilson (1927), this interval is
obtained by inverting the central limit theorem approximation to the
family of equal-tail tests of \\p = p_0\\. It is recommended by Agresti
and Coull (1998) and Brown et al. (2001). The same interval is returned
as `conf.int` by [`prop.test`](https://rdrr.io/r/stats/prop.test.html)
with `correct = FALSE`.

**Wilson with continuity correction**: A continuity-corrected
modification of the Wilson interval. This corresponds to
[`prop.test`](https://rdrr.io/r/stats/prop.test.html) with
`correct = TRUE`.

**Modified Wilson**: An adjustment of the Wilson interval for extreme
counts (i.e., \\x\\ close to 0 or \\n\\), as proposed by Brown et al.
(2001).

**Agresti-Coull**: A simplified modification of the Wilson interval
(Agresti and Coull, 1998). These intervals are never shorter than the
Wilson intervals (Brown et al., 2001). The internally used adjusted
estimator \\\tilde{p}\\ is returned as an attribute.

**Jeffreys**: The equal-tailed Bayesian interval based on the Jeffreys
prior, as described in Brown et al. (2001).

**Modified Jeffreys**: A modification of the Jeffreys interval for
boundary cases (e.g., \\x = 0\\, \\x = n\\, or near-boundary values),
following Brown et al. (2001).

**Clopper-Pearson**: The so-called exact interval, based on quantiles of
the corresponding beta distribution.

**Arcsine**: Based on the variance-stabilizing arcsine transformation
for the binomial distribution.

**Logit**: Obtained by constructing a Wald-type interval on the log-odds
scale and transforming back to the probability scale.

**Witting**: A randomized procedure (Witting, 1985) providing uniformly
optimal lower and upper confidence bounds for binomial proportions.
Repeated calls may yield slightly different results unless the random
number generator seed is fixed.

**Pratt**: Based on a highly accurate normal approximation (Pratt,
1968).

**Mid-p**: Designed to reduce the conservatism of the Clopper-Pearson
interval. The lower bound \\p_l\\ solves \$\$\frac{1}{2} f(x; n, p_l) +
(1 - F(x; n, p_l)) = \frac{\alpha}{2}\$\$ and the upper bound \\p_u\\
solves \$\$\frac{1}{2} f(x; n, p_u) + F(x - 1; n, p_u) =
\frac{\alpha}{2}\$\$ where \\f\\ and \\F\\ denote the binomial
probability mass and cumulative distribution functions. For \\x = 0\\
the lower bound is set to 0; for \\x = n\\ the upper bound is set to 1.

**Likelihood-based**: Confidence intervals obtained by profiling the
binomial deviance in the neighbourhood of the maximum likelihood
estimator.

**Blaker**: An exact interval based on the method proposed by Blaker
(2000).

**Khouadji**: A transformation-based approximation for binomial
confidence intervals. It applies a variance-stabilizing transformation
to the sample proportion, constructs a normal-based interval, and
back-transforms to the original scale. Compared to the Wald interval it
is more stable for small samples, but it is rarely used in practice
compared to methods like Wilson or exact intervals.

Some methods may produce limits outside the admissible range \\\[0,
1\]\\. In such cases, the bounds are truncated to remain within the
valid parameter space.

For the methods `"wilson"`, `"wilson-cc"`, `"wilson-mod"`,
`"agresti-coull"`, `"witting"`, and `"arcsine"`, the internally used
adjusted point estimator can be returned by setting `stdEst = FALSE`.
These estimators are typically slightly shrunk toward 0.5 compared to
the usual estimator \\x/n\\. See the cited literature for further
details.

**Required Samplesize** (by `binomCIn()`):  
The required sample size for a given confidence interval width depends
on the assumed population proportion. Since this proportion is often
unknown at the planning stage of a study, a conservative approach is to
use the worst-case scenario of \\p = 0.5\\, which maximizes the variance
and therefore yields the largest required sample size. If a more
accurate estimate of the population proportion is available, it can be
used to obtain a smaller required sample size for the same level of
precision.

## **Which interval should be used?**

The Wald interval is known to have poor coverage properties,
particularly for small sample sizes or proportions near 0 or 1. In
contrast, the Clopper-Pearson interval is conservative and often
unnecessarily wide. Brown et al. (2001) recommend the Wilson or Jeffreys
intervals for small sample sizes, and the Agresti-Coull, Wilson, or
Jeffreys intervals for larger samples, as providing more reliable
coverage than most alternatives.

## Contributors

The function is based on earlier work by Matthias Kohl in package
SLmisc, whose original implementation provided the methodological
foundation The implementations of the Pratt, Mid-p, Blaker and Khouadji
methods are based on contributions by Rand R. Wilcox, Michael Hoehle,
Ralph Scherer, and Carl Pearson, respectively.

The current implementation was written and is maintained by Andri
Signorell.

## References

Agresti A. and Coull B.A. (1998) Approximate is better than "exact" for
interval estimation of binomial proportions. *American Statistician*,
**52**, pp. 119-126.

Brown L.D., Cai T.T. and Dasgupta A. (2001) Interval estimation for a
binomial proportion *Statistical Science*, **16**(2), pp. 101-133.

Witting H. (1985) *Mathematische Statistik I*. Stuttgart: Teubner.

Pratt J. W. (1968) A normal approximation for binomial, F, Beta, and
other common, related tail probabilities *Journal of the American
Statistical Association*, 63, 1457- 1483.

Wilcox, R. R. (2005) *Introduction to robust estimation and hypothesis
testing*. Elsevier Academic Press

Newcombe, R. G. (1998) Two-sided confidence intervals for the single
proportion: comparison of seven methods, *Statistics in Medicine*,
17:857-872 https://pubmed.ncbi.nlm.nih.gov/16206245/

Blaker, H. (2000) Confidence curves and improved exact confidence
intervals for discrete distributions, *Canadian Journal of Statistics*
28 (4), 783-798

A. Khouadji (1999) Sur une méthode d’approximation des intervalles de
confiance pour une proportion binomiale.

## See also

[`binom.test`](https://rdrr.io/r/stats/binom.test.html), `binconf`

Other ci.proportion:
[`binomDiffCI()`](https://andrisignorell.github.io/lumen/reference/binomDiffCI.md),
[`binomRatioCI()`](https://andrisignorell.github.io/lumen/reference/binomRatioCI.md),
[`multinomCI()`](https://andrisignorell.github.io/lumen/reference/multinomCI.md)

## Examples

``` r

binomCI(x=37, n=43, 
        method=eval(formals(binomCI)$method))   # return all methods
#>          est       lci       uci  x  n conf.level     sides          method
#> 1  0.8604651 0.7273641 0.9344428 37 43       0.95 two.sided          wilson
#> 2  0.8604651 0.7137335 0.9419725 37 43       0.95 two.sided       wilson-cc
#> 3  0.8604651 0.7273641 0.9344428 37 43       0.95 two.sided      wilson-mod
#> 4  0.8604651 0.7568980 0.9640322 37 43       0.95 two.sided            wald
#> 5  0.8604651 0.7452701 0.9756601 37 43       0.95 two.sided         wald-cc
#> 6  0.8604651 0.7348110 0.9395927 37 43       0.95 two.sided        jeffreys
#> 7  0.8604651 0.7348110 0.9395927 37 43       0.95 two.sided    jeffreys-mod
#> 8  0.8604651 0.7206752 0.9470234 37 43       0.95 two.sided clopper-pearson
#> 9  0.8604651 0.7235600 0.9382469 37 43       0.95 two.sided   agresti-coull
#> 10 0.8604651 0.7661306 0.9472522 37 43       0.95 two.sided           pratt
#> 11 0.8604651 0.7346862 0.9424696 37 43       0.95 two.sided         arcsine
#> 12 0.8604651 0.7224337 0.9359412 37 43       0.95 two.sided           logit
#> 13 0.8604651 0.7532381 0.9301239 37 43       0.95 two.sided         witting
#> 14 0.8604651 0.7321815 0.9414281 37 43       0.95 two.sided           mid-p
#> 15 0.8604651 0.7255219 0.9374444 37 43       0.95 two.sided          blaker
#> 16 0.8604651 0.7372546 0.9420472 37 43       0.95 two.sided      likelihood
#> 17 0.8604651 0.7223441 0.9372304 37 43       0.95 two.sided        khouadji
#>    stdEst
#> 1    TRUE
#> 2    TRUE
#> 3    TRUE
#> 4    TRUE
#> 5    TRUE
#> 6    TRUE
#> 7    TRUE
#> 8    TRUE
#> 9    TRUE
#> 10   TRUE
#> 11   TRUE
#> 12   TRUE
#> 13   TRUE
#> 14   TRUE
#> 15   TRUE
#> 16   TRUE
#> 17   TRUE

prop.test(x=37, n=43, correct=FALSE) # same as method wilson
#> 
#>  1-sample proportions test without continuity correction
#> 
#> data:  37 out of 43, null probability 0.5
#> X-squared = 22.349, df = 1, p-value = 2.274e-06
#> alternative hypothesis: true p is not equal to 0.5
#> 95 percent confidence interval:
#>  0.7273641 0.9344428
#> sample estimates:
#>         p 
#> 0.8604651 
#> 
prop.test(x=37, n=43, correct=TRUE)  # same as method wilsoncc
#> 
#>  1-sample proportions test with continuity correction
#> 
#> data:  37 out of 43, null probability 0.5
#> X-squared = 20.93, df = 1, p-value = 4.763e-06
#> alternative hypothesis: true p is not equal to 0.5
#> 95 percent confidence interval:
#>  0.7137335 0.9419725
#> sample estimates:
#>         p 
#> 0.8604651 
#> 


# the confidence interval computed by binom.test
#   corresponds to the Clopper-Pearson interval
binomCI(x=42, n=43, method="clopper-pearson")
#>       est       lci       uci 
#> 0.9767442 0.8771095 0.9994114 
binom.test(x=42, n=43)$conf.int
#> [1] 0.8771095 0.9994114
#> attr(,"conf.level")
#> [1] 0.95


# all arguments are being recycled:
binomCI(x=c(42, 35, 23, 22), n=43, method="wilson")
#>         est       lci       uci  x  n conf.level     sides method stdEst
#> 1 0.9767442 0.8794101 0.9958829 42 43       0.95 two.sided wilson   TRUE
#> 2 0.8139535 0.6738300 0.9025825 35 43       0.95 two.sided wilson   TRUE
#> 3 0.5348837 0.3891564 0.6748894 23 43       0.95 two.sided wilson   TRUE
#> 4 0.5116279 0.3675231 0.6538255 22 43       0.95 two.sided wilson   TRUE
binomCI(x=c(42, 35, 23, 22), n=c(50, 60, 70, 80), method="jeffreys")
#>         est       lci       uci  x  n conf.level     sides   method stdEst
#> 1 0.8400000 0.7206737 0.9213325 42 50       0.95 two.sided jeffreys   TRUE
#> 2 0.5833333 0.4571040 0.7017365 35 60       0.95 two.sided jeffreys   TRUE
#> 3 0.3285714 0.2272016 0.4437899 23 70       0.95 two.sided jeffreys   TRUE
#> 4 0.2750000 0.1863875 0.3795587 22 80       0.95 two.sided jeffreys   TRUE

# example Table I in Newcombe (1998)
meths <- c("wald", "wald-cc", "wilson", "wilson-cc",
           "clopper-pearson","mid-p", "lik")
bedrock::setNamesX(cbind(round(cbind(
    binomCI(81, 263, m=meths)[, c("lci","uci")],
    binomCI(15, 148, m=meths)[,  c("lci","uci")],
    binomCI(0, 20, m=meths)[, c("lci","uci")],
    binomCI(1, 29, m=meths)[, c("lci","uci")]), 4)), 
  rownames=meths)
#>                    lci    uci    lci    uci lci    uci    lci    uci
#> wald            0.2522 0.3638 0.0527 0.1500   0 0.0000 0.0000 0.1009
#> wald-cc         0.2503 0.3657 0.0494 0.1534   0 0.0250 0.0000 0.1181
#> wilson          0.2553 0.3662 0.0624 0.1605   0 0.1611 0.0061 0.1718
#> wilson-cc       0.2535 0.3682 0.0598 0.1644   0 0.2005 0.0018 0.1963
#> clopper-pearson 0.2527 0.3676 0.0578 0.1617   0 0.1684 0.0009 0.1776
#> mid-p           0.2544 0.3658 0.0601 0.1581   0 0.1391 0.0017 0.1585
#> lik             0.2543 0.3655 0.0596 0.1567   0 0.0916 0.0020 0.1432

# returning p.tilde for agresti-coull ci
binomCI(x=81, n=263, meth="agresti-coull", stdEst = c(TRUE, FALSE))
#>         est       lci       uci  x   n conf.level     sides        method
#> 1 0.3079848 0.2552207 0.3662774 81 263       0.95 two.sided agresti-coull
#> 2 0.3107490 0.2552207 0.3662774 81 263       0.95 two.sided agresti-coull
#>   stdEst
#> 1   TRUE
#> 2  FALSE

# return all implemented methods
binomCI(4, 19, conf.level =0.95, 
        method = c(".all"))[, c("est","lci","uci","method")]
#>          est          lci       uci          method
#> 1  0.2105263 0.0850767663 0.4333428          wilson
#> 2  0.2105263 0.0697069249 0.4609797       wilson-cc
#> 3  0.2105263 0.0850767663 0.4333428      wilson-mod
#> 4  0.2105263 0.0272132947 0.3938393            wald
#> 5  0.2105263 0.0008975053 0.4201551         wald-cc
#> 6  0.2105263 0.0755318529 0.4262059        jeffreys
#> 7  0.2105263 0.0755318529 0.4262059    jeffreys-mod
#> 8  0.2105263 0.0605245377 0.4556531 clopper-pearson
#> 9  0.2105263 0.0795050565 0.4389145   agresti-coull
#> 10 0.2105263 0.0650967399 0.4556617           pratt
#> 11 0.2105263 0.0687042574 0.4296953         arcsine
#> 12 0.2105263 0.0813092940 0.4455116           logit
#> 13 0.2105263 0.0958686719 0.4072225         witting
#> 14 0.2105263 0.0707064177 0.4331720           mid-p
#> 15 0.2105263 0.0752938166 0.4448898          blaker
#> 16 0.2105263 0.0706525937 0.4235382      likelihood
#> 17 0.2105263 0.0814142100 0.4403249        khouadji


binomCIn(p=0.1, width=0.05, method="pratt")
#> [1] 586.9031
```
