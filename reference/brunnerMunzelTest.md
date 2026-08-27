# Brunner-Munzel Test for Comparing Stochastic Dominance Between Two Independent Samples

Tests the nonparametric Behrens-Fisher hypothesis that a randomly drawn
observation from one sample is as likely to be smaller as it is to be
larger than a randomly drawn observation from the other, without
assuming equal shapes or equal variances in the two groups.

## Usage

``` r
brunnerMunzelTest(x, ...)

# Default S3 method
brunnerMunzelTest(
  x,
  y,
  p0 = 0.5,
  alternative = c("two.sided", "less", "greater"),
  conf.level = 0.95,
  method = c("t", "permutation", "normal"),
  exact = NULL,
  nPerm = 10000L,
  ...
)

# S3 method for class 'formula'
brunnerMunzelTest(formula, data, subset, na.action = na.pass, ...)
```

## Arguments

- x:

  a numeric vector or ordered factor

- ...:

  further arguments, passed to the default method

- y:

  a numeric vector or ordered factor

- p0:

  the value of \\p\\ under the null hypothesis; must be `0.5` for
  `method = "permutation"`

- alternative:

  a character string specifying the alternative hypothesis for \\p\\,
  one of `"two.sided"` (default), `"less"` or `"greater"`

- conf.level:

  confidence level of the interval

- method:

  the inference method, one of `"t"` (default, Satterthwaite t
  approximation), `"permutation"` (studentized permutation test) or
  `"normal"` (asymptotic normal approximation)

- exact:

  logical, whether to enumerate all splits instead of sampling them;
  `NULL` (default) decides by the number of splits. Ignored unless
  `method = "permutation"`

- nPerm:

  number of Monte-Carlo resamples used when the permutation distribution
  is not enumerated

- formula:

  a formula of the form `lhs ~ rhs` where `lhs` is numeric and `rhs` a
  factor with two levels

- data:

  an optional data frame containing the model variables

- subset:

  an optional vector specifying a subset of observations

- na.action:

  a function indicating what should happen when the data contain `NA`s

## Value

An object of class `"htest"` with components

- statistic:

  the studentized Brunner-Munzel statistic

- parameter:

  the Satterthwaite degrees of freedom, `NULL` for
  `method = "permutation"`

- p.value:

  the p-value

- conf.int:

  confidence interval for \\p\\

- estimate:

  the estimated relative effect \\\hat{p}\\

- null.value:

  the value of \\p\\ under the null hypothesis

- stderr:

  the standard error of \\\hat{p}\\

- alternative:

  the alternative hypothesis

- method:

  a character string describing the test

- data.name:

  a character string giving the names of the data

## Details

The estimated quantity is the relative effect \$\$p = P(X \< Y) +
\tfrac{1}{2} P(X = Y),\$\$ the probability that a random draw from `y`
exceeds a random draw from `x`, ties counted half. It is identical to
the Mann-Whitney statistic scaled to \\\[0, 1\]\\, i.e. \\U / (n_1
n_2)\\, and is often reported as the common language effect size.

The Wilcoxon-Mann-Whitney test is a valid test of \\p = 1/2\\ only under
the hypothesis that both samples come from the same distribution. When
the two distributions differ in shape or spread, its null distribution
is wrong and the test does not keep its level. Brunner and Munzel (2000)
studentize the rank statistic with a separate variance estimate per
group and refer it to a t distribution with Satterthwaite degrees of
freedom, which is the rank analogue of the Welch correction.
[`yuenTTest`](yuenTTest.md) plays the same role among the parametric
location tests.

**Direction of `alternative`.** The alternative is stated in terms of
\\p\\, not in terms of `x` against `y`. `"greater"` therefore means \\p
\> p_0\\, that is, `y` tends to produce the larger values. This follows
the published definition of the statistic and the reported estimate, but
it is the reverse of [`t.test`](https://rdrr.io/r/stats/t.test.html) and
[`wilcox.test`](https://rdrr.io/r/stats/wilcox.test.html), where
`"greater"` refers to `x`.

**Choice of `method`.** The t approximation is liberal in small samples;
below roughly ten observations per group the studentized permutation
test of Neubert and Brunner (2007) should be preferred. Permuting the
studentized statistic rather than the raw rank statistic is what keeps
the permutation test asymptotically valid when the distributions differ,
so this is not an ordinary permutation test on ranks. It targets \\H_0:
p = 1/2\\: it is exact in finite samples under exchangeability of the
group labels, and remains asymptotically valid under the weaker
nonparametric Behrens-Fisher null \\p = 1/2\\ alone. Note that
`exact = TRUE` means exact enumeration of the permutation distribution,
which is not the same as a finite-sample exact test under the general
null. Only \\p_0 = 1/2\\ is implemented; other values of `p0` require an
approximate method.

The two-sided permutation p-value is the proportion of splits with
\\\|T^\*\| \ge \|T\|\\, which is the convention used by the
brunnermunzel package implementation based on Neubert and Brunner
(2007). It is exact under exchangeability, but it does not in general
equal twice the smaller one-sided tail, because the permutation
distribution need not be symmetric when the group sizes differ or ties
are present.

With `exact = NULL` all \\\binom{n_1 + n_2}{n_1}\\ splits are enumerated
when there are at most `1e6` of them, and `nPerm` Monte-Carlo resamples
are drawn otherwise. Monte-Carlo p-values use the \\(1 + k) / (1 + B)\\
correction and are therefore never zero.

**Confidence interval.** The interval is the studentized Wald interval
for \\p\\ and is reported for every `method`, including the permutation
methods, which supply a p-value but no interval of their own. It is
clipped to \\\[0, 1\]\\, and for a one-sided `alternative` the open end
is reported at the range limit rather than as \\\pm\infty\\.

**Non-overlapping samples.** If every observation in one sample is
smaller than every observation in the other, both variance components
are zero and neither approximation is defined. The permutation test is
used instead, with a warning. For tie-free data small enough to
enumerate, this returns the smallest attainable two-sided p-value, \\2 /
\binom{n_1 + n_2}{n_1}\\; with ties the two mirrored fully separated
splits need not both exist, so the attainable minimum can be smaller. No
Wald interval exists in that case and `conf.int` is `NA`: a zero
standard error reflects an empty variance estimate, not certainty about
\\p\\.

The formula method accepts `lhs ~ rhs` with exactly two groups on the
right-hand side; more than two levels are rejected by
[bedrock::resolveFormula](https://andrisignorell.github.io/bedrock/reference/resolveFormula.html).

Missing values are removed.

## References

Brunner, E., Munzel, U. (2000). The nonparametric Behrens-Fisher
problem: asymptotic theory and a small-sample approximation.
*Biometrical Journal*, **42**(1), 17-25.

Neubert, K., Brunner, E. (2007). A studentized permutation test for the
non-parametric Behrens-Fisher problem. *Computational Statistics and
Data Analysis*, **51**(10), 5192-5204.

## See also

[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)

Other test.location: [`hotellingsT2Test()`](hotellingsT2Test.md),
[`moodMedianTest()`](moodMedianTest.md), [`signTest()`](signTest.md),
[`tTestA()`](tTestA.md), [`vanWaerdenTest()`](vanWaerdenTest.md),
[`yuenTTest()`](yuenTTest.md), [`zTest()`](zTest.md)

## Examples

``` r
## Brunner & Munzel (2000), pain scores
x <- c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1)
y <- c(3, 3, 4, 3, 1, 2, 3, 1, 1, 5, 4)

brunnerMunzelTest(x, y)
#> 
#>  Brunner-Munzel test (t approximation)
#> 
#> data:  x and y
#> W = 3.1375, df = 17.683, p-value = 0.005786
#> alternative hypothesis: true P(X < Y) + 0.5 * P(X = Y) is not equal to 0.5
#> 95 percent confidence interval:
#>  0.5952169 0.9827052
#> sample estimates:
#> P(X < Y) + 0.5 * P(X = Y) 
#>                  0.788961 
#> 

## the estimate is the common language effect size
brunnerMunzelTest(x, y)$estimate
#> P(X < Y) + 0.5 * P(X = Y) 
#>                  0.788961 
## [1] 0.7889610

## small groups: prefer the studentized permutation test, which is
## enumerated exactly whenever the number of splits allows it
a <- c(1, 1, 2, 3, 3)
b <- c(4, 2, 5, 1, 6, 3, 7)
brunnerMunzelTest(a, b, method = "permutation")$p.value
#> [1] 0.09469697
## [1] 0.09469697

## formula interface
d <- data.frame(score = c(x, y),
                grp = rep(c("a", "b"), c(length(x), length(y))))
brunnerMunzelTest(score ~ grp, data = d)
#> 
#>  Brunner-Munzel test (t approximation)
#> 
#> data:  score ~ grp
#> W = 3.1375, df = 17.683, p-value = 0.005786
#> alternative hypothesis: true P(X < Y) + 0.5 * P(X = Y) is not equal to 0.5
#> 95 percent confidence interval:
#>  0.5952169 0.9827052
#> sample estimates:
#> P(X < Y) + 0.5 * P(X = Y) 
#>                  0.788961 
#> 
```
