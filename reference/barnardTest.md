# Barnard's Unconditional Test for Comparing Two Independent Proportions in a 2 x 2 Table

An exact unconditional test for \\2 \times 2\\ contingency tables,
offering a more powerful alternative to Fisher's exact test by avoiding
conditioning on both marginal totals.

## Usage

``` r
barnardTest(
  x,
  y = NULL,
  alternative = c("two.sided", "less", "greater"),
  method = c("csm", "z-pooled", "z-unpooled", "boschloo", "santner-snell"),
  fixed = 1,
  tsMethod = c("square", "central"),
  delta = 0,
  conf.level = NA,
  npNumbers = 100,
  npInterval = FALSE,
  beta = 0.001
)
```

## Arguments

- x:

  a numeric vector or a two-dimensional contingency table in matrix
  form. `x` and `y` can also both be factors.

- y:

  a factor object; ignored if `x` is a matrix.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"greater"` or `"less"`, referring to
  \\\pi_1 - \pi_2\\ relative to `delta`. You can specify just the
  initial letter.

- method:

  the ordering of the tables, one of `"csm"` (default), `"z-pooled"`,
  `"z-unpooled"`, `"boschloo"` or `"santner-snell"`. See Details.

- fixed:

  indicates which margin is fixed: `1` for rows (default), `2` for
  columns, or `NA` for none of both (multinomial model).

- tsMethod:

  the two-sided method, `"square"` (default) or `"central"`. See
  Details.

- delta:

  the difference \\\pi_1 - \pi_2\\ under the null hypothesis, strictly
  between -1 and 1; non-zero values give tests of non-inferiority or
  superiority. Binomial model only.

- conf.level:

  confidence level of the interval for \\\pi_1 - \pi_2\\ (binomial model
  only). The default `NA` computes no interval.

- npNumbers:

  the number of grid points for the nuisance parameter, at least 2,
  default 100.

- npInterval:

  logical, restrict the nuisance parameter to a confidence interval
  (Berger and Boos, 1994), default `FALSE`.

- beta:

  the error probability of that interval, default 0.001; only used with
  `npInterval = TRUE`.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the ordering statistic of the observed table: `z` (z-pooled,
  z-unpooled), `d` (Santner-Snell, the difference minus `delta`) or
  `p.Fisher` (Boschloo); absent for CSM, whose ordering has no
  statistic.

- `parameter`:

  the sizes `n1`, `n2` of the two samples (binomial model) or the total
  sample size `N` (multinomial model).

- `p.value`:

  the p-value of the test.

- `conf.int`:

  the confidence interval, if `conf.level` is not `NA`.

- `estimate`:

  the observed difference in proportions (binomial) or in products of
  proportions (multinomial).

- `null.value`:

  the difference under the null hypothesis.

- `alternative`:

  a character string describing the alternative hypothesis.

- `np`:

  the value(s) of the nuisance parameter(s) that maximize the p-value.

- `np.range`:

  the range of nuisance parameters considered.

- `model`, `tsmethod`:

  the sampling model and the two-sided method.

- `method`:

  a character string describing the test.

- `data.name`:

  a character string giving the name of the data.

## Details

Two fundamentally different exact tests exist for the equality of two
binomial probabilities: Fisher's exact test (Fisher, 1925), which
conditions on both margins, and Barnard's unconditional test (Barnard,
1945), which does not. Fisher's test is by far the more widely used, and
Fisher himself objected strongly to Barnard's proposal. For 2 x 2
tables, however, Barnard's test is generally more powerful, as Barnard
already pointed out; its limited use is mainly due to its computational
cost (Mehta and Senchaudhuri, 2003).

**Choice of test.** In most studies only one margin is fixed by design,
e.g. the group sizes of a randomised trial, while the number of events
is random. Conditioning on both margins, as Fisher's exact test does,
then makes the test markedly conservative: because of the discreteness
of the conditional distribution, the actual level often lies far below
the nominal one, at the cost of power, in particular for the small
samples for which exact tests are used. Unconditional tests avoid this.
Boschloo's test, which orders the tables by Fisher's p-value, is
uniformly more powerful than Fisher's test: its p-value never exceeds
Fisher's, while it still keeps the nominal level (Boschloo, 1970).
Lydersen et al. (2009) therefore recommend unconditional tests (or
Fisher's mid-p) instead of the traditional Fisher test, and the score
statistic (`"z-pooled"`) is to be preferred to the Wald statistic
(`"z-unpooled"`), which can behave poorly (Mehrotra et al., 2003).

For two independent proportions with fixed group sizes, `barnardTest()`
is thus the recommended choice; CSM and Boschloo are the most powerful
orderings, z-pooled is the most widely used one.
[`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html) remains
appropriate for r x c tables, for exact conditional inference on the
odds ratio, when both margins are fixed by design, and whenever
conditional inference is wanted on principle: whether the event margin
should be conditioned on is a long-standing debate (Yates, 1984; Little,
1989).

**Sampling model.** Unconditional exact tests can be performed for
binomial or multinomial models. The binomial model assumes that the row
or the column margins (but not both) are fixed by design, the
multinomial model only the total sample size. `fixed` selects the model:
`1` (default) for fixed rows, `2` for fixed columns, `NA` for the
multinomial model. Conditioning on both margins leads to Fisher's exact
test, see [`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html).

**Binomial model.** The fixed margin defines the two groups, of sizes
\\n_1\\ and \\n_2\\; the first category of the other classification
counts the successes \\x_1\\ and \\x_2\\. The tested difference is
\\\pi_1 - \pi_2\\ with the estimate \\\hat{p}\_1-\hat{p}\_2\\,
\\\hat{p}\_i = x_i/n_i\\, against the null value `delta`. Under the null
hypothesis the probability of a table is \$\$P(x_1, x_2) =
\binom{n_1}{x_1}(\pi+\delta)^{x_1}(1-\pi-\delta)^{n_1-x_1}
\binom{n_2}{x_2}\pi^{x_2}(1-\pi)^{n_2-x_2},\$\$ where \\\pi\\ is the
unknown nuisance parameter. The tables are ordered by how extreme they
are, the p-value for a given \\\pi\\ is the probability of all tables at
least as extreme as the observed one, and the p-value of the test is its
maximum over \\\pi\\. The maximum is searched on a grid of `npNumbers`
equally spaced values over the admissible range \\\[\max(0, -\delta),
\min(1, 1-\delta)\]\\, bounds included, and then refined locally with
[`optimize()`](https://rdrr.io/r/stats/optimize.html).

`method` defines the ordering:

- `"csm"`:

  Barnard's original convexity, symmetry and minimum ordering (CSM;
  Barnard, 1947). Starting from the most extreme table, the table that
  increases the maximal p-value least is added to the rejection region,
  step by step. Generally the most powerful choice; binomial model only.

- `"z-pooled"`:

  the score statistic \$\$Z =
  \frac{\hat{p}\_1-\hat{p}\_2-\delta}{\sqrt{\tilde{p}\_1(1-\tilde{p}\_1)/n_1 +
  \tilde{p}\_2(1-\tilde{p}\_2)/n_2}},\$\$ with the pooled proportion for
  \\\tilde{p}\_1 = \tilde{p}\_2\\ if \\\delta = 0\\, and the restricted
  maximum likelihood estimates (Miettinen and Nurminen, 1985) otherwise
  (Suissa and Shuster, 1985).

- `"z-unpooled"`:

  the Wald statistic, with the unrestricted estimates \\\hat{p}\_i\\ in
  the variance (and the same numerator
  \\\hat{p}\_1-\hat{p}\_2-\delta\\).

- `"boschloo"`:

  the p-value of Fisher's exact test (Boschloo, 1970).

- `"santner-snell"`:

  the plain difference \\\hat{p}\_1 - \hat{p}\_2\\ (Santner and Snell,
  1980).

Two-sided tests are by default based on the absolute value of the
statistic, i.e. on a region that is symmetric about the null value
(`tsMethod = "square"`). With `tsMethod = "central"` the two-sided
p-value is twice the smaller one-sided p-value.

With `npInterval = TRUE` the maximum is only searched within a
\\100(1-\beta)\\% Clopper-Pearson interval of the nuisance parameter and
\\\beta\\ is added to the p-value (Berger and Boos, 1994). The interval
is based on the pooled successes, which are binomial only for \\\delta =
0\\; it is therefore not available with a non-zero `delta`, a confidence
interval or CSM.

The confidence interval for \\\pi_1 - \pi_2\\ (binomial model only) is
computed if `conf.level` is not `NA`, and is obtained by inverting the
test over `delta`. It is computationally demanding, in particular for
CSM, where the ordering has to be rebuilt for every candidate value.

**Multinomial model.** Only the total \\N\\ is fixed. The nuisance
parameters are the probabilities of the first row and of the first
column, and the tested quantity is the difference in the products of the
cell proportions, \\(x\_{11}x\_{22} - x\_{12}x\_{21})/N^2\\. Available
for the z-pooled, z-unpooled and Boschloo orderings, with `delta = 0`
and without confidence interval.

If `x` is a matrix, it is taken as a two-dimensional contingency table,
and hence its entries should be nonnegative integers. Otherwise, both
`x` and `y` must be vectors of the same length. Incomplete cases are
removed, the vectors are coerced into factor objects, and the
contingency table is computed from these.

## Note

The algorithms follow the Exact package by Peter Calhoun
([`Exact::exact.test()`](https://rdrr.io/pkg/Exact/man/exact.test.html)),
with these deliberate differences:

- the rejection region of the ordering statistics consists of exactly
  the tables that are at least as extreme as the observed one. Exact
  builds it by walking the boundary and filling rows and columns by
  convexity, which in two-sided tests can admit tables that are less
  extreme (e.g. `matrix(c(3, 2, 1, 0), 2)`, z-pooled: 0.783 instead of
  0.568);

- the nuisance grid includes the bounds of the admissible range (Exact
  uses 0.00001 and 0.99999), where the maximum can lie;

- the z-unpooled statistic subtracts `delta` in the numerator, as the
  z-pooled one does;

- the stored CSM orderings of the companion package ExactData are not
  used: the compiled ordering is fast enough to be computed on the fly,
  for every `npNumbers`;

- when the grid maximum is attained at several nuisance values, the
  refined p-value is the largest of the refined maxima (Exact takes the
  smallest).

## References

Barnard, G.A. (1945) A new test for 2x2 tables. *Nature*, 156:177.

Barnard, G.A. (1947) Significance tests for 2x2 tables. *Biometrika*,
34:123-138.

Berger, R.L. and Boos, D.D. (1994) P values maximised over a confidence
set for the nuisance parameter. *Journal of the American Statistical
Association*, 89:1012-1016.

Boschloo, R.D. (1970) Raised conditional level of significance for the
2x2-table when testing the equality of two probabilities. *Statistica
Neerlandica*, 24:1-35.

Calhoun, P. (2024) Exact: Unconditional Exact Test. R package version
3.3.

Little, R.J.A. (1989) Testing the equality of two independent binomial
proportions. *The American Statistician*, 43:283-288.

Lydersen, S., Fagerland, M.W. and Laake, P. (2009) Recommended tests for
association in 2x2 tables. *Statistics in Medicine*, 28:1159-1175.

Mehrotra, D.V., Chan, I.S.F. and Berger, R.L. (2003) A cautionary note
on exact unconditional inference for a difference between two
independent binomial proportions. *Biometrics*, 59:441-450.

Mehta, C.R., Senchaudhuri, P. (2003) Conditional versus unconditional
exact tests for comparing two binomials. Cytel Software Corporation,
technical report.

Miettinen, O. and Nurminen, M. (1985) Comparative analysis of two rates.
*Statistics in Medicine*, 4:213-226.

Santner, T.J. and Snell, M.K. (1980) Small-sample confidence intervals
for p1 - p2 and p1/p2 in 2x2 contingency tables. *Journal of the
American Statistical Association*, 75:386-394.

Suissa, S. and Shuster, J. J. (1985) Exact unconditional sample sizes
for the 2x2 binomial trial, *Journal of the Royal Statistical Society*,
Ser. A, 148, 317-327.

Yates, F. (1984) Tests of significance for 2 x 2 contingency tables.
*Journal of the Royal Statistical Society*, Ser. A, 147:426-463.

## See also

[`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html),
[`prop.test()`](https://rdrr.io/r/stats/prop.test.html)

Other test.categorical: [`bhapkarTest()`](bhapkarTest.md),
[`breslowDayTest()`](breslowDayTest.md),
[`cochranQTest()`](cochranQTest.md), [`gTest()`](gTest.md),
[`lehmacherTest()`](lehmacherTest.md),
[`mantelTrendTest()`](mantelTrendTest.md),
[`stuartMaxwellTest()`](stuartMaxwellTest.md),
[`woolfTest()`](woolfTest.md)

## Examples

``` r
tab <- as.table(matrix(c(8, 14, 1, 3), nrow=2,
                dimnames=list(treat=c("I","II"), out=c("I","II"))))
barnardTest(tab, method="z-pooled")
#> 
#>  Barnard's unconditional exact test (z-pooled ordering)
#> 
#> data:  tab
#> z = 0.43944, n1 = 9, n2 = 17, p-value = 0.7858
#> alternative hypothesis: true difference in proportion is not equal to 0
#> sample estimates:
#> difference in proportion 
#>               0.06535948 
#> 

# the default CSM ordering
barnardTest(tab)
#> 
#>  Barnard's unconditional exact test (CSM ordering)
#> 
#> data:  tab
#> n1 = 9, n2 = 17, p-value = 0.7181
#> alternative hypothesis: true difference in proportion is not equal to 0
#> sample estimates:
#> difference in proportion 
#>               0.06535948 
#> 

# Mehta and Senchaudhuri (2003), vaccine example
tab <- as.table(matrix(c(7, 12, 8, 3), nrow=2,
                       dimnames=list(treat=c("vaccine","placebo"),
                                     infection=c("yes","no"))))
barnardTest(tab, alternative="less")
#> 
#>  Barnard's unconditional exact test (CSM ordering)
#> 
#> data:  tab
#> n1 = 15, n2 = 15, p-value = 0.03669
#> alternative hypothesis: true difference in proportion is less than 0
#> sample estimates:
#> difference in proportion 
#>               -0.3333333 
#> 

# non-inferiority margin of -0.1, with confidence interval
barnardTest(tab, method="z-pooled", alternative="greater", delta=-0.1,
            conf.level=0.95)
#> 
#>  Barnard's unconditional exact test (z-pooled ordering)
#> 
#> data:  tab
#> z = -1.339, n1 = 15, n2 = 15, p-value = 0.9444
#> alternative hypothesis: true difference in proportion is greater than -0.1
#> 95 percent confidence interval:
#>  -0.5963737  1.0000000
#> sample estimates:
#> difference in proportion 
#>               -0.3333333 
#> 

# multinomial model: only the total is fixed
barnardTest(tab, method="boschloo", fixed=NA)
#> 
#>  Barnard's unconditional exact test (Boschloo ordering)
#> 
#> data:  tab
#> p.Fisher = 0.12814, N = 30, p-value = 0.08367
#> alternative hypothesis: true difference in product of proportion is not equal to 0
#> sample estimates:
#> difference in product of proportion 
#>                         -0.08333333 
#> 
```
