# Barnard's Unconditional Test for Comparing Two Independent Proportions in a 2×2 Table

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
  useStoredCSM = FALSE,
  ...
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
  of `"two.sided"` (default), `"greater"` or `"less"`. You can specify
  just the initial letter.

- method:

  the method for finding the more extreme tables, one of `"csm"`
  (default), `"z-pooled"`, `"z-unpooled"`, `"boschloo"` or
  `"santner-snell"`. The CSM test cannot be calculated for multinomial
  models and is computationally the most demanding method (see the
  Details and the `useStoredCSM` argument).

- fixed:

  indicates which margin is fixed: `1` for rows (default), `2` for
  columns, or `NA` for none of both (multinomial model).

- useStoredCSM:

  logical, use a stored ordering matrix for the CSM test to greatly
  reduce the computation time (default is `FALSE`).

- ...:

  further arguments passed on to
  [`Exact::exact.test()`](https://rdrr.io/pkg/Exact/man/exact.test.html),
  e.g. `npNumbers` or `conf.int`.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of the test statistic used to order the tables.

- `parameter`:

  the sizes of the two samples.

- `p.value`:

  the p-value of the test.

- `estimate`:

  the observed difference in proportions.

- `null.value`:

  the difference in proportions under the null hypothesis.

- `alternative`:

  a character string describing the alternative hypothesis.

- `np`:

  the value of the nuisance parameter that maximizes the p-value.

- `np.range`:

  the range of nuisance parameters considered.

- `model`, `method`:

  character strings describing the sampling model and the method used to
  order the tables.

- `data.name`:

  a character string giving the name of the data.

## Details

There are two fundamentally different exact tests for comparing the
equality of two binomial probabilities - Fisher's exact test (Fisher,
1925), and Barnard's exact test (Barnard, 1945). Fisher's exact test
(Fisher, 1925) is the more popular of the two. In fact, Fisher was
bitterly critical of Barnard's proposal for esoteric reasons that we
will not go into here. For 2 x 2 tables, Barnard's test is more powerful
than Fisher's, as Barnard noted in his 1945 paper, much to Fisher's
chagrin. Anyway, perhaps due to its computational difficulty the
Barnard's is not widely used. (Mehta and Senchaudhuri, 2003)

Unconditional exact tests can be performed for binomial or multinomial
models. The binomial model assumes the row or column margins (but not
both) are known in advance, while the multinomial model assumes only the
total sample size is known beforehand. For the binomial model, the user
needs to specify which margin is fixed (default is rows). Conditional
tests (e.g., Fisher's exact test) have both row and column margins
fixed, but this is a very uncommon design. (See Calhoun (2019) for more
details.)

If `x` is a matrix, it is taken as a two-dimensional contingency table,
and hence its entries should be nonnegative integers. Otherwise, both
`x` and `y` must be vectors of the same length. Incomplete cases are
removed, the vectors are coerced into factor objects, and the
contingency table is computed from these.

For a 2x2 contingency table, such as \\X=\[n_1,n_2;n_3,n_4\]\\, the
normalized difference in proportions between the two categories, given
in each column, can be written with pooled variance (Score statistic) as
\$\$T(X)=\frac{\hat{p}\_2-\hat{p}\_1}{\sqrt{\hat{p}(1-\hat{p})(\frac{1}{c_1}+\frac{1}{c_2})}},\$\$
where \\\hat{p}=(n_1+n_3)/(n_1+n_2+n_3+n_4)\\,
\\\hat{p}\_2=n_2/(n_2+n_4)\\, \\\hat{p}\_1=n_1/(n_1+n_3)\\,
\\c_1=n_1+n_3\\ and \\c_2=n_2+n_4\\. Alternatively, with unpooled
variance (Wald statistic), the difference in proportions can be written
as
\$\$T(X)=\frac{\hat{p}\_2-\hat{p}\_1}{\sqrt{\frac{\hat{p}\_1(1-\hat{p}\_1)}{c_1}+\frac{\hat{p}\_2(1-\hat{p}\_2)}{c_2}}}.\$\$
The probability of observing \\X\\ is
\$\$P(X)=\frac{c_1!c_2!}{n_1!n_2!n_3!n_4!}p^{n_1+n_2}(1-p)^{n_3+n_4},\$\$
where \\p\\ is the unknown nuisance parameter.

Barnard's test considers all tables with category sizes \\c_1\\ and
\\c_2\\ for a given \\p\\. The p-value is the sum of probabilities of
the tables having a score in the rejection region, e.g. having
significantly large difference in proportions for a two-sided test. The
p-value of the test is the maximum p-value calculated over all \\p\\
between 0 and 1.

If `useStoredCSM` is set to `TRUE` a companion data package called
ExactData must be installed from GitHub.

The author states: *"The CSM test is computationally intensive due to
iteratively maximizing the p-value calculation to order the tables. The
CSM ordering matrix has been stored for all possible sample sizes less
than or equal to 100 (i.e., max(n1,n2)\<=100). Thus, using the
useStoredCSM = TRUE can greatly improve computation time. However, the
stored ordering matrix was computed with npNumbers=100 and it is
possible that the ordering matrix was not optimal for larger npNumbers.
Increasing npNumbers and setting useStoredCSM = FALSE ensures the
p-value is correctly calculated at the expense of significantly greater
computation time. The stored ordering matrix is not used in the
calculation of confidence intervals or non-inferiority tests, so CSM can
still be very computationally intensive."*

## Note

`barnardTest()` is an interface to
[`Exact::exact.test()`](https://rdrr.io/pkg/Exact/man/exact.test.html)
by Peter Calhoun; the Exact package (available on CRAN) must be
installed.

## References

Barnard, G.A. (1945) A new test for 2x2 tables. *Nature*, 156:177.

Barnard, G.A. (1947) Significance tests for 2x2 tables. *Biometrika*,
34:123-138.

Suissa, S. and Shuster, J. J. (1985) Exact unconditional sample sizes
for the 2x2 binomial trial, *Journal of the Royal Statistical Society*,
Ser. A, 148, 317-327.

Lin C.Y., Yang M.C. (2009) Improved p-value tests for comparing two
independent binomial proportions. *Communications in
Statistics-Simulation and Computation*, 38(1):78-91.

Mehta, C.R., Senchaudhuri, P. (2003) Conditional versus unconditional
exact tests for comparing two binomials. Cytel Software Corporation,
technical report.

Calhoun, P. (2019) Exact: Unconditional Exact Test. R package version
2.0.

## See also

[`fisher.test`](https://rdrr.io/r/stats/fisher.test.html)

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
#>  Z-pooled Exact Test
#> 
#> data:  tab
#> test statistic = 0.43944, first sample size = 9, second sample size =
#> 17, p-value = 0.7858
#> alternative hypothesis: true difference in proportion is not equal to 0
#> sample estimates:
#> difference in proportion 
#>               0.06535948 
#> 

# \donttest{
# the default CSM method is more powerful,
# but computationally expensive
barnardTest(tab)
#> 
#>  CSM Exact Test
#> 
#> data:  tab
#> test statistic = NA, first sample size = 9, second sample size = 17,
#> p-value = 0.7181
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
#>  CSM Exact Test
#> 
#> data:  tab
#> test statistic = NA, first sample size = 15, second sample size = 15,
#> p-value = 0.03669
#> alternative hypothesis: true difference in proportion is less than 0
#> sample estimates:
#> difference in proportion 
#>               -0.3333333 
#> 
# }
```
