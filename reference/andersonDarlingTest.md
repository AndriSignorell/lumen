# Anderson-Darling Test for Assessing Distributional Fit, Especially in the Tails

A goodness-of-fit test assessing whether a sample follows a specified
distribution. Compared to the Kolmogorov-Smirnov test, it places greater
weight on discrepancies in the tails of the distribution.

## Usage

``` r
andersonDarlingTest(x, null = "punif", ..., estimated = FALSE, nullname)
```

## Arguments

- x:

  numeric vector of data values.

- null:

  a function, or a character string giving the name of a function, to
  compute the cumulative distribution function for the null
  distribution.

- ...:

  additional arguments for the cumulative distribution function.

- estimated:

  logical value indicating whether the parameters of the distribution
  were estimated using the data `x` (composite null hypothesis), or were
  fixed in advance (simple null hypothesis, the default).

- nullname:

  optional character string describing the null distribution. By default
  the name is derived from `null`, e.g. `"uniform distribution"` for the
  default `null="punif"`.

## Value

An object of class `"htest"` representing the result of the hypothesis
test.

## Details

This command performs the Anderson-Darling test of goodness-of-fit to
the distribution specified by the argument `null`. It is assumed that
the values in `x` are independent and identically distributed random
values, with some cumulative distribution function \\F\\. The null
hypothesis is that \\F\\ is the function specified by the argument
`null`, while the alternative hypothesis is that \\F\\ is some other
function.

By default, the test assumes that all the parameters of the null
distribution are known in advance (a *simple* null hypothesis). This
test does not account for the effect of estimating the parameters.

If the parameters of the distribution were estimated (that is, if they
were calculated from the same data `x`), then this should be indicated
by setting the argument `estimated=TRUE`. The test will then use the
method of Braun (1980) to adjust for the effect of parameter estimation.

Note that Braun's method involves randomly dividing the data into \\m
\approx \sqrt{n}\\ groups, so the \\p\\-value is not exactly the same if
the test is repeated. This technique is expected to work well when the
number of observations in `x` is large. If there are too few
observations for the adjustment (\\n \le 4\\), the unadjusted test is
performed and a warning is issued.

Missing values are silently removed.

## Note

Original C code by George Marsaglia and John Marsaglia; R interface by
Adrian Baddeley, previously published in the goftest package. Rewritten
in C++ with an adapted R interface to conform to package standards.

## References

Anderson, T.W. and Darling, D.A. (1952) Asymptotic theory of certain
'goodness-of-fit' criteria based on stochastic processes. *Annals of
Mathematical Statistics* **23**, 193–212.

Anderson, T.W. and Darling, D.A. (1954) A test of goodness of fit.
*Journal of the American Statistical Association* **49**, 765–769.

Braun, H. (1980) A simple method for testing goodness-of-fit in the
presence of nuisance parameters. *Journal of the Royal Statistical
Society, Series B* **42**, 53–63.

Marsaglia, G. and Marsaglia, J. (2004) Evaluating the Anderson-Darling
distribution. *Journal of Statistical Software* **9** (2), 1–5.
[doi:10.18637/jss.v009.i02](https://doi.org/10.18637/jss.v009.i02)

## See also

[`pAD`](https://andrisignorell.github.io/lumen/reference/pAD.md) for the
null distribution of the test statistic

Other test.normality:
[`cramerVonMisesTest()`](https://andrisignorell.github.io/lumen/reference/cramerVonMisesTest.md),
[`jarqueBeraTest()`](https://andrisignorell.github.io/lumen/reference/jarqueBeraTest.md),
[`lillieTest()`](https://andrisignorell.github.io/lumen/reference/lillieTest.md),
[`pearsonTest()`](https://andrisignorell.github.io/lumen/reference/pearsonTest.md),
[`shapiroFranciaTest()`](https://andrisignorell.github.io/lumen/reference/shapiroFranciaTest.md)

## Examples

``` r
x <- rnorm(10, mean=2, sd=1)
andersonDarlingTest(x, "pnorm", mean=2, sd=1)
#> 
#>  Anderson-Darling test of goodness-of-fit; null hypothesis: Normal
#>  distribution; with parameters mean = 2, sd = 1; parameters fixed
#> 
#> data:  x
#> An = 0.95102, p-value = 0.3815
#> 
andersonDarlingTest(x, "pnorm", mean=mean(x), sd=sd(x), estimated=TRUE)
#> 
#>  Anderson-Darling test of goodness-of-fit; Braun's adjustment using 3
#>  groups; null hypothesis: Normal distribution; with parameters mean =
#>  2.24125265128704, sd = 0.603331309611038; parameters estimated from
#>  data
#> 
#> data:  x
#> Anmax = 2.2295, p-value = 0.204
#> 
```
