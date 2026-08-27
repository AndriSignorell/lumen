# Cochran's Q Test for Comparing Matched Binary Responses Across Multiple Conditions

A nonparametric test for dependent samples with dichotomous data,
assessing whether proportions differ across multiple conditions or time
points. It generalizes the McNemar test to more than two groups.

## Usage

``` r
cochranQTest(y, ...)

# Default S3 method
cochranQTest(
  y,
  groups,
  blocks,
  method = c("asymptotic", "approximate"),
  nresample = 10000,
  na.action = na.omit,
  ...
)

# S3 method for class 'formula'
cochranQTest(
  formula,
  data,
  subset,
  na.action = na.pass,
  method = c("asymptotic", "approximate"),
  nresample = 10000,
  ...
)
```

## Arguments

- y:

  either a numeric vector of data values, or a matrix with the blocks in
  the rows and the groups in the columns.

- ...:

  further arguments to be passed to or from methods.

- groups:

  a vector giving the group for the corresponding elements of `y` if
  this is a vector; ignored if `y` is a matrix. If not a factor object,
  it is coerced to one.

- blocks:

  a vector giving the block for the corresponding elements of `y` if
  this is a vector; ignored if `y` is a matrix. If not a factor object,
  it is coerced to one.

- method:

  a character string specifying how the p-value is computed, one of
  `"asymptotic"` (default, chi-squared approximation) or `"approximate"`
  (Monte Carlo permutation via the coin package).

- nresample:

  the number of Monte Carlo replicates used for `method = "approximate"`
  (default is `1e4`).

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

- formula:

  a formula of the form `y ~ groups | blocks`.

- data:

  an optional matrix or data frame (or similar: see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in the formula. By default the variables are taken from
  `environment(formula)`.

- subset:

  an optional vector specifying a subset of observations to be used.

## Value

A list with class `"htest"` containing the following components:

- `statistic`:

  the value of Cochran's chi-squared statistic.

- `parameter`:

  the degrees of freedom of the approximate chi-squared distribution of
  the test statistic (asymptotic method only).

- `p.value`:

  the p-value of the test.

- `method`:

  a character string indicating the test performed and the method used
  to compute the p-value.

- `data.name`:

  a character string giving the names of the data.

## Details

Performs Cochran's Q test for related samples with a binary response.
The test is appropriate for unreplicated complete block designs (i.e.,
matched or paired data), where each block contains exactly one
observation for each group.

Cochran's Q test is a nonparametric method for testing whether the
proportions of a binary outcome differ across multiple related groups,
while accounting for block effects. It can be regarded as an extension
of McNemar's test to more than two groups. The null hypothesis is that,
apart from block effects, the probability of a "success" (coded as 1) is
the same in all groups.

If `y` is a matrix, groups and blocks are inferred from the columns and
rows, respectively. Missing values in `y` lead to the removal of the
corresponding blocks (for both methods). Missing values are not allowed
in `groups` or `blocks`. If all blocks give identical responses, the
statistic is 0 by convention and the p-value 1.

With `method = "asymptotic"` the test statistic is referred to its
approximate chi-squared distribution. With `method = "approximate"` a
Monte Carlo permutation p-value is obtained via
[`coin::symmetry_test()`](https://rdrr.io/pkg/coin/man/SymmetryTest.html),
whose quadratic-form statistic coincides with Cochran's Q for this
design.

Cochran's Q test is closely related to the Friedman test, but is
specifically designed for binary (0/1) responses.

## References

Cochran, W.G. (1950) The Comparison of Percentages in Matched Samples.
*Biometrika*, 37 (3/4): 256-266.

## See also

[`friedman.test()`](https://rdrr.io/r/stats/friedman.test.html)

Other test.categorical: [`barnardTest()`](barnardTest.md),
[`bhapkarTest()`](bhapkarTest.md),
[`breslowDayTest()`](breslowDayTest.md), [`gTest()`](gTest.md),
[`lehmacherTest()`](lehmacherTest.md),
[`mantelTrendTest()`](mantelTrendTest.md),
[`stuartMaxwellTest()`](stuartMaxwellTest.md),
[`woolfTest()`](woolfTest.md)

## Examples

``` r
# example in:
# http://support.sas.com/documentation/cdl/en/statugfreq/63124/PDF/default/statugfreq.pdf
# pp. S. 1824

# create the dataset
d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
            rep(1:8, c(6,2,2,6,16,4,4,6)), ]
row.names(d.frm) <- NULL

# rearrange to long shape
d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[c(1:3)],
                  v.names="resp", direction="long")

# after having done the hard work of data organisation,
# performing the test is a piece of cake....
cochranQTest(resp ~ time | id, data=d.long)
#> Warning: coercing factor to binary (0/1): using 'U' as '1'
#> 
#>  Cochran's Q test (asymptotic)
#> 
#> data:  resp ~ time | id
#> Cochran's Q = 8.4706, df = 2, p-value = 0.01448
#> 

# and let's perform a post hoc analysis using mcnemar's test
z <- split(d.long, f=d.long$time)
pairwise.table(function(i, j) {
    mcnemar.test(z[[i]]$resp, z[[j]]$resp, correct=FALSE)$p.value
  },
  level.names = names(z),
  p.adjust.method = "fdr"
)
#>           A         B
#> B 1.0000000        NA
#> C 0.0350133 0.0350133
```
