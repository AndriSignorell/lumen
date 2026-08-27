# (Exact) Page Test for Detecting Ordered Treatment Effects in Blocked Designs

A nonparametric test for ordered alternatives in a two-way layout with
one observation per cell, assessing whether a monotonic trend exists
across ordered treatment conditions in a randomized complete block
design.

## Usage

``` r
pageTest(y, ...)

# Default S3 method
pageTest(y, groups, blocks, ...)

# S3 method for class 'formula'
pageTest(formula, data, subset, na.action = na.pass, ...)
```

## Arguments

- y:

  either a numeric vector of data values, or a data matrix.

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

- formula:

  a formula of the form `y ~ groups | blocks`.

- data:

  an optional data frame containing the variables in `formula`.

- subset:

  an optional vector specifying a subset of observations to be used.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

## Value

A list with class `"htest"` containing:

- statistic:

  the L-statistic with names attribute “L”.

- p.value:

  the p-value of the test.

- parameter:

  named vector with `k` (groups) and `n` (blocks).

- method:

  the character string `"Page test for ordered alternatives"` with
  `"(exact)"` or `"(asymptotic)"` appended.

- data.name:

  a character string giving the names of the data.

## Details

Performs a Page test for ordered alternatives using an exact algorithm
by Stefan Wellek (1989) with unreplicated blocked data.

`pageTest` can be used for analyzing unreplicated complete block designs
(i.e., there is exactly one observation in `y` for each combination of
levels of `groups` and `blocks`) where the normality assumption may be
violated.

The null hypothesis is that apart from an effect of `blocks`, the
location parameter of `y` is the same in each of the `groups`.

The alternative hypothesis is that the location parameter increases
monotonically across groups: \\H_A: \theta_1 \le \theta_2 \le \theta_3
\ldots\\ (where at least one inequality is strict).

If the decreasing direction is of interest, reverse the order of the
groups before calling `pageTest`.

The Page test for ordered alternatives is slightly more powerful than
the Friedman analysis of variance by ranks.

If `y` is a matrix, `groups` and `blocks` are obtained from the column
and row indices, respectively. `NA`'s are not allowed in `groups` or
`blocks`; if `y` contains `NA`'s, the corresponding blocks are removed.

For \\k \le 15\\ (number of groups), exact p-values are computed from
the pre-tabulated null distribution of Wellek (1989). For \\k \> 15\\, a
normal approximation is used.

## Note

Exact p-values are based on pre-tabulated distributions by Stefan Wellek
(1989), valid for \\k = 3, \ldots, 15\\ groups and any number of blocks.

## References

Page, E. (1963): Ordered hypotheses for multiple treatments: A
significance test for linear ranks. *Journal of the American Statistical
Association*, 58, 216–230.

Siegel, S. & Castellan, N. J. Jr. (1988): *Nonparametric statistics for
the behavioral sciences*. Boston, MA: McGraw-Hill.

Wellek, S. (1989): Computing exact p-values in Page's nonparametric test
against trend. *Biometrie und Informatik in Medizin und Biologie 20*,
163–170.

## See also

[`friedman.test`](https://rdrr.io/r/stats/friedman.test.html)

Other test.trend:
[`cochranArmitageTest()`](https://andrisignorell.github.io/lumen/reference/cochranArmitageTest.md),
[`coxStuartTest()`](https://andrisignorell.github.io/lumen/reference/coxStuartTest.md),
[`jonckheereTerpstraTest()`](https://andrisignorell.github.io/lumen/reference/jonckheereTerpstraTest.md)

## Examples

``` r

# Craig's data from Siegel & Castellan, p. 186
soa.mat <- matrix(c(.797,.873,.888,.923,.942,.956,
                    .794,.772,.908,.982,.946,.913,
                    .838,.801,.853,.951,.883,.837,
                    .815,.801,.747,.859,.887,.902),
                  nrow = 4, byrow = TRUE)
pageTest(soa.mat)
#> 
#>  Page test for ordered alternatives (exact)
#> 
#> data:  soa.mat
#> L = 342, k = 6, n = 4, p-value = 0.0005661
#> 

# Duller, pg. 236
pers <- matrix(c(
  1, 72, 72, 71.5, 69, 70, 69.5, 68, 68, 67, 68,
  2, 83, 81, 81, 82, 82.5, 81, 79, 80.5, 80, 81,
  3, 95, 92, 91.5, 89, 89, 90.5, 89, 89, 88, 88,
  4, 71, 72, 71, 70.5, 70, 71, 71, 70, 69.5, 69,
  5, 79, 79, 78.5, 77, 77.5, 78, 77.5, 76, 76.5, 76,
  6, 80, 78.5, 78, 77, 77.5, 77, 76, 76, 75.5, 75.5
), nrow = 6, byrow = TRUE)
colnames(pers) <- c("person", paste("week", 1:10))

# Alternative: week10 < week9 < ... < week1
pageTest(pers[, 11:2])
#> 
#>  Page test for ordered alternatives (exact)
#> 
#> data:  pers[, 11:2]
#> L = 2226, k = 10, n = 6, p-value = 9.037e-14
#> 

# long format and formula interface
plng <- data.frame(
  expand.grid(block = 1:9, group = c("B","C","D","A")),
  x = as.vector(
    matrix(c(3,2,1,4, 4,2,3,1, 4,1,2,3, 4,2,3,1,
             3,2,1,4, 4,1,2,3, 4,3,2,1, 3,1,2,4,
             3,1,4,2),
           nrow = 9, byrow = TRUE,
           dimnames = list(1:9, LETTERS[1:4]))[, c("B","C","D","A")]
  )
)
pageTest(x ~ group | block, data = plng)
#> 
#>  Page test for ordered alternatives (exact)
#> 
#> data:  x ~ group | block
#> L = 252, k = 4, n = 9, p-value = 0.0007053
#> 


# Sachs, pg. 464: L = 252
pers <- matrix(c(
  3,2,1,4,
  4,2,3,1,
  4,1,2,3,
  4,2,3,1,
  3,2,1,4,
  4,1,2,3,
  4,3,2,1,
  3,1,2,4,
  3,1,4,2), 
  nrow=9, byrow=TRUE, dimnames=list(1:9, LETTERS[1:4]))  

# Alternative: B < C < D < A
pageTest(pers[, c("B","C","D","A")])
#> 
#>  Page test for ordered alternatives (exact)
#> 
#> data:  pers[, c("B", "C", "D", "A")]
#> L = 252, k = 4, n = 9, p-value = 0.0007053
#> 

```
