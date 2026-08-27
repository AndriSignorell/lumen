# Runs Test for Detecting Non-Random Order in Sequences

A nonparametric test for randomness of a sequence, based on the number
of runs (consecutive sequences of identical values or values above/below
a threshold).

## Usage

``` r
runsTest(x, ...)

# S3 method for class 'formula'
runsTest(formula, data, subset, na.action = na.pass, ...)

# Default S3 method
runsTest(
  x,
  y = NULL,
  alternative = c("two.sided", "less", "greater"),
  exact = NULL,
  correct = TRUE,
  na.rm = TRUE,
  ...
)
```

## Arguments

- x:

  a dichotomous vector of data values or a (non-empty) numeric vector of
  data values.

- ...:

  further arguments passed to methods.

- formula:

  a formula of the form `lhs ~ rhs` where `lhs` gives the data values
  and `rhs` the corresponding groups (exactly two groups required).

- data:

  an optional data frame containing the variables in `formula`.

- subset:

  an optional vector specifying a subset of observations to be used.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

- y:

  an optional (non-empty) numeric vector of data values for the
  Wald-Wolfowitz two-sample test.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"less"` (fewer runs, clustering) or
  `"greater"` (more runs, alternation).

- exact:

  a logical indicating whether an exact p-value should be computed. By
  default exact values are calculated for \\n_0 + n_1 \le 30\\ and the
  normal approximation otherwise.

- correct:

  a logical indicating whether to apply a continuity correction when
  computing the test statistic. Default is `TRUE`. Ignored when
  `exact = TRUE`.

- na.rm:

  logical. If `TRUE` (default), `NA`s are removed before the test.

## Value

A list of class `"htest"` containing:

- statistic:

  the standardized z-statistic (normal approximation only; `NULL` for
  exact tests).

- parameter:

  named vector with the number of runs, `m` (\\n_0\\: observations below
  threshold) and `n` (\\n_1\\: observations above threshold).

- p.value:

  the p-value for the test.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string indicating the test performed.

- data.name:

  a character string giving the name of the data.

## Details

Performs a test whether the elements of `x` are serially independent by
counting how many runs there are above and below a threshold. If `y` is
supplied, a two-sample Wald-Wolfowitz test for the equality of two
distributions is computed.

**The runs test for randomness** requires a dichotomous sequence. For a
numeric variable `x` with more than two distinct values, the sequence is
dichotomised by comparing each observation to the median. Observations
exactly equal to the median are removed before the test, as is standard
in the runs test literature. To use a different threshold, pass a
logical vector directly: `runsTest(x > mean(x))`.

The normal approximation uses the expected number of runs under the null
\$\$\mu_r = \frac{2 n_0 n_1}{n_0 + n_1} + 1\$\$ and its variance
\$\$\sigma_r^2 = \frac{2 n_0 n_1 (2 n_0 n_1 - n_0 - n_1)}{ (n_0 + n_1)^2
(n_0 + n_1 - 1)}\$\$ as \\\hat{z} = (r - \mu_r + c) / \sigma_r\\, where
\\n_0, n_1\\ are the number of values below/above the threshold and
\\r\\ is the number of runs.

Setting `correct = TRUE` applies a continuity correction as SAS (and
SPSS for \\n \< 50\\) does: \\c = 0.5\\ if \\r \< \mu_r\\ and \\c =
-0.5\\ if \\r \> \mu_r\\.

The exact p-value is computed from the conditional distribution of the
number of runs given \\n_0\\ and \\n_1\\, implemented in C++.

**Interpretation of alternatives:**

- `"less"`: fewer runs than expected — indicates clustering or positive
  serial correlation.

- `"greater"`: more runs than expected — indicates alternation or
  negative serial correlation.

- `"two.sided"`: departure in either direction.

**The Wald-Wolfowitz test** is a two-sample nonparametric test for the
equality of two continuous distributions against general alternatives.
Exact p-values are not valid in the presence of inter-group ties; a
warning is issued when such ties are detected.

## References

Wackerly, D., Mendenhall, W., Scheaffer, R. L. (1986) *Mathematical
Statistics with Applications*, 3rd Ed., Duxbury Press, CA.

Wald, A. and Wolfowitz, J. (1940): On a test whether two samples are
from the same population. *Annals of Mathematical Statistics* **11**,
147–162.

Siegel, S. (1956) *Nonparametric Statistics for the Behavioural
Sciences*, McGraw-Hill Kogakusha, Tokyo.

## See also

[`rle`](https://rdrr.io/r/base/rle.html)

Other test.timeseries:
[`adfTest()`](https://andrisignorell.github.io/lumen/reference/adfTest.md),
[`bartelsRankTest()`](https://andrisignorell.github.io/lumen/reference/BartelsRankTest.md),
[`kpssTest()`](https://andrisignorell.github.io/lumen/reference/kpssTest.md),
[`vonNeumannTest()`](https://andrisignorell.github.io/lumen/reference/vonNeumannTest.md)

## Examples

``` r
# dichotomous character vector
x <- c("S","S","T","S","T","T","T","S","T")
runsTest(x)
#> 
#>  Runs Test for Randomness
#> 
#> data:  x
#> runs = 6, m = 4, n = 5, p-value = 0.7619
#> alternative hypothesis: true number of runs is not equal to the expected number
#> 

# numeric vector: compared against median; ties at median removed
x <- c(13,3,14,14,1,14,3,8,14,17,9,14,13,2,16,1,3,12,13,14)
runsTest(x)
#> 
#>  Runs Test for Randomness
#> 
#> data:  x
#> runs = 12, m = 9, n = 8, p-value = 0.226
#> alternative hypothesis: true number of runs is not equal to the expected number
#> 

# SPSS example
x <- c(31,23,36,43,51,44,12,26,43,75,2,3,15,18,78,24,13,27,86,61,13,7,6,8)
runsTest(x, exact = TRUE)
#> 
#>  Runs Test for Randomness
#> 
#> data:  x
#> runs = 10, m = 12, n = 12, p-value = 0.3009
#> alternative hypothesis: true number of runs is not equal to the expected number
#> 
runsTest(x, exact = FALSE)
#> 
#>  Runs Test for Randomness
#> 
#> data:  x
#> z = -1.0436, runs = 10, m = 12, n = 12, p-value = 0.2967
#> alternative hypothesis: true number of runs is not equal to the expected number
#> 

# Wald-Wolfowitz two-sample test
A <- c(35,44,39,50,48,29,60,75,49,66)
B <- c(17,23,13,24,33,21,18,16,32)
runsTest(A, B, exact = TRUE)
#> 
#>  Wald-Wolfowitz Runs Test
#> 
#> data:  A and B
#> runs = 4, m = 10, n = 9, p-value = 0.003139
#> alternative hypothesis: true number of runs is not equal to the expected number
#> 
```
