# Gompertz Distribution

The Gompertz distribution is a continuous distribution with a
non-negative real support, commonly used to model human mortality and
customer lifetime value. It is parameterized by a shape and a scale
parameter, and is characterized by an exponentially increasing hazard
rate.

## Usage

``` r
dgompertz(x, shape, rate = 1, log = FALSE)

pgompertz(q, shape, rate = 1, lower.tail = TRUE, log.p = FALSE)

qgompertz(p, shape, rate = 1, lower.tail = TRUE, log.p = FALSE)

rgompertz(n, shape, rate = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- shape, rate:

  vector of shape and rate parameters.

- log, log.p:

  logical; if TRUE, probabilities p are given as log(p).

- lower.tail:

  logical; if TRUE (default), probabilities are \\P(X \\\\\le x)\\,
  otherwise, \\P(X \> x)\\.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dgompertz()` gives the density, `pgompertz()` gives the distribution
function, `qgompertz()` gives the quantile function, and `rgompertz()`
generates random deviates.

## Details

The Gompertz distribution with `shape` parameter \\a\\ and `rate`
parameter \\b\\ has probability density function

\$\$f(x \| a, b) = be^{ax}\exp(-b/a (e^{ax} - 1))\$\$

For \\a=0\\ the Gompertz is equivalent to the exponential distribution
with constant hazard and rate \\b\\.

The probability distribution function is \$\$F(x \| a, b) = 1 -
\exp(-b/a (e^{ax} - 1))\$\$

Thus if \\a\\ is negative, letting \\x\\ tend to infinity shows that
there is a non-zero probability \\1 - \exp(b/a)\\ of living forever. On
these occasions `qgompertz()` and `rgompertz()` will return `Inf`, and
`pgompertz()` approaches \\1 - \exp(b/a)\\ rather than one.

A non-positive `rate` gives `NaN` with a warning, as in the base R
distribution functions.

**Note:** Some implementations of the Gompertz restrict \\a\\ to be
strictly positive, which ensures that the probability of survival
decreases to zero as \\x\\ increases to infinity. The more flexible
implementation given here is consistent with `streg` in Stata.

The functions `dgompertz()` and similar available in the package eha
label the parameters the other way round, so that what is called the
`shape` there is called the `rate` here, and what is called `1 / scale`
there is called the `shape` here. The terminology here is consistent
with the exponential
[`dexp()`](https://rdrr.io/r/stats/Exponential.html) and Weibull
[`dweibull()`](https://rdrr.io/r/stats/Weibull.html) distributions in R.

## Note

Based on code by Christopher Jackson previously published in the
flexsurv package, adapted to conform to package standards.

## References

Gompertz, B. (1825) On the nature of the function expressive of the law
of human mortality. *Philosophical Transactions of the Royal Society*,
**115**, 513–583.

Stata Press (2007) *Stata Release 10 Manual: Survival Analysis and
Epidemiological Tables*. Stata Press.

## See also

[distributions-overview](distributions-overview.md);
[`dexp()`](https://rdrr.io/r/stats/Exponential.html)

## Examples

``` r

dgompertz(1:3, shape = 0.1, rate = 0.2)
#> [1] 0.1791056 0.1568848 0.1341019
pgompertz(1:3, shape = 0.1, rate = 0.2)
#> [1] 0.1896928 0.3577679 0.5032744
qgompertz(seq(0.9, 0.6, -0.1), shape = 0.1, rate = 0.2)
#> [1] 7.660688 5.904049 4.712444 3.771653
rgompertz(6, shape = 0.1, rate = 0.2)
#> [1] 1.32919735 0.03655073 9.11096610 4.91622943 1.51676067 5.19308535

## for shape = 0 the Gompertz reduces to the exponential distribution
all.equal(pgompertz(1:3, shape = 0, rate = 0.2), pexp(1:3, rate = 0.2))
#> [1] TRUE

## a negative shape leaves a non-zero probability of living forever,
## for which the quantile function returns Inf
qgompertz(0.9, shape = -0.5, rate = 0.2)
#> [1] Inf

mgompertz(shape = 0.1, rate = 0.2)
#>     mean variance 
#> 3.613286 8.024897 
 
```
