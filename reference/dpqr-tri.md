# Triangular Distribution

The triangular distribution is a continuous distribution with a lower
bound, an upper bound, and a mode, producing a piecewise linear,
triangular-shaped density function. It is commonly used in risk
assessment and simulation when only the minimum, maximum, and most
likely value of a quantity are known.

## Usage

``` r
dtri(x, min = 0, max = 1, mode = 1/2)

ptri(q, min = 0, max = 1, mode = 1/2)

qtri(p, min = 0, max = 1, mode = 1/2)

rtri(n, min = 0, max = 1, mode = 1/2)
```

## Arguments

- x:

  vector of quantiles. Missing values (`NA`s) are allowed.

- min:

  vector of minimum values of the distribution of the random variable.
  The default value is `min=0`.

- max:

  vector of maximum values of the random variable. The default value is
  `max=1`.

- mode:

  vector of modes of the random variable. The default value is
  `mode=1/2`.

- q:

  vector of quantiles. Missing values (`NA`s) are allowed.

- p:

  vector of probabilities between 0 and 1. Missing values (`NA`s) are
  allowed.

- n:

  sample size. If `length(n)` is larger than 1, then `length(n)` random
  values are returned.

## Value

`dtri()` gives the density, `ptri()` gives the distribution function,
`qtri()` gives the quantile function, and `rtri()` generates random
deviates.

The triangular distribution is sometimes used as an input distribution
in probability risk assessment.

## Details

Density, distribution function, quantile function, and random generation
for the triangular distribution with parameters `min`, `max`, and
`mode`.

Let \\X\\ be a triangular random variable with parameters `min=`\\a\\,
`max=`\\b\\, and `mode=`\\c\\.

*Probability Density and Cumulative Distribution Function*  
The density function of \\X\\ is given by:

|                     |                               |                       |
|---------------------|-------------------------------|-----------------------|
| \\f(x; a, b, c) =\\ | \\\frac{2(x-a)}{(b-a)(c-a)}\\ | for \\a \le x \le c\\ |
|                     | \\\frac{2(b-x)}{(b-a)(b-c)}\\ | for \\c \le x \le b\\ |

where \\a \< c \< b\\.

The cumulative distribution function of \\X\\ is given by:

|  |  |  |
|----|----|----|
| \\F(x; a, b, c) =\\ | \\\frac{(x-a)^2}{(b-a)(c-a)}\\ | for \\a \le x \le c\\ |
|  | \\1 - \frac{(b-x)^2}{(b-a)(b-c)}\\ | for \\c \le x \le b\\ |

where \\a \< c \< b\\.

*Quantiles*  
The \\p^th\\ quantile of \\X\\ is given by:

|           |                               |                          |
|-----------|-------------------------------|--------------------------|
| \\x_p =\\ | \\a + \sqrt{(b-a)(c-a)p}\\    | for \\0 \le p \le F(c)\\ |
|           | \\b - \sqrt{(b-a)(b-c)(1-p}\\ | for \\F(c) \le p \le 1\\ |

where \\0 \le p \le 1\\.

*Random Numbers*  
Random numbers are generated using the inverse transformation method:
\$\$x = F^{-1}(u)\$\$ where \\u\\ is a random deviate from a uniform
\\\[0, 1\]\\ distribution.

*Mean and Variance*  
The mean and variance of \\X\\ are given by: \$\$E(X) = \frac{a + b +
c}{3}\$\$ \$\$Var(X) = \frac{a^2 + b^2 + c^2 - ab - ac - bc}{18}\$\$

The triangular distribution is so named because of the shape of its
probability density function. The average of two independent identically
distributed uniform random variables with parameters `min=`\\\alpha\\
and `max=`\\\beta\\ has a triangular distribution with parameters
`min=`\\\alpha\\, `max=`\\\beta\\, and `mode=`\\(\beta-\alpha)/2\\.

## Note

Based on code by Steven P. Millard previously published in the EnvStats
package, adapted to conform to package standards.

## References

Forbes, C., M. Evans, N. Hastings, and B. Peacock. (2011). Statistical
Distributions. Fourth Edition. John Wiley and Sons, Hoboken, NJ.

Johnson, N. L., S. Kotz, and N. Balakrishnan. (1995). *Continuous
Univariate Distributions, Volume 2*. Second Edition. John Wiley and
Sons, New York.

## See also

[distributions-overview](distributions-overview.md);
[Uniform](https://rdrr.io/r/stats/Uniform.html)

## Examples

``` r

# Density of a triangular distribution with parameters 
# min=10, max=15, and mode=12, evaluated at 12, 13 and 14: 

dtri(12:14, 10, 15, 12) 
#> [1] 0.4000000 0.2666667 0.1333333
## [1] 0.4000000 0.2666667 0.1333333

# The cdf of a triangular distribution with parameters 
# min=2, max=7, and mode=5, evaluated at 3, 4, and 5: 

ptri(3:5, 2, 7, 5) 
#> [1] 0.06666667 0.26666667 0.60000000
## [1] 0.06666667 0.26666667 0.60000000

# The 25'th percentile of a triangular distribution with parameters 
# min=1, max=4, and mode=3: 

qtri(0.25, 1, 4, 3) 
#> [1] 2.224745
## [1] 2.224745

# A random sample of 4 numbers from a triangular distribution with 
# parameters min=3 , max=20, and mode=12. 
# (Note: the call to set.seed simply allows you to reproduce this example.)

set.seed(10) 
rtri(4, 3, 20, 12) 
#> [1] 11.811593  9.850955 11.081885 13.539496
## [1] 11.811593  9.850955 11.081885 13.539496
```
