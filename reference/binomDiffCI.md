# Confidence Intervals for the Difference of Two Binomial Proportions

Computes confidence intervals for the difference between two independent
binomial proportions. A variety of classical and modern methods are
available, which may yield substantially different results, particularly
for small sample sizes or extreme proportions.

## Usage

``` r
binomDiffCI(
  x1,
  n1,
  x2,
  n2,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("miettinen-nurminen", "newcombe-score", "newcombe-score-cc",
    "mee-farrington-manning", "agresti-caffo", "wald", "wald-cc", "exact",
    "brown-li-jeffreys", "hauck-anderson", "beal", "haldane", "jeffreys-perks")
)
```

## Arguments

- x1:

  number of successes in the first group.

- n1:

  number of trials in the first group.

- x2:

  number of successes in the second group.

- n2:

  number of trials in the second group.

- conf.level:

  confidence level, default is 0.95.

- sides:

  a character string specifying the type of confidence interval:
  `"two.sided"` (default), `"left"`, or `"right"`. Partial matching is
  allowed.

- method:

  one of: `"miettinen-nurminen"`, `"newcombe-score"`,
  `"newcombe-score-cc"`, `"mee-farrington-manning"`, `"agresti-caffo"`,
  `"wald"`, `"wald-cc"`, `"exact"`, `"brown-li-jeffreys"`,
  `"hauck-anderson"`, `"beal"`, `"haldane"`, `"jeffreys-perks"`.

## Value

If recycling yields a single case, a named numeric vector with elements:

- `est`:

  point estimate of the difference in binomial proportions,
  `x1/n1 - x2/n2`.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

If recycling yields multiple cases, a data frame with one row per case
is returned. Its first three columns are `est`, `lci`, and `uci`; the
remaining columns contain the recycled argument values.

## Details

All arguments are vectorized and recycled according to standard R rules.

The difference in proportions is estimated by \$\$ \hat{\delta} =
\hat{p}\_1 - \hat{p}\_2 = \frac{x_1}{n_1} - \frac{x_2}{n_2}. \$\$

**Wald**: The traditional large-sample normal approximation interval
based on the asymptotic distribution of \\\hat{\delta}\\.

**Wald with continuity correction**: A continuity-corrected version of
the Wald interval. The correction term \\(1/n_1 + 1/n_2)/2\\ is added or
subtracted from the test statistic depending on its sign.

**Agresti-Caffo**: A simple adjustment of the Wald interval (Agresti and
Caffo, 2000) obtained by adding one success and one failure to each
group. This approach performs well in many practical situations.

**Newcombe score**: Based on inverting the Wilson score interval for
each proportion and combining them to obtain an interval for the
difference (Newcombe, 1998).

**Newcombe score with continuity correction**: A continuity-corrected
variant of the Newcombe score method.

**Miettinen-Nurminen**: Based on restricted maximum likelihood
estimation obtained by solving a cubic equation (Miettinen and Nurminen,
1985). Often recommended for small to moderate sample sizes.

**Mee-Farrington-Manning**: Uses the same maximum likelihood estimators
as the Miettinen-Nurminen method but applies a different correction
factor (Mee, 1984; Farrington and Manning, 1990).

**Brown-Li-Jeffreys**: A method proposed by Brown and Li (2005).

**Hauck-Anderson**: A large-sample method described by Hauck and
Anderson (1986).

**Beal**: An asymptotic method intended for use with small samples
(Beal, 1987).

**Haldane**: Described in Newcombe (1998), based on adding 0.5 to all
cells.

**Jeffreys-Perks**: Also described in Newcombe (1998), based on
Bayesian-type adjustments.

Some methods may produce limits outside the admissible parameter space
\\\[-1, 1\]\\. In such cases, interval bounds are truncated to remain
within the valid range.

**Which interval should be used?**  
The choice of method remains an active topic of discussion. The Wald
interval is known to perform poorly in many practical situations.
Reviews such as Fagerland et al. (2011) provide comparative evaluations
and recommendations.

Miettinen-Nurminen might be a sensible default. Newcombe-Score is almost
equivalent.

## References

Agresti A, Caffo B (2000). Simple and effective confidence intervals for
proportions and difference of proportions result from adding two
successes and two failures. *The American Statistician*, 54(4), 280-288.

Beal SL (1987). Asymptotic confidence intervals for the difference
between two binomial parameters for use with small samples.
*Biometrics*, 43, 941-950.

Brown L, Li X (2005). Confidence intervals for two sample binomial
distribution. *Journal of Statistical Planning and Inference*, 130(1),
359-375.

Fagerland MW, Lydersen S, Laake P (2011). Recommended confidence
intervals for two independent binomial proportions. *Statistical Methods
in Medical Research*.

Farrington CP, Manning G (1990). Test statistics and sample size
formulae for comparative binomial trials. *Statistics in Medicine*, 9,
1447-1454.

Hauck WW, Anderson S (1986). A comparison of large-sample confidence
interval methods for the difference of two binomial probabilities. *The
American Statistician*, 40(4), 318-322.

Mee RW (1984). Confidence bounds for the difference between two
probabilities. *Biometrics*, 40, 1175-1176.

Miettinen OS, Nurminen M (1985). Comparative analysis of two rates.
*Statistics in Medicine*, 4, 213-226.

Newcombe RG (1998). Interval estimation for the difference between
independent proportions. *Statistics in Medicine*, 17, 873-890.

## See also

[`binom.test`](https://rdrr.io/r/stats/binom.test.html),
[`prop.test`](https://rdrr.io/r/stats/prop.test.html)

Other ci.proportion:
[`binomCI()`](https://andrisignorell.github.io/lumen/reference/binomCI.md),
[`binomRatioCI()`](https://andrisignorell.github.io/lumen/reference/binomRatioCI.md),
[`multinomCI()`](https://andrisignorell.github.io/lumen/reference/multinomCI.md)

## Examples

``` r

x1 <- 56; n1 <- 70; x2 <- 48; n2 <- 80
meths <- c("wald", "wald-cc", "agresti-caffo", 
                    "newcombe-score", "newcombe-score-cc", 
                    "miettinen-nurminen", "mee-farrington-manning", 
                    "brown-li-jeffreys", "hauck-anderson")
                    
xci <- binomDiffCI(x1, n1, x2, n2, method=meths)
pharos::fm(xci[,-1], digits=4)
#>      lci    uci      x1      n1      x2      n2 conf.level     sides
#> 1 0.0575 0.3425 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#> 2 0.0441 0.3559 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#> 3 0.0525 0.3358 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#> 4 0.0524 0.3339 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#> 5 0.0428 0.3422 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#> 6 0.0528 0.3382 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#> 7 0.0533 0.3377 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#> 8 0.0540 0.3400 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#> 9 0.0494 0.3506 56.0000 70.0000 48.0000 80.0000     0.9500 two.sided
#>                   method
#> 1                   wald
#> 2                wald-cc
#> 3          agresti-caffo
#> 4         newcombe-score
#> 5      newcombe-score-cc
#> 6     miettinen-nurminen
#> 7 mee-farrington-manning
#> 8      brown-li-jeffreys
#> 9         hauck-anderson

x1 <- 9; n1 <- 10; x2 <- 3; n2 <- 10
yci <- binomDiffCI(x1, n1, x2, n2, method=meths)
pharos::fm(yci[, -1], digits=4)
#>      lci    uci     x1      n1     x2      n2 conf.level     sides
#> 1 0.2605 0.9395 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#> 2 0.1605 1.0000 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#> 3 0.1600 0.8400 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#> 4 0.1705 0.8090 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#> 5 0.1013 0.8387 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#> 6 0.1700 0.8406 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#> 7 0.1821 0.8370 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#> 8 0.1869 0.9040 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#> 9 0.1922 1.0000 9.0000 10.0000 3.0000 10.0000     0.9500 two.sided
#>                   method
#> 1                   wald
#> 2                wald-cc
#> 3          agresti-caffo
#> 4         newcombe-score
#> 5      newcombe-score-cc
#> 6     miettinen-nurminen
#> 7 mee-farrington-manning
#> 8      brown-li-jeffreys
#> 9         hauck-anderson

# https://www.lexjansen.com/wuss/2016/127_Final_Paper_PDF.pdf, page 9
bedrock::setNamesX(round(
  binomDiffCI(56, 70, 48, 80, 
              method=c("wald", "wald-cc", "haldane", 
                       "jeffreys-perks", "mee-farrington-manning",
                       "miettinen-nurminen", "newcombe-score", 
                       "newcombe-score-cc", 
                       "hauck-anderson", "agresti-caffo" ,
                       "brown-li-jeffreys")
  )[,c(2,3)], 4),
  rownames=c("1. Wald, no CC", "2. Wald, CC", "3. Haldane", "4. Jeffreys-Perks",
             "5. Mee", "6. Miettinen-Nurminen", "10. Score, no CC", "11. Score, CC",
             "12. Hauck-Andersen", "13. Agresti-Caffo", "16. Brown-Li"))
#>                          lci    uci
#> 1. Wald, no CC        0.0575 0.3425
#> 2. Wald, CC           0.0441 0.3559
#> 3. Haldane            0.0535 0.3351
#> 4. Jeffreys-Perks     0.0531 0.3355
#> 5. Mee                0.0533 0.3377
#> 6. Miettinen-Nurminen 0.0528 0.3382
#> 10. Score, no CC      0.0524 0.3339
#> 11. Score, CC         0.0428 0.3422
#> 12. Hauck-Andersen    0.0494 0.3506
#> 13. Agresti-Caffo     0.0525 0.3358
#> 16. Brown-Li          0.0540 0.3400
```
