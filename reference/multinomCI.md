# Confidence Intervals for Multinomial Proportions

Confidence intervals for multinomial proportions are often approximated
by independent binomial confidence intervals per class, which can work
well in practice but is strictly speaking not correct. This function
computes simultaneous confidence intervals for multinomial proportions,
using one of several methods, e.g. Sison-Glaz, Wald or Wilson.

## Usage

``` r
multinomCI(
  x,
  conf.level = 0.95,
  sides = c("two.sided", "left", "right"),
  method = c("sison-glaz", "cplus1", "goodman", "wald", "waldcc", "wilson", "qh", "fs")
)
```

## Arguments

- x:

  vector of positive integers, the number of occurrences observed in
  each class; the total number of samples is `sum(x)`.

- conf.level:

  confidence level; defaults to `0.95`.

- sides:

  character string specifying the side of the confidence interval, one
  of `"two.sided"` (default), `"left"` or `"right"`; can be abbreviated.
  `"left"` corresponds to a hypothesis of `"greater"` in
  [`t.test()`](https://rdrr.io/r/stats/t.test.html).

- method:

  character string specifying which method to use, one of `"sison-glaz"`
  (default), `"cplus1"`, `"goodman"`, `"wald"`, `"waldcc"`, `"wilson"`,
  `"qh"` or `"fs"`; can be abbreviated. See ‘Details’ for the individual
  methods.

## Value

A numeric matrix with one row per class and columns:

- `est`:

  estimated difference

- `lci`:

  lower confidence limit

- `uci`:

  upper confidence limit

The number of rows correspond to the dimension of x.

## Details

Given a vector of observations with the number of samples falling in
each class of a multinomial distribution, builds the simultaneous
confidence intervals for the multinomial probabilities according to the
method passed in `method` (see `method` below for the full list). The R
code for Sison-Glaz (1995) has been translated from the SAS code written
by May and Johnson (2000). See the references for the other methods
(`qh` = Quesenberry-Hurst, `fs` = Fitzpatrick-Scott).  
Some of the methods can yield confidence limits below 0 or above 1;
these are truncated to `[0, 1]`.

## Note

Based on code by Pablo J. Villacorta Iglesias (Sison-Glaz), Andri
Signorell (Goodman, Wald, Wilson, Fitzpatrick-Scott, Quesenberry-Hurst),
adapted to coform to package standards.

## References

Fitzpatrick, S. and Scott, A. (1987). Quick simultaneous confidence
interval for multinomial proportions. *Journal of American Statistical
Association* 82(399): 875-878.

Glaz, J., Sison, C.P. (1999) Simultaneous confidence intervals for
multinomial proportions. *Journal of Statistical Planning and Inference*
82:251-262.

Goodman, L. A. (1965) On Simultaneous Confidence Intervals for
Multinomial Proportions *Technometrics*, 7, 247-254.

May, W.L., Johnson, W.D.(2000) Constructing two-sided simultaneous
confidence intervals for multinomial proportions for small counts in a
large number of cells. *Journal of Statistical Software* 5(6) . Paper
and code available at <https://www.jstatsoft.org/v05/i06>.

Quesenberry, C.P. and Hurst, D.C. (1964). Large Sample Simultaneous
Confidence Intervals for Multinational Proportions. *Technometrics*, 6:
191-195.

Sangeetha, U., Subbiah, M., Srinivasan, M. R. (2013) Mathematical
Analysis of propensity of aberration on the methods for interval
estimation of the multinomial proportions. *IOSR Journal of
Mathematics*, e-ISSN: 2278-5728,p-ISSN: 2319-765X, Volume 7, Issue 4
(Jul. - Aug. 2013), PP 23-28

Sison, C.P and Glaz, J. (1995) Simultaneous confidence intervals and
sample size determination for multinomial proportions. *Journal of the
American Statistical Association*, 90:366-369.

Wald, A. Tests of statistical hypotheses concerning several parameters
when the number of observations is large, *Trans. Am. Math. Soc.* 54
(1943) 426-482.

Wilson, E. B. Probable inference, the law of succession and statistical
inference, *J.Am. Stat. Assoc.* 22 (1927) 209-212.

## See also

Other ci.proportion: [`binomCI()`](binomCI.md),
[`binomDiffCI()`](binomDiffCI.md), [`binomRatioCI()`](binomRatioCI.md)

## Examples

``` r

# Multinomial distribution with 3 classes, from which a sample of 79 elements
# were drawn: 23 of them belong to the first class, 12 to the
# second class and 44 to the third class. Punctual estimations
# of the probabilities from this sample would be 23/79, 12/79
# and 44/79 but we want to build 95% simultaneous confidence intervals
# for the true probabilities

multinomCI(c(23, 12, 44), conf.level=0.95)
#>            est        lci       uci
#> [1,] 0.2911392 0.18987342 0.4106703
#> [2,] 0.1518987 0.05063291 0.2714298
#> [3,] 0.5569620 0.45569620 0.6764931

# single sided
multinomCI(c(23, 12, 44), conf.level=0.95, sides="left")
#>            est        lci uci
#> [1,] 0.2911392 0.20253165   1
#> [2,] 0.1518987 0.06329114   1
#> [3,] 0.5569620 0.46835443   1
multinomCI(c(23, 12, 44), conf.level=0.95, sides="right")
#>            est lci       uci
#> [1,] 0.2911392   0 0.3938117
#> [2,] 0.1518987   0 0.2545712
#> [3,] 0.5569620   0 0.6596345


x <- c(35, 74, 22, 69)

multinomCI(x, method="goodman")
#>        est        lci       uci
#> [1,] 0.175 0.11801881 0.2516431
#> [2,] 0.370 0.28986972 0.4579951
#> [3,] 0.110 0.06611447 0.1774798
#> [4,] 0.345 0.26687843 0.4324988
multinomCI(x, method="sison-glaz")
#>        est   lci       uci
#> [1,] 0.175 0.105 0.2513437
#> [2,] 0.370 0.300 0.4463437
#> [3,] 0.110 0.040 0.1863437
#> [4,] 0.345 0.275 0.4213437
multinomCI(x, method="cplus1")
#>        est   lci   uci
#> [1,] 0.175 0.100 0.250
#> [2,] 0.370 0.295 0.445
#> [3,] 0.110 0.035 0.185
#> [4,] 0.345 0.270 0.420
multinomCI(x, method="wald")
#>        est        lci       uci
#> [1,] 0.175 0.12234021 0.2276598
#> [2,] 0.370 0.30308797 0.4369120
#> [3,] 0.110 0.06663649 0.1533635
#> [4,] 0.345 0.27911853 0.4108815
multinomCI(x, method="waldcc")
#>        est        lci       uci
#> [1,] 0.175 0.11984021 0.2301598
#> [2,] 0.370 0.30058797 0.4394120
#> [3,] 0.110 0.06413649 0.1558635
#> [4,] 0.345 0.27661853 0.4133815
multinomCI(x, method="wilson")
#>        est        lci       uci
#> [1,] 0.175 0.12860515 0.2336443
#> [2,] 0.370 0.30612608 0.4387737
#> [3,] 0.110 0.07377244 0.1609269
#> [4,] 0.345 0.28259794 0.4132441

# compare to
binomCI(x, n=sum(x))
#>     est        lci       uci  x   n conf.level     sides method stdEst
#> 1 0.175 0.12860515 0.2336443 35 200       0.95 two.sided wilson   TRUE
#> 2 0.370 0.30612608 0.4387737 74 200       0.95 two.sided wilson   TRUE
#> 3 0.110 0.07377244 0.1609269 22 200       0.95 two.sided wilson   TRUE
#> 4 0.345 0.28259794 0.4132441 69 200       0.95 two.sided wilson   TRUE

# example in Goodman (1965)
multinomCI(x=c(91, 49, 37, 43), conf.level=0.95, method="goodman")
#>            est       lci       uci
#> [1,] 0.4136364 0.3342025 0.4978332
#> [2,] 0.2227273 0.1608587 0.2998874
#> [3,] 0.1681818 0.1145513 0.2401121
#> [4,] 0.1954545 0.1374690 0.2702358

# example from Sison, Glaz (1999) in Sangeetha (2013) - Table 2
#
#    Wald          Wald_CC       Wilson        Quesnberry-Hurst  
#    LL     UL     LL     UL     LL     UL     LL     UL
# 1   0.090  0.149  0.089  0.150   0.094  0.153  0.076  0.183
# 2  0.121  0.187  0.120  0.188   0.124  0.190  0.104  0.222
# 3   0.123  0.189  0.122  0.190   0.126  0.192  0.106  0.225
# 4   0.096  0.156  0.095  0.158   0.099  0.160  0.081  0.191
# 5   0.102  0.164  0.101  0.165   0.105  0.167  0.087  0.198
# 6   0.151  0.222  0.150  0.223   0.154  0.224  0.131  0.258
# 7   0.094  0.154  0.093  0.155   0.097  0.157  0.080  0.188
 
#    Goodman        Fitzpatrick-Scott  Sison-Glaz  
#    LL     UL      LL     UL          LL    UL
# 1   0.085  0.166   0.075  0.165       0.079  0.164
# 2   0.115  0.204   0.109  0.200       0.114  0.199
# 3   0.116  0.207   0.111  0.202       0.116  0.201
# 4   0.091  0.173   0.081  0.172       0.086  0.171
# 5   0.096  0.181   0.087  0.178       0.092  0.177
# 6   0.143  0.239   0.141  0.232       0.146  0.231
# 7   0.089  0.171   0.079  0.170       0.084  0.169

x <- c(56, 72, 73, 59, 62, 87, 58)
do.call(cbind, lapply(c("wald", "waldcc", "wilson", 
                        "qh", "goodman", "fs", "sison-glaz"),
                      function(m) round(multinomCI(x, method=m)[,-1], 3)))
#>        lci   uci   lci   uci   lci   uci   lci   uci   lci   uci   lci   uci
#> [1,] 0.090 0.149 0.089 0.150 0.094 0.153 0.076 0.183 0.085 0.166 0.075 0.165
#> [2,] 0.121 0.187 0.120 0.188 0.124 0.190 0.104 0.222 0.115 0.204 0.109 0.200
#> [3,] 0.123 0.189 0.122 0.190 0.126 0.192 0.106 0.225 0.116 0.207 0.111 0.202
#> [4,] 0.096 0.156 0.095 0.158 0.099 0.160 0.081 0.191 0.091 0.173 0.081 0.172
#> [5,] 0.102 0.164 0.101 0.165 0.105 0.167 0.087 0.198 0.096 0.181 0.087 0.178
#> [6,] 0.151 0.222 0.150 0.223 0.154 0.224 0.131 0.258 0.143 0.239 0.141 0.232
#> [7,] 0.094 0.154 0.093 0.155 0.097 0.157 0.080 0.188 0.089 0.171 0.079 0.170
#>        lci   uci
#> [1,] 0.079 0.164
#> [2,] 0.113 0.199
#> [3,] 0.116 0.201
#> [4,] 0.086 0.171
#> [5,] 0.092 0.177
#> [6,] 0.146 0.231
#> [7,] 0.084 0.169
       
       
```
