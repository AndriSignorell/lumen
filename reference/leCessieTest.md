# Le Cessie-Van Houwelingen-Copas-Hosmer Global Goodness of Fit Test for Assessing Global Logistic Regression Fit

Computes the le Cessie-van Houwelingen-Copas-Hosmer unweighted sum of
squares test for global goodness of fit of a logistic regression model.

## Usage

``` r
leCessieTest(x, ...)

# S3 method for class 'glm'
leCessieTest(x, ...)

# Default S3 method
leCessieTest(x, obs, X, ...)

# S3 method for class 'LeCessieTest'
print(x, digits = 4, ...)
```

## Arguments

- x:

  a fitted binomial [`glm`](https://rdrr.io/r/stats/glm.html) object
  (`glm` method), or a numeric vector of fitted probabilities, each in
  \\\[0, 1\]\\, without missing values (default method).

- ...:

  further arguments passed to methods.

- obs:

  numeric vector of observed binary outcomes (0 or 1), of the same
  length as `x`, without missing values; unused for the `glm` method.

- X:

  the full numeric design matrix used to fit the model, including the
  intercept column, with `nrow(X) == length(x)` and without missing
  values; unused for the `glm` method. See the Details.

- digits:

  number of significant digits to display.

## Value

An object of class `c("LeCessieTest", "htest")`, which is a list with
components:

- statistic:

  the standardised Z statistic, named `"Z"`.

- p.value:

  two-sided p-value from the standard normal distribution.

- method:

  a character string describing the test.

- sse:

  observed sum of squared errors.

- expected:

  expected sum of squared errors under H0.

- sd:

  standard deviation of the sum of squared errors under H0.

- data.name:

  a character string with the name(s) of the data.

## Details

The test compares the observed sum of squared residuals to its expected
value under the null hypothesis of correct model specification. The
standardised difference follows an approximate standard normal
distribution.

Unlike the Hosmer-Lemeshow tests, this test does not rely on grouping
and is therefore sensitive to a different class of model
misspecification.

The default method requires the *full* design matrix `X` used to fit the
model, i.e. including the intercept column, as produced by
`model.matrix(fit)`. Omitting the intercept column changes the
projection used to estimate the standard deviation of the test statistic
under the null and gives a materially different, incorrect result. The
`glm` method extracts the design matrix, fitted probabilities and
observed outcomes directly from a fitted model object and is therefore
the safer interface.

## Note

Adapted from code by Matthias Kohl previously published as
`HLgof.test()` in the MKmisc package, adapted to conform to package
standards.

## References

le Cessie, S. and van Houwelingen, J.C. (1991) A goodness-of-fit test
for binary regression models based on smoothing methods. *Biometrics*,
47, 1267-1282.

Hosmer, D.W., Hosmer, T., le Cessie, S. and Lemeshow, S. (1997) A
comparison of goodness-of-fit tests for the logistic regression model.
*Statistics in Medicine*, 16, 965-980.

## See also

[`glm()`](https://rdrr.io/r/stats/glm.html)

Other test.regression:
[`bpTest()`](https://andrisignorell.github.io/lumen/reference/bpTest.md),
[`breuschGodfreyTest()`](https://andrisignorell.github.io/lumen/reference/breuschGodfreyTest.md),
[`durbinWatsonTest()`](https://andrisignorell.github.io/lumen/reference/durbinWatsonTest.md),
[`hosmerLemeshowTest()`](https://andrisignorell.github.io/lumen/reference/hosmerLemeshowTest.md)

## Examples

``` r
set.seed(111)
x1  <- factor(sample(1:3, 50, replace = TRUE))
x2  <- rnorm(50)
obs <- sample(c(0, 1), 50, replace = TRUE)

model <- glm(obs ~ x1 + x2, family = binomial)

# glm method: design matrix, fitted values and outcomes are extracted
# from the model, including the intercept column
leCessieTest(model)
#> 
#>   le Cessie-van Houwelingen-Copas-Hosmer global goodness of fit test 
#> 
#> data:  obs ~ x1 + x2 
#> Z = -0.8651 , p-value = 0.387 
#> Sum of squared errors: 12 (expected: 12.01 )
#> 

# equivalent call with explicit arguments (note: full design matrix!)
leCessieTest(fitted(model), obs, model.matrix(model))
#> 
#>   le Cessie-van Houwelingen-Copas-Hosmer global goodness of fit test 
#> 
#> data:  fitted(model) and obs 
#> Z = -0.8651 , p-value = 0.387 
#> Sum of squared errors: 12 (expected: 12.01 )
#> 
```
