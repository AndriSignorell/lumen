# Breusch-Pagan Test for Detecting Heteroscedasticity in Regression Models

Tests the null hypothesis of homoscedasticity (constant error variance)
against heteroscedasticity using the Koenker variant of the
Breusch-Pagan test. The test statistic is \\BP = n \cdot R^2\\ from an
auxiliary regression of squared residuals on fitted values,
asymptotically distributed as \\\chi^2\\ with \\k\\ degrees of freedom,
where \\k\\ is the number of predictors.

## Usage

``` r
bpTest(fit)
```

## Arguments

- fit:

  a fitted [`lm`](https://rdrr.io/r/stats/lm.html) object.

## Value

An object of class `"htest"` with the following components:

- `statistic`:

  the BP test statistic.

- `parameter`:

  degrees of freedom.

- `p.value`:

  p-value based on the \\\chi^2\\ distribution.

- `method`:

  character string describing the test.

- `data.name`:

  the formula of the fitted model.

## References

Breusch, T.S. and Pagan, A.R. (1979). A simple test for
heteroscedasticity and random coefficient variation. *Econometrica*, 47,
1287–1294.

Koenker, R. (1981). A note on studentizing a test for
heteroscedasticity. *Journal of Econometrics*, 17, 107–112.

## See also

[`lm`](https://rdrr.io/r/stats/lm.html)

Other test.regression:
[`breuschGodfreyTest()`](https://andrisignorell.github.io/lumen/reference/breuschGodfreyTest.md),
[`durbinWatsonTest()`](https://andrisignorell.github.io/lumen/reference/durbinWatsonTest.md),
[`hosmerLemeshowTest()`](https://andrisignorell.github.io/lumen/reference/hosmerLemeshowTest.md),
[`leCessieTest()`](https://andrisignorell.github.io/lumen/reference/leCessieTest.md)

## Examples

``` r
fit <- lm(Sepal.Length ~ Sepal.Width, data = iris)
bpTest(fit)
#> 
#>  Breusch-Pagan test (Koenker)
#> 
#> data:  Sepal.Length ~ Sepal.Width
#> BP = 0.78243, df = 1, p-value = 0.3764
#> 
```
