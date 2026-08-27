# Add Up Partial Confidence Intervals to a Total CI

Starting with a response variable that obtains different confidence
intervals (CI) when calculated with different explanatory variables, all
the values of the response variable should be added up. This function
returns the CI for the sum.

## Usage

``` r
sumCI(x)
```

## Arguments

- x:

  a matrix with 3 columns, containing the estimate in the first column
  followed by the lower and the upper confidence interval .

## Value

A named numeric vector with elements:

- `est`:

  point estimate, `sum(x)`.

- `lci`:

  lower confidence interval bound.

- `uci`:

  upper confidence interval bound.

## References

<https://stats.stackexchange.com/questions/223924/how-to-add-up-partial-confidence-intervals-to-create-a-total-confidence-interval>

## See also

[`binomCI`](binomCI.md)

Other ci.location: [`meanCI()`](meanCI.md), [`meanCIn()`](meanCIn.md),
[`meanDiffCI()`](meanDiffCI.md), [`medianCI()`](medianCI.md),
[`quantileCI()`](quantileCI.md)

## Examples

``` r
x <- do.call(rbind, 
             tapply(bedrock::Pizza$delivery_min, 
                    bedrock::Pizza$area, meanCI))
sumCI(x)
#>      est      lci      uci 
#> 78.47279 76.79179 80.15379 
```
