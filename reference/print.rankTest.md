# Print Method for rankTest Objects

Prints pairwise comparison results produced by
[`dunnTest`](dunnTest.md), [`conoverTest`](conoverTest.md), or
[`nemenyiTest`](nemenyiTest.md).

## Usage

``` r
# S3 method for class 'rankTest'
print(x, digits = getOption("digits", 3), ...)
```

## Arguments

- x:

  an object of class `"rankTest"`.

- digits:

  number of significant digits used for printing numeric values. Passed
  to [`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html).
  Defaults to `getOption("digits", 3)`.

- ...:

  further arguments passed to
  [`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html) or
  [`print.default`](https://rdrr.io/r/base/print.default.html).

## Value

`x`, invisibly.

## See also

[`dunnTest`](dunnTest.md), [`conoverTest`](conoverTest.md),
[`nemenyiTest`](nemenyiTest.md)
