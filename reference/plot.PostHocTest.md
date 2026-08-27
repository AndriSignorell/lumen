# Plot Routine for Posthoc Tests Results

The plot method visualizes pairwise comparisons from a `PostHocTest`
object using dot plots. Each component of `x` is plotted separately,
showing estimated differences along with confidence intervals. A
vertical reference line at zero is added.

## Usage

``` r
# S3 method for class 'PostHocTest'
plot(x, ...)
```

## Arguments

- x:

  an object of class `"PostHocTest"`, typically returned by
  [`postHocTest`](postHoc.md).

- ...:

  additional graphical parameters passed to
  [`plotDot`](https://andrisignorell.github.io/pharos/reference/plotDot.html)
  and base plotting functions.

## Value

Invisibly returns `NULL`.

## Details

For each factor in `x`, a dot plot is produced displaying pairwise
differences in means and their confidence intervals. The confidence
level is taken from the `"conf.level"` attribute of `x`.

## See also

[`pharos::plotDot()`](https://andrisignorell.github.io/pharos/reference/plotDot.html)

Other test.posthoc: [`conoverTest()`](conoverTest.md),
[`dscfTest()`](dscfTest.md), [`dunnTest()`](dunnTest.md),
[`dunnettTest()`](dunnettTest.md),
[`gamesHowellTest()`](gamesHowellTest.md),
[`nemenyiTest()`](nemenyiTest.md), [`postHoc`](postHoc.md),
[`scheffeTest()`](scheffeTest.md), [`signifDiff()`](signifDiff.md),
[`steelTest()`](steelTest.md)
