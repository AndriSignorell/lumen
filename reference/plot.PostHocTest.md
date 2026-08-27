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
  [`postHocTest`](https://andrisignorell.github.io/lumen/reference/postHoc.md).

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

Other test.posthoc:
[`conoverTest()`](https://andrisignorell.github.io/lumen/reference/conoverTest.md),
[`dscfTest()`](https://andrisignorell.github.io/lumen/reference/dscfTest.md),
[`dunnTest()`](https://andrisignorell.github.io/lumen/reference/dunnTest.md),
[`dunnettTest()`](https://andrisignorell.github.io/lumen/reference/dunnettTest.md),
[`gamesHowellTest()`](https://andrisignorell.github.io/lumen/reference/gamesHowellTest.md),
[`nemenyiTest()`](https://andrisignorell.github.io/lumen/reference/nemenyiTest.md),
[`postHoc`](https://andrisignorell.github.io/lumen/reference/postHoc.md),
[`scheffeTest()`](https://andrisignorell.github.io/lumen/reference/scheffeTest.md),
[`signifDiff()`](https://andrisignorell.github.io/lumen/reference/signifDiff.md),
[`steelTest()`](https://andrisignorell.github.io/lumen/reference/steelTest.md)
