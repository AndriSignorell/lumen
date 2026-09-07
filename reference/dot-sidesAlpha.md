# Two-Sided Level Behind a One-Sided Bound

Internal helper shared by the `binom*CI()` family. A one-sided bound at
level \\\gamma\\ is the corresponding end of the two-sided interval at
level \\2\gamma - 1\\ (design_rules 4.1), so the tail probability the
methods have to work with is \\2(1 - \gamma)\\ rather than \\1 -
\gamma\\.

## Usage

``` r
.sidesAlpha(conf.level, sides)
```

## Arguments

- conf.level:

  the confidence level as requested by the user.

- sides:

  one of `"two.sided"`, `"left"` or `"right"`, already matched.

## Value

a single numeric giving the alpha the interval methods use.
