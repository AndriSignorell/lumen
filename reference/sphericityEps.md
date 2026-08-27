# Greenhouse-Geisser And Huynh-Feldt Epsilons

Computes the Greenhouse-Geisser and Huynh-Feldt epsilon correction
factors for assessing and correcting for violations of the sphericity
assumption in repeated measures ANOVA.

## Usage

``` r
sphericityEps(S, p, nGroups, n, method = c("both", "gg", "hf"))
```

## Arguments

- S:

  pxp covariance matrix.

- p:

  dimension of observation vectors.

- nGroups:

  number of groups.

- n:

  number of subjects.

- method:

  a character string specifying which epsilon to return, must be one of
  `"both"` (default), `"gg"` for Greenhouse-Geisser, or `"hf"` for
  Huynh-Feldt.

## Value

A numeric value.

## Note

Based on code by Hans Rudolf Roth, adapted to conform to package
standards.

## References

Vonesh, E.F., Chinchilli, V.M. (1997) *Linear and Nonlinear Models for
the Analysis of Repeated Measurements* Marcel Dekker, New York, p.84-86

Crowder, M.J., Hand, D.J. (1990) *Analysis of Repeated Measures*.
Chapman & Hall, London, p.54-55

## See also

[`aov`](https://rdrr.io/r/stats/aov.html)

## Examples

``` r

# a 4x4 covariance matrix among 4 repeated measurements, 20 subjects,
# one between-subject group
set.seed(1)
A <- matrix(rnorm(16), 4, 4)
S <- A %*% t(A) + diag(4)

sphericityEps(S, p = 4, nGroups = 1, n = 20)
#>        gg        hf 
#> 0.6922218 0.7786745 
sphericityEps(S, p = 4, nGroups = 1, n = 20, method = "gg")
#>        gg 
#> 0.6922218 
```
