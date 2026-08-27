# Compute Scores for Ordinal Contingency Tables

A utility function computing score transformations of raw data,
including normal scores, exponential scores, and Savage scores,
typically used as a preprocessing step for nonparametric tests.

## Usage

``` r
scores(x, MARGIN = 1, method = c("table", "ranks", "ridit", "modridit"))
```

## Arguments

- x:

  a contingency table (matrix or array of counts).

- MARGIN:

  an integer indicating the margin over which to compute the scores.
  Defaults to `1` (rows). Use `2` for columns.

- method:

  a character string specifying the scoring method. One of:

  - `"table"`: uses numeric dimnames if available, otherwise assigns
    sequential integers.

  - `"ranks"`: mid-ranks based on cumulative frequencies.

  - `"ridit"`: ridit scores (ranks divided by total count).

  - `"modridit"`: modified ridit scores (ranks divided by total count +
    1).

## Value

A numeric vector of scores corresponding to the levels of the selected
margin.

## Details

Computes score values for the levels of a contingency table margin.
These scores are used in several statistical procedures such as the
Cochran-Armitage test and correlation measures for ordinal data.

The function supports different scoring methods, including simple
table-based scores, ranks, and ridit-type transformations.

For `method = "table"`, numeric dimension names are used as scores if
available. Otherwise, consecutive integers starting from 1 are assigned.

For rank-based methods, scores are computed as midpoints of cumulative
frequencies along the selected margin.

Ridit and modified ridit scores are normalized versions of these ranks.

## References

Lecoutre, E. (2005). R-help mailing list discussion.
<https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html>

## See also

[`cochranArmitageTest`](https://andrisignorell.github.io/lumen/reference/cochranArmitageTest.md),
[`cor`](https://rdrr.io/r/stats/cor.html)
