# Check and resolve verbose level

Resolves verbosity level using the following priority:

- function argument

- global option `DescTools.verbose`

- default (2)

## Usage

``` r
.checkVerbose(verbose = NULL)
```

## Arguments

- verbose:

  Optional integer (1-3).

## Value

Integer in {1, 2, 3}.
