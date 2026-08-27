# lumen (development version)

## New features

- Inferential layer of the DescToolsX package suite: classical, robust,
  nonparametric and specialised hypothesis tests, confidence intervals
  for the common estimands, and probability distributions not available
  in base R.
- Tests return ordinary `"htest"` objects, so they print and combine
  like the procedures in base R.
- Confidence intervals follow one argument convention throughout
  (`conf.level`, `sides`, `method`), across proportions, counts, location,
  scale and correlation.
- `bootCI()` provides bootstrap intervals for arbitrary statistics; the
  resampling is implemented in C++ via Rcpp, RcppArmadillo and
  RcppParallel.
- `postHoc()` collects the post-hoc procedures for multiple comparisons
  behind one interface, with a plot method for the resulting intervals.

## Acknowledgements

Parts of the code and documentation were reviewed with the help of large
language models (OpenAI Codex, Anthropic Claude). Every suggestion was
assessed, edited and verified by the maintainer, who remains solely
responsible for the content of this package.
