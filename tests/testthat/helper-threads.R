# Evaluate `expr` with RcppParallel limited to `n` threads and restore the
# caller's setting afterwards. setThreadOptions() records its value in
# RCPP_PARALLEL_NUM_THREADS; an unset variable means "auto".
withThreads <- function(n, expr) {
  old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA)
  on.exit(RcppParallel::setThreadOptions(
    numThreads = if (is.na(old)) "auto" else as.integer(old)))
  RcppParallel::setThreadOptions(numThreads = n)
  force(expr)
}
