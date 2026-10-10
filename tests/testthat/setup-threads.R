# CRAN allows at most two cores per check; RcppParallel would otherwise use
# all of them and trigger the "CPU time > elapsed time" NOTE. R CMD check
# runs the tests in a fresh process, so the setting is not restored.
RcppParallel::setThreadOptions(numThreads = 2)
