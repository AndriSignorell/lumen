
#pragma once

// ============================================================
// boot_framework.h
//
// Generic parallel bootstrap framework using RcppParallel.
//
// To add a new bootstrapped statistic:
//
//   1. Create a new .cpp file
//   2. #include "boot_framework.h"
//   3. Define a StatFn struct with:
//        double compute(const arma::mat& X, const arma::vec& y) const
//      or for scalar-input functions use the helpers below.
//   4. Call run_boot() or run_boot_matrix()
//
// ============================================================

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <random>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;


// ============================================================
// scalar bootstrap worker
//
// StatFn must implement:
//   double compute(const arma::mat& X, const arma::vec& y) const
//
// ============================================================

template <typename StatFn>
struct ScalarBootWorker : public Worker {

  const RMatrix<double> X;
  const RVector<double> y;
  const std::size_t     n;
  const std::size_t     p;
  unsigned int          base_seed;
  RVector<double>       stats;
  StatFn                fn;

  ScalarBootWorker(NumericMatrix X,
                   NumericVector y,
                   NumericVector stats,
                   unsigned int  base_seed,
                   StatFn        fn)
    : X(X), y(y),
      n(X.nrow()),
      p(X.ncol()),
      base_seed(base_seed),
      stats(stats),
      fn(fn) {}

  void operator()(std::size_t begin, std::size_t end) {

    std::vector<std::size_t> idx(n);

    for (std::size_t b = begin; b < end; b++) {

      // thread-local RNG seeded per iteration
      std::mt19937 rng(base_seed + static_cast<unsigned>(b));
      std::uniform_int_distribution<std::size_t> dist(0, n - 1);

      for (std::size_t i = 0; i < n; i++)
        idx[i] = dist(rng);

      // build bootstrap sample
      arma::mat Xb(n, p);
      arma::vec yb(n);

      for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < p; j++)
          Xb(i, j) = X(idx[i], j);
        yb[i] = y[idx[i]];
      }

      try {
        stats[b] = fn.compute(Xb, yb);
      } catch (...) {
        stats[b] = NA_REAL;
      }
    }
  }
};


// ============================================================
// matrix bootstrap worker
//
// StatFn must implement:
//   arma::vec compute(const arma::mat& X, const arma::vec& y) const
//
// Used for vector-valued statistics (e.g. regression coefficients).
//
// ============================================================

template <typename StatFn>
struct MatrixBootWorker : public Worker {

  const RMatrix<double> X;
  const RVector<double> y;
  const std::size_t     n;
  const std::size_t     p;
  const std::size_t     k;     // number of output statistics
  unsigned int          base_seed;
  RMatrix<double>       stats; // B x k
  StatFn                fn;

  MatrixBootWorker(NumericMatrix X,
                   NumericVector y,
                   NumericMatrix stats,
                   unsigned int  base_seed,
                   StatFn        fn)
    : X(X), y(y),
      n(X.nrow()),
      p(X.ncol()),
      k(stats.ncol()),
      base_seed(base_seed),
      stats(stats),
      fn(fn) {}

  void operator()(std::size_t begin, std::size_t end) {

    std::vector<std::size_t> idx(n);

    for (std::size_t b = begin; b < end; b++) {

      std::mt19937 rng(base_seed + static_cast<unsigned>(b));
      std::uniform_int_distribution<std::size_t> dist(0, n - 1);

      for (std::size_t i = 0; i < n; i++)
        idx[i] = dist(rng);

      arma::mat Xb(n, p);
      arma::vec yb(n);

      for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < p; j++)
          Xb(i, j) = X(idx[i], j);
        yb[i] = y[idx[i]];
      }

      try {
        arma::vec result = fn.compute(Xb, yb);
        for (std::size_t j = 0; j < k; j++)
          stats(b, j) = result[j];
      } catch (...) {
        for (std::size_t j = 0; j < k; j++)
          stats(b, j) = NA_REAL;
      }
    }
  }
};


// same quantile type=7 implementation as in R
inline double quantile_type7(const NumericVector& x,
                             double p){
  
  int n = x.size();
  
  if(n == 1)
    return x[0];
  
  double h = (n - 1) * p;
  int    i = std::floor(h);
  double frac = h - i;
  
  if(i >= n - 1)
    return x[n - 1];
  
  return x[i] +
    frac * (x[i + 1] - x[i]);
}


// ============================================================
// percentile CI from sorted bootstrap distribution
// ============================================================


inline NumericVector boot_percentile_ci(NumericVector stats,
                                        double        est,
                                        double        alpha) {
  
  LogicalVector ok = !is_na(stats);
  NumericVector s  = stats[ok];
  
  if (s.size() == 0)
    Rcpp::stop("All bootstrap samples failed.");
  
  std::sort(s.begin(), s.end());
  
  double lo = quantile_type7(s, alpha / 2.0);
  double hi = quantile_type7(s, 1.0 - alpha / 2.0);
  
  NumericVector out(3);
  
  out[0] = est;
  out[1] = lo;
  out[2] = hi;
  
  out.attr("names") =
    CharacterVector::create("est", "lci", "uci");
  
  return out;
}


// ============================================================
// run_boot: scalar statistic -> named vector (est, lci, uci)
// ============================================================

template <typename StatFn>
NumericVector run_boot(NumericMatrix X,
                       NumericVector y,
                       int           R,
                       double        alpha,
                       int           seed,
                       StatFn        fn) {

  unsigned int base_seed = (seed < 0)
    ? std::random_device{}()
    : static_cast<unsigned int>(seed);

  NumericVector stats(R);

  ScalarBootWorker<StatFn> worker(X, y, stats, base_seed, fn);
  parallelFor(0, R, worker);

  // original estimate on full data
  arma::mat Xa = Rcpp::as<arma::mat>(X);
  arma::vec ya = Rcpp::as<arma::vec>(y);
  double    est = fn.compute(Xa, ya);

  return boot_percentile_ci(stats, est, alpha);
}


// ============================================================
// run_boot_matrix: vector statistic -> matrix (k rows x 3 cols)
//
// Each row: one parameter, columns: est / lci / uci
// ============================================================


template <typename StatFn>
NumericMatrix run_boot_matrix(NumericMatrix X,
                              NumericVector y,
                              int           R,
                              double        alpha,
                              int           seed,
                              StatFn        fn) {
  
  unsigned int base_seed = (seed < 0)
  ? std::random_device{}()
    : static_cast<unsigned int>(seed);
  
  // determine k from a trial run on full data
  arma::mat Xa  = Rcpp::as<arma::mat>(X);
  arma::vec ya  = Rcpp::as<arma::vec>(y);
  arma::vec est = fn.compute(Xa, ya);
  
  int k = est.n_elem;
  
  NumericMatrix stats(R, k);   // R x k
  
  MatrixBootWorker<StatFn> worker(
      X, y, stats, base_seed, fn
  );
  
  parallelFor(0, R, worker);
  
  // output: k rows x 3 cols
  NumericMatrix out(k, 3);
  
  colnames(out) =
    CharacterVector::create(
      "est", "lci", "uci"
    );
  
  for (int j = 0; j < k; j++) {
    
    NumericVector s = stats(_, j);
    
    LogicalVector ok = !is_na(s);
      
      NumericVector sv = s[ok];
      out(j, 0) = est[j];
      
      if (sv.size() == 0) {
        out(j, 1) = NA_REAL;
        out(j, 2) = NA_REAL;
        
        continue;
      }
      
      std::sort(sv.begin(), sv.end());
      double lci = quantile_type7( sv, alpha / 2.0 );
      double uci = quantile_type7( sv, 1.0 - alpha / 2.0 );
      
      out(j, 1) = lci;
      out(j, 2) = uci;
  }
  
  return out;
}



// ============================================================
// convenience helper: wrap a scalar vector x into a 1-column
// matrix so scalar-only functions fit the X/y interface.
//
// Usage:
//   run_boot(vec_to_matrix(x), dummy_vec(x.size()), ...)
//
// ============================================================

inline NumericMatrix vec_to_matrix(NumericVector x) {
  NumericMatrix m(x.size(), 1);
  m(_, 0) = x;
  return m;
}

inline NumericVector dummy_vec(int n) {
  return NumericVector(n, 0.0);
}


