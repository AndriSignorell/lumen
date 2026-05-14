
// ============================================================
// median_boot.cpp
//
// Bootstrap confidence interval for the median.
// Uses the generic parallel bootstrap framework.
//
// R interface:
//   median_boot_cpp(x, B = 1000, alpha = 0.05, seed = -1)
//
// ============================================================

#include "boot_framework.h"

// ============================================================
// statistic: median of x
// ============================================================

double median_cpp(arma::vec x){
  
  size_t n = x.n_elem;
  
  if(n == 0)
    return NA_REAL;
  
  double* begin = x.memptr();
  
  size_t mid = n / 2;
  
  std::nth_element(begin, begin + mid, begin + n);
  double med = begin[mid];
  
  if(n % 2 == 0){
    double lower = *std::max_element(begin, begin + mid);
    med = 0.5 * (med + lower);
  }
  
  return med;
}


struct MedianFn {
  
  double compute(const arma::mat& X, 
                 const arma::vec&) const {
    
    arma::vec tmp = X.col(0);
    
    return median_cpp(tmp);
  }
};



// [[Rcpp::export]]
NumericVector median_boot_cpp(NumericVector x,
                              int    R     = 1000,
                              double alpha = 0.05,
                              int    seed  = -1) {

  return run_boot(
    vec_to_matrix(x),   // X: n x 1 (unused by MedianFn)
    dummy_vec(x.size()),// y: zeros  (unused by MedianFn)
    R,
    alpha,
    seed,
    MedianFn()
  );
}


