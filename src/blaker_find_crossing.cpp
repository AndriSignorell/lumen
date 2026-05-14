
#include <Rcpp.h>
using namespace Rcpp;


#include <Rcpp.h>
using namespace Rcpp;

// ======================================================
// internal helper: acceptability function for Blaker CI

inline double accept_bin(int x, int n, double p) {
  
  double p1 = 1.0 - R::pbinom(x - 1, n, p, 1, 0);
  double p2 = R::pbinom(x, n, p, 1, 0);
  double q1 = R::qbinom(p1, n, p, 1, 0);
  double q2 = R::qbinom(1.0 - p2, n, p, 1, 0);
  double a1 = p1 + R::pbinom(q1 - 1, n, p, 1, 0);
  double a2 = p2 + 1.0 - R::pbinom(q2, n, p, 1, 0);
  
  return std::min(a1, a2);
  
}

// ======================================================
// [[Rcpp::export]]
double blaker_find_crossing(
    int    x,
    int    n,
    double alpha,
    double lo,
    double hi,
    bool   from_left,
    double tol        = 1e-10,
    int    safe_steps = 50
) {
  
  constexpr double eps = 1e-14;
  
  double left  = lo;
  double right = hi;
  
  // evaluate endpoints once
  double lo_val = accept_bin(x, n, lo);
  double hi_val = accept_bin(x, n, hi);
  
  bool lo_above = lo_val >= alpha - eps;
  bool hi_above = hi_val >= alpha - eps;
  
  if (from_left) {
    
    // looking for smallest p where accept_bin >= alpha
    if (lo_above)
      return lo;
    
    if (!hi_above) {
      // no bracket found – coarse scan left to right
      double step = (hi - lo) / safe_steps;
      double p    = lo;
      double val  = lo_val;
      
      while (p <= hi) {
        p   += step;
        val  = accept_bin(x, n, p);
        if (val >= alpha - eps) break;
      }
      
      if (p > hi)
        // no crossing found – conservative fallback
        return hi;
      
      left  = p - step;
      right = p;
    }
    
  } else {
    
    // looking for largest p where accept_bin >= alpha
    if (hi_above)
      return hi;
    
    if (!lo_above) {
      // no bracket found – coarse scan right to left
      double step = (hi - lo) / safe_steps;
      double p    = hi;
      double val  = hi_val;
      
      while (p >= lo) {
        p   -= step;
        val  = accept_bin(x, n, p);
        if (val >= alpha - eps) break;
      }
      
      if (p < lo)
        // no crossing found – conservative fallback
        return lo;
      
      left  = p;
      right = p + step;
    }
    
  }
  
  // bisection refinement
  for (int i = 0; i < 80; i++) {
    
    double mid = 0.5 * (left + right);
    double val = accept_bin(x, n, mid);
    
    if (from_left) {
      if (val >= alpha - eps) right = mid; else left = mid;
    } else {
      if (val >= alpha - eps) left  = mid; else right = mid;
    }
    
    if (std::abs(right - left) < tol) break;
    
  }
  
  return from_left ? right : left;
  
}

