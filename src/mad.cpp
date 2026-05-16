
// ============================================================
  // mad_boot.cpp
//
  // Bootstrap confidence interval for the Median Absolute
// Deviation (MAD).
//
  // The MAD is defined as
//
  //   MAD(x) = constant * median( |x - median(x)| )
  //
    // which matches base-R mad() with the same default constant
  // (1.4826, making the estimator consistent for sigma under
      // normality).
  //
    // Uses the generic parallel bootstrap framework.
  //
    // R interface:
    //   mad_boot_cpp(x, B = 1000, alpha = 0.05,
                      //                constant = 1.4826, seed = -1,
                      //                method = "perc")
  //
    // ============================================================
    
    #include "boot_framework.h"
    
    // ============================================================
    // mad_cpp: replicates base-R mad(x, constant = constant)
  // ============================================================
    
    double mad_cpp(const arma::vec& x, double constant = 1.4826) {
      
      size_t n = x.n_elem;
      
      if (n == 0)
        return NA_REAL;
      
      // --- median of x -----------------------------------------
        arma::vec tmp = x;           // working copy (nth_element mutates)
        double*   beg = tmp.memptr();
        size_t    mid = n / 2;
        
        std::nth_element(beg, beg + mid, beg + n);
        double med = beg[mid];
        
        if (n % 2 == 0) {
          double lower = *std::max_element(beg, beg + mid);
          med = 0.5 * (med + lower);
        }
        
        // --- absolute deviations from median ----------------------
          arma::vec dev(n);
        for (size_t i = 0; i < n; i++)
          dev[i] = std::abs(x[i] - med);
        
        // --- median of deviations ---------------------------------
          double* dbeg = dev.memptr();
        std::nth_element(dbeg, dbeg + mid, dbeg + n);
        double mdev = dbeg[mid];
        
        if (n % 2 == 0) {
          double lower = *std::max_element(dbeg, dbeg + mid);
          mdev = 0.5 * (mdev + lower);
        }
        
        return constant * mdev;
    }
  
  
  // ============================================================
    // StatFn struct — captures the constant at construction time
  // ============================================================
    
    struct MadFn {
      
      double constant;
      
      explicit MadFn(double c = 1.4826) : constant(c) {}
      
      double compute(const arma::mat& X,
                     const arma::vec& /* y */) const {
                       
                       arma::vec tmp = X.col(0);
                       return mad_cpp(tmp, constant);
                     }
    };
  
  
  // ============================================================
    // Exported R function
  // ============================================================
    
    // [[Rcpp::export]]
  NumericVector mad_boot_cpp(NumericVector x,
                             int    R        = 1000,
                             double alpha    = 0.05,
                             double constant = 1.4826,
                             int    seed     = -1,
                             String method   = "perc") {
    
    return run_boot(
      vec_to_matrix(x),
      dummy_vec(x.size()),
      R,
      alpha,
      seed,
      MadFn(constant),
      std::string(method.get_cstring())
    );
  }
  
  
  
  
  // ============================================================
  // mad_two_sample_boot.cpp
  //
  // Bootstrap confidence intervals for two-sample MAD statistics:
  //
  //   mad_diff_boot_cpp   — MAD(x) - MAD(y)
  //   mad_ratio_boot_cpp  — (MAD(x) / MAD(y))^2
  //
  // Both statistics operate on a two-column matrix X:
  //   col 0 → x sample
  //   col 1 → y sample
  //
  // Resampling is done independently within each column (same
  // bootstrap index is applied to both, which is correct because
  // the two samples are independent and of potentially different
  // length — see two_col_matrix() below).
  //
  // Uses the generic parallel bootstrap framework.
  //
  // R interfaces:
  //   mad_diff_boot_cpp(x, y, R = 1000, alpha = 0.05,
  //                    constant = 1.4826, seed = -1,
  //                    method = "perc")
  //
  //   mad_ratio_boot_cpp(x, y, R = 1000, alpha = 0.05,
  //                     constant = 1.4826, seed = -1,
  //                     method = "perc")
  //
  // ============================================================
  
#include "boot_framework.h"
  
  // ============================================================
  // mad_cpp (local copy — same as in mad_boot.cpp)
  // ============================================================
  
  static double mad_cpp2(const arma::vec& x, double constant) {
    
    size_t n = x.n_elem;
    if (n == 0) return NA_REAL;
    
    arma::vec tmp = x;
    double*   beg = tmp.memptr();
    size_t    mid = n / 2;
    
    std::nth_element(beg, beg + mid, beg + n);
    double med = beg[mid];
    if (n % 2 == 0) {
      double lower = *std::max_element(beg, beg + mid);
      med = 0.5 * (med + lower);
    }
    
    arma::vec dev(n);
    for (size_t i = 0; i < n; i++)
      dev[i] = std::abs(x[i] - med);
    
    double* dbeg = dev.memptr();
    std::nth_element(dbeg, dbeg + mid, dbeg + n);
    double mdev = dbeg[mid];
    if (n % 2 == 0) {
      double lower = *std::max_element(dbeg, dbeg + mid);
      mdev = 0.5 * (mdev + lower);
    }
    
    return constant * mdev;
  }
  
  
  // ============================================================
  // Two-sample bootstrap layout
  //
  // The framework resamples rows of a single matrix X.  For two
  // independent samples of (potentially) unequal length we stack
  // them into an n_max x 2 matrix padded with NA, and let each
  // StatFn extract col 0 / col 1 and drop NAs before computing.
  //
  // two_col_matrix(): build the padded matrix from x and y.
  // extract_col():    pull a column and strip NA padding.
  // ============================================================
  
  inline NumericMatrix two_col_matrix(NumericVector x, NumericVector y) {
    
    int nx = x.size();
    int ny = y.size();
    int n  = std::max(nx, ny);
    
    NumericMatrix m(n, 2);
    std::fill(m.begin(), m.end(), NA_REAL);
    
    for (int i = 0; i < nx; i++) m(i, 0) = x[i];
    for (int i = 0; i < ny; i++) m(i, 1) = y[i];
    
    return m;
  }
  
  // extract non-NA values from column j of an arma::mat
  inline arma::vec extract_col(const arma::mat& X, arma::uword j) {
    
    arma::vec col = X.col(j);
    // count non-NA
    arma::uword k = 0;
    for (arma::uword i = 0; i < col.n_elem; i++){
      if (!std::isnan(col[i])) k++;
    }
      
      arma::vec out(k);
      arma::uword idx = 0;
      for (arma::uword i = 0; i < col.n_elem; i++){
        if (!std::isnan(col[i])) out[idx++] = col[i];
      }
        
        return out;
  }
  
  
  // ============================================================
  // StatFn: MAD(x) - MAD(y)
  // ============================================================
  
  struct MadDiffFn {
    
    double constant;
    explicit MadDiffFn(double c = 1.4826) : constant(c) {}
    
    double compute(const arma::mat& X,
                   const arma::vec& /* y */) const {
      
      arma::vec xv = extract_col(X, 0);
      arma::vec yv = extract_col(X, 1);
      
      double mx = mad_cpp2(xv, constant);
      double my = mad_cpp2(yv, constant);
      
      if (std::isnan(mx) || std::isnan(my)) return NA_REAL;
      
      return mx - my;
    }
  };
  
  
  // ============================================================
  // StatFn: (MAD(x) / MAD(y))^2
  // ============================================================
  
  struct MadRatioFn {
    
    double constant;
    explicit MadRatioFn(double c = 1.4826) : constant(c) {}
    
    double compute(const arma::mat& X,
                   const arma::vec& /* y */) const {
      
      arma::vec xv = extract_col(X, 0);
      arma::vec yv = extract_col(X, 1);
      
      double mx = mad_cpp2(xv, constant);
      double my = mad_cpp2(yv, constant);
      
      if (std::isnan(mx) || std::isnan(my) || my == 0.0) return NA_REAL;
      
      return (mx / my) * (mx / my);
    }
  };
  
  
  // ============================================================
  // Exported R functions
  // ============================================================
  
  // [[Rcpp::export]]
  NumericVector mad_diff_boot_cpp(NumericVector x,
                                  NumericVector y,
                                  int    R        = 1000,
                                  double alpha    = 0.05,
                                  double constant = 1.4826,
                                  int    seed     = -1,
                                  String method   = "perc") {
    
    return run_boot(
      two_col_matrix(x, y),
      dummy_vec(std::max(x.size(), y.size())),
      R,
      alpha,
      seed,
      MadDiffFn(constant),
      std::string(method.get_cstring())
    );
  }
  
  
  // [[Rcpp::export]]
  NumericVector mad_ratio_boot_cpp(NumericVector x,
                                   NumericVector y,
                                   int    R        = 1000,
                                   double alpha    = 0.05,
                                   double constant = 1.4826,
                                   int    seed     = -1,
                                   String method   = "perc") {
    
    return run_boot(
      two_col_matrix(x, y),
      dummy_vec(std::max(x.size(), y.size())),
      R,
      alpha,
      seed,
      MadRatioFn(constant),
      std::string(method.get_cstring())
    );
  }