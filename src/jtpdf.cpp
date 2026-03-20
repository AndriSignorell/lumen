

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector jtpdf_cpp(int mxsum,
                    NumericVector cgsize,
                    NumericVector pdf0,
                    NumericVector pdf1) {
  
  int ng = cgsize.size();
  NumericVector pdf(mxsum);
  
  int i, j, g;
  int m, n, mn0, mn1;
  double di, dm, dn;
  
  // -----------------------------
    // initial: last two groups
  // -----------------------------
    m = cgsize[ng-2] - cgsize[ng-1];
    n = cgsize[ng-1];
    mn1 = m * n;
    
    dm = (double) m;
    dn = (double) n;
    
    for (i = 0; i <= mn1; i++) {
      di = (double) i;
      pdf[i] = R::dwilcox(di, dm, dn, false); // density
    }
    
    // -----------------------------
      // iterate backwards over groups
    // -----------------------------
      for (g = ng-3; g >= 0; g--) {
        
        // copy pdf -> pdf1 and reset pdf
        for (i = 0; i <= mn1; i++) {
          pdf1[i] = pdf[i];
          pdf[i] = 0.0;
        }
        
        // current MW for group g vs rest
        m = cgsize[g] - cgsize[g+1];
        n = cgsize[g+1];
        mn0 = m * n;
        
        dm = (double) m;
        dn = (double) n;
        
        for (i = 0; i <= mn0; i++) {
          di = (double) i;
          pdf0[i] = R::dwilcox(di, dm, dn, false);
        }
        
        // convolution
        for (i = 0; i <= mn0; i++) {
          for (j = 0; j <= mn1; j++) {
            pdf[i + j] += pdf0[i] * pdf1[j];
          }
        }
        
        mn1 = mn0 + mn1;
      }
    
    return pdf;
}


