
#include <Rcpp.h>
using namespace Rcpp;

// log version of druns_nom for a single r
double druns_nom_log(int r, int n1, int n2) {
  
  if (r % 2 == 0) {
    // even: r = 2k
    int k = r / 2;
    return std::log(2.0) +
      R::lchoose(n1 - 1, k - 1) +
      R::lchoose(n2 - 1, k - 1);
    
  } else {
    // odd: r = 2k+1
    int k = (r - 1) / 2;
    
    double a = R::lchoose(n1 - 1, k - 1) + R::lchoose(n2 - 1, k);
    double b = R::lchoose(n1 - 1, k)     + R::lchoose(n2 - 1, k - 1);
    
    // log(exp(a)+exp(b)) trick
    double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
  }
}


// [[Rcpp::export]]
double pruns_rcpp(int r, int n1, int n2, std::string alternative = "two.sided") {
  
  int n = n1 + n2;
  
  if (r <= 1)
    stop("Number of runs must be > 1");
  if (r > n)
    stop("Number of runs must be < n1+n2");
  if (n1 < 1 || n2 < 1)
    return 0.0;
  
  double E = 1.0 + 2.0 * n1 * n2 / n;
  
  double log_denom = R::lchoose(n, n1);
  
  int rmax = (n1 == n2) ? 2 * n1 : 2 * std::min(n1, n2) + 1;
  
  double pL = 0.0;
  double pU = 0.0;
  double p2 = 0.0;
  
  // accumulate in log-space (via exp differences)
  for (int rv = 2; rv <= rmax; ++rv) {
    
    double logp = druns_nom_log(rv, n1, n2) - log_denom;
    double p = std::exp(logp);
    
    if (rv <= r) pL += p;
    if (rv <= r - 1) pU += p;
    if (std::abs(rv - E) >= std::abs(r - E)) p2 += p;
  }
  
  pU = 1.0 - pU;
  
  if (alternative == "less") return pL;
  if (alternative == "greater") return pU;
  
  return p2;
}