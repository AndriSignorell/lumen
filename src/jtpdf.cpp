
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector jtpdf_cpp(IntegerVector gsize) {
  
  int ng = gsize.size();
  int N  = sum(gsize);
  
  IntegerVector cs(ng);
  cs[0] = gsize[0];
  for (int i = 1; i < ng; i++) cs[i] = cs[i-1] + gsize[i];
  
  int maxJ = 0;
  for (int i = 0; i < ng - 1; i++)
    maxJ += gsize[i] * (N - cs[i]);
  
  bool even  = (maxJ % 2 == 0);
  int  upper = even ? maxJ / 2 : (maxJ - 1) / 2;
  
  double num_comb = 1.0;
  int rem = N;
  for (int i = 0; i < ng; i++) {
    for (int j = 1; j <= gsize[i]; j++) {
      num_comb *= (double)rem-- / (double)j;
    }
  }

  std::vector<double> freq(upper + 1, 0.0);
  freq[0] = 1.0;
  
  auto update = [&](int m, int n) {
    
    if ((n + 1) <= upper) {
      int p = std::min(m + n, upper);
      for (int t = n + 1; t <= p; t++)
        for (int u = upper; u >= t; u--)   // downwards
          freq[u] -= freq[u - t];
    }
    
    int q = std::min(m, upper);
    for (int s = 1; s <= q; s++)
      for (int u = s; u <= upper; u++)     // upwards
        freq[u] += freq[u - s];
  };
  
  for (int i = 0; i < ng - 1; i++)
    update(gsize[i], N - cs[i]);
  
  NumericVector prob(maxJ + 1);

  for (int i = 0; i <= upper; i++)
    prob[i] = freq[i] / num_comb;
  
  for (int i = 0; i <= upper; i++)          
    prob[maxJ - i] = prob[i];
  
  if (even)
    prob[upper] = freq[upper] / num_comb; // mean only with even maxJ
  
  
  return prob;
}

