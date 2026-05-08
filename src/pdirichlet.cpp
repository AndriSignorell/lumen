

#include <Rcpp.h>
#include <RcppParallel.h>
#include <random>

using namespace Rcpp;
using namespace RcppParallel;

struct DirichletCDFWorker : public Worker {
  
  const RVector<double> alpha;
  const RVector<double> q;
  int k;
  
  RVector<int> results;
  
  DirichletCDFWorker(const NumericVector& alpha,
                     const NumericVector& q,
                     IntegerVector& results)
    : alpha(alpha), q(q), k(alpha.size()), results(results) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    // thread-local RNG
    std::mt19937 rng(std::random_device{}());
    
    for (std::size_t i = begin; i < end; i++) {
      
      std::vector<double> x(k);
      double sum = 0.0;
      
      for (int j = 0; j < k; j++) {
        std::gamma_distribution<double> g(alpha[j], 1.0);
        x[j] = g(rng);
        sum += x[j];
      }
      
      bool ok = true;
      for (int j = 0; j < k; j++) {
        x[j] /= sum;
        if (x[j] > q[j]) {
          ok = false;
          break;
        }
      }
      
      results[i] = ok;
    }
  }
};


// [[Rcpp::export]]
double pdirichlet_cpp(NumericVector q,
                               NumericVector alpha,
                               int n_sim) {
  
  if (alpha.size() != q.size()) {
    stop("q and alpha must have same length");
  }
  
  IntegerVector results(n_sim);
  
  DirichletCDFWorker worker(alpha, q, results);
  
  parallelFor(0, n_sim, worker);
  
  return mean(results);
}
