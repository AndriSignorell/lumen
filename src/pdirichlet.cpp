
#include <Rcpp.h>
#include <RcppParallel.h>
#include <random>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

struct DirichletCDFWorker : public Worker {

  const RVector<double> alpha;
  const RVector<double> q;
  const int k;
  const std::uint32_t seed;

  RVector<int> results;

  DirichletCDFWorker(const NumericVector& alpha,
                     const NumericVector& q,
                     const std::uint32_t seed,
                     IntegerVector& results)
    : alpha(alpha), q(q), k(alpha.size()), seed(seed), results(results) {}

  void operator()(std::size_t begin, std::size_t end) {

    // the R RNG is not thread-safe, so each range runs its own generator,
    // seeded from R's stream so that set.seed() governs the result
    std::mt19937 rng(seed + static_cast<std::uint32_t>(begin));

    // the gamma distributions carry state and are worth keeping
    std::vector<std::gamma_distribution<double> > gammas;
    gammas.reserve(k);
    for (int j = 0; j < k; j++) {
      gammas.push_back(std::gamma_distribution<double>(alpha[j], 1.0));
    }

    std::vector<double> x(k);

    for (std::size_t i = begin; i < end; i++) {

      double sum = 0.0;

      for (int j = 0; j < k; j++) {
        x[j] = gammas[j](rng);
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
double pdirichlet_cpp(const NumericVector& q,
                      const NumericVector& alpha,
                      const int n_sim) {

  // 'q' and 'alpha' are checked in pdirichlet(); this guards direct callers
  if (alpha.size() != q.size()) {
    stop("q and alpha must have same length");
  }

  const std::uint32_t seed =
    static_cast<std::uint32_t>(R::unif_rand() * 4294967295.0);

  IntegerVector results(n_sim);

  DirichletCDFWorker worker(alpha, q, seed, results);

  parallelFor(0, n_sim, worker);

  return mean(results);
}
