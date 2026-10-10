
#include <Rcpp.h>
#include <RcppParallel.h>
#include <random>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

// draws are cut into fixed blocks of kBlock; block b always gets the
// generator seeded with seed + b, so the result depends only on the seed
// and n_sim, never on how the scheduler splits the work across threads
static const std::size_t kBlock = 4096;

struct DirichletCDFWorker : public Worker {

  const RVector<double> alpha;
  const RVector<double> q;
  const int k;
  const std::uint32_t seed;
  const std::size_t n_sim;

  RVector<int> hits;   // one count per block

  DirichletCDFWorker(const NumericVector& alpha,
                     const NumericVector& q,
                     const std::uint32_t seed,
                     const std::size_t n_sim,
                     IntegerVector& hits)
    : alpha(alpha), q(q), k(alpha.size()), seed(seed), n_sim(n_sim),
      hits(hits) {}

  void operator()(std::size_t begin, std::size_t end) {

    std::vector<double> x(k);

    for (std::size_t b = begin; b < end; b++) {

      // the R RNG is not thread-safe, so each block runs its own generator,
      // seeded from R's stream so that set.seed() governs the result
      std::mt19937 rng(seed + static_cast<std::uint32_t>(b));

      // fresh distributions per block: they carry state between draws
      std::vector<std::gamma_distribution<double> > gammas;
      gammas.reserve(k);
      for (int j = 0; j < k; j++) {
        gammas.push_back(std::gamma_distribution<double>(alpha[j], 1.0));
      }

      const std::size_t lo = b * kBlock;
      const std::size_t hi = std::min(lo + kBlock, n_sim);
      int cnt = 0;

      for (std::size_t i = lo; i < hi; i++) {

        double sum = 0.0;
        for (int j = 0; j < k; j++) {
          x[j] = gammas[j](rng);
          sum += x[j];
        }

        bool ok = true;
        for (int j = 0; j < k; j++) {
          if (x[j] / sum > q[j]) {
            ok = false;
            break;
          }
        }
        cnt += ok;
      }

      hits[b] = cnt;
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

  const std::size_t n = static_cast<std::size_t>(n_sim);
  const std::size_t n_block = (n + kBlock - 1) / kBlock;

  IntegerVector hits(n_block);

  DirichletCDFWorker worker(alpha, q, seed, n, hits);

  parallelFor(0, n_block, worker);

  double total = 0.0;
  for (std::size_t b = 0; b < n_block; b++) total += hits[b];

  return total / static_cast<double>(n);
}
