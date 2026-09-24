// Admissible exact two-sided confidence interval for the number M of
// successes in a finite population, X ~ Hyper(M, N - M, n), after
// Wang, W. (2015). Exact optimal confidence intervals for hypergeometric
// parameters. Journal of the American Statistical Association 110,
// 1491-1499.
//
// Starting from the two-sided Clopper-Pearson type interval (Konijn 1973),
// the integer limits are shrunk pairwise, U[k] = N - L[n - k], from the middle
// of the sample space {0, ..., n} outwards, each as far towards its one-sided
// Clopper-Pearson counterpart as the coverage allows. The interval for x
// depends only on the sample points between x and n - x.
//
// The parameter is discrete, so the construction needs no tolerances: raising
// L[k] from a - 1 to a (the inner lower limits follow to keep L monotone)
// removes exactly the value M = a - 1 from every interval that loses it, and
// the mirror value N - a + 1 on the other side. The coverage changes at these
// two points only, so one step is accepted iff the coverage there keeps the
// level. A failed step cannot be passed by a larger one, which would remove
// the same point from the same intervals, so the scan stops at the first
// failure and the result is exact. The same holds for the order L[j] <= U[j]
// and the monotonicity of L, which a larger step can only violate further.
//
// L and U are nondecreasing, the intervals covering M form a contiguous index
// range, and the coverage is a difference of two phyper() values.

#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <utility>
#include <vector>

namespace {

class WangHyper {
public:
  WangHyper(int n, int N, double alpha)
    : n_(n), N_(N), level_(1.0 - alpha), L_(n + 1), L1_(n + 1) {
    for (int k = 0; k <= n; ++k) {
      L_[k]  = cpLower(k, alpha / 2.0);   // two-sided CP type
      L1_[k] = cpLower(k, alpha);         // one-sided CP type
    }
  }

  std::pair<int, int> solve(int x) {
    if (n_ % 2 == 0) shrinkLower(n_ / 2, n_ / 2 + 1, n_);   // middle point
    const int kmin = std::min(x, n_ - x);
    for (int k = (n_ - 1) / 2; k >= kmin; --k) {
      shrinkLower(k, k + 1,      n_ - k - 1);
      shrinkUpper(k, n_ - k + 1, n_);
      Rcpp::checkUserInterrupt();
    }
    return std::make_pair(L_[x], U(x));
  }

private:
  const int n_, N_;
  const double level_;
  std::vector<int> L_, L1_;
  std::vector<std::pair<int, int> > undo_;

  int U(int k) const { return N_ - L_[n_ - k]; }

  // P(X >= x | M)
  double upperTail(int x, int M) const {
    return x <= 0 ? 1.0 : R::phyper(x - 1, M, N_ - M, n_, 0, 0);
  }

  // smallest M with P(X >= x | M) > a; the tail is nondecreasing in M and
  // equals 1 at M = N - n + x, so the search range is [x, N - n + x]
  int cpLower(int x, double a) const {
    if (x == 0) return 0;
    int lo = x, hi = N_ - n_ + x;
    while (lo < hi) {
      const int mid = lo + (hi - lo) / 2;
      if (upperTail(x, mid) > a) hi = mid; else lo = mid + 1;
    }
    return lo;
  }

  double coverage(int M) const {
    // jhi: last j with L[j] <= M;  jlo: first j with U[j] >= M, i.e. with
    // L[n - j] <= N - M
    const int jhi = int(std::upper_bound(L_.begin(), L_.end(), M) - L_.begin()) - 1;
    const int m   = int(std::upper_bound(L_.begin(), L_.end(), N_ - M) - L_.begin());
    const int jlo = n_ - m + 1;
    if (jlo > jhi) return 0.0;
    const double hi = R::phyper(jhi, M, N_ - M, n_, 1, 0);
    const double lo = jlo > 0 ? R::phyper(jlo - 1, M, N_ - M, n_, 1, 0) : 0.0;
    return hi - lo;
  }

  // a step removes M and its mirror N - M. The coverage is a rational number
  // and often hits the level exactly (e.g. 9/10 at n = 1, N = 10); as a
  // difference of two cdf values it carries rounding errors of order 1e-16,
  // hence the tolerance
  bool stepOk(int M) const {
    const double lev = level_ - 1e-12;
    return coverage(M) >= lev && coverage(N_ - M) >= lev;
  }

  // set L[k] = v and raise L[from..to] to at least v (keeps L monotone)
  void setLower(int k, int v, int from, int to) {
    undo_.clear();
    undo_.push_back(std::make_pair(k, L_[k]));
    L_[k] = v;
    for (int j = from; j <= to && L_[j] < v; ++j) {
      undo_.push_back(std::make_pair(j, L_[j]));
      L_[j] = v;
    }
  }
  // the step keeps the family valid: L[j] <= U[j] and L nondecreasing
  // around every changed index (the mirror interval gives the same
  // conditions). At levels of the usual size both hold by themselves; at low
  // levels the coverage alone admits crossed limits (hyper: n = 6, N = 7,
  // x = 3, 20 %) or a middle limit passing its neighbour (n = 13, N = 36, 20 %)
  bool valid() const {
    for (std::size_t i = 0; i < undo_.size(); ++i) {
      const int j = undo_[i].first;
      if (L_[j] > U(j)) return false;
      if (j > 0  && L_[j - 1] > L_[j]) return false;
      if (j < n_ && L_[j] > L_[j + 1]) return false;
    }
    return true;
  }

  void rollback() {
    for (std::size_t i = undo_.size(); i-- > 0; ) L_[undo_[i].first] = undo_[i].second;
    undo_.clear();
  }

  // raise the lower limit of point k towards the one-sided CP limit
  void shrinkLower(int k, int from, int to) {
    for (int a = L_[k] + 1; a <= L1_[k]; ++a) {
      setLower(k, a, from, to);
      if (!valid() || !stepOk(a - 1)) { rollback(); break; }
    }
  }

  // lower the upper limit of point k, i.e. raise L[n-k] to N - b
  void shrinkUpper(int k, int from, int to) {
    for (int b = U(k) - 1; b >= N_ - L1_[n_ - k]; --b) {
      setLower(n_ - k, N_ - b, from, to);
      if (!valid() || !stepOk(b + 1)) { rollback(); break; }
    }
  }
};

} // namespace

// x, n and N are taken as double, so that a non-integer count is rejected
// instead of being truncated silently by the conversion to int
// [[Rcpp::export(.wangHyperCI)]]
Rcpp::IntegerVector wangHyperCI(double x, double n, double N, double alpha) {
  if (!(N >= 1 && N < INT_MAX && N == std::floor(N)))
    Rcpp::stop("'N' must be a positive integer");
  if (!(n >= 1 && n <= N && n == std::floor(n)))
    Rcpp::stop("'n' must be an integer in [1, N]");
  if (!(x >= 0 && x <= n && x == std::floor(x)))
    Rcpp::stop("'x' must be an integer in [0, n]");
  if (!(alpha > 0 && alpha < 1))
    Rcpp::stop("'alpha' must lie in (0, 1)");
  WangHyper w(int(n), int(N), alpha);
  const std::pair<int, int> ci = w.solve(int(x));
  return Rcpp::IntegerVector::create(ci.first, ci.second);
}
