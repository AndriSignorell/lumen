// Admissible exact two-sided confidence interval for a binomial proportion
// after Wang, W. (2014). An iterative construction of confidence intervals
// for a proportion. Statistica Sinica 24, 1389-1410.
//
// Starting from the two-sided Clopper-Pearson interval, the limits are shrunk
// pairwise from the middle of the sample space {0, ..., n} outwards; each
// limit is moved by bisection as far towards its one-sided Clopper-Pearson
// counterpart as the infimum coverage allows. The interval for x depends only
// on the sample points between x and n - x, so the construction stops there.
//
// Storage: only the lower limits L[0..n] are kept; the upper limits follow from
// the symmetry U[k] = 1 - L[n - k]. Both sequences are nondecreasing, so the
// set of intervals covering p is a contiguous index range and the coverage is
// a difference of two binomial cdf values (O(log n) instead of O(n)).
//
// A tentative move is accepted only if the one-sided limits of the coverage
// function hold the level at every limit inside the changed region; they are
// evaluated at the breakpoint itself, never at a shifted p, so arbitrarily
// narrow pieces between two limits are covered as well.

#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <utility>
#include <vector>

namespace {

class WangBinom {
public:
  WangBinom(int n, double alpha, double tol)
    : n_(n), level_(1.0 - alpha), tol_(tol), L_(n + 1), L1_(n + 1) {
    for (int k = 1; k <= n; ++k) {
      L_[k]  = R::qbeta(alpha / 2.0, k, n - k + 1.0, 1, 0);  // two-sided CP
      L1_[k] = R::qbeta(alpha,       k, n - k + 1.0, 1, 0);  // one-sided CP
    }
  }

  // limits for x, computing only the sample points needed
  std::pair<double, double> solve(int x) {
    if (n_ % 2 == 0) shrinkLower(n_ / 2, n_ / 2 + 1, n_);   // middle point
    const int kmin = std::min(x, n_ - x);
    for (int k = (n_ - 1) / 2; k >= kmin; --k) {
      shrinkLower(k, k + 1,      n_ - k - 1);
      shrinkUpper(k, n_ - k + 1, n_);
    }
    return std::make_pair(L_[x], U(x));
  }

private:
  const int n_;
  const double level_, tol_;
  std::vector<double> L_, L1_;
  std::vector<std::pair<int, double> > undo_;

  double U(int k) const { return 1.0 - L_[n_ - k]; }

  bool converged(double lo, double hi) const {
    return hi - lo <= tol_ * std::max(std::fabs(hi), 1e-300);
  }

  // first j with U(j) > p (strict) or U(j) >= p; U is nondecreasing and is
  // compared as computed, so the limits are consistent with the breakpoints
  int firstU(double p, bool strict) const {
    int lo = 0, hi = n_ + 1;
    while (lo < hi) {
      const int mid = lo + (hi - lo) / 2;
      const double u = U(mid);
      if (strict ? u > p : u >= p) hi = mid; else lo = mid + 1;
    }
    return lo;
  }

  // one-sided limits of the coverage function at p, evaluated at p itself.
  // left:  intervals with L[j] <  p <= U[j]  (coverage at p - 0)
  // right: intervals with L[j] <= p <  U[j]  (coverage at p + 0)
  // The covering indices form a contiguous range because L and U are
  // nondecreasing; the probability of a contiguous range is unimodal in p,
  // so its minimum on a piece between breakpoints lies at a piece end.
  double coverage(double p, bool left) const {
    if (p <= 0.0 || p >= 1.0) return 1.0;
    const int jhi = int((left ? std::lower_bound(L_.begin(), L_.end(), p)
                              : std::upper_bound(L_.begin(), L_.end(), p))
                        - L_.begin()) - 1;
    const int jlo = firstU(p, !left);
    if (jlo > jhi) return 0.0;
    const double hi = R::pbinom(jhi, n_, p, 1, 0);
    const double lo = jlo > 0 ? R::pbinom(jlo - 1, n_, p, 1, 0) : 0.0;
    return hi - lo;
  }

  // infimum of the coverage on (lo, hi): right limit at lo, left limit at hi,
  // both limits at every limit L[j], U[j] inside. A tentative move changes the
  // coverage on one region and its mirror image; symmetry holds only up to the
  // rounding of 1 - L, so both regions are checked
  bool coverageOk(double lo, double hi) const {
    if (coverage(lo, false) < level_ || coverage(hi, true) < level_)
      return false;
    for (int j = int(std::upper_bound(L_.begin(), L_.end(), lo) - L_.begin());
         j <= n_ && L_[j] < hi; ++j)
      if (coverage(L_[j], true) < level_ || coverage(L_[j], false) < level_)
        return false;
    for (int j = firstU(lo, true); j <= n_ && U(j) < hi; ++j) {
      const double u = U(j);
      if (coverage(u, true) < level_ || coverage(u, false) < level_)
        return false;
    }
    return true;
  }

  // raise L[from..to] to at least v (keeps L monotone), remember old values
  void raise(int from, int to, double v) {
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

  // raise the lower limit of point k (and lower U[n-k] symmetrically) by
  // bisection between the current, verified value a0 and the one-sided CP
  // limit a1; the coverage changes on (a0, a) only
  void shrinkLower(int k, int from, int to) {
    double a0 = L_[k], a1 = L1_[k];
    while (a1 > a0 && !converged(a0, a1)) {
      const double a = 0.5 * (a0 + a1);
      if (a <= a0 || a >= a1) break;                    // no representable midpoint
      undo_.clear();
      undo_.push_back(std::make_pair(k, L_[k]));
      L_[k] = a;
      raise(from, to, a);
      if (valid() && coverageOk(a0, a) && coverageOk(1.0 - a, 1.0 - a0)) { a0 = a; undo_.clear(); } else { a1 = a; rollback(); }
    }
  }

  // lower the upper limit of point k, i.e. raise L[n-k] to 1 - b, between the
  // one-sided CP limit b0 and the current, verified value b1; the coverage
  // changes on (b, b1) only
  void shrinkUpper(int k, int from, int to) {
    double b0 = 1.0 - L1_[n_ - k], b1 = U(k);
    while (b1 > b0 && !converged(b0, b1)) {
      const double b = 0.5 * (b0 + b1);
      if (b <= b0 || b >= b1) break;
      undo_.clear();
      undo_.push_back(std::make_pair(n_ - k, L_[n_ - k]));
      L_[n_ - k] = 1.0 - b;
      raise(from, to, 1.0 - b);
      if (valid() && coverageOk(U(k), b1) && coverageOk(1.0 - b1, L_[n_ - k])) { b1 = U(k); undo_.clear(); } else { b0 = b; rollback(); }
    }
  }
};

} // namespace

// x and n are taken as double, so that a non-integer count is rejected
// instead of being truncated silently by the conversion to int
// [[Rcpp::export(.wangBinomCI)]]
Rcpp::NumericVector wangBinomCI(double x, double n, double alpha,
                                double tol = 1e-10) {
  if (!(n >= 1 && n < INT_MAX && n == std::floor(n)))   // n + 1 must fit
    Rcpp::stop("'n' must be a positive integer");
  if (!(x >= 0 && x <= n && x == std::floor(x)))
    Rcpp::stop("'x' must be an integer in [0, n]");
  if (!(alpha > 0 && alpha < 1))
    Rcpp::stop("'alpha' must lie in (0, 1)");
  if (!(tol > 0 && tol < 1))
    Rcpp::stop("'tol' must lie in (0, 1)");
  WangBinom w(int(n), alpha, tol);
  const std::pair<double, double> ci = w.solve(int(x));
  return Rcpp::NumericVector::create(ci.first, ci.second);
}
