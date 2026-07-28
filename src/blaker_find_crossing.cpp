#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;


// ======================================================
// internal helper: acceptability function for Blaker CI
//
// every tail probability is evaluated on the tail it belongs to
// (lower_tail = 0 where needed). Computing an upper tail as
// 1 - pbinom(.) cancels catastrophically for large n: as soon as the
// complement rounds to 1, the whole tail is lost, and qbinom(1 - p2)
// then silently returns n.

inline double accept_bin(int x, int n, double p) {

  if (p <= 0.0) return (x == 0) ? 1.0 : 0.0;
  if (p >= 1.0) return (x == n) ? 1.0 : 0.0;

  const double p1 = R::pbinom(x - 1, n, p, 0, 0);   // P(X >= x)
  const double p2 = R::pbinom(x,     n, p, 1, 0);   // P(X <= x)

  const double q1 = R::qbinom(p1, n, p, 1, 0);      // F^-1(p1)
  const double q2 = R::qbinom(p2, n, p, 0, 0);      // F^-1(1 - p2)

  const double a1 = p1 + R::pbinom(q1 - 1, n, p, 1, 0);
  const double a2 = p2 + R::pbinom(q2,     n, p, 0, 0);

  return std::min(std::min(a1, a2), 1.0);

}


// ======================================================
// smallest (from_left = true) resp. largest (from_left = false) p in
// [lo, hi] with accept_bin(x, n, p) >= alpha.
//
// The caller is expected to supply a bracket that actually contains
// the crossing, i.e. [clopper-pearson limit, x/n] resp. [x/n,
// clopper-pearson limit]. The coarse scan below is only a safety net
// for degenerate brackets, never the working path.

// [[Rcpp::export]]
double blaker_find_crossing(
    int    x,
    int    n,
    double alpha,
    double lo,
    double hi,
    bool   from_left,
    double tol        = 1e-12,   // relative
    int    safe_steps = 2000
) {

  constexpr double eps = 1e-14;

  if (!(hi > lo))                       // degenerate / inverted bracket
    return from_left ? hi : lo;

  double left  = lo;
  double right = hi;

  // evaluate endpoints once
  const double lo_val = accept_bin(x, n, lo);
  const double hi_val = accept_bin(x, n, hi);

  const bool lo_above = lo_val >= alpha - eps;
  const bool hi_above = hi_val >= alpha - eps;

  const double step = (hi - lo) / safe_steps;
  int i = 1;

  if (from_left) {

    if (lo_above)
      return lo;

    if (!hi_above) {
      // no bracket - coarse scan left to right, stays inside [lo, hi]
      for (; i <= safe_steps; i++)
        if (accept_bin(x, n, lo + i * step) >= alpha - eps) break;

      if (i > safe_steps)
        return hi;                      // conservative fallback

      left  = lo + (i - 1) * step;
      right = lo + i * step;
    }

  } else {

    if (hi_above)
      return hi;

    if (!lo_above) {
      // no bracket - coarse scan right to left
      for (; i <= safe_steps; i++)
        if (accept_bin(x, n, hi - i * step) >= alpha - eps) break;

      if (i > safe_steps)
        return lo;                      // conservative fallback

      left  = hi - i * step;
      right = hi - (i - 1) * step;
    }

  }

  // bisection refinement, relative stopping rule
  for (int k = 0; k < 200; k++) {

    const double mid = 0.5 * (left + right);

    if (mid <= left || mid >= right)    // machine precision reached
      break;

    if (accept_bin(x, n, mid) >= alpha - eps) {
      if (from_left) right = mid; else left = mid;
    } else {
      if (from_left) left = mid; else right = mid;
    }

    if (right - left <= tol * std::max(right, 1e-12))
      break;

  }

  return from_left ? right : left;

}
