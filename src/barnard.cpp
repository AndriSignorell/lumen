
// Barnard's unconditional exact test for 2 x 2 tables - compiled kernels.
//
// The algorithms follow the R package Exact (Peter Calhoun, GPL-2):
// the statistics, the CSM ordering and the multinomial model are reproduced
// step by step. Deliberate departures from Exact are documented in
// ?barnardTest; apart from those, results agree with Exact::exact.test().
// What changes is the cost:
//
//  - the CSM ordering keeps the probability of the current rejection region
//    on the nuisance grid and adds a candidate in O(npNumbers) instead of
//    re-summing the whole region, so that no stored ordering matrices
//    (ExactData) are needed;
//  - the multinomial model aggregates the region by its two margins, which
//    turns the double sum over tables and grid points into a small matrix
//    product.
//
// All values that Exact passes through signif(, 12) or round(, 12) are
// treated identically here (R::fprec, R::fround), as ties are decided on
// those rounded values.

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>

using namespace Rcpp;

namespace {

inline double sig12(double v) { return R::fprec(v, 12.0); }

// success probability of the first group; the grid reaches the admissible
// bounds exactly, where p + delta may miss 0 or 1 by a rounding error
inline double shifted(double p, double delta) {
  return std::min(1.0, std::max(0.0, p + delta));
}

enum Stat { Z_POOLED, Z_UNPOOLED, SANTNER, BOSCHLOO };

Stat statCode(const std::string& m) {
  if (m == "z-pooled")   return Z_POOLED;
  if (m == "z-unpooled") return Z_UNPOOLED;
  if (m == "santner-snell") return SANTNER;
  if (m == "boschloo")   return BOSCHLOO;
  Rcpp::stop("unknown method '%s'", m);
}


// score statistic; for delta != 0 with the restricted ML estimates of
// Miettinen and Nurminen / Farrington and Manning (closed-form cubic root)
double zPooled(int x, int y, int n1, int n2, double delta) {
  const double N = n1 + n2;
  const double p1 = double(x) / n1, p2 = double(y) / n2;
  const double num = p1 - p2 - delta;
  double den;

  if (delta == 0) {
    const double ph = (x + y) / N;
    den = std::sqrt(ph * (1 - ph) * (1.0 / n1 + 1.0 / n2));
  } else {
    const double theta = double(n2) / n1;
    const double a = 1 + theta;
    const double b = -(1 + theta + p1 + theta * p2 + delta * (theta + 2));
    const double c = delta * delta + delta * (2 * p1 + theta + 1) + p1 +
      theta * p2;
    const double d = -p1 * delta * (1 + delta);
    const double v = b * b * b / ((3 * a) * (3 * a) * (3 * a)) -
      b * c / (6 * a * a) + d / (2 * a);
    const double sgn = (v > 0) - (v < 0);
    const double u = sgn * std::sqrt(b * b / ((3 * a) * (3 * a)) - c / (3 * a));
    const double ratio = v / (u * u * u);
    double w = (M_PI + std::acos(ratio)) / 3;          // NaN outside [-1, 1]
    if (v == 0 && u == 0) w = (M_PI + std::acos(0.0)) / 3;
    if (ratio >= 1)       w = (M_PI + std::acos(1.0)) / 3;
    const double p1D = 2 * u * std::cos(w) - b / (3 * a);
    const double p2D = p1D - delta;
    den = std::sqrt(p2D * (1 - p2D) / n2 + p1D * (1 - p1D) / n1);
  }
  if (num == 0 && den == 0) return 0.0;
  return num / den;
}

double zUnpooled(int x, int y, int n1, int n2, double delta) {
  const double p1 = double(x) / n1, p2 = double(y) / n2;
  const double num = p1 - p2 - delta;
  const double den = std::sqrt(p2 * (1 - p2) / n2 + p1 * (1 - p1) / n1);
  if (num == 0 && den == 0) return 0.0;
  return num / den;
}

// Fisher's p-value as ordering statistic (Boschloo), as in stats::fisher.test
double fisherP(int x, int y, int n1, int n2, const std::string& alt) {
  const int k = x + y;
  if (alt == "less")
    return R::phyper(x, n1, n2, k, true, false);
  if (alt == "greater")
    return R::phyper(x - 1, n1, n2, k, false, false);

  const double relErr = 1 + 1e-7;
  const int lo = std::max(0, k - n2), hi = std::min(k, n1);
  std::vector<double> d(hi - lo + 1);
  double dmax = -std::numeric_limits<double>::infinity();
  for (int i = lo; i <= hi; ++i) {
    d[i - lo] = R::dhyper(i, n1, n2, k, true);
    dmax = std::max(dmax, d[i - lo]);
  }
  long double s = 0;
  for (double& di : d) { di = std::exp(di - dmax); s += di; }
  const double sd = double(s);
  for (double& di : d) di /= sd;
  const double ref = d[x - lo] * relErr;
  long double p = 0;
  for (double di : d) if (di <= ref) p += di;
  return double(p);
}

double statValue(Stat st, int x, int y, int n1, int n2, double delta,
                 const std::string& alt) {
  switch (st) {
  case Z_POOLED:   return zPooled(x, y, n1, n2, delta);
  case Z_UNPOOLED: return zUnpooled(x, y, n1, n2, delta);
  case SANTNER:    return double(x) / n1 - double(y) / n2 - delta;
  case BOSCHLOO:   return fisherP(x, y, n1, n2, alt);
  }
  return NA_REAL;
}

// binomial probabilities on the nuisance grid, A[i * K + k], B[j * K + k]
void binomLookup(int n1, int n2, const NumericVector& grid, double delta,
                 std::vector<double>& A, std::vector<double>& B) {
  const int K = grid.size();
  A.assign((n1 + 1) * K, 0.0);
  B.assign((n2 + 1) * K, 0.0);
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i <= n1; ++i)
      A[i * K + k] = R::dbinom(i, n1, shifted(grid[k], delta), false);
    for (int j = 0; j <= n2; ++j)
      B[j * K + k] = R::dbinom(j, n2, grid[k], false);
  }
}

} // namespace


// Statistic for every table with margins n1, n2 (rows x = 0..n1, columns
// y = 0..n2), rounded to 12 significant digits.
// [[Rcpp::export(.barnardStat)]]
NumericMatrix barnardStat(int n1, int n2, std::string method,
                          std::string alternative, double delta) {
  const Stat st = statCode(method);
  NumericMatrix TX(n1 + 1, n2 + 1);
  for (int i = 0; i <= n1; ++i)
    for (int j = 0; j <= n2; ++j)
      TX(i, j) = sig12(statValue(st, i, j, n1, n2, delta, alternative));
  return TX;
}


// Probability of the rejection region for every nuisance value in 'grid'
// (Exact's maxPvalue(), before rounding): sum over the region of
// dbinom(x, n1, p + delta) * dbinom(y, n2, p).
// [[Rcpp::export(.barnardProb)]]
NumericVector barnardProb(IntegerMatrix M, NumericVector grid, double delta) {
  const int n1 = M.nrow() - 1, n2 = M.ncol() - 1, K = grid.size();
  int maxX1 = -1, minX2 = n2 + 1;
  for (int i = 0; i <= n1; ++i)
    for (int j = 0; j <= n2; ++j)
      if (M(i, j) == 1) { maxX1 = std::max(maxX1, i); minX2 = std::min(minX2, j); }

  NumericVector prob(K);
  if (maxX1 < 0) return prob;

  std::vector<double> a(maxX1 + 1), b(n2 + 1);
  for (int k = 0; k < K; ++k) {
    const double p = grid[k];
    for (int i = 0; i <= maxX1; ++i) a[i] = R::dbinom(i, n1, shifted(p, delta), false);
    for (int j = minX2; j <= n2; ++j) b[j] = R::dbinom(j, n2, p, false);
    long double s = 0;
    for (int i = 0; i <= maxX1; ++i) {
      double row = 0;
      for (int j = minX2; j <= n2; ++j) if (M(i, j) == 1) row += b[j];
      s += a[i] * row;
    }
    prob[k] = double(s);
  }
  return prob;
}


// CSM ordering (Barnard 1947), Exact's moreExtremeCSM()/csmTemp() and, for
// a two-sided test with delta != 0, csmTemp2sidedDelta(). Starting from the
// most extreme table, the candidate on the boundary of the current region
// whose addition raises the maximal p-value least is added, until the
// observed table (x, y) has been reached. Returns the rejection region.
// [[Rcpp::export(.barnardCSM)]]
IntegerMatrix barnardCSM(int n1, int n2, int x, int y, NumericVector grid,
                         double delta, bool twoSided, bool doubleP) {
  const int K = grid.size();
  const bool deltaVariant = twoSided && delta != 0;
  const bool mirror = twoSided && !deltaVariant;
  const double mult = doubleP ? 2.0 : 1.0;

  std::vector<double> A, B;
  binomLookup(n1, n2, grid, delta, A, B);

  // tie-breaking statistic: z-pooled, -|z| for two-sided tests
  std::vector<double> tieTX((n1 + 1) * (n2 + 1));
  for (int i = 0; i <= n1; ++i)
    for (int j = 0; j <= n2; ++j) {
      double t = sig12(zPooled(i, j, n1, n2, delta));
      tieTX[i * (n2 + 1) + j] = twoSided ? -std::fabs(t) : t;
    }

  IntegerMatrix M(n1 + 1, n2 + 1);
  std::fill(M.begin(), M.end(), -1);
  std::vector<long double> cur(K, 0.0L);

  auto addToSum = [&](int i, int j) {
    for (int k = 0; k < K; ++k) cur[k] += A[i * K + k] * B[j * K + k];
  };
  auto finish = [&]() {
    for (auto& v : M) if (v == -1) v = 0;
    return M;
  };

  M(0, n2) = 1;
  if (n2 >= 1) M(0, n2 - 1) = 0;
  if (n1 >= 1) M(1, n2) = 0;
  if (twoSided) M(n1, 0) = (delta == 0) ? 1 : 0;

  for (int i = 0; i <= n1; ++i)
    for (int j = 0; j <= n2; ++j)
      if (M(i, j) == 1) {
        if (i == x && j == y) return finish();
        addToSum(i, j);
      }

  // maximal p-value of the current region plus the given tables, rounded
  // as in Exact (signif 12 per grid point, then round 12)
  auto candValue = [&](int i, int j, bool withMirror) {
    const int mi = n1 - i, mj = n2 - j;
    const bool self = (mi == i && mj == j);
    double best = R_NegInf;
    for (int k = 0; k < K; ++k) {
      long double s = cur[k] + (long double)(A[i * K + k] * B[j * K + k]);
      if (withMirror && !self) s += A[mi * K + k] * B[mj * K + k];
      const double v = sig12(mult * double(s));
      if (!ISNAN(v) && v > best) best = v;
    }
    return R::fround(best, 12.0);
  };

  std::vector<std::pair<int, int>> cand;
  std::vector<int> rowBest(n1 + 1), colBest(n2 + 1);

  // collect boundary candidates among the 0-cells: per row the extreme
  // column, then per column the extreme row
  auto collect = [&](bool lowerTri, bool useFlag, bool flagValue) {
    std::fill(rowBest.begin(), rowBest.end(), -1);
    for (int i = 0; i <= n1; ++i)
      for (int j = 0; j <= n2; ++j) {
        if (M(i, j) != 0) continue;
        if (useFlag && ((double(i) / n1 - double(j) / n2 < delta) != flagValue))
          continue;
        if (rowBest[i] < 0) rowBest[i] = j;
        else rowBest[i] = lowerTri ? std::min(rowBest[i], j)
                                   : std::max(rowBest[i], j);
      }
    std::fill(colBest.begin(), colBest.end(), -1);
    for (int i = 0; i <= n1; ++i) {
      const int j = rowBest[i];
      if (j < 0) continue;
      if (colBest[j] < 0) colBest[j] = i;
      else colBest[j] = lowerTri ? std::max(colBest[j], i)
                                 : std::min(colBest[j], i);
    }
    for (int j = 0; j <= n2; ++j)
      if (colBest[j] >= 0) cand.emplace_back(colBest[j], j);
  };

  for (long iter = 1; ; ++iter) {
    if ((iter & 1023) == 0) Rcpp::checkUserInterrupt();

    cand.clear();
    if (deltaVariant) {
      collect(true,  true, false);   // below the line p1 - p2 = delta
      collect(false, true, true);    // above
    } else {
      collect(false, false, false);
    }
    if (cand.empty())
      Rcpp::stop("CSM ordering ran out of candidate tables");

    std::vector<double> val(cand.size());
    double smallest = R_PosInf;
    for (size_t c = 0; c < cand.size(); ++c) {
      val[c] = candValue(cand[c].first, cand[c].second, mirror);
      smallest = std::min(smallest, val[c]);
    }

    std::vector<std::pair<int, int>> add;
    for (size_t c = 0; c < cand.size(); ++c)
      if (val[c] == smallest) add.push_back(cand[c]);

    if (add.size() > 1) {
      double tmin = R_PosInf;
      for (const auto& c : add) {
        const double t = tieTX[c.first * (n2 + 1) + c.second];
        if (!ISNAN(t)) tmin = std::min(tmin, t);
      }
      if (R_FINITE(tmin) || tmin == R_NegInf) {
        std::vector<std::pair<int, int>> keep;
        for (const auto& c : add)
          if (tieTX[c.first * (n2 + 1) + c.second] <= tmin) keep.push_back(c);
        add.swap(keep);
      }
    }

    for (const auto& c : add) {
      const int i = c.first, j = c.second;
      if (M(i, j) != 1) { M(i, j) = 1; addToSum(i, j); }
      if (mirror && M(n1 - i, n2 - j) != 1) {
        M(n1 - i, n2 - j) = 1;
        addToSum(n1 - i, n2 - j);
      }
      if (j >= 1 && M(i, j - 1) == -1) M(i, j - 1) = 0;
      if (i < n1 && M(i + 1, j) == -1) M(i + 1, j) = 0;
      if (deltaVariant) {
        if (j < n2 && M(i, j + 1) == -1) M(i, j + 1) = 0;
        if (i >= 1 && M(i - 1, j) == -1) M(i - 1, j) = 0;
      }
    }

    for (const auto& c : add)
      if ((c.first == x && c.second == y) ||
          (mirror && n1 - c.first == x && n2 - c.second == y))
        return finish();
  }
}


// Multinomial model: statistic for the observed table and the rejection
// region aggregated by its margins, W[r, c] = sum over the region tables
// with row-1 total r and column-1 total c of the multinomial coefficient
// divided by choose(N, r) * choose(N, c), so that
//   P(region | p1, p2) = sum W[r, c] dbinom(r, N, p1) dbinom(c, N, p2).
// [[Rcpp::export(.barnardMultinom)]]
List barnardMultinom(int N, int x11, int x12, int x21, std::string method,
                     std::string alternative) {
  const Stat st = statCode(method);

  auto stat = [&](int a, int b, int c) -> double {
    const int N1 = a + c;                 // column-1 total = first group
    if (N1 == 0 || N1 == N) return (st == BOSCHLOO) ? 1.0 : 0.0;
    return sig12(statValue(st, a, b, N1, N - N1, 0.0, alternative));
  };

  const double TXO = stat(x11, x12, x21);
  if (ISNAN(TXO)) Rcpp::stop("test statistic is NA");

  auto inRegion = [&](double t) -> bool {
    if (ISNAN(t)) return false;
    if (alternative == "less" || st == BOSCHLOO) return t <= TXO;
    if (alternative == "greater") return t >= TXO;
    return std::fabs(t) >= std::fabs(TXO);
  };

  std::vector<double> lf(N + 1);
  for (int i = 0; i <= N; ++i) lf[i] = std::lgamma(i + 1.0);

  std::vector<long double> W((N + 1) * (N + 1), 0.0L);
  for (int a = 0; a <= N; ++a) {
    for (int b = 0; a + b <= N; ++b) {
      for (int c = 0; a + b + c <= N; ++c) {
        if (!inRegion(stat(a, b, c))) continue;
        const int d = N - a - b - c, r = a + b, col = a + c;
        const double lw = lf[N] - lf[a] - lf[b] - lf[c] - lf[d] -
          (lf[N] - lf[r] - lf[N - r]) - (lf[N] - lf[col] - lf[N - col]);
        W[r * (N + 1) + col] += std::exp(lw);
      }
    }
    Rcpp::checkUserInterrupt();
  }

  NumericMatrix Wm(N + 1, N + 1);
  for (int r = 0; r <= N; ++r)
    for (int c = 0; c <= N; ++c) Wm(r, c) = double(W[r * (N + 1) + c]);

  return List::create(Named("statistic") = TXO, Named("W") = Wm);
}


// P(region | p1, p2) for all combinations of p1 (rows) and p2 (columns)
// [[Rcpp::export(.barnardMultinomProb)]]
NumericMatrix barnardMultinomProb(NumericMatrix W, NumericVector p1,
                                  NumericVector p2) {
  const int N = W.nrow() - 1;
  const int K1 = p1.size(), K2 = p2.size();

  // V[r, k2] = sum_c W[r, c] dbinom(c, N, p2[k2])
  std::vector<double> bc(N + 1);
  NumericMatrix V(N + 1, K2);
  for (int k = 0; k < K2; ++k) {
    for (int c = 0; c <= N; ++c) bc[c] = R::dbinom(c, N, p2[k], false);
    for (int r = 0; r <= N; ++r) {
      long double s = 0;
      for (int c = 0; c <= N; ++c) s += W(r, c) * bc[c];
      V(r, k) = double(s);
    }
  }

  NumericMatrix P(K1, K2);
  std::vector<double> ar(N + 1);
  for (int k1 = 0; k1 < K1; ++k1) {
    for (int r = 0; r <= N; ++r) ar[r] = R::dbinom(r, N, p1[k1], false);
    for (int k2 = 0; k2 < K2; ++k2) {
      long double s = 0;
      for (int r = 0; r <= N; ++r) s += ar[r] * V(r, k2);
      P(k1, k2) = double(s);
    }
  }
  return P;
}
