// Brunner-Munzel test: core statistic and studentized permutation test.
//
// Speed notes
// -----------
// The pooled midranks are invariant under permutation of the group labels:
// permuting labels reshuffles which pooled ranks land in which group, but never
// changes the pooled ranks themselves. Enumerating subsets of the *sorted*
// pooled-rank vector therefore yields both groups already in ascending order,
// so the within-group midranks needed by the studentized statistic can be read
// off in a single O(n) tie-run scan. No sorting and no allocation happen inside
// the permutation loop; the cost per replicate is O(N) rather than O(N log N).

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// ---------------------------------------------------------------- helpers --

// Midranks of an ascending-sorted vector, written into `out`.
static inline void midranksSorted(const double* v, int n, double* out) {
  int i = 0;
  while (i < n) {
    int j = i;
    while (j + 1 < n && v[j + 1] == v[i]) ++j;
    const double mr = 0.5 * (i + j) + 1.0;      // mean of (i+1) .. (j+1)
    for (int k = i; k <= j; ++k) out[k] = mr;
    i = j + 1;
  }
}

// Group means and the two Brunner-Munzel variance components. `rx`/`ry` are the
// pooled midranks of each group, ascending; `wx`/`wy` are scratch buffers.
static inline void bmParts(const double* rx, int n1, const double* ry, int n2,
                           double* wx, double* wy,
                           double& m1, double& m2, double& v1, double& v2) {
  midranksSorted(rx, n1, wx);
  midranksSorted(ry, n2, wy);

  m1 = 0.0; m2 = 0.0;
  for (int i = 0; i < n1; ++i) m1 += rx[i];
  for (int j = 0; j < n2; ++j) m2 += ry[j];
  m1 /= n1;
  m2 /= n2;

  const double c1 = m1 - 0.5 * (n1 + 1.0);
  const double c2 = m2 - 0.5 * (n2 + 1.0);

  double d;
  v1 = 0.0; v2 = 0.0;
  for (int i = 0; i < n1; ++i) { d = rx[i] - wx[i] - c1; v1 += d * d; }
  for (int j = 0; j < n2; ++j) { d = ry[j] - wy[j] - c2; v2 += d * d; }
  v1 /= (n1 - 1.0);
  v2 /= (n2 - 1.0);
}

// Studentized statistic only (permutation inner loop).
static inline double bmStatistic(const double* rx, int n1,
                                 const double* ry, int n2,
                                 double* wx, double* wy) {
  double m1, m2, v1, v2;
  bmParts(rx, n1, ry, n2, wx, wy, m1, m2, v1, v2);
  const double den = n1 * v1 + n2 * v2;
  if (!(den > 0.0)) {
    // Degenerate split (e.g. complete separation): +/-Inf is the correct
    // ordering for the permutation count, and never an arbitrary tie.
    if (m2 > m1) return R_PosInf;
    if (m2 < m1) return R_NegInf;
    return 0.0;                                  // every observation tied
  }
  return (double) n1 * n2 * (m2 - m1) / (n1 + n2) / std::sqrt(den);
}

// Split the sorted pooled ranks according to `in`, preserving ascending order.
static inline void gather(const double* rs, const char* in, int N,
                          double* rx, double* ry) {
  int a = 0, b = 0;
  for (int p = 0; p < N; ++p) {
    if (in[p]) rx[a++] = rs[p]; else ry[b++] = rs[p];
  }
}

// ------------------------------------------------------------------ core ---

//' @noRd
// [[Rcpp::export(.bmCore)]]
NumericVector bmCore(NumericVector rx, NumericVector ry) {
  const int n1 = rx.size(), n2 = ry.size();
  std::vector<double> sx(rx.begin(), rx.end()), sy(ry.begin(), ry.end());
  std::sort(sx.begin(), sx.end());
  std::sort(sy.begin(), sy.end());
  std::vector<double> wx(n1), wy(n2);

  double m1, m2, v1, v2;
  bmParts(sx.data(), n1, sy.data(), n2, wx.data(), wy.data(), m1, m2, v1, v2);

  const double phat = (m2 - 0.5 * (n2 + 1.0)) / n1;
  const double den  = n1 * v1 + n2 * v2;
  const double se   = std::sqrt(v1 / (n1 * (double) n2 * n2) +
                                v2 / (n2 * (double) n1 * n1));

  double stat, df;
  if (!(den > 0.0)) {
    stat = (m2 > m1) ? R_PosInf : ((m2 < m1) ? R_NegInf : 0.0);
    df   = NA_REAL;
  } else {
    stat = (double) n1 * n2 * (m2 - m1) / (n1 + n2) / std::sqrt(den);
    df   = den * den / ((n1 * v1) * (n1 * v1) / (n1 - 1.0) +
                        (n2 * v2) * (n2 * v2) / (n2 - 1.0));
  }
  return NumericVector::create(
    _["phat"] = phat, _["statistic"] = stat, _["df"] = df,
    _["se"] = se, _["v1"] = v1, _["v2"] = v2);
}

// ----------------------------------------------------------- permutation ---

//' @noRd
// [[Rcpp::export(.bmPermExact)]]
NumericVector bmPermExact(NumericVector rSorted, int n1, double tObs) {
  const int N = rSorted.size(), n2 = N - n1;
  const double* rs = rSorted.begin();
  // tObs is infinite under complete separation; a relative tolerance would
  // then be Inf and tObs - tol would be NaN, silencing every comparison.
  const double tol = std::isfinite(tObs) ? 1e-9 * (1.0 + std::fabs(tObs)) : 0.0;

  std::vector<int>    idx(n1);
  std::vector<double> rx(n1), ry(n2), wx(n1), wy(n2);
  std::vector<char>   in(N, 0);
  for (int i = 0; i < n1; ++i) { idx[i] = i; in[i] = 1; }

  double cLess = 0.0, cGreater = 0.0, cAbs = 0.0, total = 0.0;
  R_xlen_t sinceCheck = 0;

  while (true) {
    gather(rs, in.data(), N, rx.data(), ry.data());
    const double t = bmStatistic(rx.data(), n1, ry.data(), n2,
                                 wx.data(), wy.data());
    if (t <=  tObs + tol)                            cLess    += 1.0;
    if (t >=  tObs - tol)                            cGreater += 1.0;
    if (std::fabs(t) >= std::fabs(tObs) - tol)       cAbs     += 1.0;
    total += 1.0;

    if (++sinceCheck >= 1048576) { Rcpp::checkUserInterrupt(); sinceCheck = 0; }

    // next combination in lexicographic order
    int i = n1 - 1;
    while (i >= 0 && idx[i] == i + N - n1) --i;
    if (i < 0) break;
    std::fill(in.begin(), in.end(), 0);
    ++idx[i];
    for (int j = i + 1; j < n1; ++j) idx[j] = idx[j - 1] + 1;
    for (int j = 0; j < n1; ++j) in[idx[j]] = 1;
  }
  return NumericVector::create(_["less"] = cLess, _["greater"] = cGreater,
                               _["two.sided"] = cAbs, _["n"] = total);
}

//' @noRd
// [[Rcpp::export(.bmPermMC)]]
NumericVector bmPermMC(NumericVector rSorted, int n1, double tObs, int nPerm) {
  const int N = rSorted.size(), n2 = N - n1;
  const double* rs = rSorted.begin();
  // tObs is infinite under complete separation; a relative tolerance would
  // then be Inf and tObs - tol would be NaN, silencing every comparison.
  const double tol = std::isfinite(tObs) ? 1e-9 * (1.0 + std::fabs(tObs)) : 0.0;

  std::vector<int>    pool(N);
  std::vector<double> rx(n1), ry(n2), wx(n1), wy(n2);
  std::vector<char>   in(N);
  for (int i = 0; i < N; ++i) pool[i] = i;

  double cLess = 0.0, cGreater = 0.0, cAbs = 0.0;

  for (int b = 0; b < nPerm; ++b) {
    // partial Fisher-Yates: draw n1 distinct positions, O(n1)
    for (int i = 0; i < n1; ++i) {
      const int j = i + (int) (unif_rand() * (N - i));
      std::swap(pool[i], pool[j < N ? j : N - 1]);
    }
    std::fill(in.begin(), in.end(), 0);
    for (int i = 0; i < n1; ++i) in[pool[i]] = 1;

    gather(rs, in.data(), N, rx.data(), ry.data());
    const double t = bmStatistic(rx.data(), n1, ry.data(), n2,
                                 wx.data(), wy.data());
    if (t <=  tObs + tol)                        cLess    += 1.0;
    if (t >=  tObs - tol)                        cGreater += 1.0;
    if (std::fabs(t) >= std::fabs(tObs) - tol)   cAbs     += 1.0;

    if ((b & 0xFFFF) == 0xFFFF) Rcpp::checkUserInterrupt();
  }
  return NumericVector::create(_["less"] = cLess, _["greater"] = cGreater,
                               _["two.sided"] = cAbs, _["n"] = (double) nPerm);
}
