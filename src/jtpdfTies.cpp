#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <functional>
#include <cmath>
using namespace Rcpp;


// Exact null distribution of the Jonckheere-Terpstra statistic with ties.
//
// With ties JT depends on the data only through the k x m table of counts of
// group by distinct value, and under the null hypothesis that table follows
// the multiple hypergeometric distribution with the group sizes and the tie
// counts as its margins. The distribution is therefore built by walking
// through the groups and splitting off one row of the table at a time, the
// state being the vector of values not yet assigned. The contribution of a
// group is its Mann-Whitney statistic against the whole remaining pool,
// which is what the sum over the pairs k < l amounts to, so that the last
// group adds nothing.
//
// Everything is counted twice, so that the half weights of the ties stay
// integral: the returned vector holds P(2 * JT = i), i = 0, ..., 2M, with M
// the largest attainable value of JT.
//
// [[Rcpp::export]]
NumericVector jtpdfTies_cpp(IntegerVector gsize, IntegerVector cnt) {

  const int ng = gsize.size();
  const int m  = cnt.size();

  if (ng < 2) stop("at least two groups are needed");
  if (m < 1)  stop("'cnt' must not be empty");

  for (int i = 0; i < ng; i++)
    if (gsize[i] < 1) stop("group sizes must be positive");

  for (int a = 0; a < m; a++)
    if (cnt[a] < 1) stop("tie counts must be positive");

  const int N = sum(gsize);

  if (sum(cnt) != N)
    stop("the tie counts must sum to the number of observations");

  // largest attainable value of 2 * JT
  int maxJ2 = 0;

  for (int i = 0, ahead = N; i < ng; i++) {
    ahead -= gsize[i];
    maxJ2 += 2 * gsize[i] * ahead;
  }

  // mixed radix encoding of the counts not yet assigned. The states times
  // the support are what the recursion has to hold; the R wrapper keeps the
  // regular calls well below this limit, which is here for direct calls
  const long long maxCells = 1000000000LL;
  const long long support = (long long)maxJ2 + 1;

  std::vector<long long> radix(m);
  long long code0 = 0;
  long long f = 1;

  for (int a = 0; a < m; a++) {

    radix[a] = f;
    code0 += f * (long long)cnt[a];

    if (f > maxCells / support / (long long)(cnt[a] + 1))
      stop("the table of group by distinct value is too large to enumerate");

    f *= (long long)(cnt[a] + 1);
  }

  // log binomial coefficients for the hypergeometric weights
  std::vector< std::vector<double> > lch(N + 1);

  for (int i = 0; i <= N; i++) {
    lch[i].resize(i + 1);
    for (int j = 0; j <= i; j++) lch[i][j] = R::lchoose((double)i, (double)j);
  }

  typedef std::unordered_map<long long, std::vector<double> > Layer;

  Layer layer;
  layer[code0] = std::vector<double>(maxJ2 + 1, 0.0);
  layer[code0][0] = 1.0;

  std::vector<int> rem(m), prefix(m + 1);

  for (int grp = 0; grp < ng - 1; grp++) {

    Layer next;
    const int nGrp = gsize[grp];

    for (Layer::const_iterator it = layer.begin(); it != layer.end(); ++it) {

      const std::vector<double>& from = it->second;
      long long code = it->first;

      int nRem = 0;

      for (int i = 0; i < m; i++) {
        rem[i] = (int)((code / radix[i]) % (long long)(cnt[i] + 1));
        nRem += rem[i];
      }

      // prefix[i] = number of values still available at positions 0, ..., i-1
      prefix[0] = 0;
      for (int i = 0; i < m; i++) prefix[i + 1] = prefix[i] + rem[i];

      const double lnorm = lch[nRem][nGrp];

      // walk the compositions of nGrp bounded by rem, from the largest value
      // downwards, so that the pool lying above the current value is known:
      // above = number of values already left in the pool at higher positions
      std::function<void(int, int, double, int, int, long long)> walk =
        [&](int i, int todo, double lw, int gain, int above, long long code2) {

          if (i < 0) {
            if (todo != 0) return;

            const double w = std::exp(lw - lnorm);

            Layer::iterator dst = next.find(code2);

            if (dst == next.end())
              dst = next.insert(std::make_pair(
                code2, std::vector<double>(maxJ2 + 1, 0.0))).first;

            std::vector<double>& to = dst->second;

            for (int j = 0; j + gain <= maxJ2; j++)
              if (from[j] != 0.0) to[j + gain] += from[j] * w;

            return;
          }

          // nothing left below to take the rest of the group from
          if (todo > prefix[i + 1]) return;

          const int hi = std::min(rem[i], todo);

          for (int a = 0; a <= hi; a++) {

            const int left = rem[i] - a;

            walk(i - 1,
                 todo - a,
                 lw + lch[rem[i]][a],
                 gain + a * (2 * above + left),
                 above + left,
                 code2 + radix[i] * (long long)left);
          }
        };

      walk(m - 1, nGrp, 0.0, 0, 0, 0);
    }

    layer.swap(next);

    Rcpp::checkUserInterrupt();
  }

  NumericVector prob(maxJ2 + 1);

  for (Layer::const_iterator it = layer.begin(); it != layer.end(); ++it) {
    const std::vector<double>& d = it->second;
    for (int j = 0; j <= maxJ2; j++) prob[j] += d[j];
  }

  return prob;
}
