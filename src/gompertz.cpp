#include <algorithm>
#include <cmath>

#define STRICT_R_HEADERS
#include <Rcpp.h>


// --- distribution helpers ---
inline double below_distribution(bool lower_tail, bool give_log) {
  if (lower_tail) {
    return give_log ? R_NegInf : 0;
  } else {
    return give_log ? 0 : R_NegInf;
  }
}

// NA and NaN are distinct on the way out, as they are in the base R
// distribution functions: NA propagates NA, NaN propagates NaN
inline bool missing(const double a, const double b, const double c) {
  return ISNAN(a) || ISNAN(b) || ISNAN(c);
}

inline double missing_value(const double a, const double b, const double c) {
  return (R_IsNA(a) || R_IsNA(b) || R_IsNA(c)) ? NA_REAL : R_NaN;
}


// --- gompertz implementation ---
namespace {

namespace gompertz {

// below this the shape is treated as zero and the exponential limit is used;
// .checkGompertz() and qgompertz() on the R side use the same value
constexpr double shape_tol = 1e-12;

inline double exprel(const double x) {
  if (x != 0.0) {
    return expm1(x) / x;
  } else {
    return 1.0;
  }
}

// out-of-range parameters give NaN, as in the base R distribution functions;
// the accompanying warning is issued once per call, not once per element
inline bool bad(const double shape, const double rate) {
  return rate <= 0;
}

inline double safe_coeff(const double q,
                         const double shape,
                         const double rate) {
  if (!std::isinf(q)) {
    const double scale_q = shape * q;
    return -rate * q * exprel(scale_q);
  }
  // q is +Inf here, negative q having been dealt with by the caller:
  // exp(shape * q) diverges for a positive shape and vanishes for a
  // negative one, where a mass of exp(rate/shape) escapes to infinity
  return (shape > 0) ? R_NegInf : rate / shape;
}

// ---------------- density ----------------
class density {
public:
  typedef double result_type;

  inline double operator()(const double x,
                         const double shape,
                         const double rate) const {

    if (missing(x, shape, rate)) {
      return missing_value(x, shape, rate);
    }

    if (bad(shape, rate)) {
      return R_NaN;
    }

    if (x < 0 || std::isinf(x)) {
      return R_NegInf;
    }

    const double scale_x = shape * x;
    const double shift   = x * exprel(scale_x);

    return std::log(rate) + scale_x - rate * shift;
  }
};

// ---------------- cdf ----------------
class cdf {
public:
  typedef double result_type;

  cdf(bool lower_tail_, bool give_log_) :
    lower_tail(lower_tail_),
    give_log(give_log_) {}

  inline double operator()(const double q,
                         const double shape,
                         const double rate) const {

    if (missing(q, shape, rate)) {
      return missing_value(q, shape, rate);
    }

    if (bad(shape, rate)) {
      return R_NaN;
    }

    if (q < 0) {
      return below_distribution(lower_tail, give_log);
    }

    // numerically stable check instead of shape != 0
    if (std::abs(shape) > shape_tol) {

      const double coeff = safe_coeff(q, shape, rate);

      if ((!give_log) && (lower_tail)) {
        return -expm1(coeff);
      }

      if ((!give_log) && (!lower_tail)) {
        return std::exp(coeff);
      }

      if (give_log && lower_tail) {
        return log1p(-std::exp(coeff));
      }

      return coeff;

    } else {
      return R::pexp(q * rate, 1.0, lower_tail, give_log);
    }
  }

private:
  bool lower_tail;
  bool give_log;
};

} // namespace gompertz
} // unnamed namespace


// ---------------- exported functions ----------------

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
dgompertz_cpp(const Rcpp::NumericVector& x,
              const Rcpp::NumericVector& shape,
              const Rcpp::NumericVector& rate,
              const bool log) {

  if (x.size() == 0) return x;

  const R_xlen_t size = std::max({x.size(), shape.size(), rate.size()});

  Rcpp::NumericVector out = Rcpp::mapply(
    Rcpp::rep_len(x, size),
    Rcpp::rep_len(shape, size),
    Rcpp::rep_len(rate, size),
    gompertz::density()
  );

  if (!log) out = Rcpp::exp(out);
  return out;
}


// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
pgompertz_cpp(const Rcpp::NumericVector& q,
               const Rcpp::NumericVector& shape,
               const Rcpp::NumericVector& rate,
               const bool lower_tail,
               const bool give_log) {

  if (q.size() == 0) return q;

  const R_xlen_t size = std::max({q.size(), shape.size(), rate.size()});

  return Rcpp::mapply(
    Rcpp::rep_len(q, size),
    Rcpp::rep_len(shape, size),
    Rcpp::rep_len(rate, size),
    gompertz::cdf(lower_tail, give_log)
  );
}


// the single definition of parameter validity, shared by d, p, q and r
// through .checkGompertz() on the R side
//
// [[Rcpp::export(name="checkGompertz_cpp", rng=false)]]
Rcpp::LogicalVector
check_gompertz_cpp(const Rcpp::NumericVector& shape,
               const Rcpp::NumericVector& rate) {

  if (shape.size() == 0 && rate.size() == 0) {
    return Rcpp::LogicalVector(0);
  }

  const R_xlen_t size = std::max(shape.size(), rate.size());

  return !Rcpp::mapply(
      Rcpp::rep_len(shape, size),
      Rcpp::rep_len(rate, size),
      gompertz::bad
  );
}
