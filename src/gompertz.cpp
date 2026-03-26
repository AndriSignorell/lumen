
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

template <typename T1>
inline Rcpp::NumericVector perhaps_exp(const T1& y, bool log) {
  return log ? y : Rcpp::exp(y);
}

// --- rep_len helper (like flexsurv) ---
namespace flexsurv {

template <int RTYPE, bool NA, typename T>
inline Rcpp::sugar::Rep_len<RTYPE, NA, T>
rep_len(const Rcpp::VectorBase<RTYPE, NA, T>& t, R_xlen_t len) {
  if (t.size() == 0) {
    Rcpp::stop("zero length vector provided");
  } else {
    return Rcpp::rep_len(t, len);
  }
}

}

// --- gompertz implementation ---
namespace {

namespace gompertz {

inline double exprel(const double x) {
  if (x != 0.0) {
    return expm1(x) / x;
  } else {
    return 1.0;
  }
}

inline double safe_coeff(const double q,
                         const double shape,
                         const double rate) {
  if (!std::isinf(q)) {
    const double scale_q = shape * q;
    return -rate * q * exprel(scale_q);
  } else {
    return R_NegInf;
  }
}

inline bool bad(const double shape, const double rate) {
  if (rate <= 0) {
    Rcpp::warning("Non-positive rate parameter");
    return true;
  }
  return false;
}

// ---------------- density ----------------
class density {
public:
  typedef double result_type;
  
  inline double operator()(const double x,
                         const double shape,
                         const double rate) const {
    
    if (Rcpp::NumericVector::is_na(x) ||
        Rcpp::NumericVector::is_na(shape) ||
        Rcpp::NumericVector::is_na(rate)) {
      return NA_REAL;
    }
    
    if (bad(shape, rate)) {
      return NA_REAL;
    }
    
    if (x < 0) {
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
    
    if (Rcpp::NumericVector::is_na(q) ||
        Rcpp::NumericVector::is_na(shape) ||
        Rcpp::NumericVector::is_na(rate)) {
      return NA_REAL;
    }
    
    if (bad(shape, rate)) {
      return NA_REAL;
    }
    
    if (q < 0) {
      return below_distribution(lower_tail, give_log);
    }
    
    // numerically stable check instead of shape != 0
    if (std::abs(shape) > 1e-12) {
      
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
  
  return perhaps_exp(
    Rcpp::mapply(
      flexsurv::rep_len(x, size),
      flexsurv::rep_len(shape, size),
      flexsurv::rep_len(rate, size),
      gompertz::density()
    ),
    log
  );
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
    flexsurv::rep_len(q, size),
    flexsurv::rep_len(shape, size),
    flexsurv::rep_len(rate, size),
    gompertz::cdf(lower_tail, give_log)
  );
}

// [[Rcpp::export(name="check.gompertz", rng=false)]]
Rcpp::LogicalVector
check_gompertz(const Rcpp::NumericVector& shape,
               const Rcpp::NumericVector& rate) {
  
  if (shape.size() == 0 && rate.size() == 0) {
    return Rcpp::LogicalVector(0);
  }
  
  const R_xlen_t size = std::max(shape.size(), rate.size());
  
  return !Rcpp::mapply(
      flexsurv::rep_len(shape, size),
      flexsurv::rep_len(rate, size),
      gompertz::bad
  );
}

