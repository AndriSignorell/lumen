// ============================================================
// mad_boot.cpp
//
// Bootstrap-Intervall fuer die Median Absolute Deviation (MAD),
// definiert wie base-R mad(x, constant).
//
// Gegenueber der bisherigen Fassung nur zwei Aenderungen:
//
//   - die Statistik kommt aus mad_impl.h statt als lokale Kopie
//     (sie stand zweimal in derselben Datei, als mad_cpp und mad_cpp2)
//   - Argumentpruefung, die bisher ganz fehlte
//
// Die Zweistichprobenfassungen sind in mad_two_sample_boot.cpp und
// benutzen NICHT run_boot - das Ziehschema laesst sich mit einem
// Zeilen-Resampler nicht ausdruecken. Siehe die Begruendung dort.
// ============================================================

#include "boot_framework.h"
#include "mad_impl.h"

using namespace Rcpp;


namespace {

struct MadFn {

  double constant;

  explicit MadFn(double c = 1.4826) : constant(c) {}

  double compute(const arma::mat& X, const arma::vec& /* y */) const {
    const arma::vec col = X.col(0);
    std::vector<double> work(col.begin(), col.end());
    std::vector<double> dev(work.size());
    return desctoolsx::madImpl(work, dev, constant);
  }
};

}  // namespace


// [[Rcpp::export]]
NumericVector mad_boot_cpp(NumericVector x,
                           int    R        = 1000,
                           double alpha    = 0.05,
                           double constant = 1.4826,
                           int    seed     = -1,
                           String method   = "perc") {

  if (R < 2)
    stop("'R' must be at least 2.");
  if (alpha <= 0.0 || alpha >= 1.0)
    stop("'alpha' must lie in (0, 1).");
  if (!(constant > 0.0))
    stop("'constant' must be a positive number.");
  if (x.size() < 2)
    stop("'x' must contain at least 2 observations.");

  return run_boot(
    vec_to_matrix(x),
    dummy_vec(x.size()),
    R,
    alpha,
    seed,
    MadFn(constant),
    std::string(method.get_cstring())
  );
}
