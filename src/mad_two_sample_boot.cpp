// ============================================================
// mad_two_sample_boot.cpp
//
// Bootstrap-Intervalle fuer Zweistichproben-MAD-Statistiken:
//
//   mad_diff_boot_cpp    MAD(x) - MAD(y)
//   mad_ratio_boot_cpp   (MAD(x) / MAD(y))^2
//
// ------------------------------------------------------------
// WAS HIER FALSCH WAR
// ------------------------------------------------------------
// Die alte Fassung stapelte x und y in eine n_max x 2 Matrix, mit NA
// aufgefuellt, und liess run_boot ZEILEN ziehen. Der Kommentar dazu -
// "resampling independently within each column (same bootstrap index is
// applied to both)" - widerspricht sich selbst: derselbe Zeilenindex auf
// beide Spalten IST der gepaarte Bootstrap. Drei Folgen:
//
//  1. Ueber den gemeinsamen Multiplizitaetsvektor korrelieren MAD(x*)
//     und MAD(y*). Die Varianz der Differenz stimmt nicht, und zwar
//     auch bei gleichen Laengen.
//
//  2. Bei UNGLEICHEN Laengen wird auf n_max aufgefuellt und ueber n_max
//     Zeilen gezogen. Die Zahl der nicht-NA-Werte pro Replikation ist
//     damit binomialverteilt statt fest: bei nx=50, ny=200 schwankt die
//     effektive x-Stichprobe um 50 herum. Die Bootstrapverteilung mischt
//     Stichprobenvariabilitaet mit zufaelliger Stichprobengroesse.
//
//  3. extract_col() verwarf NA und unterschied dabei nicht zwischen
//     Auffuell-NA und echtem fehlendem Wert in den Daten.
//
// Ein eigener Worker mit ZWEI unabhaengigen Indexziehungen loest alle
// drei auf einmal: keine Matrix, keine Auffuellung, keine NA-Frage.
//
// ------------------------------------------------------------
// AENDERUNG AM VERHALTEN, die du kennen solltest
// ------------------------------------------------------------
// NA in x oder y ist jetzt ein FEHLER statt stiller Ausschluss. Wer
// na.rm will, entscheidet das auf der R-Seite - hier war es bisher eine
// Nebenwirkung der Auffuellmechanik und nicht als Zusage gemeint.
//
// Ebenfalls geaendert: seed < 0 zieht den Basis-Seed jetzt aus R's RNG
// statt aus std::random_device. Damit steht auch dieser Pfad unter
// set.seed(); vorher war er die einzige Stelle, an der set.seed()
// nachweislich wirkungslos blieb.
// ============================================================

#include "boot_framework.h"   // quantile_type7, boot_valid_sorted, BCa
#include "mad_impl.h"

#include <random>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;


namespace {

enum MadTwoStat { MAD_DIFF = 0, MAD_RATIO = 1 };


// Puffer werden vom Aufrufer gestellt und hier verbraucht.
inline double madTwoStat(std::vector<double>& xw, std::vector<double>& xd,
                         std::vector<double>& yw, std::vector<double>& yd,
                         double constant, int stat) {

  const double mx = desctoolsx::madImpl(xw, xd, constant);
  const double my = desctoolsx::madImpl(yw, yd, constant);

  if (ISNAN(mx) || ISNAN(my))
    return NA_REAL;

  if (stat == MAD_DIFF)
    return mx - my;

  if (my == 0.0)
    return NA_REAL;

  const double r = mx / my;
  return r * r;
}


// ------------------------------------------------------------
// Worker: zwei unabhaengige Ziehungen, jede in ihrer eigenen Laenge
// ------------------------------------------------------------
struct MadTwoWorker : public Worker {

  const RVector<double> x;
  const RVector<double> y;
  const std::size_t     nx;
  const std::size_t     ny;
  const double          constant;
  const unsigned int    seed;
  const int             stat;
  RVector<double>       out;

  MadTwoWorker(NumericVector x_, NumericVector y_, double constant_,
               unsigned int seed_, int stat_, NumericVector out_)
    : x(x_), y(y_), nx(x_.size()), ny(y_.size()),
      constant(constant_), seed(seed_), stat(stat_), out(out_) {}

  void operator()(std::size_t begin, std::size_t end) {

    // je Chunk einmal
    std::vector<double> xw(nx), xd(nx), yw(ny), yd(ny);

    for (std::size_t b = begin; b < end; b++) {

      // seed + b wie im uebrigen Framework: das Ergebnis haengt am
      // Replikationsindex, nicht an der Threadzahl (F14)
      std::mt19937 rng(seed + static_cast<unsigned int>(b));

      std::uniform_int_distribution<std::size_t> dx(0, nx - 1);
      std::uniform_int_distribution<std::size_t> dy(0, ny - 1);

      // ZWEI Ziehungen. Das ist der ganze Unterschied zur alten
      // Fassung - und der Grund, warum die Laengen nicht mehr
      // uebereinstimmen muessen.
      for (std::size_t i = 0; i < nx; i++) xw[i] = x[dx(rng)];
      for (std::size_t i = 0; i < ny; i++) yw[i] = y[dy(rng)];

      out[b] = madTwoStat(xw, xd, yw, yd, constant, stat);
    }
  }
};


// ------------------------------------------------------------
// Jackknife ueber BEIDE Stichproben, fuer die BCa-Beschleunigung
// ------------------------------------------------------------
// nx + ny Werte: je einmal eine Beobachtung aus x weglassen, je einmal
// eine aus y. Das ist die uebliche Konstruktion fuer eine Statistik aus
// zwei unabhaengigen Stichproben.
NumericVector madTwoJackknife(const NumericVector& x,
                              const NumericVector& y,
                              double constant, int stat) {

  const std::size_t nx = x.size();
  const std::size_t ny = y.size();

  NumericVector jack(nx + ny);

  std::vector<double> xs(nx), ys(ny);
  for (std::size_t i = 0; i < nx; i++) xs[i] = x[i];
  for (std::size_t j = 0; j < ny; j++) ys[j] = y[j];

  {
    std::vector<double> xw(nx - 1), xd(nx - 1), yw(ny), yd(ny);

    for (std::size_t i = 0; i < nx; i++) {
      std::size_t k = 0;
      for (std::size_t t = 0; t < nx; t++)
        if (t != i) xw[k++] = xs[t];
      yw = ys;
      jack[i] = madTwoStat(xw, xd, yw, yd, constant, stat);
    }
  }

  {
    std::vector<double> xw(nx), xd(nx), yw(ny - 1), yd(ny - 1);

    for (std::size_t j = 0; j < ny; j++) {
      std::size_t k = 0;
      for (std::size_t t = 0; t < ny; t++)
        if (t != j) yw[k++] = ys[t];
      xw = xs;
      jack[nx + j] = madTwoStat(xw, xd, yw, yd, constant, stat);
    }
  }

  return jack;
}


void checkSample(const NumericVector& v, const char* name) {

  if (v.size() < 2)
    stop("'%s' must contain at least 2 observations.", name);

  for (R_xlen_t i = 0; i < v.size(); i++)
    if (ISNAN(v[i]))
      stop("'%s' contains missing values; remove them before calling.",
           name);
}


NumericVector madTwoBoot(NumericVector x, NumericVector y,
                         int R, double alpha, double constant,
                         int seed, std::string method, int stat) {

  if (R < 2)
    stop("'R' must be at least 2.");
  if (!(alpha > 0.0) || !(alpha < 1.0))
    stop("'alpha' must lie in (0, 1).");
  if (!(constant > 0.0))
    stop("'constant' must be a positive number.");
  if (method != "perc" && method != "bca")
    stop("Unknown bootstrap method '%s'. Use 'perc' or 'bca'.",
         method.c_str());

  checkSample(x, "x");
  checkSample(y, "y");

  // Basis-Seed aus R's RNG statt aus std::random_device: so bleibt
  // set.seed() die einzige Steuerung. Das hier laeuft auf dem
  // Hauptthread, Rcpp haelt den RNG-Zustand.
  const unsigned int base_seed =
    (seed < 0) ? static_cast<unsigned int>(R::unif_rand() * 4294967295.0)
               : static_cast<unsigned int>(seed);

  NumericVector stats(R);

  MadTwoWorker worker(x, y, constant, base_seed, stat, stats);
  parallelFor(0, R, worker);

  // Punktschaetzer aus den ORIGINALstichproben
  std::vector<double> xw(x.begin(), x.end()), xd(x.size());
  std::vector<double> yw(y.begin(), y.end()), yd(y.size());
  const double est = madTwoStat(xw, xd, yw, yd, constant, stat);

  NumericVector s = boot_valid_sorted(stats);

  if (s.size() == 0)
    stop("All bootstrap samples failed.");

  // Wie viele Replikationen ausgefallen sind, gehoert gesagt - stumm
  // aus 99 statt 999 Werten ein Intervall zu bilden ist keine Auskunft.
  if (s.size() < static_cast<R_xlen_t>(R))
    warning("%d of %d bootstrap replicates were not usable.",
            (int) (R - s.size()), R);

  double pl = alpha / 2.0;
  double pu = 1.0 - alpha / 2.0;

  if (method == "bca") {
    const double z0  = bca_z0(s, est);
    const double acc = bca_acceleration(
                         madTwoJackknife(x, y, constant, stat));
    std::pair<double, double> lv = bca_quantile_levels(z0, acc, alpha);
    pl = lv.first;
    pu = lv.second;
  }

  NumericVector out(3);
  out[0] = est;
  out[1] = quantile_type7(s, pl);
  out[2] = quantile_type7(s, pu);
  out.attr("names") = CharacterVector::create("est", "lci", "uci");

  return out;
}

}  // namespace


// ============================================================
// Exportierte Funktionen - Signaturen unveraendert
// ============================================================

// [[Rcpp::export]]
NumericVector mad_diff_boot_cpp(NumericVector x,
                                NumericVector y,
                                int    R        = 1000,
                                double alpha    = 0.05,
                                double constant = 1.4826,
                                int    seed     = -1,
                                String method   = "perc") {

  return madTwoBoot(x, y, R, alpha, constant, seed,
                    std::string(method.get_cstring()), MAD_DIFF);
}


// [[Rcpp::export]]
NumericVector mad_ratio_boot_cpp(NumericVector x,
                                 NumericVector y,
                                 int    R        = 1000,
                                 double alpha    = 0.05,
                                 double constant = 1.4826,
                                 int    seed     = -1,
                                 String method   = "perc") {

  return madTwoBoot(x, y, R, alpha, constant, seed,
                    std::string(method.get_cstring()), MAD_RATIO);
}
