#ifndef DESCTOOLSX_MAD_IMPL_H
#define DESCTOOLSX_MAD_IMPL_H

// Median absolute deviation, matching base-R mad(x, constant).
//
// Eine Implementierung statt bisher zwei (mad_cpp und mad_cpp2, Kopien
// voneinander in derselben Datei).
//
// madImpl() ARBEITET AUF DEM UEBERGEBENEN PUFFER und veraendert dessen
// Reihenfolge - nth_element permutiert. Das ist Absicht: der Bootstrap
// legt den Resample-Puffer ohnehin einmal pro Chunk an und fuellt ihn
// pro Replikation neu, eine zusaetzliche Kopie waere reine Verschwendung.

#include <R_ext/Arith.h>   // NA_REAL, ISNAN
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace desctoolsx {

// median of v, permutes v
inline double medianInplaceMad(std::vector<double>& v) {

  const std::size_t n    = v.size();
  const std::size_t mid  = n / 2;

  std::nth_element(v.begin(), v.begin() + mid, v.end());
  const double upper = v[mid];

  if (n % 2 == 1)
    return upper;

  const double lower = *std::max_element(v.begin(), v.begin() + mid);
  return 0.5 * (upper + lower);
}


// work is consumed; dev must have the same length as work
inline double madImpl(std::vector<double>& work,
                      std::vector<double>& dev,
                      double constant) {

  const std::size_t n = work.size();

  if (n == 0)
    return NA_REAL;

  // nth_element() braucht eine strikte schwache Ordnung; ein NaN macht
  // die Vergleiche inkonsistent, was undefiniertes Verhalten ist - nicht
  // bloss ein NaN im Ergebnis.
  for (std::size_t i = 0; i < n; i++)
    if (ISNAN(work[i]))
      return NA_REAL;

  // dev VOR der Permutation fuellen waere falsch; medianInplaceMad
  // permutiert work, aber die Abweichungen haengen nur am Median, nicht
  // an der Reihenfolge - deshalb erst Median, dann dev aus dem
  // (permutierten) work.
  const double med = medianInplaceMad(work);

  for (std::size_t i = 0; i < n; i++)
    dev[i] = std::fabs(work[i] - med);

  return constant * medianInplaceMad(dev);
}


// bequeme Fassung fuer Aufrufer ausserhalb einer Schleife
inline double madImpl(const std::vector<double>& x, double constant) {
  std::vector<double> work = x;
  std::vector<double> dev(x.size());
  return madImpl(work, dev, constant);
}

}  // namespace desctoolsx

#endif  // DESCTOOLSX_MAD_IMPL_H
