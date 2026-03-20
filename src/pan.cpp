

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double pan(NumericVector A, int M, double C, int N) {
  
  const double ZERO = 0.0;
  const double ONE  = 1.0;
  const double HALF = 0.5;
  const double TWO  = 2.0;
  const double EPS  = 1e-14;
  
  int D, H, I, J1, J2, J3, J4, K, L1, L2, NU, N2;
  double NUM, PIN, SGN, SUM, SUM1, U, V, Y, X;
  
  // ordering
  if (A[1] > A[M]) {
    H = M; K = -1; I = 1;
  } else {
    H = 1; K = 1; I = M;
  }
  
  X = A[0];
  
  bool found = false;
  for (NU = H; (K == 1 ? NU <= I : NU >= I); NU += K) {
    if (A[NU] >= X) {
      found = true;
      break;
    }
  }
  
  if (!found) {
    if (C >= ZERO) return ONE;
  }
  
  if (NU == H && C <= ZERO) return ZERO;
  
  if (K == 1) NU = NU - 1;
  
  H = M - NU;
  
  if (C == ZERO) {
    Y = H - NU;
  } else {
    Y = C * (A[1] - A[M]);
  }
  
  if (Y >= ZERO) {
    D = 2;
    H = NU;
    K = -K;
    J1 = 0;
    J2 = 2;
    J3 = 3;
    J4 = 1;
  } else {
    D = -2;
    NU = NU + 1;
    J1 = M - 2;
    J2 = M - 1;
    J3 = M + 1;
    J4 = M;
  }
  
  PIN = TWO * std::atan(ONE) / N;
  SUM = HALF * (K + 1);
  SGN = K / (double) N;
  N2 = 2 * N - 1;
  
  for (L1 = H - 2 * (H / 2); L1 >= 0; L1--) {
    
    for (L2 = J2; (D == 2 ? L2 <= NU : L2 >= NU); L2 += D) {
      
      double A_j4 = A[J4];
      double A_l2 = A[L2];
      
      U = HALF * (A_j4 + A_l2);
      V = HALF * (A_j4 - A_l2);
      
      SUM1 = ZERO;
      
      for (I = 1; I <= N2; I += 2) {
        
        Y = U - V * std::cos((double) I * PIN);
        NUM = Y - X;
        
        if (std::abs(NUM) < EPS) continue;
        
        // -------- LOG DOMAIN --------
          double logProd = -C / NUM;
          
          // first product
          for (int k = 1; k <= J1; k++) {
            double denom = Y - A[k];
            if (std::abs(denom) < EPS) {
              logProd = -INFINITY;
              break;
            }
            logProd += std::log(std::abs(NUM)) - std::log(std::abs(denom));
          }
          
          // second product
          for (int k = J3; k <= M; k++) {
            double denom = Y - A[k];
            if (std::abs(denom) < EPS) {
              logProd = -INFINITY;
              break;
            }
            logProd += std::log(std::abs(NUM)) - std::log(std::abs(denom));
          }
          
          if (!std::isfinite(logProd)) continue;
          
          double term = std::exp(0.5 * logProd);
          SUM1 += term;
      }
      
      SGN = -SGN;
      SUM += SGN * SUM1;
      
      J1 += D;
      J3 += D;
      J4 += D;
    }
    
    if (D == 2) {
      J3 = J3 - 1;
    } else {
      J1 = J1 + 1;
    }
    
    J2 = 0;
    NU = 0;
  }
  
  return SUM;
}


