#include "linear_algebra.h"
#include <assert.h>

#define E(i, j) eigvals.data[(i)*N + (j)]

std::pair<Vector, Matrix> Matrix::eigs() const {
  assert(is_symmetric());
  Matrix eigvecs = Matrix::Eye(N);
  Matrix eigvals = *this;
  size_t iter_recalc = 0;
  for (size_t _ = 0; _ < N * 30; _++) {
    size_t maxi = 0, maxj = 0;
    double val = 0;
    double norm = 0;
    for (size_t i = 0; i < N; i++) {
      for (size_t j = i + 1; j < N; j++) {
        if (std::abs(E(i, j)) > val) {
          maxi = i;
          maxj = j;
          val = std::abs(E(i, j));
        }
        norm += E(i, j) * E(i, j);
      }
    }
    if (norm < 1e-10)
      break;
    double w = (E(maxj, maxj) - E(maxi, maxi)) / (2 * E(maxi, maxj));
    double t;
    if (w > 0) {
      t = 1.0 / (std::sqrt(w * w + 1) + w);
    } else {
      t = -1.0 / (std::sqrt(w * w + 1) - w);
    }
    double sqrt1t2 = std::sqrt(1 + t * t);
    double s = t / sqrt1t2;
    double c = 1 / sqrt1t2;
    double tau = (sqrt1t2 + 1) / (sqrt1t2 + 1 + t * t);
    // Update eigvals
    E(maxi, maxi) += t * val;
    E(maxj, maxj) -= t * val;
    E(maxi, maxj) = 0;
    E(maxj, maxi) = 0;
    for (size_t k = 0; k < N; k++) {
      if (k == maxi || k == maxj)
        continue;
      double aik = E(maxi, k);
      double ajk = E(maxj, k);
      E(k, maxi) = E(maxi, k) = tau * aik - s * ajk;
      E(k, maxj) = E(maxj, k) = tau * ajk + s * aik;
    }
    // Update eigvecs
    Vector coli = eigvecs.col(maxi);
    Vector colj = eigvecs.col(maxj);
    eigvecs.col(maxi) = c * coli - s * colj;
    eigvecs.col(maxj) = s * coli + c * colj;
    // Hackish hack to fix numerical instability issues.
    if (iter_recalc++ > 25) {
      iter_recalc = 0;
      eigvals = eigvecs.transpose() * (*this) * eigvecs;
    }
  }
  eigvals = eigvecs.transpose() * (*this) * eigvecs;
  return {eigvals.diag(), std::move(eigvecs)};
}
