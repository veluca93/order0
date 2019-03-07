#include "linear_algebra.h"

namespace {
constexpr double kEps = 1e-6;
}

std::pair<Matrix, Matrix> Matrix::trid() const {
  assert(is_symmetric());
  return hessemberg();
}

std::pair<Matrix, Matrix> Matrix::hessemberg() const {
  Matrix hess = *this;
  Matrix base = Matrix::Eye(N);
  for (size_t k = 0; k < N - 1; k++) {
    double s = 0;
    for (size_t i = k + 1; i < N; i++) {
      s += hess(i, k) * hess(i, k);
    }
    s = std::sqrt(s);
    if (s < kEps) {
      continue;
    }
    double sign = hess(k + 1, k) < 0 ? -1 : 1;
    double z = 0.5 * (1 + sign * hess(k + 1, k) / s);
    Matrix v(N, 1);
    v(k + 1, 0) = std::sqrt(z);
    for (size_t i = k + 2; i < N; i++) {
      v(i, 0) = sign * hess(i, k) / (2 * v(k + 1, 0) * s);
    }
    base.rhouseholder(v);
    hess.rhouseholder(v);
    hess.lhouseholder(v);
  }
  return {hess, base};
}

std::pair<Matrix, Matrix> Matrix::qr() const {
  Matrix R = *this;
  Matrix Q = Matrix::Eye(N);
  for (size_t k = 0; k < N - 1; k++) {
    double xnorm = 0;
    for (size_t i = k; i < N; i++) {
      xnorm += R(i, k) * R(i, k);
    }
    xnorm = std::sqrt(xnorm);
    double sign = R(k, k) < 0 ? -1 : 1;
    Matrix v(N, 1);
    v(k, 0) = R(k, k) + sign * xnorm;
    for (size_t i = k + 1; i < N; i++) {
      v(i, 0) = R(i, k);
    }
    v /= Norm(v.col(0));
    Q.rhouseholder(v);
    R.lhouseholder(v);
  }
  return {Q, R};
}

void Matrix::hess_qr(std::vector<double> *A, std::vector<double> *B) {
  A->resize(N - 1);
  B->resize(N - 1);
  for (size_t k = 0; k < N - 1; k++) {
    double xnorm = 0;
    xnorm += (*this)(k, k) * (*this)(k, k);
    xnorm += (*this)(k + 1, k) * (*this)(k + 1, k);
    xnorm = std::sqrt(xnorm);
    double sign = (*this)(k, k) < 0 ? -1 : 1;
    double a = (*this)(k, k) + sign * xnorm;
    double b = (*this)(k + 1, k);
    double norm = std::sqrt(a * a + b * b);
    (*A)[k] = a / norm;
    (*B)[k] = b / norm;
    this->lhouseholder2(k, (*A)[k], (*B)[k]);
  }
}

std::pair<Vector, Matrix> Matrix::eigs() const {
  assert(is_square());
  bool sym = is_symmetric();
  auto hs = hessemberg();
  Matrix eigvecs = hs.second;
  Matrix eigvals = hs.first;
  // Shifted QR method.
  std::vector<double> A, B;
  while (true) {
    double offd_norm = 0;
    if (sym) {
      for (size_t i = 0; i < N - 1; i++) {
        offd_norm = std::max(offd_norm, std::abs(eigvals(i, i + 1)));
      }
    } else {
      for (size_t i = 0; i < N; i++) {
        if (i != 0) {
          offd_norm = std::max(offd_norm, std::abs(eigvals(i, i - 1)));
        }
        for (size_t j = i + 1; j < N; j++) {
          offd_norm = std::max(offd_norm, std::abs(eigvals(i, j)));
        }
      }
    }
    if (offd_norm < kEps)
      break;
    size_t last = N - 1;
    while (last > 0 && std::abs(eigvals(last, last - 1)) < kEps) {
      last--;
    }
    assert(N != 1);
    double mu = eigvals(last, last);
    for (size_t i = 0; i < N; i++) {
      eigvals(i, i) -= mu;
    }
    eigvals.hess_qr(&A, &B);
    for (size_t i = 0; i < N - 1; i++) {
      eigvecs.rhouseholder2(i, A[i], B[i]);
      eigvals.rhouseholder2(i, A[i], B[i]);
    }
    for (size_t i = 0; i < N; i++) {
      eigvals(i, i) += mu;
    }
  }
  return {eigvals.diag(), std::move(eigvecs)};
}
