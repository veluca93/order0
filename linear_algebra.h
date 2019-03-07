#pragma once
#include <assert.h>
#include <random>
#include <valarray>

using Vector = std::valarray<double>;

inline double Norm(const Vector &vec) { return std::sqrt((vec * vec).sum()); }

class Matrix {
public:
  explicit Matrix(size_t N, size_t M) : data(N * M), N(N), M(M) {}

  static Matrix Eye(size_t N) { return Diag(Vector(1.0, N)); }

  static Matrix Diag(const Vector &vec) {
    Matrix ret(vec.size(), vec.size());
    ret.diag() = vec;
    return ret;
  }

  template <typename URNG>
  static Matrix Rand(URNG &&rng, size_t N, size_t M = 0) {
    if (M == 0)
      M = N;
    Matrix ret(N, M);
    Vector rnd_temp(M);
    std::uniform_real_distribution dist(0.0, 1.0);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < M; j++) {
        rnd_temp[j] = dist(rng);
      }
      ret[i] = rnd_temp;
    }
    return ret;
  }

  std::pair<size_t, size_t> size() const { return {N, M}; }

  std::slice_array<double> diag() {
    assert(N == M);
    return data[std::slice(0, N, N + 1)];
  }

  double &operator()(size_t i, size_t j) { return data[i * M + j]; }
  double operator()(size_t i, size_t j) const { return data[i * M + j]; }

  // i-th row
  std::slice_array<double> operator[](size_t i) {
    assert(i < N);
    return data[std::slice(i * M, M, 1)];
  }
  std::slice_array<double> operator[](size_t i) const {
    return (*const_cast<Matrix *>(this))[i];
  }

  std::slice_array<double> col(size_t i) {
    assert(i < M);
    return data[std::slice(i, N, M)];
  }
  std::slice_array<double> col(size_t i) const {
    return const_cast<Matrix *>(this)->col(i);
  }

  std::pair<Vector, Matrix> eigs() const;

  Matrix transpose() const {
    Matrix ret(M, N);
    for (size_t i = 0; i < M; i++) {
      ret[i] = col(i);
    }
    return ret;
  }

  Matrix make_symmetric() const {
    assert(M == N);
    return (*this + transpose()) * 0.5;
  }

  bool is_square() const { return M == N; }

  bool is_symmetric() const {
    assert(M == N);
    std::valarray<bool> vec(true, N);
    for (size_t i = 0; i < N; i++) {
      vec &= (Vector(col(i)) == Vector((*this)[i]));
    }
    return vec.min();
  }

  Matrix operator*(double x) const {
    Matrix ret = *this;
    ret *= x;
    return ret;
  }

  Matrix &operator*=(double x) {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] *= Vector(x, M);
    }
    return *this;
  }

  Matrix operator/(double x) const { return (*this) * (1.0 / x); }

  Matrix &operator/=(double x) { return *this *= 1.0 / x; }

  Vector operator*(const Vector &x) const {
    assert(x.size() == M);
    Vector ret(N);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < M; j++) {
        ret[i] += (*this)(i, j) * x[j];
      }
    }
    return ret;
  }

  Matrix operator*(const Matrix &other) const {
    assert(M == other.N);
    Matrix ret(N, other.M);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < other.M; j++) {
        for (size_t k = 0; k < M; k++) {
          ret(i, j) += (*this)(i, k) * other(k, j);
        }
      }
    }
    return ret;
  }

  Matrix &operator-=(const Matrix &other) {
    assert(size() == other.size());
    for (size_t i = 0; i < N; i++) {
      (*this)[i] -= other[i];
    }
    return *this;
  }

  Matrix operator-(const Matrix &other) const {
    Matrix ret = *this;
    ret -= other;
    return ret;
  }

  Matrix &operator+=(const Matrix &other) {
    assert(size() == other.size());
    for (size_t i = 0; i < N; i++) {
      (*this)[i] += other[i];
    }
    return *this;
  }

  Matrix operator+(const Matrix &other) const {
    Matrix ret = *this;
    ret += other;
    return ret;
  }

private:
  Vector data;
  size_t N;
  size_t M;
};
