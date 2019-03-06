#include <random>
#include <valarray>

using Vector = std::valarray<double>;

class Matrix {
public:
  explicit Matrix(size_t N) : data(N * N), N(N) {}

  static Matrix Eye(size_t N) { return Diag(Vector(1.0, N)); }

  static Matrix Diag(const Vector &vec) {
    Matrix ret(vec.size());
    ret.diag() = vec;
    return ret;
  }

  template <typename URNG> static Matrix Rand(URNG &&rng, size_t N) {
    Matrix ret(N);
    Vector rnd_temp(N);
    std::uniform_real_distribution dist(0.0, 1.0);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        rnd_temp[j] = dist(rng);
      }
      ret[i] = rnd_temp;
    }
    return ret;
  }

  size_t size() const { return N; }

  std::slice_array<double> diag() { return data[std::slice(0, N, N + 1)]; }

  double &operator()(size_t i, size_t j) { return data[i * N + j]; }
  double operator()(size_t i, size_t j) const { return data[i * N + j]; }

  // i-th row
  std::slice_array<double> operator[](size_t i) {
    return data[std::slice(i * N, N, 1)];
  }
  std::slice_array<double> operator[](size_t i) const {
    return (*const_cast<Matrix *>(this))[i];
  }

  std::slice_array<double> col(size_t i) { return data[std::slice(i, N, N)]; }
  std::slice_array<double> col(size_t i) const {
    return const_cast<Matrix *>(this)->col(i);
  }

  std::pair<Vector, Matrix> eigs() const;

  Matrix transpose() const {
    Matrix ret(N);
    for (size_t i = 0; i < N; i++) {
      ret[i] = col(i);
    }
    return ret;
  }

  Matrix make_symmetric() const { return (*this + transpose()) * 0.5; }

  bool is_symmetric() const {
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
      (*this)[i] *= Vector(x, N);
    }
    return *this;
  }

  Matrix operator/(double x) const { return (*this) * (1.0 / x); }

  Matrix &operator/=(double x) { return *this *= 1.0 / x; }

  Vector operator*(const Vector &x) const {
    Vector ret(N);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        ret[i] += (*this)(i, j) * x[j];
      }
    }
    return ret;
  }

  Matrix operator*(const Matrix &other) {
    Matrix ret(N);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        for (size_t k = 0; k < N; k++) {
          ret(i, j) += (*this)(i, k) * other(k, j);
        }
      }
    }
    return ret;
  }

  Matrix &operator-=(const Matrix &other) {
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
};
