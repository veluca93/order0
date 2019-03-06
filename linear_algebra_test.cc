#include "linear_algebra.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>

namespace {

void CheckMatrixNear(const Matrix &a, const Matrix &b, double tol = 1e-5) {
  EXPECT_EQ(a.size(), b.size());
  if (a.size() != b.size())
    return;
  size_t N = a.size();
  double max_err = 0;
  for (size_t i = 0; i < N; i++) {
    max_err = std::max(max_err, std::abs(Vector(a[i]) - Vector(b[i])).max());
  }
  EXPECT_LT(max_err, tol);
}

void CheckVectorNear(const Vector &a, const Vector &b, double tol = 1e-5) {
  EXPECT_EQ(a.size(), b.size());
  if (a.size() != b.size())
    return;
  double max_err = std::abs(a - b).max();
  EXPECT_LT(max_err, tol);
}

TEST(LinearAlgebraTest, MatrixSumDiffInverse) {
  std::mt19937 rng;
  Matrix m1 = Matrix::Rand(rng, 20);
  Matrix m2 = Matrix::Rand(rng, 20);
  Matrix sum = m1 + m2;
  Matrix orig = sum - m2;
  CheckMatrixNear(orig, m1);
}

TEST(LinearAlgebraTest, MatrixVectorProduct) {
  Matrix m(3);
  m[0] = {1, 2, 3};
  m[1] = {3, 4, 5};
  m[2] = {5, 6, 7};
  Vector v{0.1, 0.2, 0.3};
  Vector result{1.4, 2.6, 3.8};
  CheckVectorNear(m * v, result);
}

TEST(LinearAlgebraTest, MatrixMatrixProduct) {
  std::mt19937 rng;
  Matrix m1 = Matrix::Rand(rng, 20);
  Matrix m2 = Matrix::Rand(rng, 20);
  Vector vec(0.1, 20);
  Vector p1 = m1 * (m2 * vec);
  Vector p2 = (m1 * m2) * vec;
  CheckVectorNear(p1, p2);
}

TEST(LinearAlgebraTest, MatrixEyeProduct) {
  std::mt19937 rng;
  Matrix m = Matrix::Rand(rng, 20);
  Matrix eye = Matrix::Eye(20);
  CheckMatrixNear(m, m * eye);
  CheckMatrixNear(m, eye * m);
}

TEST(LinearAlgebraTest, EigsTestSmall) {
  std::mt19937 rng;
  Matrix m(4);
  m[0] = {2, -1, -1, 0};
  m[1] = {-1, 3, -1, -1};
  m[2] = {-1, -1, 3, -1};
  m[3] = {0, -1, -1, 2};
  auto [eigvals, eigvecs] = m.eigs();
  Matrix diag = Matrix::Diag(eigvals);
  CheckMatrixNear(Matrix::Eye(4), eigvecs * eigvecs.transpose());
  CheckMatrixNear(m, eigvecs * diag * eigvecs.transpose());
}

TEST(LinearAlgebraTest, EigsTest) {
  for (size_t i = 1; i < 24; i++) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, i).make_symmetric();
    auto [eigvals, eigvecs] = m.eigs();
    Matrix diag = Matrix::Diag(eigvals);
    CheckMatrixNear(Matrix::Eye(i), eigvecs * eigvecs.transpose());
    CheckMatrixNear(m, eigvecs * diag * eigvecs.transpose());
  }
}

TEST(LinearAlgebraTest, TestMakeSymmetric) {
  std::mt19937 rng;
  Matrix m1 = Matrix::Rand(rng, 20);
  Matrix sym = m1.make_symmetric();
  EXPECT_TRUE(sym.is_symmetric());
}

} // namespace
