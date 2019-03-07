#include "linear_algebra.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>

namespace {

void CheckMatrixNear(const Matrix &a, const Matrix &b, double tol = 1e-5) {
  EXPECT_EQ(a.size(), b.size());
  if (a.size() != b.size())
    return;
  size_t N = a.size().first;
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
  Matrix m(3, 3);
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

TEST(LinearAlgebraTest, HouseholderMul) {
  for (size_t N = 1; N < 300; N += 30) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, N);
    Matrix v = Matrix::Rand(rng, N, 1);
    Matrix H = Matrix::Eye(N) - v * v.transpose() * 2;
    Matrix mH = m;
    mH.rhouseholder(v);
    Matrix Hm = m;
    Hm.lhouseholder(v);
    CheckMatrixNear(H * m, Hm);
    CheckMatrixNear(m * H, mH);
  }
}

TEST(LinearAlgebraTest, Householder2MulR) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (size_t N = 1; N < 300; N += 30) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, N);
    double a = dist(rng);
    double b = dist(rng);
    for (size_t i = 0; i < N - 1; i++) {
      Matrix v(N, 1);
      v[i] = a;
      v[i + 1] = b;
      Matrix mH = m;
      Matrix mH2 = m;
      mH.rhouseholder(v);
      mH2.rhouseholder2(i, a, b);
      CheckMatrixNear(mH2, mH);
    }
  }
}

TEST(LinearAlgebraTest, Householder2MulL) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (size_t N = 1; N < 300; N += 30) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, N);
    double a = dist(rng);
    double b = dist(rng);
    for (size_t i = 0; i < N - 1; i++) {
      Matrix v(N, 1);
      v[i] = a;
      v[i + 1] = b;
      Matrix Hm = m;
      Matrix Hm2 = m;
      Hm.lhouseholder(v);
      Hm2.lhouseholder2(i, a, b);
      CheckMatrixNear(Hm2, Hm);
    }
  }
}

TEST(LinearAlgebraTest, HessembergTest) {
  for (size_t N = 1; N < 300; N += 30) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, N);
    auto [hess, base] = m.hessemberg();
    // Check orthogonality of the base.
    CheckMatrixNear(Matrix::Eye(N), base * base.transpose());
    // Check similarity with original matrix.
    CheckMatrixNear(m, base * hess * base.transpose());
    // Check that hess is hessemberg.
    for (size_t i = 0; i < N; i++) {
      for (size_t j = i + 2; j < N; j++) {
        EXPECT_NEAR(hess(j, i), 0.0, 1e-6);
      }
    }
  }
}

TEST(LinearAlgebraTest, TridiagTest) {
  for (size_t N = 1; N < 300; N += 30) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, N).make_symmetric();
    auto [trid, base] = m.trid();
    // Check orthogonality of the base.
    CheckMatrixNear(Matrix::Eye(N), base * base.transpose());
    // Check similarity with original matrix.
    CheckMatrixNear(m, base * trid * base.transpose());
    // Check that trid is still symmetric.
    CheckMatrixNear(trid, trid.make_symmetric());
    // Check that trid is tridiagonal.
    for (size_t i = 0; i < N; i++) {
      for (size_t j = i + 2; j < N; j++) {
        EXPECT_NEAR(trid(i, j), 0.0, 1e-6);
      }
    }
  }
}

TEST(LinearAlgebraTest, QRTest) {
  for (size_t N = 1; N < 300; N += 30) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, N);
    auto [Q, R] = m.qr();
    // Check orthogonality of Q.
    CheckMatrixNear(Matrix::Eye(N), Q * Q.transpose());
    // Check Q*R ~ m.
    CheckMatrixNear(m, Q * R);
    // Check that R is upper triangular.
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < i; j++) {
        EXPECT_NEAR(R(i, j), 0.0, 1e-6);
      }
    }
  }
}

TEST(LinearAlgebraTest, HessQRTest) {
  for (size_t N = 1; N < 300; N += 30) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, N).hessemberg().first;
    auto [A, B, R] = m.hess_qr();
    ASSERT_EQ(A.size(), N - 1);
    ASSERT_EQ(B.size(), N - 1);
    // Check that R is upper triangular.
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < i; j++) {
        EXPECT_NEAR(R(i, j), 0.0, 1e-6);
      }
    }
    // Check that left-multiplying all the householder vectors given by A and B
    // to R2 gives back m.
    for (size_t i = N - 1; i > 0; i--) {
      R.lhouseholder2(i - 1, A[i - 1], B[i - 1]);
    }
    CheckMatrixNear(m, R);
  }
}

TEST(LinearAlgebraTest, EigsTestSmall) {
  std::mt19937 rng;
  Matrix m(4, 4);
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
  for (size_t N = 1; N < 200; N += 40) {
    std::mt19937 rng;
    Matrix m = Matrix::Rand(rng, N).make_symmetric();
    auto [eigvals, eigvecs] = m.eigs();
    Matrix diag = Matrix::Diag(eigvals);
    CheckMatrixNear(Matrix::Eye(N), eigvecs * eigvecs.transpose());
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
