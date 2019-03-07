#pragma once
#include "linear_algebra.h"
#include <random>

class MultivariateNormalDistribution {
public:
  MultivariateNormalDistribution(Vector mean, const Matrix &covariance)
      : mean(std::move(mean)), mul(0, 0) {
    auto [eigenval, eigenvec] = covariance.make_symmetric().eigs();
    mul = eigenvec * Matrix::Diag(std::sqrt(eigenval));
  }

  MultivariateNormalDistribution(Vector mean, const Vector &stddev,
                                 const Matrix &basis)
      : mean(std::move(mean)), mul(basis * Matrix::Diag(stddev)) {}

  template <typename URNG> Vector operator()(URNG &&gen) {
    Vector z(mean.size());
    for (size_t i = 0; i < mean.size(); i++) {
      z[i] = dist(gen);
    }
    return mean + mul * z;
  }

private:
  std::normal_distribution<double> dist{0.0, 1.0};
  Vector mean;
  Matrix mul;
};
