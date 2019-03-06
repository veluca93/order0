#include "multivariate_normal.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

TEST(MultivariateNormalTest, Check1d) {
  Vector mean_v{5};
  Matrix variance_m(1);
  variance_m(0, 0) = 3;
  std::mt19937 rng;
  std::vector<double> samples;
  MultivariateNormalDistribution dist(mean_v, variance_m);
  size_t kNumSamples = 1 << 22;
  for (size_t i = 0; i < kNumSamples; i++) {
    samples.push_back(dist(rng)[0]);
  }
  double mean = 0;
  for (size_t i = 0; i < kNumSamples; i++) {
    mean += samples[i];
  }
  mean /= kNumSamples;
  EXPECT_NEAR(mean, mean_v[0], 1e-3);
  double variance = 0;
  for (size_t i = 0; i < kNumSamples; i++) {
    double v = samples[i] - mean;
    variance += v * v;
  }
  variance /= kNumSamples;
  EXPECT_NEAR(variance, variance_m(0, 0), 1e-2);
}

TEST(MultivariateNormalTest, CheckCovariance) {
  Vector mean_v{-2, -1, 0, 1, 2};
  std::mt19937 rng(123);
  size_t size = mean_v.size();
  Matrix variance_m = Matrix::Rand(rng, size);
  variance_m = variance_m * variance_m.transpose();
  variance_m = variance_m.make_symmetric();
  std::vector<Vector> samples;
  MultivariateNormalDistribution dist(mean_v, variance_m);
  size_t kNumSamples = 1 << 22;
  for (size_t i = 0; i < kNumSamples; i++) {
    samples.push_back(dist(rng));
  }
  Vector mean(0.0, size);
  for (size_t i = 0; i < size; i++) {
    for (size_t k = 0; k < kNumSamples; k++) {
      mean[i] += samples[k][i];
    }
    mean[i] /= kNumSamples;
    EXPECT_NEAR(mean[i], mean_v[i], 1e-3);
  }
  for (size_t i = 0; i < size; i++) {
    for (size_t j = i; j < size; j++) {
      double covariance = 0;
      for (size_t k = 0; k < kNumSamples; k++) {
        double v = samples[k][i] - mean[i];
        double w = samples[k][j] - mean[j];
        covariance += v * w;
      }
      covariance /= kNumSamples;
      EXPECT_NEAR(covariance, variance_m(i, j), 1e-3);
    }
  }
}

} // namespace
