#include "multivariate_normal.h"
#include "optimizers/base.h"

class CMAES : public Optimizer {
public:
  using Optimizer::Optimizer;
  using PointWithScore = std::pair<double, Point>;
  static constexpr double kInitialSigma = 2.0;

  void Run() override {
    std::vector<PointWithScore> points;
    size_t N = dim;
    Point xmean(0.0, N);
    double sigma = kInitialSigma;

    // Points to evaluate
    size_t lambda = 4 + std::floor(3 * std::log(N));

    // Points to use for updates
    double mu = lambda / 2.0;

    // Weight array
    Vector weights(0.0, size_t(mu));
    for (size_t i = 0; i < size_t(mu); i++) {
      weights[i] = log(mu + 0.5) - log(i + 1);
    }
    weights /= weights.sum();

    // Variance-effectiveness of \sum w_i x_i
    double mueff = weights.sum() * weights.sum() / ((weights * weights).sum());

    // C decay constant
    double cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N);

    // Sigma decay constant
    double cs = (mueff + 2) / (N + mueff + 5);

    // Learning rate for rank-one C updates
    double c1 = 2 / ((N + 1.3) * (N + 1.3) + mueff);

    // and rank-mu C updates
    double cmu = std::min(1 - c1, 2 * (mueff - 2 + 1 / mueff) /
                                      ((N + 2) * (N + 2) + mueff));

    // Sigma dampening
    double damps =
        1 + 2 * std::max(0., std::sqrt((mueff - 1) / (N + 1)) - 1) + cs;

    // Evolution paths for C and sigma
    Vector pc(0.0, N);
    Vector ps(0.0, N);

    // Coordinate system
    Matrix B = Matrix::Eye(N);

    // Dimension scaling
    Vector D(1.0, N);

    // Covariance matrix
    Matrix C = B * Matrix::Diag(D * D) * B.transpose();
    Matrix invsqrtC = B * Matrix::Diag(1.0 / D) * B.transpose();

    // Approx. expected value of ||N(0, I)||
    const double chiN =
        std::sqrt(1.0 * N) * (1 - 1 / (4. * N) + 1 / (21. * N * N));

    size_t counteval = 0;

    std::mt19937 rng;
    auto evaluate_all_points = [&]() {
      points.resize(lambda);
      MultivariateNormalDistribution dist(xmean, C * (sigma * sigma));
      for (size_t i = 0; i < lambda; i++) {
        points[i].second = dist(rng);
        counteval++;
      }
      std::vector<Point> pts;
      for (const PointWithScore &ps : points) {
        pts.push_back(ps.second);
      }
      auto res = runner->Run(pts);
      for (size_t i = 0; i < pts.size(); i++) {
        points[i].first = res[i];
      }
      std::sort(points.begin(), points.end(),
                [](const auto &a, const auto &b) { return a.first < b.first; });
    };

    while (true) {
      evaluate_all_points();

      // Update mean
      Vector xold = xmean;
      xmean = Vector(0.0, N);
      for (size_t i = 0; i < size_t(mu); i++) {
        xmean += weights[i] * points[i].second;
      }

      // Update evolution paths
      ps = (1 - cs) * ps +
           invsqrtC * (xmean - xold) / sigma * std::sqrt(cs * (2 - cs) * mueff);
      bool hsig = Norm(ps) /
                      std::sqrt(1 - std::pow(1 - cs, 2. * counteval / lambda)) /
                      chiN <
                  1.4 + 2. / (N + 1);
      pc = (1 - cc) * pc +
           hsig * std::sqrt(cc * (2 - cc) * mueff) * (xmean - xold) / sigma;

      // Adapt C
      // This should be a N * mu matrix, but is N*N here to simplify
      // implementation.
      Matrix artmp(N, mu);
      for (size_t i = 0; i < size_t(mu); i++) {
        artmp.col(i) = (points[i].second - xold) / sigma;
      }

      // N*N matrix version of pc.
      Matrix pcm(N, 1);
      pcm.col(0) = pc;

      C = C * (1 - c1 - cmu) +
          (pcm * pcm.transpose() + C * (1 - hsig) * cc * (2 - cc)) * c1 +
          artmp * Matrix::Diag(weights) * artmp.transpose() * cmu;

      // Adapt sigma
      sigma = sigma * std::exp((cs / damps) * (Norm(ps) / chiN - 1));

      // Decomposition of C.
      // TODO: this makes the algorithm have higher complexity.
      C = C.make_symmetric();
      auto [eigvals, eigvecs] = C.eigs();
      D = std::sqrt(eigvals);
      invsqrtC = B * Matrix::Diag(1.0 / D) * B.transpose();
    }
  }
};

namespace {
RegisterOptimizer<CMAES> r("cma_es");
}
