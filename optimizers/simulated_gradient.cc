#include "optimizers/base.h"
#include <math.h>

class SimulatedGradient : public Optimizer {
public:
  using Optimizer::Optimizer;
  using PointWithScore = std::pair<double, Point>;
  static constexpr double kGradientStep = 1e-5;
  static constexpr double kLearnRate = 1e-2;
  static constexpr double kLearnRateDecayExp = 0.3;
  static constexpr double kGradientClip = 1.0;
  static constexpr double kInvGradientStep = 1.0 / kGradientStep;

  static constexpr double kB1 = 0.9;
  static constexpr double kB2 = 0.999;
  static constexpr double kEps = 1e-8;

  void Run() override {
    Point p(dim, 0.0);
    Point momentum(dim, 0.0);
    Point variance(dim, 0.0);
    double b1pow = 1;
    double b2pow = 1;

    size_t iter = 0;
    while (b1pow *= kB1, b2pow *= kB2, ++iter) {
      double val = runner->Run({p})[0];
      PrintBest(p, val);
      std::vector<Point> gradient_helpers(dim, p);
      for (size_t i = 0; i < dim; i++) {
        gradient_helpers[i][i] += kGradientStep;
      }
      std::vector<double> vals = runner->Run(gradient_helpers);
      Point estimated_gradient(dim, val);
      for (size_t i = 0; i < dim; i++) {
        estimated_gradient[i] = (vals[i] - val) * kInvGradientStep;
        if (estimated_gradient[i] > kGradientClip)
          estimated_gradient[i] = kGradientClip;
        if (estimated_gradient[i] < -kGradientClip)
          estimated_gradient[i] = -kGradientClip;
      }

      // ADAM
      for (size_t i = 0; i < dim; i++) {
        momentum[i] = kB1 * momentum[i] + (1 - kB1) * estimated_gradient[i];
        variance[i] = kB2 * variance[i] +
                      (1 - kB2) * estimated_gradient[i] * estimated_gradient[i];
      }
      double m_scale = 1.0 / (1.0 - b1pow);
      double v_scale = 1.0 / (1.0 - b2pow);

      const double learn_rate = kLearnRate / pow(iter, kLearnRateDecayExp);
      for (size_t i = 0; i < dim; i++) {
        p[i] -= learn_rate * momentum[i] * m_scale /
                (std::sqrt(variance[i] * v_scale) + kEps);
      }
    }
  }
};

namespace {
RegisterOptimizer<SimulatedGradient> r("simulated_gradient");
}
