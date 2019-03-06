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
    Point p(0.0, dim);
    Point momentum(0.0, dim);
    Point variance(0.0, dim);
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
      Point estimated_gradient(val, dim);
      for (size_t i = 0; i < dim; i++) {
        estimated_gradient[i] = (vals[i] - val) * kInvGradientStep;
        if (estimated_gradient[i] > kGradientClip)
          estimated_gradient[i] = kGradientClip;
        if (estimated_gradient[i] < -kGradientClip)
          estimated_gradient[i] = -kGradientClip;
      }

      // ADAM
      momentum = kB1 * momentum + (1 - kB1) * estimated_gradient;
      variance =
          kB2 * variance + (1 - kB2) * estimated_gradient * estimated_gradient;
      double m_scale = 1.0 / (1.0 - b1pow);
      double v_scale = 1.0 / (1.0 - b2pow);

      const double learn_rate = kLearnRate / pow(iter, kLearnRateDecayExp);
      p -= learn_rate * momentum * m_scale /
           (std::sqrt(variance * v_scale) + kEps);
    }
  }
};

namespace {
RegisterOptimizer<SimulatedGradient> r("simulated_gradient");
}
