#include "optimizers/base.h"

class Annealing : public Optimizer {
public:
  using Optimizer::Optimizer;
  using PointWithScore = std::pair<double, Point>;
  static constexpr double kNeighSigma = 1.0;
  static constexpr double kTemperatureDecay = 0.95;
  static constexpr double kInitialTemperature = 100;
  static constexpr size_t kNumPoints = 4;

  void Run() override {
    auto evaluate_all_points = [this](std::vector<PointWithScore> &points) {
      std::vector<Point> pts;
      for (const PointWithScore &ps : points) {
        pts.push_back(ps.second);
      }
      auto res = runner->Run(pts);
      for (size_t i = 0; i < pts.size(); i++) {
        points[i].first = res[i];
      }
    };
    std::mt19937 rng;
    std::normal_distribution dist(0.0, kNeighSigma);
    auto add_noise = [&](Point p) {
      for (double &d : p) {
        d += dist(rng);
      }
      return p;
    };
    std::vector<PointWithScore> points;
    for (size_t i = 0; i < kNumPoints; i++) {
      points.emplace_back(0.0, add_noise(Point(0.0, dim)));
    }
    evaluate_all_points(points);
    double temperature = kInitialTemperature;
    while (true) {
      auto new_points = points;
      for (auto &point : new_points) {
        point.second = add_noise(point.second);
      }
      evaluate_all_points(new_points);
      for (size_t i = 0; i < kNumPoints; i++) {
        double prob =
            std::exp(-(new_points[i].first - points[i].first) / temperature);
        if (prob > 1 || std::bernoulli_distribution(prob)(rng)) {
          points[i] = new_points[i];
        }
      }
      temperature *= kTemperatureDecay;
    }
  }
};

namespace {
RegisterOptimizer<Annealing> r("annealing");
}
