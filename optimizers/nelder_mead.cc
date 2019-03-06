#include "optimizers/base.h"

class NelderMead : public Optimizer {
public:
  using Optimizer::Optimizer;
  using PointWithScore = std::pair<double, Point>;
  static constexpr double kAlpha = 1;
  static constexpr double kGamma = 2;
  static constexpr double kRho = 0.5;
  static constexpr double kSigma = 0.5;

  void Run() override {
    std::vector<PointWithScore> points;
    points.push_back({0, Point(0.0, dim)});
    for (size_t i = 0; i < dim; i++) {
      Point p(0.0, dim);
      p[i] = 1.0;
      points.push_back({0, p});
    }
    auto evaluate_all_points = [&]() {
      std::vector<Point> pts;
      for (const PointWithScore &ps : points) {
        pts.push_back(ps.second);
      }
      auto res = runner->Run(pts);
      for (size_t i = 0; i < pts.size(); i++) {
        points[i].first = res[i];
      }
    };
    evaluate_all_points();
    double inv_dim = 1.0 / dim;
    while (true) {
      std::sort(points.begin(), points.end(),
                [](const auto &a, const auto &b) { return a.first < b.first; });
      Point centroid(0.0, dim);
      // Skip worst point.
      for (size_t j = 0; j < dim; j++) {
        centroid += points[j].second;
      }
      centroid *= inv_dim;
      // Reflection
      Point reflected = centroid + kAlpha * (centroid - points.back().second);
      double rval = runner->Run({reflected})[0];
      if (points.front().first <= rval && rval < points[dim - 1].first) {
        points.back().first = rval;
        points.back().second = reflected;
        continue;
      }
      // Expansion
      if (rval < points.front().first) {
        Point expanded = centroid + kGamma * (reflected - centroid);
        double eval = runner->Run({expanded})[0];
        if (eval < rval) {
          points.back().first = eval;
          points.back().second = expanded;

        } else {
          points.back().first = rval;
          points.back().second = reflected;
        }
        continue;
      }
      // Contraction
      Point contracted = centroid + kRho * (points.back().second - centroid);
      double cval = runner->Run({contracted})[0];
      if (cval < points.back().first) {
        points.back().first = cval;
        points.back().second = contracted;
        continue;
      }
      // Shrink
      for (size_t j = 1; j < points.size(); j++) {
        points[j].second = points.front().second +
                           kSigma * (points[j].second - points.front().second);
      }
      evaluate_all_points();
    }
  }
};

namespace {
RegisterOptimizer<NelderMead> r("nelder_mead");
}
