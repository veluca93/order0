#include "optimizers/base.h"

class LineSearch : public Optimizer {
public:
  using Optimizer::Optimizer;
  void Run() override {
    Point best(0.0, dim);
    double best_val = runner->Run({best})[0];
    double last_best = 0;
    double lim1 = 1e-2;
    double lim2 = 1e-3;
    for (size_t i = 0;; i++) {
      if (last_best < best_val * (1 + 1e-6)) {
        lim1 *= 0.9;
        lim2 *= 0.9;
      }
      size_t current_dimension = i % dim;
      Point next;
      double next_val = 1e90;
      for (int direction : {1, -1}) {
        Point current = best;
        double delta = direction * lim1;
        double current_val = best_val;
        while (true) {
          current[current_dimension] += delta;
          delta *= 2;
          double val = runner->Run({current})[0];
          if (val >= current_val) {
            delta *= 0.5;
            current[current_dimension] -= delta;
            break;
          }
          current_val = val;
        }
        delta *= 0.5;
        while (delta > lim2) {
          Point candidate = current;
          candidate[current_dimension] += delta;
          double val = runner->Run({candidate})[0];
          if (val < current_val) {
            current = candidate;
            current_val = val;
          }
          delta *= 0.5;
        }
        if (current_val < next_val) {
          next = current;
          next_val = current_val;
        }
      }
      last_best = best_val;
      best_val = next_val;
      best = next;
      PrintBest(best, best_val);
    }
  }
};

namespace {
RegisterOptimizer<LineSearch> r("line_search");
}
