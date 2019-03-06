#include "runner.h"
#include <stdio.h>

namespace {
double EvaluatePoint(const std::string &path, const Point &point) {
  for (size_t i = 0; i < point.size(); i++) {
    setenv(("VAR" + std::to_string(i)).c_str(),
           std::to_string(point[i]).c_str(), 1);
  }
  FILE *proc = popen(path.c_str(), "r");
  if (proc == nullptr) {
    perror("popen");
    exit(1);
  }
  double val = 0;
  if (fscanf(proc, "%lf", &val) != 1) {
    return 1e100;
  }
  int ret = pclose(proc);
  if (ret == -1) {
    perror("pclose");
    exit(1);
  }
  if (ret != 0) {
    return 1e99;
  }
  return val;
}
} // namespace

std::vector<double> Runner::Run(const std::vector<Point> &points) {
  std::vector<double> ret(points.size(), 0);

  //#pragma omp parallel for
  for (size_t i = 0; i < points.size(); i++) {
    ret[i] = EvaluatePoint(path, points[i]);
  }

  num_runs += points.size();

  for (size_t i = 0; i < points.size(); i++) {
    if (ret[i] < best_val) {
      best_val = ret[i];
      best = points[i];
    }
  }

  fprintf(stderr, "%12lu evals, best: %8.4f at ", num_runs, best_val);
  for (size_t i = 0; i < best.size(); i++) {
    fprintf(stderr, " %12.8f", best[i]);
  }
  fprintf(stderr, "\n");

  return ret;
}
