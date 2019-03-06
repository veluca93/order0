#pragma once
#include "linear_algebra.h"
#include <string>
#include <vector>

using Point = Vector;

class Runner {
public:
  Runner(std::string path) : path(path) {}
  std::vector<double> Run(const std::vector<Point> &points);

private:
  std::string path;
  size_t num_runs = 0;
  Point best;
  double best_val = 1e101;
};
