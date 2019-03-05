#pragma once
#include <string>
#include <vector>

using Point = std::vector<double>;

class Runner {
public:
  Runner(std::string path) : path(path) {}
  std::vector<double> Run(const std::vector<Point> &points);
  size_t NumRuns() const { return num_runs; }

private:
  std::string path;
  size_t num_runs = 0;
};
