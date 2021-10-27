#include <array>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#ifndef N
#define N 2
#endif

// Shifted to have a minimum in (M_PI/1, M_PI/2, .... M_PI/N).

int main() {
  std::array<double, N> vars;
  for (size_t i = 0; i < N; i++) {
    vars[i] = atof(getenv(("VAR" + std::to_string(i)).c_str()));
  }
  constexpr double A = 10;
  double value = A * N;
  for (size_t i = 0; i < N; i++) {
    double x = vars[i] - M_PI / (i + 1);
    value += x * x - A * std::cos(2 * M_PI * x);
  }
  printf("%f\n", value);
}
