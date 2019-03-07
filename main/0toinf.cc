#include <array>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef N
#define N 2
#endif

int main() {
  std::array<double, N> vars;
  for (size_t i = 0; i < N; i++) {
    vars[i] = atof(getenv(("VAR" + std::to_string(i)).c_str()));
  }
  double x = 0;
  for (size_t i = 0; i < N; i++) {
    x += vars[i] * vars[i];
  }
  x = std::sqrt(x);
  double value = 1 / (std::abs(x) + 1);
  value += 5 * (std::sin(x * x) + 1) / (x * x + 2);
  printf("%f\n", value);
}
