#include <array>
#include <stdio.h>
#include <stdlib.h>

#ifndef N
#define N 2
#endif

int main() {
  static_assert(N % 2 == 0, "N must be even!");
  std::array<double, N> vars;
  for (size_t i = 0; i < N; i++) {
    vars[i] = atof(getenv(("VAR" + std::to_string(i)).c_str()));
  }
  double value = 0;
  for (size_t i = 0; i < N / 2; i++) {
    double x = vars[2 * i];
    double y = vars[2 * i + 1];
    double t1 = x * x - y;
    double t2 = x - 1;
    value += 100 * t1 * t1 + t2 * t2;
  }
  printf("%f\n", value);
}
