#include <array>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

// Shifted to have the minimum in M_E/4,M_PI.

int main() {
  std::array<double, 2> vars;
  for (size_t i = 0; i < 2; i++) {
    vars[i] = atof(getenv(("VAR" + std::to_string(i)).c_str()));
  }
  double value = 0;
  double x = vars[0] - M_E / 4;
  double y = vars[1] - M_PI;
  value = -20 * std::exp(-0.2 * std::sqrt(0.5 * (x * x + y * y)));
  value += M_E + 20;
  value -= std::exp(0.5 * (std::cos(2 * M_PI * x) + std::cos(2 * M_PI * y)));
  printf("%f\n", value);
}
