#include <array>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

int main() {
  std::array<double, 2> vars;
  for (size_t i = 0; i < 2; i++) {
    vars[i] = atof(getenv(("VAR" + std::to_string(i)).c_str()));
  }
  double value = 0;
  double x = vars[0];
  double y = vars[1];
  double penalty = 0;
  if (x < -10) {
    penalty += -10 - x;
    x = -10;
  }
  if (x > 10) {
    penalty += x - 10;
    x = 10;
  }
  if (y < -10) {
    penalty += -10 - y;
    y = -10;
  }
  if (y > 10) {
    penalty += y - 10;
    y = 10;
  }
  value = std::sin(x) * std::cos(y);
  value *= exp(std::abs(1 - std::sqrt(x * x + y * y) / M_PI));
  value = -std::abs(value);
  printf("%f\n", value + penalty);
}
