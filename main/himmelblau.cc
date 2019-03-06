#include <array>
#include <stdio.h>
#include <stdlib.h>

int main() {
  std::array<double, 2> vars;
  for (size_t i = 0; i < 2; i++) {
    vars[i] = atof(getenv(("VAR" + std::to_string(i)).c_str()));
  }
  double value = 0;
  double x = vars[0];
  double y = vars[1];
  double t1 = x * x + y - 11;
  double t2 = x + y * y - 7;
  value += t1 * t1 + t2 * t2;
  printf("%f\n", value);
}
