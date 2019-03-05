#include "optimizers/base.h"

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s N target optimizer\n", argv[0]);
    exit(1);
  }
  Runner runner(argv[2]);
  auto optimizer = Optimizer::Create(argv[3], &runner, atoll(argv[1]));
  optimizer->Run();
}
