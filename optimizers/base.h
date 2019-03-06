#pragma once
#include "runner.h"
#include <functional>
#include <memory>
#include <unordered_map>

class Optimizer {
public:
  virtual void Run() = 0;
  Optimizer(Runner *runner, size_t dim) : runner(runner), dim(dim) {}
  virtual ~Optimizer() = default;

  static std::unique_ptr<Optimizer> Create(const std::string &name,
                                           Runner *runner, size_t dim) {
    if (!Registry().count(name)) {
      fprintf(stderr, "Unknown optimizer!\n");
      exit(1);
    }
    return Registry()[name](runner, dim);
  }

  static void
  Register(const std::string &name,
           std::function<std::unique_ptr<Optimizer>(Runner *, size_t)> fun) {
    Registry().emplace(name, std::move(fun));
  }

protected:
  Runner *runner;
  size_t dim;

private:
  using reg_t = std::unordered_map<
      std::string, std::function<std::unique_ptr<Optimizer>(Runner *, size_t)>>;
  static reg_t &Registry() {
    static reg_t r;
    return r;
  }
};

template <typename T> class RegisterOptimizer {
public:
  RegisterOptimizer(const std::string &name) {
    Optimizer::Register(
        name, [](Runner *r, size_t d) { return std::make_unique<T>(r, d); });
  }
};
