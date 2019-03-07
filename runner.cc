#include "runner.h"
#include <spawn.h>
#include <stdio.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

double EvaluatePoint(const std::string &path, const Point &point) {
  int pipe_fds[2];
  if (pipe(pipe_fds) == -1) {
    perror("pipe");
    exit(1);
  }
  posix_spawn_file_actions_t actions;
  int ret = posix_spawn_file_actions_init(&actions);
  if (ret != 0)
    exit(1);
  ret = posix_spawn_file_actions_addclose(&actions, pipe_fds[0]);
  if (ret != 0)
    exit(1);
  ret = posix_spawn_file_actions_addclose(&actions, STDIN_FILENO);
  if (ret != 0)
    exit(1);
  ret = posix_spawn_file_actions_adddup2(&actions, pipe_fds[1], STDOUT_FILENO);
  if (ret != 0)
    exit(1);
  ret = posix_spawn_file_actions_addclose(&actions, pipe_fds[1]);
  if (ret != 0)
    exit(1);
  std::vector<std::vector<char>> args;
  auto add_arg = [&args](std::string s) {
    std::vector<char> arg(s.size() + 1);
    std::copy(s.begin(), s.end(), arg.begin());
    arg.back() = '\0';
    args.push_back(std::move(arg));
  };
  add_arg(path);

  std::vector<std::vector<char>> envs;
  auto add_env = [&envs](std::string k, std::string v) {
    std::string s = k + "=" + v;
    std::vector<char> env(s.size() + 1);
    std::copy(s.begin(), s.end(), env.begin());
    env.back() = '\0';
    envs.push_back(std::move(env));
  };
  for (size_t i = 0; i < point.size(); i++) {
    add_env("VAR" + std::to_string(i), std::to_string(point[i]));
  }

  std::vector<char *> args_list(args.size() + 1);
  for (size_t i = 0; i < args.size(); i++)
    args_list[i] = args[i].data();
  args_list.back() = nullptr;

  std::vector<char *> environ(envs.size() + 1);
  for (size_t i = 0; i < envs.size(); i++)
    environ[i] = envs[i].data();
  environ.back() = nullptr;

  int child_pid = 0;
  ret = posix_spawn(&child_pid, args_list[0], &actions, nullptr,
                    args_list.data(), environ.data());
  close(pipe_fds[1]);
  if (ret != 0) {
    close(pipe_fds[0]);
    exit(1);
  }
  int child_status = 0;
  if (waitpid(child_pid, &child_status, 0) == -1) {
    close(pipe_fds[0]);
    exit(1);
  }
  if (child_status != 0) {
    return 1e99;
  }
  char val_buf[1024] = {};
  if (read(pipe_fds[0], val_buf, 1024) == -1) {
    close(pipe_fds[0]);
    return 1e100;
  }
  close(pipe_fds[0]);
  double val;
  int sscanf_ret = sscanf(val_buf, "%lf", &val);
  if (sscanf_ret != 1) {
    return 1e99;
  }
  return val;
}
} // namespace

std::vector<double> Runner::Run(const std::vector<Point> &points) {
  std::vector<double> ret(points.size(), 0);

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t i = 0; i < points.size(); i++) {
    ret[i] = EvaluatePoint(path, points[i]);
  }

  num_runs += points.size();

  for (size_t i = 0; i < points.size(); i++) {
    if (ret[i] < best_val) {
      best_val = ret[i];
      best = points[i];
    }
  }

  fprintf(stderr, "%12lu evals, best: %8.4f at ", num_runs, best_val);
  for (size_t i = 0; i < best.size(); i++) {
    fprintf(stderr, " %12.8f", best[i]);
  }
  fprintf(stderr, "\n");

  return ret;
}
