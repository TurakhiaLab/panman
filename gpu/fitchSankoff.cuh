
#include "panmanUtils.cuh"
#ifndef UTILS
#include "utils.hpp"
#endif
#include <chrono>
#define DEBUG_GPU 0

#define INF 2000000000

void fitch_sankoff_on_gpu(panmanUtils::Tree* T, std::unordered_map<std::string, std::string>& seqs, utility::util* u);