#include "constructionBackend.hpp"
#include "inferenceColumn.hpp"
#include <tbb/parallel_for.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>
namespace panmanUtils {
#ifdef PANMAN_HAS_CUDA
bool constructionCudaAvailable(std::string& reason);
std::vector<int> constructionCudaInfer(const InferenceTree&, const std::vector<int>&,
                                       const std::vector<int>&, int);
#endif
ConstructionOptions& constructionOptions() { static ConstructionOptions options; return options; }
ConstructionInference::ConstructionInference(InferenceTree input): tree(std::move(input)) {
    const auto& requested = constructionOptions().backend;
    if(requested != "auto" && requested != "cpu" && requested != "cuda")
        throw std::runtime_error("--backend must be auto, cpu, or cuda");
    if(tree.size() == 0) throw std::runtime_error("Empty inference tree");
    automatic = requested == "auto";
    std::string reason = "this build has no CUDA support";
#ifdef PANMAN_HAS_CUDA
    if(requested != "cpu") useCuda = constructionCudaAvailable(reason);
#endif
    if(requested == "cuda" && !useCuda) throw std::runtime_error("CUDA requested but unavailable: " + reason);
    std::cerr << "Construction backend: " << (useCuda ? "cuda" : "cpu");
    if(automatic && !useCuda) std::cerr << " (auto: " << reason << ")";
    std::cerr << '\n';
}
size_t ConstructionInference::capacity(int states) const {
    // Bounded workspace, independent of total alignment/genome size.
    return std::max(size_t(1), std::min(size_t(256), size_t(64*1024*1024) /
                                    (tree.size() * (states+2) * sizeof(int))));
}
std::vector<int> ConstructionInference::infer(const std::vector<int>& observed,
        const std::vector<int>& hints, int states) {
    if(states < 1 || states > 32 || observed.size() != hints.size() * tree.size())
        throw std::runtime_error("Invalid inference batch");
    if(hints.empty()) return {};
    if(useCuda) {
#ifdef PANMAN_HAS_CUDA
        try { return constructionCudaInfer(tree, observed, hints, states); }
        catch(const std::exception& error) {
            if(!automatic) throw;
            useCuda = false;
            std::cerr << "CUDA inference failed; auto falling back to CPU: " << error.what() << '\n';
        }
#endif
    }
    std::vector<int> output(observed.size());
    std::vector<int> costs(observed.size() * states);
    tbb::parallel_for(size_t(0), hints.size(), [&](size_t col) {
        inferColumn((int)tree.size(), tree.offsets.data(), tree.children.data(), tree.parents.data(),
                    observed.data()+col*tree.size(), hints[col], states,
                    costs.data()+col*tree.size()*states, output.data()+col*tree.size());
    });
    return output;
}
}
