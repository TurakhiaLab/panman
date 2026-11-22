#pragma once
#include <cstddef>
#include <string>
#include <vector>
namespace panmanUtils {
struct ConstructionOptions {
    std::string backend = "auto";
    size_t matchBudget = 1000000;
    size_t refinementBudget = 2000000;
    unsigned refinementPasses = 2;
};
ConstructionOptions& constructionOptions();
struct InferenceTree {
    // Postorder indices: every child precedes its parent; root is last.
    std::vector<int> offsets, children, parents;
    size_t size() const { return parents.size(); }
};
class ConstructionInference {
    InferenceTree tree;
    bool useCuda = false, automatic = true;
public:
    explicit ConstructionInference(InferenceTree input);
    size_t capacity(int states) const;
    std::vector<int> infer(const std::vector<int>& observed,
                           const std::vector<int>& rootHints, int states);
};
}
