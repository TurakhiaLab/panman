#include "../src/constructionBackend.hpp"
#include "../src/inferenceColumn.hpp"
#include <cuda_runtime.h>
#include <stdexcept>
namespace panmanUtils {
namespace {
void checked(cudaError_t code) {
    if(code != cudaSuccess) throw std::runtime_error(cudaGetErrorString(code));
}
struct DeviceInts {
    int* data = nullptr;
    explicit DeviceInts(size_t count) { if(count) checked(cudaMalloc((void**)&data, count*sizeof(int))); }
    ~DeviceInts() { if(data) cudaFree(data); }
    DeviceInts(const DeviceInts&) = delete;
    void upload(const std::vector<int>& values) {
        if(!values.empty()) checked(cudaMemcpy(data, values.data(), values.size()*sizeof(int), cudaMemcpyHostToDevice));
    }
};
__global__ void inferenceKernel(int nodes, int columns, const int* offsets, const int* children,
        const int* parents, const int* observed, const int* hints, int states,
        int* costs, int* assigned) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    if(col < columns) inferColumn(nodes, offsets, children, parents,
        observed+(size_t)col*nodes, hints[col], states,
        costs+(size_t)col*nodes*states, assigned+(size_t)col*nodes);
}
}
bool constructionCudaAvailable(std::string& reason) {
    int count = 0;
    cudaError_t status = cudaGetDeviceCount(&count);
    if(status != cudaSuccess) { reason = cudaGetErrorString(status); return false; }
    if(count == 0) { reason = "no CUDA device found"; return false; }
    return true;
}
std::vector<int> constructionCudaInfer(const InferenceTree& tree,
        const std::vector<int>& observed, const std::vector<int>& hints, int states) {
    DeviceInts offsets(tree.offsets.size()), children(tree.children.size()), parents(tree.parents.size());
    DeviceInts input(observed.size()), rootHints(hints.size()), costs(observed.size()*states), output(observed.size());
    offsets.upload(tree.offsets); children.upload(tree.children); parents.upload(tree.parents);
    input.upload(observed); rootHints.upload(hints);
    inferenceKernel<<<(hints.size()+63)/64,64>>>((int)tree.size(), (int)hints.size(), offsets.data,
        children.data, parents.data, input.data, rootHints.data, states, costs.data, output.data);
    checked(cudaGetLastError()); checked(cudaDeviceSynchronize());
    std::vector<int> result(observed.size());
    if(!result.empty()) checked(cudaMemcpy(result.data(), output.data, result.size()*sizeof(int), cudaMemcpyDeviceToHost));
    return result;
}
}
