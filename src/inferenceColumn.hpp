#pragma once
#ifdef __CUDACC__
#define PANMAN_HD __host__ __device__
#else
#define PANMAN_HD
#endif
namespace panmanUtils {
// Unit-cost Sankoff, O(nodes * states). Missing observations (-1) are unknown,
// not gaps. Ties always choose the lowest state. Used unchanged on CPU/CUDA.
PANMAN_HD inline void inferColumn(int nodes, const int* offsets, const int* children,
        const int* parents, const int* observed, int hint, int states,
        int* costs, int* assigned) {
    constexpr int inf = 100000000;
    for(int node = 0; node < nodes; ++node) {
        int* row = costs + node * states;
        if(offsets[node] == offsets[node+1]) {
            for(int s = 0; s < states; ++s)
                row[s] = observed[node] < 0 || observed[node] == s ? 0 : inf;
        } else {
            for(int s = 0; s < states; ++s) row[s] = 0;
            for(int j = offsets[node]; j < offsets[node+1]; ++j) {
                const int* child = costs + children[j] * states;
                int minimum = inf;
                for(int s = 0; s < states; ++s) if(child[s] < minimum) minimum = child[s];
                for(int s = 0; s < states; ++s)
                    row[s] += child[s] < minimum + 1 ? child[s] : minimum + 1;
            }
        }
    }
    for(int node = nodes-1; node >= 0; --node) {
        int best = 0, score = inf;
        for(int s = 0; s < states; ++s) {
            int value = costs[node * states + s];
            if(parents[node] >= 0) value += s != assigned[parents[node]];
            if(value < score) { score = value; best = s; }
        }
        assigned[node] = parents[node] < 0 && hint >= 0 ? hint : best;
    }
}
}
#undef PANMAN_HD
