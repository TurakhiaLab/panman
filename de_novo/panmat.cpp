#include "panmat.hpp"
#include <unordered_map>
#include <bits/stdc++.h>    

void panmat::sample::find_splitter(vector<int8_t> s, compression_param_t param)
{
    map<vector<int8_t>, int8_t> splitter_count;

    uint32_t kmer_length = param.kmer_length;
    uint32_t s_size = s.size();

    for (auto i = 0; i < s_size; i += kmer_length/2)
    {
        vector<int8_t> s_curr;
        for (auto j = i; j < i + kmer_length/2; j++)
        {
            s_curr.push_back(s[j]);
        }
        splitter_count[s_curr]++;
    } 

    for (auto s_curr: splitter_count)
    {
        if (s_curr.second == 1)
        {
            splitters.push_back(s_curr.first);
        }
    }
}
