#include <fstream>
#include <unordered_map>
#include <iostream>
// #include "file_io.hpp"

#define KMER_SIZE  17

using namespace std;

int main(int argc, char* argv[])
{
    ifstream f(argv[1]);

    string name;
    string seq;
    string d;
    getline(f, name);
    while (true)
    {
        getline(f, d);
        if (f.eof()) break;
        seq.append(d);
    }

    unordered_map<string, int> splitters;
    unordered_map<string, int> splitters_pos;
    for (auto i = 0; i < seq.size(); i+=KMER_SIZE)
    {
        string splitter = seq.substr(i, KMER_SIZE);
        if (splitters.find(splitter) == splitters.end())
        {
            splitters[splitter] = 1;
            splitters_pos[splitter] = i;
        }
        else
        {   
            splitters[splitter]++;
        }
    }

    // Printing
    cout << splitters.size() << std::endl;

    // filtering 
    int c = 0;
    for (auto e: splitters)
    {
        if (e.second == 1)
        {
            c++;
        }
    }
    cout << "total splitters: " << c << endl;

    // // Printing
    // cout << splitters.size() << std::endl;

}



// int main(int argc, char *argv[])
// {
//     ifstream f(argv[1]);
//     read_t r (f);
    
//     for (auto i = 0; i < r.id.size(); i++)
//     {
//         std::cout << r.id[i] << std::endl;
//         for (auto j = 0; j< r.sequence[i].size(); j++)
//         {
//             pair<char, char> p = panmat::int82char(r.sequence[i][j]);
//             // std::cout << "int8_t: " << (r.sequence[i][j]&0xFF) << std::endl; 
//             std::cout << p.first <<p.second;
//             // break;
//         }
//         std::cout << std::endl;
//     }

//     return 0;
// }