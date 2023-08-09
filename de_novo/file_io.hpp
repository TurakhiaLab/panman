#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "panmat.hpp"

using namespace std;

class read_t
{
    public:
        vector<string> id;
        vector<vector<int8_t>> sequence;
        read_t(ifstream& f);    
    
};

read_t::read_t(ifstream& f)
    {   
        string d;
        vector<int8_t> seq;
        string s;
        while (true)
        {
            getline(f, d);
            if (f.eof()) break;
            if (d[0] == '>')
            {
                if (id.size() > 0)
                {
                    // Move pointer back
                    // f.seekg(-d.size(), ios::cur);
                    // break;
                
                    for (auto i = 0; i < s.size(); i+=2)
                    {
                        if (s[i] == '\n') continue;
                        int8_t c = panmat::nuc2int8(s[i]) << 4;
                        if (i <= s.size() - 2)
                            c |= panmat::nuc2int8(s[i+1]);
                        seq.push_back(c);
                    }
                    sequence.push_back(seq); 
                    seq.clear();
                }
                
                id.push_back(d.substr(1));
            }
            else 
                s.append(d);
                
            d.clear();
            
        } 
        for (auto i = 0; i < s.size(); i+=2)
        {
            if (s[i] == '\n') continue;
            int8_t c = panmat::nuc2int8(s[i]) << 4;
            if (i <= s.size() - 2)
                c |= panmat::nuc2int8(s[i+1]);
            seq.push_back(c);
        }
        sequence.push_back(seq); 
        seq.clear();
        
        
    };

