#include <fstream>
#include <vector>

using namespace std;

namespace panmat {

    struct compression_param_t
    {
        uint32_t kmer_length;
    };

    class sample
    {
        private:
            vector<vector<int8_t>> splitters;    
        public:
            void find_splitter(vector<int8_t>, compression_param_t);
    };
        
    
    
    enum FILE_TYPE {
        PANMAT = 0,
        GFA = 1,
        PANGRAPH=2,
        MSA = 3,
        MSA_OPTIMIZE = 4,
        FASTA = 5,
        NWK = 6
    };

    enum BlockMutationType {
        BI = 1,
        BD = 0,
        BIn = 2,
        NONE = 1000
    };

    uint8_t nuc2int8(char nuc)
    {
        switch(nuc)
        {
            case 'A':
                return 1;
            case 'C':
                return 2;
            case 'G':
                return 4;
            case 'T':
                return 8;
            case 'R':
                return 5;
            case 'Y':
                return 10;
            case 'S':
                return 6;
            case 'W':
                return 9;
            case 'K':
                return 12;
            case 'M':
                return 3;
            case 'B':
                return 14;
            case 'D':
                return 13;
            case 'H':
                return 11;
            case 'V':
                return 7;
            case 'N':
                return 15;
            default:
                return 0;
        }
    }

    char int2char(uint8_t v)
    {
        switch (v)
        {
            case 1:
                return 'A';
            case 2:
                return 'C';
            case 3:
                return 'M';
            case 4:
                return 'G';
            case 5:
                return 'R';
            case 6:
                return 'S';
            case 7:
                return 'V';
            case 8:
                return 'T';
            case 9:
                return 'W';
            case 10:
                return 'Y';
            case 11:
                return 'H';
            case 12:
                return 'K';
            case 13:
                return 'D';
            case 14:
                return 'B';
            case 15:
                return 'N';
            default:
                return '-';
        }
    }

    pair<char, char> int82char(uint8_t v)
    {
        pair<char, char> chars;

        chars.first = int2char(v >> 4);
        chars.second = int2char(v & 0x0F);
        // std::cout << ((v)&0xFF) << chars.first <<chars.second << " ";
        return chars;
    }

}
