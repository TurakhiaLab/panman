#pragma once
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace panmanUtils {
struct MafCoordinate {
    int64_t start, sourceSize;
    size_t width;
    std::string source, lowercase;
};
namespace mafCoordinateCodec {
inline const std::string& marker() {
    static const std::string value="panman:maf-coordinates:v2\n";return value;
}
inline const std::string& digits() {
    static const std::string value="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";return value;
}
inline void invalid(){throw std::runtime_error("Invalid compact MAF coordinate annotation");}
// Printable VLQ: five payload bits and one continuation bit per character.
// No binary bytes/NULs enter Cap'n Proto Text, and no schema changes are needed.
inline void put(std::string& out,uint64_t value) {
    do {unsigned digit=value&31;value>>=5;out+=digits()[digit|(value?32:0)];}while(value);
}
inline uint64_t get(const std::string& text,size_t& position) {
    uint64_t value=0;
    for(unsigned shift=0;shift<64;shift+=5) {
        if(position>=text.size())invalid();
        size_t digit=digits().find(text[position++]);
        if(digit==std::string::npos || (digit&31)>(UINT64_MAX>>shift))invalid();
        value|=uint64_t(digit&31)<<shift;
        if(!(digit&32))return value;
    }
    invalid();return 0;
}
class Encoder {
    std::string source,payload;
    int64_t total=0,lastStart=0;
    size_t lastId=0,lastWidth=0,count=0;
public:
    bool empty()const{return count==0;}
    bool full()const{return payload.size()>=65536;}
    bool accepts(const std::string& name,int64_t length)const{return empty() || (source==name && total==length);}
    void add(size_t id,int64_t start,int64_t length,size_t width,const std::string& name,const std::string& lowercase) {
        if(!width || !accepts(name,length) || start<lastStart || start>length ||
           (count && id<=lastId) || id-lastId>(UINT64_MAX>>2))invalid();
        if(empty()){source=name;total=length;}
        bool changed=width!=lastWidth,masked=lowercase!=".";
        put(payload,(uint64_t(id-lastId)<<2)|(changed?1:0)|(masked?2:0));
        put(payload,uint64_t(start-lastStart));
        if(changed)put(payload,width);
        if(masked) {
            size_t runs=0;for(char c:lowercase)runs+=c==',';
            if(!runs)invalid();put(payload,runs);
            std::istringstream input(lowercase);size_t first,end,previous=0;char colon,comma;
            for(size_t i=0;i<runs;++i) {
                if(!(input>>first>>colon>>end>>comma) || colon!=':' || comma!=',' ||
                   first<previous || end<=first || end>width)invalid();
                put(payload,first-previous);put(payload,end-first);previous=end;
            }
            input>>std::ws;if(!input.eof())invalid();
        }
        lastId=id;lastStart=start;lastWidth=width;++count;
    }
    std::string finish() {
        if(empty())return {};
        std::string result=marker()+source+" "+std::to_string(total)+" "+std::to_string(count)+"\n"+payload;
        *this=Encoder{};return result;
    }
};
inline void decode(const std::string& annotation,std::map<size_t,MafCoordinate>& rows) {
    if(annotation.compare(0,marker().size(),marker())!=0)invalid();
    size_t headerEnd=annotation.find('\n',marker().size());if(headerEnd==std::string::npos)invalid();
    std::istringstream header(annotation.substr(marker().size(),headerEnd-marker().size()));
    std::string source;int64_t total;uint64_t count;
    if(!(header>>source>>total>>count) || total<0 || !count)invalid();
    header>>std::ws;if(!header.eof())invalid();
    size_t position=headerEnd+1;uint64_t id=0,start=0,width=0;
    if(count>(annotation.size()-position)/2)invalid();
    for(uint64_t n=0;n<count;++n) {
        uint64_t flags=get(annotation,position),delta=flags>>2;
        if((n && !delta) || delta>SIZE_MAX-id)invalid();id+=delta;
        delta=get(annotation,position);if(delta>uint64_t(total)-start)invalid();start+=delta;
        if(flags&1)width=get(annotation,position);
        if(!width || width>SIZE_MAX)invalid();
        MafCoordinate row{int64_t(start),total,size_t(width),source,"."};
        if(flags&2) {
            uint64_t runs=get(annotation,position),previous=0;
            if(!runs || runs>(annotation.size()-position)/2)invalid();row.lowercase.clear();
            for(uint64_t i=0;i<runs;++i) {
                uint64_t offset=get(annotation,position),length=get(annotation,position);
                if(offset>width-previous || !length || length>width-previous-offset)invalid();
                uint64_t first=previous+offset;previous=first+length;
                row.lowercase+=std::to_string(first)+":"+std::to_string(previous)+",";
            }
        }
        if(!rows.emplace(size_t(id),std::move(row)).second)invalid();
    }
    if(position!=annotation.size())invalid();
}
}
}
