// Original MAF positions live in existing tip annotations. Older files lack
// this table and continue using their historical concatenated coordinates.
std::map<size_t,panmanUtils::MafCoordinate> panmanUtils::Tree::mafCoordinates(const std::string& sample) const {
    std::map<size_t,MafCoordinate> result;
    auto node=allNodes.find(sample);if(node==allNodes.end())return result;
    const std::string marker="panman:maf-coordinates:v1\n";
    for(const auto& annotation:node->second->annotations) {
        if(annotation.compare(0,mafCoordinateCodec::marker().size(),mafCoordinateCodec::marker())==0) {
            mafCoordinateCodec::decode(annotation,result);continue;
        }
        if(annotation.compare(0,marker.size(),marker)!=0) {
            if(annotation.rfind("panman:maf-coordinates:",0)==0)
                throw std::runtime_error("Unsupported MAF coordinate annotation version");
            continue;
        }
        // Read the original absolute/text encoding unchanged for existing files.
        std::istringstream input(annotation.substr(marker.size()));
        size_t id;MafCoordinate row;
        while(input>>id) {
            if(!(input>>row.start>>row.sourceSize>>row.width>>row.source>>row.lowercase) || row.start<0 ||
               row.sourceSize<row.start || !result.emplace(id,row).second)
                throw std::runtime_error("Invalid original MAF coordinates for "+sample);
        }
        if(!input.eof())throw std::runtime_error("Invalid original MAF coordinate annotation");
    }
    return result;
}

// Shared additive metadata and chromosome resolver. Legacy files retain their
// existing reconstruction path when sequenceChromosomeRanges is absent.
void panmanUtils::Tree::validateChromosomeRanges() const {
    if(chromosomeRanges.empty()) {
        if(!homologyGroups.empty()) throw std::runtime_error("Homology groups require chromosome ranges");
        return;
    }
    std::set<int64_t> ids;
    for(const auto& b:blocks) if(b.secondaryBlockId==-1) ids.insert(b.primaryBlockId);
    if(ids.empty() || *ids.begin()!=0 || (uint64_t)*ids.rbegin()+1!=ids.size())
        throw std::runtime_error("Chromosome ranges require sequential primary block IDs");
    std::map<std::string,std::vector<std::pair<uint64_t,uint64_t>>> bySample;
    std::set<std::pair<std::string,std::string>> names;
    for(const auto& path:chromosomeRanges) {
        if(!allNodes.count(path.sequenceId) || path.chromosomeName.empty() ||
           !names.insert({path.sequenceId,path.chromosomeName}).second)
            throw std::runtime_error("Invalid or duplicate sequence/chromosome range owner");
        uint64_t last=0;
        for(const auto& range:path.ranges) {
            if(range.first>=range.second || range.first<last || range.second>ids.size())
                throw std::runtime_error("Invalid chromosome block range");
            bySample[path.sequenceId].push_back(range);last=range.second;
        }
    }
    for(auto& item:bySample) {
        std::sort(item.second.begin(),item.second.end());
        for(size_t i=1;i<item.second.size();++i)
            if(item.second[i].first<item.second[i-1].second)
                throw std::runtime_error("Chromosome ranges overlap within a genome");
    }
    std::set<int64_t> grouped;
    for(const auto& group:homologyGroups) {
        if(group.size()<2) throw std::runtime_error("Homology group needs at least two block instances");
        for(int64_t id:group)
            if(id<0 || (id&0xFFFFFFFFLL)!=0 || !ids.count(id>>32) || !grouped.insert(id).second)
                throw std::runtime_error("Invalid or duplicate homology group member");
    }
}
void panmanUtils::Tree::readChromosomeMetadata(const panman::Tree::Reader& input) {
    chromosomeRanges.clear();homologyGroups.clear();
    if(input.hasSequenceChromosomeRanges()) for(auto source:input.getSequenceChromosomeRanges()) {
        ChromosomeRanges path;path.sequenceId=source.getSequenceId().cStr();path.chromosomeName=source.getChromosomeName().cStr();
        for(auto range:source.getRanges()) path.ranges.push_back({range.getStart(),range.getEnd()});
        chromosomeRanges.push_back(std::move(path));
    }
    if(input.hasHomologyGroups()) for(auto source:input.getHomologyGroups()) {
        std::vector<int64_t> group;for(auto id:source.getBlockIds())group.push_back(id);
        homologyGroups.push_back(std::move(group));
    }
    validateChromosomeRanges();
}
void panmanUtils::Tree::writeChromosomeMetadata(panman::Tree::Builder output,
        const std::set<std::string>& owners) const {
    auto catalog=output.initChrLists(chrList.empty()?1:chrList.size());
    if(chrList.empty()) {
        catalog[0].setChrIdx(0);catalog[0].setChrName("chr0");
        auto ids=catalog[0].initBlockIds(blocks.size());
        for(size_t i=0;i<blocks.size();++i)ids.set(i,(int64_t)blocks[i].primaryBlockId<<32);
    } else for(size_t i=0;i<chrList.size();++i) {
        catalog[i].setChrIdx(chrList[i].chrIdx);catalog[i].setChrName(chrList[i].chrName);
        auto ids=catalog[i].initBlockIds(chrList[i].blockIds.size());
        for(size_t j=0;j<chrList[i].blockIds.size();++j)ids.set(j,chrList[i].blockIds[j]);
    }
    std::vector<std::string> keys;
    for(const auto& item:sequenceBlockOrders) if(owners.count(item.first))keys.push_back(item.first);
    std::sort(keys.begin(),keys.end());
    auto legacy=output.initSequenceBlockOrders(keys.size());
    for(size_t i=0;i<keys.size();++i) {
        legacy[i].setSequenceId(keys[i]);const auto& order=sequenceBlockOrders.at(keys[i]);
        auto ids=legacy[i].initOrderedBlockIds(order.size());
        for(size_t j=0;j<order.size();++j)ids.set(j,order[j]);
    }
    std::vector<const ChromosomeRanges*> retained;
    for(const auto& path:chromosomeRanges)if(owners.count(path.sequenceId))retained.push_back(&path);
    auto paths=output.initSequenceChromosomeRanges(retained.size());
    for(size_t i=0;i<retained.size();++i) {
        paths[i].setSequenceId(retained[i]->sequenceId);paths[i].setChromosomeName(retained[i]->chromosomeName);
        auto ranges=paths[i].initRanges(retained[i]->ranges.size());
        for(size_t j=0;j<retained[i]->ranges.size();++j) {
            ranges[j].setStart(retained[i]->ranges[j].first);ranges[j].setEnd(retained[i]->ranges[j].second);
        }
    }
    auto groups=output.initHomologyGroups(homologyGroups.size());
    for(size_t i=0;i<homologyGroups.size();++i) {
        auto ids=groups[i].initBlockIds(homologyGroups[i].size());
        for(size_t j=0;j<homologyGroups[i].size();++j)ids.set(j,homologyGroups[i][j]);
    }
}
std::vector<size_t> panmanUtils::Tree::chromosomeBlockIds(const std::string& sample,
        const std::string& chromosome) const {
    std::vector<size_t> ids;
    for(const auto& path:chromosomeRanges) if(path.sequenceId==sample &&
        (chromosome.empty() || path.chromosomeName==chromosome)) {
        for(const auto& range:path.ranges)
            for(uint64_t id=range.first;id<range.second;++id)ids.push_back(id);
    }
    return ids;
}

void panmanUtils::Tree::printRangeMAF(std::ostream& output) {
    struct Row {std::string name,text;int64_t start,size,total;bool strand;bool original=false;};
    std::vector<size_t> leader(blocks.size());std::iota(leader.begin(),leader.end(),0);
    for(const auto& group:homologyGroups)for(int64_t id:group)leader[id>>32]=group.front()>>32;
    std::map<size_t,std::vector<Row>> records;
    std::set<std::string> samples;
    for(const auto& item:allNodes)if(item.second->children.empty())samples.insert(item.first);
    for(const auto& sample:samples) {
        const auto coordinates=mafCoordinates(sample);
        sequence_t sequence;blockExists_t present;blockStrand_t strands;
        getSequenceFromReference(sequence,present,strands,sample,true);
        for(const auto& path:chromosomeRanges) if(path.sequenceId==sample) {
            std::vector<Row> rows;int64_t total=0;
            for(const auto& range:path.ranges) for(size_t id=range.first;id<range.second;++id) {
                if(id>=present.size() || !present[id].first)continue;
                std::string text;
                for(const auto& position:sequence[id].first) {
                    for(char c:position.second)text+=c=='x'?'-':c;
                    text+=position.first=='x'?'-':position.first;
                }
                int64_t size=std::count_if(text.begin(),text.end(),[](char c){return c!='-';});
                rows.push_back({sample+(chrList.size()==1 && chrList[0].chrName=="chr0"?std::string():"."+path.chromosomeName),std::move(text),total,size,0,strands[id].first});
                total+=size;
            }
            size_t index=0;
            for(const auto& range:path.ranges)for(size_t id=range.first;id<range.second;++id) {
                if(id>=present.size() || !present[id].first)continue;
                auto row=std::move(rows[index++]);row.total=total;
                auto original=coordinates.find(id);
                if(original!=coordinates.end()) {
                    const auto& source=original->second;
                    if(row.text.size()<source.width || row.size>source.sourceSize-source.start)
                        throw std::runtime_error("Original MAF coordinates no longer match reconstructed sequence");
                    row.name=source.source;row.start=source.start;row.total=source.sourceSize;
                    // Remove only the reconstruction sentinel; real all-gap
                    // columns in the original alignment must survive a round-trip.
                    row.text.resize(source.width);row.original=true;
                    if(source.lowercase!=".") {
                        std::istringstream runs(source.lowercase);size_t first,last;char colon,comma;
                        while(runs>>first>>colon>>last>>comma) {
                            if(colon!=':' || comma!=',' || first>last || last>row.text.size())
                                throw std::runtime_error("Invalid MAF lowercase annotation");
                            for(size_t k=first;k<last;++k)row.text[k]=std::tolower(static_cast<unsigned char>(row.text[k]));
                        }
                        if(!runs.eof())throw std::runtime_error("Invalid MAF lowercase annotation");
                    }
                }
                if(!row.strand)row.start=row.total-row.start-row.size;
                records[leader[id]].push_back(std::move(row));
            }
        }
    }
    output<<"##maf version=1\n";
    for(auto& record:records) {
        size_t width=0;for(const auto& row:record.second)width=std::max(width,row.text.size());
        bool original=std::all_of(record.second.begin(),record.second.end(),[](const Row& row){return row.original;});
        std::vector<bool> keep(width,original);
        for(const auto& row:record.second)for(size_t i=0;i<row.text.size();++i)if(row.text[i]!='-')keep[i]=true;
        if(std::find(keep.begin(),keep.end(),true)==keep.end())continue;
        output<<"\na\n";
        for(const auto& row:record.second) {
            output<<"s\t"<<row.name<<'\t'<<row.start<<'\t'<<row.size<<'\t'<<(row.strand?'+':'-')<<'\t'<<row.total<<'\t';
            for(size_t i=0;i<width;++i)if(keep[i])output<<(i<row.text.size()?row.text[i]:'-');
            output<<'\n';
        }
    }
    output<<'\n';
}
