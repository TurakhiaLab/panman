// Included after parseMafGappedConsensus; construction-specific implementation.
#include "mafLayout.hpp"
#include <tbb/parallel_for.h>
#include "constructionBackend.hpp"
#include "constructionProgress.hpp"

void panmanUtils::Tree::buildFromMAF(std::ifstream& input, const std::string& reference) {
    MAFBlockGraph graph;
    readMAF(input, graph);
    if(graph.consensusByGroup.empty()) throw std::runtime_error("MAF contained no sequences.");
    if(!reference.empty() && (allNodes.find(reference)==allNodes.end() || !allNodes.at(reference)->children.empty()))
        throw std::runtime_error("MAF reference must exactly match a Newick tip.");
    std::vector<std::string> samples=graph.seqOrder;
    if(!reference.empty()) {
        auto it=std::find(samples.begin(),samples.end(),reference);
        if(it!=samples.end()) std::rotate(samples.begin(),it,it+1);
    }
    ConstructionProgress progress;
    progress.update("Preparing chromosome paths",0,samples.size(),"genomes",true);
    std::vector<mafLayout::Path> paths;
    struct SourcePath {std::string sample,chromosome; std::vector<const MAFBlockGraph::Occ*> rows;};
    std::vector<SourcePath> sources;
    for(size_t s=0;s<samples.size();++s) {
        std::map<std::string,std::vector<const MAFBlockGraph::Occ*>> byChromosome;
        size_t grouped=0;
        for(const auto& occurrence:graph.seqOccs[samples[s]]) {
            byChromosome[occurrence.chromosome].push_back(&occurrence);
            if((++grouped & 65535)==0)progress.update("Preparing chromosome paths",s,samples.size(),
                "genome="+samples[s]+"; grouped blocks="+std::to_string(grouped));
        }
        for(const auto& chromosome:graph.chromosomeOrder) {
            auto found=byChromosome.find(chromosome);if(found==byChromosome.end()) continue;
            auto rows=std::move(found->second);
            std::stable_sort(rows.begin(),rows.end(),[](const auto* a,const auto* b){return a->forwardStart<b->forwardStart;});
            mafLayout::Path path;path.sample=s;
            for(const auto* row:rows) path.groups.push_back(row->group);
            paths.push_back(std::move(path));sources.push_back({samples[s],chromosome,std::move(rows)});
        }
        progress.update("Preparing chromosome paths",s+1,samples.size(),"genome="+samples[s]+"; paths="+
            std::to_string(paths.size()),s+1==samples.size());
    }
    const auto& options=constructionOptions();
    auto layout=mafLayout::construct(paths,options.matchBudget,options.refinementPasses,options.refinementBudget,
        [&](const std::string& phase,size_t done,size_t total,const std::string& detail,bool force) {
            progress.update(phase,done,total,detail,force);
        }, [](size_t count,const std::function<void(size_t)>& job) {
            tbb::parallel_for(size_t(0),count,[&](size_t i){job(i);});
        });
    std::cerr<<"MAF layout: "<<layout.stats.initial<<" initial instances, "<<layout.stats.merged
             <<" merged, "<<layout.instances.size()<<" final; "<<layout.stats.omittedMatches
             <<" repeat groups omitted from anchors; index builds "<<layout.stats.indexBuilds
             <<", reuses "<<layout.stats.indexReuses<<"; batches "<<layout.stats.batches
             <<", largest "<<layout.stats.largestBatch<<"; refinement work "<<layout.stats.refinementWork<<'\n';
    // Coordinate deltas and per-chromosome headers avoid repeating absolute
    // positions/source names. Independent chunks bound annotation size and keep
    // overlaps lossless; all data still uses the existing annotation field.
    for(size_t p=0;p<paths.size();++p) {
        mafCoordinateCodec::Encoder encoder;
        auto& annotations=allNodes.at(sources[p].sample)->annotations;
        for(size_t j=0;j<paths[p].instances.size();++j) {
            const auto& row=*sources[p].rows[j];
            if(!encoder.accepts(row.source,row.sourceSize))annotations.push_back(encoder.finish());
            encoder.add(paths[p].instances[j],row.forwardStart,row.sourceSize,row.row.size(),row.source,row.lowercase);
            if(encoder.full())annotations.push_back(encoder.finish());
        }
        if(!encoder.empty())annotations.push_back(encoder.finish());
    }
    if(layout.instances.size()>size_t(std::numeric_limits<int32_t>::max()))
        throw std::runtime_error("MAF layout exceeds the 32-bit primary block ID limit.");
    // Resolve every path before inference. Mutations never need renumbering.
    progress.update("Assigning chromosome ranges",0,paths.size(),"paths",true);
    std::map<int,std::vector<int64_t>> homology;
    std::map<std::string,size_t> catalog;
    for(const auto& name:graph.chromosomeOrder) {
        catalog[name]=chrList.size();std::vector<int64_t> ids;
        chrList.emplace_back(chrList.size(),ids);chrList.back().chrName=name;
    }
    std::vector<bool> legacyAssigned(layout.instances.size(),false);
    std::map<size_t,std::map<size_t,size_t>> ownership;
    for(size_t p=0;p<paths.size();++p) {
        ChromosomeRanges ranges;ranges.sequenceId=sources[p].sample;ranges.chromosomeName=sources[p].chromosome;
        chromosomeRanges.push_back(std::move(ranges));
        for(size_t id:paths[p].instances) {
            if(!ownership[paths[p].sample].emplace(id,p).second)
                throw std::runtime_error("A block instance occurs twice in one genome");
            sequenceBlockOrders[sources[p].sample].push_back((int64_t)id<<32);
            if(!legacyAssigned[id]) {
                chrList[catalog.at(sources[p].chromosome)].blockIds.push_back((int64_t)id<<32);
                legacyAssigned[id]=true;
            }
        }
        progress.update("Assigning chromosome ranges",p+1,paths.size(),"genome="+sources[p].sample,p+1==paths.size());
    }
    for(const auto& sample:ownership) {
        size_t previous=size_t(-1);
        for(const auto& occurrence:sample.second) {
            auto& ranges=chromosomeRanges[occurrence.second].ranges;
            if(previous==occurrence.second)ranges.back().second=occurrence.first+1;
            else ranges.push_back({occurrence.first,occurrence.first+1});
            previous=occurrence.second;
        }
    }
    for(auto& chr:chrList) {
        std::sort(chr.blockIds.begin(),chr.blockIds.end());
        chr.blockIds.erase(std::unique(chr.blockIds.begin(),chr.blockIds.end()),chr.blockIds.end());
    }
    for(size_t id=0;id<layout.instances.size();++id) homology[layout.instances[id].group].push_back((int64_t)id<<32);
    for(auto& group:homology) if(group.second.size()>1) homologyGroups.push_back(std::move(group.second));
    for(const auto& sample:samples) {sequenceInverted[sample]=false;rotationIndexes[sample]=0;}

    progress.update("Chromosome ranges ready",chromosomeRanges.size(),chromosomeRanges.size(),
        "homology copy groups="+std::to_string(homologyGroups.size()),true);
    // One deterministic flat tree is shared by CPU and CUDA inference.
    std::vector<Node*> nodes;
    std::unordered_map<std::string,int> nodeId;
    std::function<void(Node*)> postorder=[&](Node* n){for(auto* c:n->children)postorder(c);nodeId[n->identifier]=nodes.size();nodes.push_back(n);};
    postorder(root);
    InferenceTree flat;flat.parents.assign(nodes.size(),-1);flat.offsets.push_back(0);
    for(size_t n=0;n<nodes.size();++n) {
        for(auto* c:nodes[n]->children) {int id=nodeId.at(c->identifier);flat.children.push_back(id);flat.parents[id]=n;}
        flat.offsets.push_back(flat.children.size());
    }
    const auto parent=flat.parents;
    ConstructionInference inference(std::move(flat));
    struct Site {int block,position,gap,baseline;};
    auto appendMutation=[&](Node* node,const Site& site,int state,int ancestor) {
        if(state==ancestor) return;
        int type=ancestor==0?NI:(state==0?ND:NS);
        auto& mutations=node->nucMutation;
        const int capacity=mutationPayloadCapacity(alphabet),bits=getCodec(alphabet).bitsPerSymbol;
        if(!mutations.empty()) {
            auto& last=mutations.back();int length=last.mutInfo>>4;
            bool adjacent=site.gap<0 ? last.nucGapPosition<0 && site.position==last.nucPosition+length
                : last.nucGapPosition>=0 && site.position==last.nucPosition && site.gap==last.nucGapPosition+length;
            if(last.primaryBlockId==site.block && (last.mutInfo&15)==type && length<capacity && adjacent) {
                last.nucs|=(uint32_t)state<<(bits*(capacity-1-length));last.mutInfo+=16;return;
            }
        }
        mutations.emplace_back(std::make_tuple(site.block,-1,site.position,site.gap,type,state));
    };
    // Batch presence/orientation states without constructing a dense genomes x
    // blocks matrix. State 0=absent, 1=forward, 2=reverse.
    const size_t N=nodes.size();
    progress.update("Block presence/orientation inference",0,layout.instances.size(),"block instances",true);
    for(size_t first=0;first<layout.instances.size();) {
        size_t count=std::min(inference.capacity(3),layout.instances.size()-first);
        std::vector<int> observed(count*N,-1),hints(count,-1);
        for(size_t j=0;j<count;++j) {
            for(size_t n=0;n<N;++n) if(nodes[n]->children.empty()) observed[j*N+n]=0;
            for(const auto& row:layout.instances[first+j].rows) {
                const auto& source=sources[row.second.first];const auto* occurrence=source.rows[row.second.second];
                observed[j*N+nodeId.at(source.sample)]=occurrence->strand?1:2;
            }
            if(!reference.empty()) hints[j]=observed[j*N+nodeId.at(reference)];
        }
        auto states=inference.infer(observed,hints,3);
        for(size_t j=0;j<count;++j) for(size_t n=0;n<N;++n) {
            int state=states[j*N+n],ancestor=parent[n]<0?0:states[j*N+parent[n]];
            if(state==ancestor)continue;
            auto mutation=ancestor==0?std::make_pair(BI,state==2):(state==0?std::make_pair(BD,false):std::make_pair(BD,true));
            nodes[n]->blockMutation.emplace_back(first+j,mutation);
        }
        first+=count;
        progress.update("Block presence/orientation inference",first,layout.instances.size(),"block instances",first==layout.instances.size());
    }
    const int stateCount=getCodec(alphabet).stateCount;
    std::vector<int> observed,hints;std::vector<Site> sites;
    const size_t capacity=inference.capacity(stateCount);
    size_t completedColumns=0, totalColumns=0, completedBlocks=0;
    for(const auto& instance:layout.instances)totalColumns+=graph.consensusByGroup[instance.group].size();
    progress.update("Nucleotide/protein inference",0,totalColumns,"columns",true);
    auto flush=[&](){
        if(sites.empty())return;
        auto states=inference.infer(observed,hints,stateCount);
        for(size_t j=0;j<sites.size();++j) for(size_t n=0;n<N;++n)
            appendMutation(nodes[n],sites[j],states[j*N+n],parent[n]<0?sites[j].baseline:states[j*N+parent[n]]);
        completedColumns+=sites.size();
        progress.update("Nucleotide/protein inference",completedColumns,totalColumns,"columns; blocks="+
            std::to_string(completedBlocks)+"/"+std::to_string(layout.instances.size()),completedColumns==totalColumns);
        observed.clear();hints.clear();sites.clear();
    };
    for(size_t id=0;id<layout.instances.size();++id) {
        const auto& instance=layout.instances[id];
        const std::string& consensus=graph.consensusByGroup[instance.group];
        auto parsed=parseMafGappedConsensus(consensus);
        blocks.emplace_back(id,parsed.gapFree,alphabet,parsed.gapFree.size());
        GapList gap;gap.primaryBlockId=id;gap.secondaryBlockId=-1;
        for(const auto& item:parsed.gaps) {gap.nucPosition.push_back(item.first);gap.nucGapLength.push_back(item.second);}
        gaps.push_back(std::move(gap));
        for(size_t col=0;col<consensus.size();++col) {
            size_t offset=observed.size();observed.resize(offset+N,-1);
            for(const auto& row:instance.rows) {
                const auto& source=sources[row.second.first];const auto* occurrence=source.rows[row.second.second];
                observed[offset+nodeId.at(source.sample)]=getCodeFromSymbol(occurrence->row[col],alphabet);
            }
            hints.push_back(reference.empty()?-1:observed[offset+nodeId.at(reference)]);
            sites.push_back({(int)id,parsed.oldToNew[col].first,parsed.oldToNew[col].second,getCodeFromSymbol(consensus[col],alphabet)});
            if(sites.size()==capacity)flush();
        }
        ++completedBlocks;
    }
    flush();
    progress.update("Validating chromosome metadata",0,chromosomeRanges.size(),"range records",true);
    validateChromosomeRanges();
    progress.update("Construction complete",layout.instances.size(),layout.instances.size(),
        "block instances; chromosome paths="+std::to_string(chromosomeRanges.size()),true);
}
