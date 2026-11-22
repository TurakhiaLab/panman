#pragma once
#include <algorithm>
#include <cstddef>
#include <map>
#include <functional>
#include <string>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
namespace panmanUtils { namespace mafLayout {
struct Path {
    size_t sample;
    std::vector<int> groups;
    std::vector<size_t> instances;
};
struct Instance {
    int group;
    // sample -> (path index, occurrence index); a sample cannot use a node twice.
    std::map<size_t, std::pair<size_t,size_t>> rows;
};
struct Statistics { size_t initial = 0, merged = 0, omittedMatches = 0, refinementWork = 0, indexBuilds = 0, indexReuses = 0, batches = 0, largestBatch = 0; };
struct Layout { std::vector<Instance> instances; Statistics stats; };
struct Match { size_t root, incoming; };
using Index = std::unordered_map<int,std::vector<size_t>>;
using Executor = std::function<void(size_t,const std::function<void(size_t)>&)>;
using Progress = std::function<void(const std::string&,size_t,size_t,const std::string&,bool)>;
// Exact sparse LCS within a matching-pair budget. Frequent identities exceeding
// the budget are omitted from anchors; unmatched rows are retained as copies.
inline std::vector<Match> chain(const std::vector<size_t>& root,
        const std::vector<Instance>& nodes, const Path& path, size_t budget, size_t& omitted,
        const std::function<void(const std::string&)>& heartbeat = {}, const Index* shared = nullptr) {
    Index index;
    if(shared) {
        for(int group:path.groups) {
            if(index.count(group)) continue;
            auto found=shared->find(group);
            auto& positions=index[group];
            if(found!=shared->end()) for(size_t pos:found->second)
                if(!nodes[root[pos]].rows.count(path.sample)) positions.push_back(pos);
        }
    } else
    for(size_t i=0;i<root.size();++i) {
        if(heartbeat && i && (i & 65535)==0)heartbeat("indexing root blocks="+std::to_string(i)+"/"+std::to_string(root.size()));
        if(!nodes[root[i]].rows.count(path.sample)) index[nodes[root[i]].group].push_back(i);
    }
    std::map<int,size_t> frequency;
    for(int group:path.groups) ++frequency[group];
    std::set<int> allowed;
    for(const auto& f:frequency) {
        auto found=index.find(f.first);
        if(found==index.end()) continue;
        if(found->second.size() <= budget / f.second) {
            budget-=found->second.size()*f.second; allowed.insert(f.first);
        } else ++omitted;
    }
    struct Link { size_t root, incoming, previous; };
    const size_t none=size_t(-1);
    std::vector<Link> links;
    std::vector<size_t> tails, tailLinks;
    for(size_t j=0;j<path.groups.size();++j) {
        if(heartbeat && j && (j & 65535)==0)heartbeat("chaining incoming blocks="+std::to_string(j)+"/"+std::to_string(path.groups.size()));
        int group=path.groups[j];
        if(!allowed.count(group)) continue;
        const auto& positions=index.at(group);
        for(auto it=positions.rbegin();it!=positions.rend();++it) {
            size_t k=std::lower_bound(tails.begin(),tails.end(),*it)-tails.begin();
            if(k<tails.size() && tails[k]==*it) continue;
            links.push_back({*it,j,k ? tailLinks[k-1] : none});
            if(k==tails.size()) { tails.push_back(*it); tailLinks.push_back(links.size()-1); }
            else { tails[k]=*it; tailLinks[k]=links.size()-1; }
        }
    }
    std::vector<Match> matches;
    if(!tailLinks.empty()) for(size_t p=tailLinks.back();p!=none;p=links[p].previous)
        matches.push_back({links[p].root,links[p].incoming});
    std::reverse(matches.begin(),matches.end());
    return matches;
}
inline Layout construct(std::vector<Path>& paths, size_t matchBudget, unsigned passes,
                        size_t refinementBudget, const Progress& progress = {}, const Executor& execute = {}) {
    Layout result;
    auto& nodes=result.instances;
    std::vector<size_t> root;
    if(progress)progress("Progressive layout",0,paths.size(),"chromosome paths",true);
    // The index contains positions only; ownership is checked afresh by chain.
    Index index;
    bool indexValid=false;
    const size_t scratchLimit=256*1024*1024;
    // Saturating arithmetic keeps huge/repetitive paths from overflowing estimates.
    auto add=[&](size_t a,size_t b){return a>=scratchLimit || b>=scratchLimit-a?scratchLimit:a+b;};
    auto multiply=[&](size_t a,size_t b){return b && a>scratchLimit/b?scratchLimit:a*b;};
    auto estimateScratch=[&](const Path& path) {
        std::map<int,size_t> frequency;
        for(int group:path.groups)++frequency[group];
        size_t pairs=0,positions=0;
        for(const auto& f:frequency) {
            auto found=index.find(f.first);
            if(found==index.end())continue;
            size_t count=found->second.size();
            positions=add(positions,count);
            // Ignoring ownership gives an upper bound even when filtering changes
            // which groups fit the chaining budget.
            pairs=std::min(matchBudget,add(pairs,multiply(count,f.second)));
        }
        return add(4096,add(multiply(path.groups.size(),128),
            add(multiply(positions,16),multiply(pairs,64))));
    };
    // Boundaries depend on input and estimated memory, never worker count.
    for(size_t begin=0;begin<paths.size();) {
        if(!indexValid) {
            index.clear();
            for(size_t i=0;i<root.size();++i) {
                index[nodes[root[i]].group].push_back(i);
                if(progress && i && (i & 65535)==0)progress("Progressive layout",begin,paths.size(),
                    "indexing shared root="+std::to_string(i)+"/"+std::to_string(root.size()),false);
            }
            indexValid=true;++result.stats.indexBuilds;
        } else ++result.stats.indexReuses;
        size_t end=begin,scratch=0;
        // Seed a single path. Short scaffold paths can subsequently share much
        // larger batches than long chromosomes without reserving matchBudget each.
        while(end<paths.size() && end-begin<(begin?64:1)) {
            size_t estimate=estimateScratch(paths[end]);
            if(end>begin && estimate>scratchLimit-scratch)break;
            scratch=add(scratch,estimate);++end;
        }
        ++result.stats.batches;
        result.stats.largestBatch=std::max(result.stats.largestBatch,end-begin);
        std::vector<std::vector<Match>> plans(end-begin);
        std::vector<size_t> omitted(end-begin,0);
        auto align=[&](size_t job) {
            plans[job]=chain(root,nodes,paths[begin+job],matchBudget,omitted[job],{},&index);
        };
        if(execute)execute(end-begin,align);
        else for(size_t job=0;job<end-begin;++job)align(job);
        // Only this deterministic commit mutates nodes. Recheck ownership because
        // two chromosomes of the same sample may have proposed the same anchor.
        std::map<size_t,std::vector<size_t>> insertions;
        for(size_t p=begin;p<end;++p) {
            auto& path=paths[p];path.instances.resize(path.groups.size());
            result.stats.omittedMatches+=omitted[p-begin];
            size_t incoming=0;
            const auto& matches=plans[p-begin];
            for(size_t m=0;m<=matches.size();++m) {
                if(m<matches.size() && nodes[root[matches[m].root]].rows.count(path.sample))continue;
                size_t gap=m<matches.size()?matches[m].root:root.size();
                size_t stop=m<matches.size()?matches[m].incoming:path.groups.size();
                while(incoming<stop) {
                    size_t id=nodes.size();
                    nodes.push_back({path.groups[incoming],{{path.sample,{p,incoming}}}});
                    insertions[gap].push_back(id);path.instances[incoming++]=id;
                }
                if(m<matches.size()) {
                    size_t id=root[gap];
                    nodes[id].rows.emplace(path.sample,std::make_pair(p,stop));
                    path.instances[stop]=id;incoming=stop+1;
                }
            }
        }
        if(!insertions.empty()) {
            std::vector<size_t> next;next.reserve(nodes.size());
            auto insertion=insertions.begin();
            for(size_t i=0;i<=root.size();++i) {
                if(insertion!=insertions.end() && insertion->first==i) {
                    next.insert(next.end(),insertion->second.begin(),insertion->second.end());++insertion;
                }
                if(i<root.size())next.push_back(root[i]);
            }
            root=std::move(next);indexValid=false;
        }
        if(progress)progress("Progressive layout",end,paths.size(),"root instances="+std::to_string(root.size())+
            "; batch paths="+std::to_string(end-begin)+"; index builds="+std::to_string(result.stats.indexBuilds)+
            "; index reuses="+std::to_string(result.stats.indexReuses),end==paths.size());
        begin=end;
    }
    result.stats.initial=nodes.size();
    // Required edges come only from tip/chromosome paths, not arbitrary root
    // adjacencies. Every commit checks the current graph, never a stale batch.
    if(progress)progress("Building ordering graph",0,paths.size(),"root instances="+std::to_string(nodes.size()),true);
    std::vector<std::set<size_t>> out(nodes.size()), in(nodes.size());
    for(const auto& path:paths) for(size_t i=1;i<path.instances.size();++i) {
        size_t a=path.instances[i-1],b=path.instances[i];
        out[a].insert(b); in[b].insert(a);
    }
    if(progress)progress("Building ordering graph",paths.size(),paths.size(),"complete",true);
    std::vector<size_t> representative(nodes.size());
    std::iota(representative.begin(),representative.end(),0);
    auto find=[&](size_t n) {
        size_t r=n; while(representative[r]!=r) r=representative[r];
        while(representative[n]!=n) { size_t next=representative[n]; representative[n]=r;n=next; }
        return r;
    };
    size_t work=0;
    unsigned activePass=0;
    if(progress)progress("Refinement",0,passes,"work budget="+std::to_string(refinementBudget),true);
    auto reportRefinement=[&](bool force=false) {
        if(progress)progress("Refinement",activePass,passes,"work="+std::to_string(work)+"/"+
            std::to_string(refinementBudget)+"; merged instances="+std::to_string(result.stats.merged),force);
    };
    std::vector<size_t> visited(nodes.size(),0); size_t stamp=0;
    auto reachable=[&](size_t a,size_t b) {
        ++stamp; std::vector<size_t> stack{a}; visited[a]=stamp;
        while(!stack.empty()) {
            if(work>=refinementBudget) return true; // unknown => conservatively reject
            size_t u=stack.back();stack.pop_back();++work;
            if((work & 65535)==0)reportRefinement();
            for(size_t v:out[u]) {
                if(work>=refinementBudget) return true;
                ++work;
                if(v==b) return true;
                if(visited[v]!=stamp) {visited[v]=stamp;stack.push_back(v);}
            }
        }
        return false;
    };
    std::map<int,std::vector<size_t>> groups;
    for(size_t id:root) groups[nodes[id].group].push_back(id);
    for(unsigned pass=0;pass<passes && work<refinementBudget;++pass) {
        activePass=pass+1;
        bool changed=false;
        for(auto& group:groups) {
            reportRefinement();
            const auto& ids=group.second;
            for(size_t i=0;i<ids.size() && work<refinementBudget;++i) {
                size_t a=find(ids[i]); if(a!=ids[i]) continue;
                for(size_t j=i+1;j<ids.size() && work<refinementBudget;++j) {
                    ++work;size_t b=find(ids[j]);if(a==b || b!=ids[j]) continue;
                    bool conflict=false;
                    for(const auto& row:nodes[b].rows) if(nodes[a].rows.count(row.first)) {conflict=true;break;}
                    if(conflict || reachable(a,b) || reachable(b,a)) continue;
                    for(size_t v:out[b]) {in[v].erase(b);in[v].insert(a);out[a].insert(v);}
                    for(size_t u:in[b]) {out[u].erase(b);out[u].insert(a);in[a].insert(u);}
                    out[b].clear();in[b].clear();
                    nodes[a].rows.insert(nodes[b].rows.begin(),nodes[b].rows.end());
                    nodes[b].rows.clear();representative[b]=a;
                    ++result.stats.merged;changed=true;
                }
            }
            if(work>=refinementBudget) break;
        }
        if(!changed) break;
    }
    reportRefinement(true);
    if(progress)progress("Refinement complete",result.stats.merged,0,
        work>=refinementBudget?"work budget reached":"merged instances",true);
    result.stats.refinementWork=work;
    if(progress)progress("Ordering refined graph",0,2,"stable and continuity layouts",true);
    size_t orderAttempt=0;
    std::vector<size_t> rank(nodes.size());
    for(size_t i=0;i<root.size();++i) rank[root[i]]=i;
    auto order=[&](bool continuity) {
        std::vector<size_t> degree(nodes.size());
        std::set<std::pair<size_t,size_t>> ready;
        for(size_t id:root) if(find(id)==id) {
            degree[id]=in[id].size();if(!degree[id]) ready.insert({rank[id],id});
        }
        std::vector<size_t> sorted;
        size_t preferred=size_t(-1);
        while(!ready.empty()) {
            auto it=ready.begin();
            if(continuity && preferred!=size_t(-1)) {
                auto candidate=ready.find({rank[preferred],preferred});
                if(candidate!=ready.end()) it=candidate;
            }
            size_t id=it->second;ready.erase(it);sorted.push_back(id);preferred=size_t(-1);
            if(progress && (sorted.size() & 65535)==0)progress("Ordering refined graph",orderAttempt,2,
                "ordered instances="+std::to_string(sorted.size()),false);
            // Prefer a newly ready successor with most shared tip paths.
            size_t weight=0;
            for(size_t v:out[id]) {
                if(--degree[v]==0) ready.insert({rank[v],v});
                if(degree[v]) continue;
                size_t shared=0;
                for(const auto& row:nodes[id].rows) {
                    auto other=nodes[v].rows.find(row.first);
                    if(other!=nodes[v].rows.end() && row.second.first==other->second.first) ++shared;
                }
                if(shared>weight || (shared==weight && preferred!=size_t(-1) && rank[v]<rank[preferred])) {
                    weight=shared;preferred=v;
                }
            }
        }
        if(sorted.size()+result.stats.merged!=nodes.size()) throw std::runtime_error("Cyclic MAF layout");
        ++orderAttempt;
        if(progress)progress("Ordering refined graph",orderAttempt,2,"ordered instances="+std::to_string(sorted.size()),orderAttempt==2);
        return sorted;
    };
    auto fragments=[&](const std::vector<size_t>& sorted) {
        std::map<size_t,size_t> lastPath;size_t count=0;
        for(size_t id:sorted) for(const auto& row:nodes[id].rows) {
            auto old=lastPath.find(row.first);
            if(old==lastPath.end() || old->second!=row.second.first) ++count;
            lastPath[row.first]=row.second.first;
        }
        return count;
    };
    auto sorted=order(false), candidate=order(true);
    if(fragments(candidate)<fragments(sorted)) sorted=std::move(candidate);
    std::vector<size_t> finalId(nodes.size());
    std::vector<Instance> compact;compact.reserve(sorted.size());
    for(size_t id:sorted) {finalId[id]=compact.size();compact.push_back(std::move(nodes[id]));}
    for(auto& path:paths) for(size_t& id:path.instances) id=finalId[find(id)];
    nodes=std::move(compact);
    // Explicitly validate reconstruction order and ownership after refinement.
    if(progress)progress("Validating layout",0,paths.size(),"tip paths",true);
    size_t validatedPaths=0;
    for(const auto& path:paths) {
        for(size_t i=0;i<path.instances.size();++i) {
            if(nodes[path.instances[i]].group!=path.groups[i] ||
               (i && path.instances[i]<=path.instances[i-1]))
                throw std::runtime_error("MAF refinement violated a tip path");
        }
        ++validatedPaths;
        if(progress)progress("Validating layout",validatedPaths,paths.size(),"tip paths",validatedPaths==paths.size());
    }
    return result;
}
}}
