#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <atomic>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/task_scheduler_init.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <json/json.h>
#include "panman.capnp.h"
#include "common.hpp"

#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <kj/std/iostream.h>

void writeNodes(capnp::List<panman::Node>::Builder& nodesBuilder, size_t nodeIndex){
    if (nodeIndex==3) return;
    panman::Node::Builder n = nodesBuilder[nodeIndex++];
    
    capnp::MallocMessageBuilder message;
    panman::Mutation::Builder mut_ = message.initRoot<panman::Mutation>();
    capnp::List<panman::NucMut>::Builder nm = mut_.initNucMutation(2);

    std::map< std::pair< int32_t, int32_t >, std::pair< std::vector< panman::NucMut::Builder >, int > > blockToMutations;



    for(size_t i = 0; i < 2; i++) {
        nm[i].setNucPosition(i);
        nm[i].setNucGapPosition(0);
        nm[i].setMutInfo(i*nodeIndex<< 16);

        blockToMutations[std::make_pair(i, 0)].first.push_back(nm[i]);

    }
    ::capnp::List<panman::Mutation>::Builder mutationsBuilder = n.initMutations(blockToMutations.size());
    size_t blockToMutationsCount=0;
    for(auto &u: blockToMutations) {
        panman::Mutation::Builder mutation = mutationsBuilder[blockToMutationsCount++];
        mutation.setBlockMutExist((u.second.second != 2));
        mutation.setBlockMutInfo(u.second.second);

        int32_t primaryBlockId = u.first.first;
        int32_t secondaryBlockId = u.first.second;
        if(secondaryBlockId != -1) {
            mutation.setBlockId(((int64_t)primaryBlockId << 32) + secondaryBlockId);
            mutation.setBlockGapExist(true);
        } else {
            mutation.setBlockId(((int64_t)primaryBlockId << 32));
            mutation.setBlockGapExist(false);
        }
        ::capnp::List<panman::NucMut>::Builder nucMutationBuilder = mutation.initNucMutation(u.second.first.size());
        for(auto i=0; i<u.second.first.size();i++) {
            std::cout << "\tBefore " << " " << 
                        u.second.first[i].getMutInfo() << " " <<
                        u.second.first[i].getNucPosition() << " " <<
                        u.second.first[i].getNucGapPosition() << " " <<
                        u.second.first[i].getNucGapExist() << " " <<
                    std::endl;
            nucMutationBuilder[i].setMutInfo(u.second.first[i].getMutInfo());
            nucMutationBuilder[i].setNucPosition(u.second.first[i].getNucPosition());
            nucMutationBuilder[i].setNucGapPosition(u.second.first[i].getNucGapPosition());
            nucMutationBuilder[i].setNucGapExist(u.second.first[i].getNucGapExist());
        }
    }

    writeNodes(nodesBuilder, nodeIndex);
}

void writeNewick(panman::Tree::Builder& treeToWrite){
    std::string newick = "((A,B),C);";
    treeToWrite.setNewick(newick);
}

void writeToProto(){
    std::ofstream outputFile("output.capnp");
    kj::std::StdOutputStream outputStream(outputFile);
    capnp::MallocMessageBuilder message;
    panman::TreeGroup::Builder treeGroupToWrite = message.initRoot<panman::TreeGroup>();

    capnp::List<panman::Tree>::Builder treestoWriteBuilder = treeGroupToWrite.initTrees(1);
    size_t treesCount = 0;
    panman::Tree::Builder treeToWrite = treestoWriteBuilder[treesCount++];
    writeNewick(treeToWrite);

    capnp::List<panman::Node>::Builder nodesBuilder = treeToWrite.initNodes(3);
    size_t nodeIndex=0;
    writeNodes(nodesBuilder, nodeIndex);
    ::capnp::writeMessage(outputStream, message);

}

class tree {
public:
    tree(const panman::Tree::Reader &inputTree) {
        protoToTree(inputTree);
    }
    void protoToTree(const panman::Tree::Reader &inputTree) {
        std::cout << inputTree.getNewick().cStr() << std::endl;
        for (auto v: inputTree.getNodes()){
            for (auto u: v.getMutations()){
                for (auto q: u.getNucMutation()){
                    std::cout << "\tAfter " << " " << 
                        q.getMutInfo() << " " <<
                        q.getNucPosition() << " " <<
                        q.getNucGapPosition() << " " <<
                        q.getNucGapExist() << " " <<
                    std::endl;
                }

            }
        }
    }   
    
};



void readProto() {
    std::ifstream inputFile("output.capnp");
    kj::std::StdInputStream inputStream(inputFile);
    capnp::InputStreamMessageReader messageReader(inputStream);
    panman::TreeGroup::Reader TG = messageReader.getRoot<panman::TreeGroup>();

    std::vector<tree> trees;

    for (auto treeFromTG: TG.getTrees()){
        trees.emplace_back(treeFromTG);
    }
}

int main(int argc, char ** argv) {

    writeToProto();

    readProto();

    return 0;
}