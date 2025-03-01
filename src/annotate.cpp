#include "panmanUtils.hpp"

void panmanUtils::Tree::annotate(std::ifstream& fin) {
    std::string line;
    char delim = '\t';
    while(getline(fin, line)) {
        std::string word;
        std::string nodeId;

        // Extract node ID
        size_t i = 0;
        for(; i < line.length() && line[i]!=delim; i++) {
            word+=line[i];
        }

        word = stripString(word);

        if(word.length()) {
            nodeId = word;
            word = "";
        } else {
            std::cout << "File in incorrect format. Line: " << line << std::endl;
            return;
        }

        if(i >= line.length()) {
            // comma not found
            std::cout << "File in incorrect format. Line: " << line << std::endl;
            return;
        }

        if(allNodes.find(nodeId) == allNodes.end()) {
            // std::cout << "Node ID not found. Line: " << nodeId << " [" << line << "]" << std::endl;
            // for (auto a: allNodes) {
            //     std::cout << a->
            // }
            return;
        }
        Node* nodeToAnnotate = allNodes[nodeId];
        
        std::cout << "node before annotation: " << nodeToAnnotate->identifier << " " << nodeToAnnotate->annotations.size() << std::endl;

        // Extract annotations
        for(; i < line.length(); i++) {
            if(line[i] != delim) {
                word += line[i];
            } else {
                word = stripString(word);
                if(word.length()) {
                    std::string annotation = word;
                    nodeToAnnotate->annotations.push_back(annotation);
                    annotationsToNodes[annotation].push_back(nodeId);
                    word = "";
                }
            }
        }

        word = stripString(word);
        if(word.length()) {
            std::string annotation = word;
            nodeToAnnotate->annotations.push_back(annotation);
            annotationsToNodes[annotation].push_back(nodeId);
            word = "";
        }

        std::cout << "node annotated: " << nodeToAnnotate->identifier << " " << nodeToAnnotate->annotations[0] << std::endl;

    }
}

std::vector< std::string > panmanUtils::Tree::searchByAnnotation(std::string annotation) {
    if(annotationsToNodes.find(annotation) != annotationsToNodes.end()) {
        return annotationsToNodes[annotation];
    }
    return {};
}
