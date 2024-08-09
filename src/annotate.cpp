#include "panmanUtils.hpp"

void panmanUtils::Tree::annotate(std::ifstream& fin) {
    std::string line;
    while(getline(fin, line)) {
        std::string word;
        std::string nodeId;

        // Extract node ID
        size_t i = 0;
        for(; i < line.length() && line[i]!=','; i++) {
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
            std::cout << "Node ID not found. Line: " << line << std::endl;
            return;
        }

        Node* nodeToAnnotate = allNodes[nodeId];

        // Extract annotations
        for(; i < line.length(); i++) {
            if(line[i] != ',') {
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

    }
}

std::vector< std::string > panmanUtils::Tree::searchByAnnotation(std::string annotation) {
    if(annotationsToNodes.find(annotation) != annotationsToNodes.end()) {
        return annotationsToNodes[annotation];
    }
    return {};
}
