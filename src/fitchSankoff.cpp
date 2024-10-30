
#include "panmanUtils.hpp"

int panmanUtils::Tree::nucFitchForwardPassOpt(
    Node* node,
    std::unordered_map< std::string, int >& states) {
    if(node->children.size() == 0) {
        return states[node->identifier];
    }

    std::vector< int > childStates;
    for(auto child: node->children) {
        childStates.push_back(nucFitchForwardPass(child, states));
    }

    int orStates = 0, andStates = childStates[0];

    for(auto u: childStates) {
        orStates |= u;
        andStates &= u;
    }

    if(andStates) {
        return states[node->identifier] = andStates;
    }
    return states[node->identifier] = orStates;
}

int panmanUtils::Tree::nucFitchForwardPass(Node* node,
        std::unordered_map< std::string, int >& states, int refState) {
    if(node->children.size() == 0) {
        if(states.find(node->identifier) == states.end()) {
            std::cerr << "Node ID not found" << std::endl;
            return states[node->identifier] = 0;
        }
        return states[node->identifier];
    }
    std::vector< int > childStates;
    for(auto child: node->children) {
        childStates.push_back(nucFitchForwardPass(child, states, refState));
    }
    //for root
    int orStates = 0, andStates = childStates[0];
    if (node->parent==nullptr) {
        return states[node->identifier] = refState;
    } 
    for(auto u: childStates) {
        orStates |= u;
        andStates &= u;
    }
    if(andStates) {
        return states[node->identifier] = andStates;
    }
    return states[node->identifier] = orStates;
}

void panmanUtils::Tree::nucFitchBackwardPassOpt(
    Node* node,
    std::unordered_map< std::string, int >& states,
    int parentState,
    int defaultState) {
    if(node == root && defaultState != (1 << 28)) {
        states[node->identifier] = defaultState;
    } else {
        if(states[node->identifier] == 0) {
            std::cout << "Issue\n";
            return;
        }
        if(node == root) {
            // The root sequence should take any of its values and not care about the parent state
            int currentState = 1;
            while(!(states[node->identifier] & currentState)) {
                currentState <<= 1;
            }
            states[node->identifier] = currentState;
            // } else if (node->children.size()==0) {
            //     return;
        } else if(parentState & states[node->identifier]) {
            states[node->identifier] = parentState;
        } else {
            int currentState = 1;
            while(!(states[node->identifier] & currentState)) {
                currentState <<= 1;
            }
            states[node->identifier] = currentState;
        }
    }

    for(auto child: node->children) {
        nucFitchBackwardPass(child, states, states[node->identifier]);
    }
}


void panmanUtils::Tree::nucFitchBackwardPass(Node* node,
        std::unordered_map< std::string, int >& states, int parentState, int defaultState) {
    if(node == root && defaultState != (1 << 28)) {
        states[node->identifier] = defaultState;
    } else {
        if(states[node->identifier] == 0) {
            return;
        }
        if(node == root) {
            // The root sequence should take any of its values and not care about the parent state
            int currentState = 1;
            while(!(states[node->identifier] & currentState)) {
                currentState <<= 1;
            }
            states[node->identifier] = currentState;
        } else if(parentState & states[node->identifier]) {
            states[node->identifier] = parentState;
        } else {
            int currentState = 1;
            while(!(states[node->identifier] & currentState)) {
                currentState <<= 1;
            }
            states[node->identifier] = currentState;
        }
    }

    for(auto child: node->children) {
        nucFitchBackwardPass(child, states, states[node->identifier]);
    }
}

void panmanUtils::Tree::nucFitchAssignMutations(Node* node,
        std::unordered_map< std::string, int >& states,
        std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType, char > >& mutations,
        int parentState) {

    if(states[node->identifier] == 0) {
        return;
    }

    if(parentState != states[node->identifier]) {
        if(parentState == 1) {
            // insertion
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0) {
                currentState >>= 1;
                code++;
            }
            code--;

            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NI, nuc);
        } else if(states[node->identifier] == 1) {
            // deletion
            mutations[node->identifier] = std::make_pair(NucMutationType::ND, '-');
        } else {
            // substitution
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0) {
                currentState >>= 1;
                code++;
            }
            code--;

            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NS, nuc);
        }
    }
    for(auto child: node->children) {
        nucFitchAssignMutations(child, states, mutations, states[node->identifier]);
    }
}

void panmanUtils::Tree::nucFitchAssignMutationsOpt(
    Node* node,
    std::unordered_map< std::string, int >& states,
    std::unordered_map< std::string, std::pair< panmanUtils::NucMutationType,
    char > >& mutations,
    int parentState) {

    if(states[node->identifier] == 0) {
        return;
    }

    if(parentState != states[node->identifier]) {
        if(parentState == 1) {
            // insertion
            int v = 0;
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0) {
                if(currentState & 1)
                    v += std::pow(2,code);
                currentState >>= 1;
                code++;
            }
            code--;

            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NI, nuc);
        } else if(states[node->identifier] == 1) {
            // deletion
            mutations[node->identifier] = std::make_pair(NucMutationType::ND, '-');
        } else {
            // substitution
            int v = 0;
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0) {
                if(currentState & 1)
                    v += std::pow(2,code);
                currentState >>= 1;
                code++;
            }
            code--;

            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NS, nuc);
        }
    }
    for(auto child: node->children) {
        nucFitchAssignMutations(child, states, mutations, states[node->identifier]);
    }
}


int panmanUtils::Tree::blockFitchForwardPassNew(Node* node,
        std::unordered_map< std::string, int >& states) {
    if(node->children.size() == 0) {
        if(states.find(node->identifier) == states.end()) {
            return states[node->identifier] = 0;
        }
        return states[node->identifier];
    }
    std::vector< int > childStates;
    for(auto child: node->children) {
        childStates.push_back(blockFitchForwardPassNew(child, states));
    }
    int orStates = 0, andStates = childStates[0];
    for(auto u: childStates) {
        orStates |= u;
        andStates &= u;
    }
    if(andStates) {
        return states[node->identifier] = andStates;
    }
    return states[node->identifier] = orStates;
}

void panmanUtils::Tree::blockFitchBackwardPassNew(Node* node,
        std::unordered_map< std::string, int >& states, int parentState, int defaultValue) {
    if(node == root && defaultValue != (1 << 28)) {
        states[node->identifier] = defaultValue;
    } else {
        if(states[node->identifier] == 0) {
            return;
        }
        if(parentState & states[node->identifier]) {
            states[node->identifier] = parentState;
        } else {
            int currentState = 1;
            while(!(states[node->identifier] & currentState)) {
                currentState <<= 1;
            }
            states[node->identifier] = currentState;
        }
    }

    for(auto child: node->children) {
        blockFitchBackwardPassNew(child, states, states[node->identifier]);
    }

}

void panmanUtils::Tree::blockFitchAssignMutationsNew(Node* node,
        std::unordered_map< std::string, int >& states,
        std::unordered_map< std::string,
        std::pair< panmanUtils::BlockMutationType, bool > >& mutations, int parentState) {
    if(states[node->identifier] == 0) {
        return;
    }
    if(parentState != states[node->identifier]) {
        if(parentState == 1) {
            // insertion
            int code = 0, currentState = states[node->identifier];
            while(currentState > 0) {
                currentState >>= 1;
                code++;
            }
            code--;
            if(code == 2) {
                // insertion of inverted block
                mutations[node->identifier] = std::make_pair(BlockMutationType::BI, true);
            } else {
                // insertion of forward strand
                mutations[node->identifier] = std::make_pair(BlockMutationType::BI, false);
            }

        } else if(states[node->identifier] == 1) {
            // deletion
            mutations[node->identifier] = std::make_pair(BlockMutationType::BD, false);
        } else {
            // inversion

            mutations[node->identifier] = std::make_pair(BlockMutationType::BD, true);
        }
    }
    for(auto child: node->children) {
        blockFitchAssignMutationsNew(child, states, mutations, states[node->identifier]);
    }
}


std::vector< int > panmanUtils::Tree::nucSankoffForwardPassOpt(Node* node,
        std::unordered_map< std::string, std::vector< int > >& stateSets) {

    if(node->children.size() == 0) {
        if(stateSets.find(node->identifier) == stateSets.end()) {
            std::vector< int > blankState(16, SANKOFF_INF);
            stateSets[node->identifier] = blankState;
        }
        return stateSets[node->identifier];
    }


    std::vector< std::vector< int > > childStates;
    for(auto child: node->children) {
        childStates.push_back(nucSankoffForwardPassOpt(child, stateSets));
    }

    bool minExists = false;
    for(size_t j = 0; j < childStates.size(); j++) {
        for(int k = 0; k < 16; k++) {
            if(childStates[j][k] < SANKOFF_INF) {
                minExists = true;
                break;
            }
        }
    }

    if(!minExists) {
        std::vector< int > currentState(16, SANKOFF_INF);
        return stateSets[node->identifier] = currentState;
    }

    std::vector< int > currentState(16, 0);
    for(int i = 0; i < 16; i++) {
        for(size_t j = 0; j < childStates.size(); j++) {
            int minVal = SANKOFF_INF;
            for(int k = 0; k < 16; k++) {
                minVal = std::min(minVal, (i != k) + childStates[j][k]);
            }
            if(minVal < SANKOFF_INF) {
                currentState[i] += minVal;
            }
        }
    }

    return stateSets[node->identifier] = currentState;
}

std::vector< int > panmanUtils::Tree::nucSankoffForwardPass(Node* node,
        std::unordered_map< std::string, std::vector< int > >& stateSets) {

    if(node->children.size() == 0) {
        if(stateSets.find(node->identifier) == stateSets.end()) {
            std::vector< int > blankState(16, SANKOFF_INF);
            stateSets[node->identifier] = blankState;
        }
        return stateSets[node->identifier];
    }


    std::vector< std::vector< int > > childStates;
    for(auto child: node->children) {
        childStates.push_back(nucSankoffForwardPass(child, stateSets));
    }

    bool minExists = false;
    for(size_t j = 0; j < childStates.size(); j++) {
        for(int k = 0; k < 16; k++) {
            if(childStates[j][k] < SANKOFF_INF) {
                minExists = true;
                break;
            }
        }
    }

    if(!minExists) {
        std::vector< int > currentState(16, SANKOFF_INF);
        return stateSets[node->identifier] = currentState;
    }

    std::vector< int > currentState(16, 0);
    for(int i = 0; i < 16; i++) {
        for(size_t j = 0; j < childStates.size(); j++) {
            int minVal = SANKOFF_INF;
            for(int k = 0; k < 16; k++) {
                minVal = std::min(minVal, (i != k) + childStates[j][k]);
            }
            if(minVal < SANKOFF_INF) {
                currentState[i] += minVal;
            }
        }
    }

    return stateSets[node->identifier] = currentState;
}

void panmanUtils::Tree::nucSankoffBackwardPass(Node* node,
        std::unordered_map< std::string, std::vector< int > >& stateSets,
        std::unordered_map< std::string, int >& states, int parentPtr,
        int defaultValue) {

    if(node == root && defaultValue != (1 << 28)) {
        states[node->identifier] = defaultValue;
    } else {
        if(node == root) {
            int minVal = SANKOFF_INF;
            int minPtr = -1;
            for(int i = 0; i < 16; i++) {
                if(stateSets[node->identifier][i] < minVal) {
                    minVal = stateSets[node->identifier][i];
                    minPtr = i;
                }
            }
            assert(minPtr != -1);

            states[node->identifier] = minPtr;
        } else {

            states[node->identifier] = parentPtr;
        }
    }

    if(states[node->identifier] == -1) {
        for(auto child: node->children) {
            nucSankoffBackwardPass(child, stateSets, states, -1);
        }
    } else {
        for(auto child: node->children) {
            int minPtr = -1;
            int minVal = SANKOFF_INF;
            for(int i = 0; i < 16; i++) {
                if((i != states[node->identifier]) + stateSets[child->identifier][i] < minVal) {
                    minVal = (i != states[node->identifier]) + stateSets[child->identifier][i];
                    minPtr = i;
                }
            }

            nucSankoffBackwardPass(child, stateSets, states, minPtr);
        }
    }
}

void panmanUtils::Tree::nucSankoffBackwardPassOpt(Node* node,
        std::unordered_map< std::string, std::vector< int > >& stateSets,
        std::unordered_map< std::string, int >& states, int parentPtr,
        int defaultValue) {

    if(node == root && defaultValue != (1 << 28)) {
        states[node->identifier] = defaultValue;
    } else {
        if(node == root) {
            int minVal = SANKOFF_INF;
            int minPtr = -1;
            for(int i = 0; i < 16; i++) {
                if(stateSets[node->identifier][i] < minVal) {
                    minVal = stateSets[node->identifier][i];
                    minPtr = i;
                }
            }
            assert(minPtr != -1);

            states[node->identifier] = minPtr;
        } else {

            states[node->identifier] = parentPtr;
        }
    }

    if(states[node->identifier] == -1) {
        for(auto child: node->children) {
            nucSankoffBackwardPassOpt(child, stateSets, states, -1);
        }
    } else {
        for(auto child: node->children) {
            int minPtr = -1;
            int minVal = SANKOFF_INF;
            for(int i = 0; i < 16; i++) {
                if((i != states[node->identifier]) + stateSets[child->identifier][i] < minVal) {
                    minVal = (i != states[node->identifier]) + stateSets[child->identifier][i];
                    minPtr = i;
                }
            }

            nucSankoffBackwardPassOpt(child, stateSets, states, minPtr);
        }
    }
}



void panmanUtils::Tree::nucSankoffAssignMutationsOpt(Node* node,
        std::unordered_map< std::string, int >& states,
        std::unordered_map< std::string,
        std::pair< panmanUtils::NucMutationType, char > >& mutations, int parentState) {
    if(states[node->identifier] == -1) {
        return;
    }
    if(parentState != states[node->identifier]) {
        if(parentState == 0) {
            // insertion
            int code = states[node->identifier];
            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NI, nuc);
        } else if(states[node->identifier] == 0) {
            // deletion
            mutations[node->identifier] = std::make_pair(NucMutationType::ND, '-');
        } else {
            // substitution
            int code = states[node->identifier];
            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NS, nuc);
        }
    }

    for(auto child: node->children) {
        nucSankoffAssignMutations(child, states, mutations, states[node->identifier]);
    }
}


void panmanUtils::Tree::nucSankoffAssignMutations(Node* node,
        std::unordered_map< std::string, int >& states,
        std::unordered_map< std::string,
        std::pair< panmanUtils::NucMutationType, char > >& mutations, int parentState) {
    if(states[node->identifier] == -1) {
        return;
    }
    if(parentState != states[node->identifier]) {
        if(parentState == 0) {
            // insertion
            int code = states[node->identifier];
            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NI, nuc);
        } else if(states[node->identifier] == 0) {
            // deletion
            mutations[node->identifier] = std::make_pair(NucMutationType::ND, '-');
        } else {
            // substitution
            int code = states[node->identifier];
            char nuc = getNucleotideFromCode(code);
            mutations[node->identifier] = std::make_pair(NucMutationType::NS, nuc);
        }
    }

    for(auto child: node->children) {
        nucSankoffAssignMutations(child, states, mutations, states[node->identifier]);
    }
}



std::vector< int > panmanUtils::Tree::blockSankoffForwardPass(Node* node,
        std::unordered_map< std::string, std::vector< int > >& stateSets) {

    if(node->children.size() == 0) {
        if(stateSets.find(node->identifier) == stateSets.end()) {
            std::vector< int > blankState({0,SANKOFF_INF,SANKOFF_INF});
            stateSets[node->identifier] = blankState;
        }
        return stateSets[node->identifier];
    }

    std::vector< std::vector< int > > childStates;
    for(auto child: node->children) {
        childStates.push_back(blockSankoffForwardPass(child, stateSets));
    }

    std::vector< int > currentState(3);
    for(int i = 0; i < 3; i++) {
        for(size_t j = 0; j < childStates.size(); j++) {
            int minVal = SANKOFF_INF;
            for(int k = 0; k < 3; k++) {
                minVal = std::min(minVal, (i != k) + childStates[j][k]);
            }
            currentState[i] += minVal;
        }
    }

    return stateSets[node->identifier] = currentState;
}

void panmanUtils::Tree::blockSankoffBackwardPass(Node* node,
        std::unordered_map< std::string, std::vector< int > >& stateSets,
        std::unordered_map< std::string, int >& states, int parentPtr,
        int defaultValue) {

    if(node == root && defaultValue != (1 << 28)) {
        states[node->identifier] = defaultValue;
    } else {
        if(node == root) {
            int minVal = SANKOFF_INF;
            int minPtr = -1;
            for(int i = 0; i < 3; i++) {
                if(stateSets[node->identifier][i] < minVal) {
                    minVal = stateSets[node->identifier][i];
                    minPtr = i;
                }
            }
            if(minPtr == -1) {
                states[node->identifier] = -1;
                return;
            }

            states[node->identifier] = minPtr;
        } else {
            bool stateExists = false;
            for(int i = 0; i < 3; i++) {
                if(stateSets[node->identifier][i] < SANKOFF_INF) {
                    stateExists = true;
                }
            }
            if(!stateExists) {
                states[node->identifier] = -1;
                return;
            }
            states[node->identifier] = parentPtr;
        }
    }

    for(auto child: node->children) {
        int minPtr = -1;
        int minVal = SANKOFF_INF;
        for(int i = 0; i < 3; i++) {
            if((i != states[node->identifier]) + stateSets[child->identifier][i] < minVal) {
                minVal = (i != states[node->identifier]) + stateSets[child->identifier][i];
                minPtr = i;
            }
        }
        blockSankoffBackwardPass(child, stateSets, states, minPtr);
    }
}

void panmanUtils::Tree::blockSankoffAssignMutations(Node* node,
        std::unordered_map< std::string, int >& states,
        std::unordered_map< std::string,
        std::pair< panmanUtils::BlockMutationType, bool > >& mutations, int parentState) {
    if(states[node->identifier] == -1) {
        return;
    }
    if(parentState != states[node->identifier]) {
        if(parentState == 0) {
            // insertion
            int code = states[node->identifier];
            if(code == 2) {
                // insertion of inverted block
                mutations[node->identifier] = std::make_pair(BlockMutationType::BI, true);
            } else {
                // insertion of forward strand
                mutations[node->identifier] = std::make_pair(BlockMutationType::BI, false);
            }

        } else if(states[node->identifier] == 0) {
            // deletion
            mutations[node->identifier] = std::make_pair(BlockMutationType::BD, false);
        } else {
            // inversion
            mutations[node->identifier] = std::make_pair(BlockMutationType::BD, true);
        }
    }
    for(auto child: node->children) {
        blockSankoffAssignMutations(child, states, mutations, states[node->identifier]);
    }
}
