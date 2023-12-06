from treeswift import *
# from ete3 import Tree
from sys import argv
from readJson import *
import json
import time
import mutation_annotation_pb2
from copy import deepcopy
import time
import sys
import networkx as nx



f = open("SARS100.gfa", "r")
data = f.readlines()
f.close()

newick_read = open("/home/AD.UCSD.EDU/swalia/data/sars/pangraph/sars_100.nwk", "r")
tree_read = newick_read.readlines()[0].split("\n")[0]
tree_old = tree_read
newick_read.close()
tree = read_tree_newick(tree_old)

innerNode = 0
for node in tree.traverse_preorder():
    if (node.get_label() == None):
        node.set_label(innerNode)
        innerNode += 1
    else:
        continue

def data_extract_from_gfa(data):
    segments = dict()
    paths = dict()
    for each_line in data:
        each_type = each_line[0]
        fields = each_line.split("\t")
        fields[-1] = fields[-1].split("\n")[0]

        if (each_type == "H"):
            version = fields[1].split(":")[-1]
            print(version)
        
        elif (each_type == "S"):
            sequence = fields[2].split("\n")[0]
            segments[fields[1]] = sequence

        elif (each_type == "P"):
            sequence_name = fields[1]
            path = fields[2].split(",")
            for i in range (len(path)):
                # Currently does not suppots reverse complement 
                path[i] = path[i].split("+")[0]

            paths[sequence_name] = path
        else:
            continue

    return (segments, paths)

def blockAlignment(seq1, seq2):
    INF = 10000
    VER = 2
    HOR = 1
    DIAG = 0
    len1 = len (seq1)
    len2 = len (seq2)

    if (len1 == 0 and len2 == 0):
        return ([],[],[])
    elif (len1 == 0):
        alnS2 = seq2
        alnS1 = len(seq2) * ["-"]
        consSeq = alnS2
        return (consSeq, alnS1, alnS2)
    elif (len2 == 0):
        alnS1 = seq1
        alnS2 = len(seq1) * ["-"]
        consSeq = alnS1
        return (consSeq, alnS1, alnS2)
    else:
        ## Score Matrix initialisation
        S = [[None for k in range (len2 + 1)] for j in range (len1 + 1)]
        S [0][0] = (0, None)


        ## penalties
        mismatch = -INF
        match = 4
        gap = -1

        ## Filling Score Matrix axes
        for i in range (1, len1 + 1): # s1 axis
            S [i][0] = (S [i - 1][0][0] + gap, VER) 
        for j in range(1, len2 + 1): # s2 axis
            S [0][j] = (S [0][j - 1][0] + gap, HOR)
        # print (S)

        ## Filling rest of Score Matrix
        for i in range(1, len1 + 1): # s1 axis
            for j in range(1, len2 + 1): # s2 axis
                if (seq1 [i - 1] == seq2 [j - 1]):
                    reward = match
                else:
                    reward = mismatch
                verScore  = S [i - 1][j][0]     + gap
                horScore  = S [i][j - 1][0]     + gap
                diagScore = S [i - 1][j - 1][0] + reward

                listScore = [diagScore, horScore, verScore]
                maxValue  = max (listScore)
                maxIdx    = listScore.index(maxValue)

                S[i][j] = (maxValue, maxIdx)

        ## Traceback
        alnS1 = []; alnS2 = []
        i = len1
        j = len2

        while ((i > 0) or (j > 0)):
            pointer = S[i][j][1]

            if (pointer == 0):
                alnS1 = [seq1 [i - 1]] + alnS1  
                alnS2 = [seq2 [j - 1]] + alnS2  
                i -= 1; j -=1
            
            elif (pointer == 1):
                alnS1 = ["-"]          + alnS1  
                alnS2 = [seq2 [j - 1]] + alnS2
                j -= 1

            else:
                alnS1 = [seq1 [i - 1]] + alnS1  
                alnS2 = ["-"]          + alnS2
                i -= 1


        if (len (alnS1) != len (alnS2)):
            print("Error")
            exit()

        else:
            consSeq = list()
            for i in range (len (alnS1)):
                if (alnS1 [i] == "-"):
                    consSeq.append(alnS2 [i])
                else:
                    consSeq.append(alnS1 [i])

        return (consSeq, alnS1, alnS2)

segments, paths = data_extract_from_gfa(data)

## Build a DAG
graph = nx.DiGraph()
for each_sequnece in paths:
    edges = []
    length_path = len(paths[each_sequnece])
    for i in range(length_path):
        if (i >= length_path - 1):
            continue
        else:
            edges.append((paths[each_sequnece][i], paths[each_sequnece][i+1]))
    
    graph.add_edges_from(edges)


global_list = list(nx.topological_sort(graph))

aligned_paths = dict()
count = 0
for each_sequnece in paths:
    print (count)
    _, _, aligned_paths[each_sequnece] = blockAlignment(global_list, paths[each_sequnece])
    count += 1
consSeq = []
mGap = dict()
block_id = dict()

for coor in range(len(global_list)):
    print (coor)
    currentChar = dict()
    for node in tree.traverse_postorder():
        nodeLabel = node.get_label()
        if (node.is_root()):
            rootLabel = nodeLabel
        
        if (node.is_leaf()):
            currentChar[nodeLabel] = [aligned_paths[nodeLabel][coor]]

        else:
            childrenList = node.child_nodes()
            X = set(currentChar[childrenList[0].get_label()])
            Y = set(currentChar[childrenList[1].get_label()])
            if (list (X & Y) == []):
                currentChar[nodeLabel] = list (X | Y)
            else:
                currentChar[nodeLabel] = list (X & Y)

    for node in tree.traverse_preorder():
        nodeLabel = node.get_label()
        
        if node.is_root():
            if "-" in currentChar[nodeLabel]:
                if len(currentChar[nodeLabel]) > 1:
                    currentChar[nodeLabel].remove("-")
                    currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                else:
                    currentChar[nodeLabel] = ["-"]
            else:
                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]

            if currentChar[nodeLabel][0] != "-":
                if rootLabel not in mGap:
                    mGap[rootLabel] = ["I:" + str(coor)]
                else:
                    mGap[rootLabel].append("I:" + str(coor))



        else:
            
            parentChar = currentChar[node.get_parent().get_label()][0]
            if parentChar in currentChar[nodeLabel]:
                currentChar[nodeLabel] = [parentChar]
            else:
                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                ## Insertion
                if parentChar == "-":
                    if nodeLabel not in mGap:
                        mGap[nodeLabel] = ["I:" + str(coor)]                                
                    else:
                        mGap[nodeLabel].append("I:" + str(coor))                               

                ## Deletion
                else:
                    if nodeLabel not in mGap:
                        mGap[nodeLabel] = ["D:" + str(coor)]                                
                    else:
                        mGap[nodeLabel].append("D:" + str(coor))

def seq2int( seq ):
    int_list = list()
    bin_compact = ""
    for i in range( len(seq) ):
        bin = char2bin( seq[i] )
        bin_compact += bin
        
        if( len(bin_compact) == 32):
            dec = int( bin_compact, 2 )
            bin_compact = ""
            int_list.append(dec)

        elif ( i == len(seq) - 1):
            bin_compact =  bin_compact + "0" * (32 - len(bin_compact))
            dec = int( bin_compact, 2 )
            bin_compact = ""
            int_list.append(dec)
        
        else:
            continue
    return (int_list)

def char2bin(char):
    bin = "1111"
    if ( char == "A" or char == "a"):
        bin = "0001"
    elif (char == "C" or char == "c"):
        bin = "0010"
    elif (char == "G" or char == "g"):
        bin = "0100"
    elif (char == "T" or char == "t"):
        bin = "1000"
    elif (char == "R" or char == "r"):
        bin = "0101"
    elif (char == "Y" or char == "y"):
        bin = "1010"
    elif (char == "S" or char == "s"):
        bin = "0110"
    elif (char == "W" or char == "w"):
        bin = "1001"
    elif (char == "K" or char == "k"):
        bin = "1100"
    elif (char == "M" or char == "m"):
        bin = "0011"
    elif (char == "B" or char == "b"):
        bin = "1110"
    elif (char == "D" or char == "d"):
        bin = "1101"
    elif (char == "H" or char == "h"):
        bin = "1011"
    elif (char == "V" or char == "v"):
        bin = "0111"
    else:
        bin = "1111"
    return (bin)


MAT = mutation_annotation_pb2.tree()

# Read the address book.
try:
  f = open(sys.argv[1], "rb")
  MAT.ParseFromString(f.read())
  f.close()
except IOError:
  print (sys.argv[1] + ": Could not open MAT.  Creating a new one.")


# Newick
MAT.newick = tree_old

# Nodes
for node in tree.traverse_preorder():
    nodeLabel = node.get_label()
    newNode = MAT.nodes.add()

    if nodeLabel in mGap:
        for each_mut in mGap[nodeLabel]: 
            fields = each_mut.split(":")
            
            newMut = newNode.mutations.add()

            newMut.blockId = int(fields[1])
            newMut.blockGapExist = False
            if (fields[0] == "I"):
                newMut.blockMutInfo = True
            else:
                newMut.blockMutInfo = False


# blocks
for i in range (len(global_list)):
    each_block = global_list[i]
    blocks = MAT.blocks.add()
    blocks.blockId = i
    blocks.blockGapExist = False
    blocks.consensusSeq.extend(seq2int(segments[each_block]))

# No gaps and blockGaps

f = open(sys.argv[1], "wb")
f.write(MAT.SerializeToString())
f.close()