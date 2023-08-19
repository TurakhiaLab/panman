from treeswift import *
from sys import argv
from readJson import *
import time
import sys
import networkx as nx
import networkx.algorithms.dag as dag

def distribution(blocksLen):
    max_ = 0
    fil_blocks = dict()
    for i in blocksLen:
        if blocksLen[i] > max_:
            max_ = blocksLen[i]
    count = 0
    cutoff = max_*0.4
    for i in blocksLen:
        if blocksLen[i] > cutoff:
            fil_blocks[i] = blocksLen[i]
            count += 1
    return max_, count, fil_blocks

def filter_(blockList, fil_blocks, common_block):
    filter_block_id_list = []
    filter_block_index_list = []
    for i in range(len(blockList)):
        if blockList[i] in fil_blocks:
            filter_block_id_list.append(blockList[i])
            filter_block_index_list.append(i)
            common_block[blockList[i]] += 1
    return (filter_block_id_list, filter_block_index_list)

def align(seq1, seq2):
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

C1=0
C2=0

def seq2graph(graph, block_list, duplicates):
    for i in range(len(block_list) - 1):
        x = block_list[i]
        y = block_list[i+1]
        if (i > 0):
            x1=block_list[i-1]
            y1=block_list[i]
            boolean = nx.has_path(graph, x1, y1)
            if (boolean == False):
                print("Error", i)
                break
        try:
            boolean = nx.has_path(graph, y, x)
            if (boolean):
                if y not in duplicates:
                    duplicates[y] = [y+"0"]
                    y = y + "0"
                    block_list[i+1] = y
                    graph.add_edges_from([(x,y)])
                else:
                    done = False
                    for d in duplicates[y]:
                        if (done == True):
                            break
                        try:
                            boolean = nx.has_path(graph, d, x)
                            if (boolean == False):
                                done = True
                                block_list[i+1] = d
                                graph.add_edges_from([(x,d)])
                                C1+=1
                        except:
                            block_list[i+1] = y
                            graph.add_edges_from([(x,d)])
                            done = True
                            C1+=1
                    if (done == False):
                        d = d + "0"
                        duplicates[y].append(d)
                        block_list[i+1] = d
                        graph.add_edges_from([(x,d)])
                        C2+=1

            else:
                graph.add_edges_from([(x,y)])
        except:
            graph.add_edges_from([(x,y)])


fjson = open(argv[1])
data = json.load(fjson)
fjson.close()

seqPaths, seqStrands = getPaths(data)
#### Writing seqpaths to files
count = 0
for p in seqPaths:
    path = "ecoli_" + str(count) + ".csv"
    f = open(path,"w")
    for b in seqPaths[p]:
        f.write(b)
        f.write(",")
    count += 1
    f.close()
####
exit()
blocksLen = getBlocksLen (data)
max_, count, fil_blocks = distribution(blocksLen)

common_block = dict()
for i in fil_blocks:
    common_block[i] = 0

block_filter = {}
block_filter_idx = {}
for seq in seqPaths:
    block_filter[seq], block_filter_idx[seq] = filter_(seqPaths[seq], fil_blocks, common_block)

## Removing duplicates with same sequence
for seq in seqPaths:
    map_ = {}
    i = 0
    for block in seqPaths[seq]:
        if block not in map_:
            map_[block] = 0
        else:
            seqPaths[seq][i] = seqPaths[seq][i] + str(map_[block])
            map_[block] += 1
        i += 1

# Create a DAG
graph = nx.DiGraph()
print(graph.number_of_nodes())
duplicates = {}
for seq in seqPaths:
    seq2graph(graph, seqPaths[seq], duplicates)
    print(seq)
print(graph.number_of_nodes())
print(len(list(nx.topological_sort(graph))))
print(C1, C2)

# for seq in seqPaths:
#     # print(seqPaths[seq])
#     # break
#     for i in range(len(seqPaths[seq]) - 1):
#         x = seqPaths[seq][i]
#         y = seqPaths[seq][i+1]
#         boolean = nx.has_path(graph, x, y)
#         if (boolean == False):
#             print(seq, "Error", x, y)


exit()


seqPaths = dict()
seqPathsSize = dict()
total_blks = dict()


valid_block = dict()
valid_block_idx = dict()
for path in paths:
    seq = path['name']
    offset = path['offset']
    # print(offset)
    seqPaths[seq] = list()
    seqPathsSize[seq] = list()
    for block in path ["blocks"]:
        blockId = block ["id"]
        seqPaths [seq].append (blockId)

        if blockId not in total_blks:
            total_blks[blockId] = 1
        else:
            total_blks[blockId] += 1
    positions = path['position']
    max_size = 0
    
    size = 0
    
    # to find size of blocks
    SIZE=True
    if (SIZE):
        for i in range(len(positions)):
            if (i == len(positions) - 1):
                break
            curr_pos = int(positions[i])
            next_pos = int(positions[i + 1])
            if (next_pos < curr_pos):
                size = next_pos
            else:
                size = next_pos - curr_pos
            seqPathsSize[seq].append(size)
        # print(max_size)
    
    # to find max
    MAX=False
    if (MAX):
        for i in range(len(positions)):
            if (i == len(positions) - 1):
                break
            curr_pos = int(positions[i])
            next_pos = int(positions[i + 1])
            if (next_pos < curr_pos):
                size = next_pos
            else:
                size = next_pos - curr_pos
            if (size > max_size):
                max_size = size
        print(max_size)

    # filteration
    FIL = True
    if (FIL):
        vb = list()
        vb_idx = list()
        for i in range(len(positions)):
            if (i == len(positions) - 1):
                break
            curr_pos = int(positions[i])
            next_pos = int(positions[i + 1])
            if (next_pos < curr_pos):
                size = next_pos
            else:
                size = next_pos - curr_pos
            if (size > 30000):
                vb.append(seqPaths[seq][i])
                vb_idx.append(i)
        valid_block[seq] = vb
        valid_block_idx[seq] = vb_idx

    #
    # blk = "YTZILRHUHM"
    # size = 0
    # for i in range(len(positions)):
    #     if (i == len(positions) - 1):
    #         break
    #     curr_pos = int(positions[i])
    #     next_pos = int(positions[i + 1])
    #     if (next_pos < curr_pos):
    #         size += next_pos
    #     else:
    #         size += next_pos - curr_pos
    #     if (size >= abs(offset)):
    #         print (i)
    #         break


for name in valid_block:
    print(name, valid_block_idx[name])  
print("Total blocks: ", len(total_blks))


# rotateArray(seqPaths["NZ_CP007391.1"], len(seqPaths["NZ_CP007391.1"]), 514)
# rotateArray(seqPaths["NC_002695.2"], len(seqPaths["NC_002695.2"]), 142)
# rotateArray(seqPaths["NZ_CP007265.1"], len(seqPaths["NZ_CP007265.1"]), 39)
# rotateArray(seqPaths["NZ_CP007390.1"], len(seqPaths["NZ_CP007390.1"]), 1161)
# rotateArray(seqPaths["NZ_CP035751.1"], len(seqPaths["NZ_CP035751.1"]), 1331)
# rotateArray(seqPaths["NC_000913.3"], len(seqPaths["NC_000913.3"]), 733)
# rotateArray(seqPaths["NZ_CP006262.1"], len(seqPaths["NZ_CP006262.1"]), 263)
# rotateArray(seqPaths["NZ_CP023541.1"], len(seqPaths["NZ_CP023541.1"]), 1499)
# rotateArray(seqPaths["NZ_CP006632.1"], len(seqPaths["NZ_CP006632.1"]), 1406)
# rotateArray(seqPaths["NZ_CP006027.1"], len(seqPaths["NZ_CP006027.1"]), 310)

for each in seqPaths:
    print(seqPaths[each][0])


##########################

# seqPaths = dict()
# seqPaths["A"] = ["a", "b", "c", "d", "e"]
# seqPaths["B"] = ["a", "b", "c", "a"]
# seqPaths["C"] = ["a", "b", "c", "d", "a"]


##########################


# build a DAG
DAG = True
if (DAG):
    graph = nx.DiGraph()
    count = 0
    for s in seqPaths:
        print(count)
        count += 1
        for i in range(len(seqPaths[s])):
            if (i >= len(seqPaths[s]) - 1):
                continue
            else:
                x = seqPaths[s][i]
                y = seqPaths[s][i+1]
                try:
                    boolean = nx.has_path(graph, y, x)
                    if (boolean):
                        y = y + "-" + x
                        seqPaths[s][i+1] = y
                        graph.add_edges_from([(x,y)])
                    else:
                        graph.add_edges_from([(x,y)])
                except:
                    graph.add_edges_from([(x,y)])

        print (len(list(nx.topological_sort(graph))))
        # if (count == 2):
        #     break

    global_list = list(nx.topological_sort(graph))

    print("Global list generated", len(global_list))

## Build a Graph
GRAPH = False
if (GRAPH):
    graph = nx.DiGraph()
    for s in seqPaths:
        for i in range(len(seqPaths[s])):
            if (i >= len(seqPaths[s]) - 1):
                continue
            else:
                x = seqPaths[s][i]
                y = seqPaths[s][i+1]
                graph.add_edges_from([(x,y)])


    for n, ns in graph.adjacency():
        print(n,ns)
        if (n == "d"):
            ns.pop("a")
        print(n,ns)
        
        # for items in ns.items():
            # print(items)

# LVJGTQEKPR HYJVKGPEHZ OMVRHLJYIA VCKDVNACYU
# FUDTRVFLVU HKEOTXMHXD
# NC_002695.2
# 3771043
# FUDTRVFLVU
# CGTATGTGGCTTCTGATGCGCAAGCTGAAGAAAAATGAGCATGGAGAATAATATGAATTTTTTAATGCGCGCTATATTCAGTCTGCTGTTGCTTTTTACTCTCTCTATTCCTGTCATTTCTGACTGTGTTGCAA
#  GTATGTGGCTTCTGATGCGCAAGCTGAAGAAAAATGAGCATGGAGAATAATATGAATTTTTTAATGCGCGCTATATTCAGTCTGCTGTTGCTTTTTACTCTCTCTATTCCTGTCATTTCTGACTGTGTTGCAATGGCCATTGAAAGTCGCTTCAAATATATGATGCTACTTTTTTAAATGGTTTTTACCTGTCGGCATCCGCTCAAAACGGGCGGTTGTCGATAAACGCTCACTTGGTTAATCATTTCACTCTTCA

def blockAlignment(seq1, seq2):
    INF = 10000
    VER = 2
    HOR = 1
    DIAG = 0
    len1 = len (seq1)
    len2 = len (seq2)

    # # # # # # Seq 2 # # # # 
    #
    #S
    #e
    #q 
    #1
    #
    #

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

    # print (S)

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

def rotateArray(arr, n, d):
    temp = []
    i = 0
    while (i < d):
        temp.append(arr[i])
        i = i + 1
    i = 0
    while (d < n):
        arr[i] = arr[d]
        i = i + 1
        d = d + 1
    arr[:] = arr[: i] + temp
    return arr

# path_list = list()
# for paths in seqPaths:
#     path_list.append(paths)

# a,b,c = blockAlignment(seqPaths[path_list[0]], seqPaths[path_list[1]])

# print(len(a))
#                     1           2               3           4               5           6               7           8               9           10              11          12              13
# NZ_CP035751.1 ['YTZILRHUHM', 'AIFOCERTKC', 'RGNVPKSXPD', 'JYBDVUIKKR', 'OAFXTYOAIL', 'XGBLBPHUSZ', 'RWIKZYCNTS', 'FMKSIYMNQA', 'EOCFVYXJHJ', 'RIXEEERGXW', 'OHQCHMAMGF', 'OINPFWAXOM', 'DGQWMSQLCL']
#                     1           2               4           5           6               7           8               9               10          11              12          13
# NC_000913.3 ['YTZILRHUHM', 'AIFOCERTKC', 'JYBDVUIKKR', 'OAFXTYOAIL', 'XGBLBPHUSZ', 'RWIKZYCNTS', 'FMKSIYMNQA', 'EOCFVYXJHJ', 'RIXEEERGXW', 'OHQCHMAMGF', 'OINPFWAXOM', 'DGQWMSQLCL']
#                     7               8           9              10           11             12           13          1               2               4           5              6
# NZ_CP006262.1 ['RWIKZYCNTS', 'FMKSIYMNQA', 'EOCFVYXJHJ', 'RIXEEERGXW', 'OHQCHMAMGF', 'OINPFWAXOM', 'DGQWMSQLCL', 'YTZILRHUHM', 'AIFOCERTKC', 'JYBDVUIKKR', 'OAFXTYOAIL', 'XGBLBPHUSZ']
#                     10            11            12            13            1              2            4           5               6           7               8               9
# NZ_CP023541.1 ['RIXEEERGXW', 'OHQCHMAMGF', 'OINPFWAXOM', 'DGQWMSQLCL', 'YTZILRHUHM', 'AIFOCERTKC', 'JYBDVUIKKR', 'OAFXTYOAIL', 'XGBLBPHUSZ', 'RWIKZYCNTS', 'FMKSIYMNQA', 'EOCFVYXJHJ']
#                     10            11            12            13            1              2            4           5               6             6.1          6.2           7               8               9
# NZ_CP006632.1 ['RIXEEERGXW', 'OHQCHMAMGF', 'OINPFWAXOM', 'DGQWMSQLCL', 'YTZILRHUHM', 'AIFOCERTKC', 'JYBDVUIKKR', 'OAFXTYOAIL', 'XGBLBPHUSZ', 'IDTCOIQVGG', 'SODUVLBJXB', 'RWIKZYCNTS', 'FMKSIYMNQA', 'EOCFVYXJHJ']
#                     7             8             9            10            11            12            13            1              2            4           5               6 
# NZ_CP006027.1 ['RWIKZYCNTS', 'FMKSIYMNQA', 'EOCFVYXJHJ', 'RIXEEERGXW', 'OHQCHMAMGF', 'OINPFWAXOM', 'DGQWMSQLCL', 'YTZILRHUHM', 'AIFOCERTKC', 'JYBDVUIKKR', 'OAFXTYOAIL', 'XGBLBPHUSZ']