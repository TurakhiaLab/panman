from treeswift import *
from sys import argv
from readJson import *
import json

INF = 10000
VER = 2
HOR = 1
DIAG = 0

def blockAlignment(seq1, seq2):
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
    
def funLCS(S1, S2, S1idx, S2idx):
    m = len (S1); n = len (S2)
    L = [[0 for x in range(n + 1)] for x in range(m + 1)]

    # Building the mtrix in bottom-up way
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif S1[i - 1] == S2[j - 1]:
                L[i][j] = L[i - 1][j - 1] + 1
            else:
                L[i][j] = max(L[i - 1][j], L[i][j - 1])

    index = L[m][n]

    lcs_algo = [""] * (index + 1)
    lcs_algo[index] = ""
    indexSeq1 = list()
    indexSeq2 = list()

    i = m
    j = n
    while i > 0 and j > 0:

        if S1[i - 1] == S2[j - 1]:
            lcs_algo[index - 1] = S1[i - 1]
            i -= 1
            j -= 1
            index -= 1
            indexSeq1 = [i] + indexSeq1
            indexSeq2 = [j] + indexSeq2

        elif L[i-1][j] > L[i][j-1]:
            i -= 1
        else:
            j -= 1

    return (lcs_algo[:-1], indexSeq1, indexSeq2)

def globalBlocks (seqPaths):
    idxList = list()
    seqList = list()
    for eachSeq in seqPaths:
        seqList.append(eachSeq)

    globalBlock = list ()
    blockSeq0 = seqPaths [seqList [0]]
    blockSeq1 = seqPaths [seqList [1]]
    globalBlock, idxSeq0, idxSeq1 = funLCS (blockSeq0, blockSeq1, len(blockSeq0), len(blockSeq1))

    for i in range (len (idxSeq1)):
        idxList.append ( [idxSeq0[i], idxSeq1[i]] )

    # print ( len(idxList), len(globalBlock) )

    seqCount = 2
    for eachSeq in seqList[2:]:
        seqCount += 1
        blockSeqi =  seqPaths [eachSeq]
        globalBlock, idxSeqi, idxSeqGlobal = funLCS (blockSeqi, globalBlock, len(blockSeqi), len(globalBlock)) 

        for i in range (len (idxSeqGlobal)):
            idxList [idxSeqGlobal [i]].append(idxSeqi [i])

        idxListFinal = list()
        for element in idxList:
            if (len(element) == seqCount):
                idxListFinal.append (element)
        idxList = idxListFinal
    
    return (globalBlock, idxList) 

f = open("/home/AD.UCSD.EDU/swalia/data/ecoli/pangraph/ecoli_10.json","r")
data = json.load(f)
f.close()

## Import Tree
tree = read_tree_newick("/home/AD.UCSD.EDU/swalia/data/ecoli/pangraph/ecoli_10.nwk")
print(tree)
innerNode = 0
for node in tree.traverse_postorder():
    if (node.get_label() == None):
        node.set_label(innerNode)
        innerNode += 1
    else:
        # print (node.get_label())
        continue
    # print (node.get_label())

seqPaths = getPaths (data)
seqName = [i for i in seqPaths]

gBlocks, gIdx = globalBlocks (seqPaths)

### Generating Gaps #####
gapListGlobal = []
for i in range (len (gIdx)):

    gapListLocal = []
    for j in range (len (gIdx [i])):

        currIdx = gIdx[i][j]
        if (i == 0):
            prevIdx = 0
        else:
            prevIdx = gIdx[i-1][j] + 1

        gapListLocal, X, Y = blockAlignment(gapListLocal, seqPaths[seqName[j]][prevIdx: currIdx])

    gapListGlobal.append (gapListLocal)


### Appending gaps in each Sequence ###
blockList = dict()
for j in range (len (seqName)):
    
    blockListLocal = []

    for i in range (len (gIdx)):
        currIdx = gIdx[i][j]
        if (i == 0):
            prevIdx = 0
        else:
            prevIdx = gIdx[i-1][j] + 1
        
        gapListLocal, X, Y = blockAlignment(gapListGlobal[i], seqPaths[seqName[j]][prevIdx: currIdx])

        blockListLocal.append(Y)
        if (i <= (len(gIdx) - 1)):
            blockListLocal.append(gBlocks[i])

    blockList[seqName[j]] = blockListLocal

blockListIdx = dict()
tempMapBlocks = dict()
for j in range (len (seqName)):
    tempMapBlocks [seqName [j]] = dict()
    blockListIdx [seqName [j]] = list()

    for i in range (len (blockList [seqName [j]])):
        elem = blockList [seqName [j]][i]

        if (type (elem) != list):
            if elem == "-":
                blockListIdx [seqName [j]].append("-")
                continue
            elif elem not in tempMapBlocks [seqName [j]]:
                tempMapBlocks [seqName [j]][elem] = 1
            else:
                tempMapBlocks [seqName [j]][elem] += 1
            blockListIdx [seqName [j]].append(tempMapBlocks [seqName [j]][elem])

        else:
            blockListIdx [seqName [j]].append(list())
            for z in range (len (elem)):
                element = blockList [seqName [j]][i][z]
                if element == "-":
                    blockListIdx [seqName [j]][i].append("-")
                    continue
                elif element not in tempMapBlocks [seqName [j]]:
                    tempMapBlocks [seqName [j]][element] = 1
                else:
                    tempMapBlocks [seqName [j]][element] += 1
                blockListIdx [seqName [j]][i].append(tempMapBlocks [seqName [j]][element])


'''
### Check if Seq lenght is correct
for j in range (len (seqName)):
    seqLen = len (seqPaths[seqName[j]])

    seqLenWithGap = 0
    seqLenWithoutGap = 0
    for elem in blockList[seqName[j]]:
        if type(elem) == list:
            seqLenWithGap += len(elem)
            for element in elem:
                if (element != "-"):
                    seqLenWithoutGap += 1
        else:
            seqLenWithGap += 1
            seqLenWithoutGap += 1

    
    print (seqLen, seqLenWithoutGap, seqLenWithGap)
'''

### Apply Fitch on Gaps ###
gMutation = dict()
for globalIdx in range (len (blockList[seqName[0]])):
    if (globalIdx%2 == 0):
        # for i in range (len (seqName)):
        #     print (seqName[i], (blockList [seqName[i]][globalIdx]))
        gapIdx = int(globalIdx/2)
        gLen = len (blockList [seqName[0]][globalIdx])

        for i in range (gLen):

            currentChar = dict()
            gMutationLocal = dict()
            for node in tree.traverse_postorder():
                nodeLabel = node.get_label()

                if (node.is_root()):
                    rootLabel = nodeLabel
                
                if (node.is_leaf()):
                    currentChar[nodeLabel] = [blockList [nodeLabel][globalIdx][i]]
                    # print (nodeLabel, blockList [nodeLabel][gIdx][i])

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
                if (node.is_leaf() == False):
                    currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                    if (node.is_root()):
                        if (currentChar[nodeLabel][0] != "-"):
                            gMutationLocal[nodeLabel] = str(globalIdx) + ":" + str(i) + ":I:" + str(currentChar[nodeLabel][0])
                        # else:
                        #     gMutationLocal[nodeLabel] = ""

                    childrenList = node.child_nodes()
                    child0Label = childrenList[0].get_label()
                    child1Label = childrenList[1].get_label()
                    
                    Z = set(currentChar[nodeLabel])
                    X = set(currentChar[child0Label])
                    if (list (X & Z) == []):
                        currentChar[child0Label] = [currentChar[child0Label][0]]
                        if (currentChar[child0Label][0] == "-"):
                            gMutationLocal[child0Label] = str(globalIdx) + ":" + str(i) + ":D"
                        else:
                            gMutationLocal[child0Label] = str(globalIdx) + ":" + str(i) + ":I:" + str(currentChar[child0Label][0])
                    else:
                        currentChar[child0Label] = currentChar[nodeLabel]

                    Y = set(currentChar[child1Label])
                    if (list (Y & Z) == []):
                        currentChar[child1Label] = [currentChar[child1Label][0]]
                        if (currentChar[child1Label][0] == "-"):
                            gMutationLocal[child1Label] = str(globalIdx) + ":" + str(i) + ":D"
                        else:
                            gMutationLocal[child1Label] = str(globalIdx) + ":" + str(i) + ":I:" + str(currentChar[child1Label][0])
                    else:
                        currentChar[child1Label] = currentChar[nodeLabel]


            for element in gMutationLocal:
                if element not in gMutation:
                    gMutation[element] = [gMutationLocal[element]]
                else:
                    gMutation[element].append(gMutationLocal[element])
    
    else:
         for node in tree.traverse_preorder():
            if (node.is_root()):
                nodeLabel = node.get_label()
                if nodeLabel in gMutation:
                    gMutation[nodeLabel].append(str(globalIdx) + ":I:" + blockList[seqName[0]][globalIdx])
                else:
                    gMutation[nodeLabel] = [str(globalIdx) + ":I:" + blockList[seqName[0]][globalIdx]]



### Apply fitch on Substitutions ###
blockData = getBlocks (data)
blockMap = dict()
counterGlobal = 0
counterLocal = 0
gMutationLocal = dict()
for i in range (len (blockList [seqName[0]])):

    # elem - just to check if it's a list or not 
    elem = blockList[seqName [0]][i]

    if type(elem) == list:
        print ("It's a list")
        for y in range (len (elem)):
            element = gapListGlobal[int(i/2)][y]
            print (element)
            blockInfo = getBlockInfo (element, blockData)
            consSeq   = getBlockSeq (blockInfo)
            blockSubs = getBlockSubs (blockInfo)

            print (blockSubs)
            blockOcc = dict()
            for z in blockList:
                if blockListIdx [z][i][y] != "-":
                    blockOcc [z] = blockListIdx [z][i][y]
            
            # print (blockOcc)

            coorList = set()
            for z in blockOcc:
                coorList.update (blockSubs[z][blockOcc [z]][0])
            
            print (coorList, len(consSeq))
            for coor in coorList:
                ## Apply fitch for each coordinate
                currentChar = dict()
                for node in tree.traverse_postorder():
                    nodeLabel = node.get_label()
                    if (node.is_root()):
                        rootLabel = nodeLabel
                    
                    if (node.is_leaf()):
                        if nodeLabel not in blockOcc:
                            currentChar[nodeLabel] = [consSeq[coor]]
                        else:
                            idxList  = blockSubs [nodeLabel][blockOcc[nodeLabel]][0]
                            charList = blockSubs [nodeLabel][blockOcc[nodeLabel]][1]
                            if coor in idxList:
                                charIdx = idxList.index(coor)
                                currentChar[nodeLabel] = [charList[charIdx]]
                            else:
                                currentChar[nodeLabel] = [consSeq[coor]]
                        # print (nodeLabel, blockList [nodeLabel][gIdx][i])

                    else:
                        childrenList = node.child_nodes()
                        X = set(currentChar[childrenList[0].get_label()])
                        Y = set(currentChar[childrenList[1].get_label()])
                        if (list (X & Y) == []):
                            currentChar[nodeLabel] = list (X | Y)
                        else:
                            currentChar[nodeLabel] = list (X & Y)

                # print (currentChar)
                for node in tree.traverse_preorder():
                    nodeLabel = node.get_label()
                    
                    if node.is_root():
                        if consSeq[coor] in currentChar:
                            currentChar[nodeLabel] = [consSeq[coor]]
                        else:
                            currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                            if nodeLabel not in gMutationLocal:
                                gMutationLocal[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor) + ":" + currentChar[nodeLabel][0]]
                            else:
                                gMutationLocal[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor) + ":" + currentChar[nodeLabel][0])

                    else:
                        if currentChar[node.get_parent().get_label()][0] in currentChar[nodeLabel]:
                            currentChar[nodeLabel] = [currentChar[node.get_parent().get_label()][0]]
                        else:
                            currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                            if nodeLabel not in gMutationLocal:
                                gMutationLocal[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor) + ":" + currentChar[nodeLabel][0]]
                            else:
                                gMutationLocal[nodeLabel].append (str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor) + ":" + currentChar[nodeLabel][0])
            counterLocal += 1

    else:
        print ("It's not a list")

        blockInfo = getBlockInfo (elem, blockData)
        consSeq   = getBlockSeq (blockInfo)
        blockSubs = getBlockSubs (blockInfo)
        
        blockOcc = dict()
        for z in blockList:
            blockOcc [z] = blockListIdx [z][i]

        coorList = set()
        for z in blockOcc:
            coorList.update (blockSubs[z][blockOcc [z]][0])

        for coor in coorList:
            ## Apply fitch for each coordinate
            currentChar = dict()
            for node in tree.traverse_postorder():
                nodeLabel = node.get_label()
                if (node.is_root()):
                    rootLabel = nodeLabel
                
                if (node.is_leaf()):
                    idxList = blockSubs [nodeLabel][blockOcc[nodeLabel]][0]
                    charList = blockSubs [nodeLabel][blockOcc[nodeLabel]][1]
                    if coor in idxList:
                        charIdx = idxList.index(coor)
                        currentChar[nodeLabel] = [charList[charIdx]]
                    else:
                        currentChar[nodeLabel] = [consSeq[coor]]
                    # print (nodeLabel, blockList [nodeLabel][gIdx][i])

                else:
                    childrenList = node.child_nodes()
                    X = set(currentChar[childrenList[0].get_label()])
                    Y = set(currentChar[childrenList[1].get_label()])
                    if (list (X & Y) == []):
                        currentChar[nodeLabel] = list (X | Y)
                    else:
                        currentChar[nodeLabel] = list (X & Y)

            # print (currentChar)
            for node in tree.traverse_preorder():
                nodeLabel = node.get_label()
                
                if node.is_root():
                    if consSeq[coor] in currentChar:
                        currentChar[nodeLabel] = [consSeq[coor]]
                    else:
                        currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                        if nodeLabel not in gMutationLocal:
                            gMutationLocal[nodeLabel] = [str(counterGlobal) + ":-1:S:" + str(coor) + ":" + currentChar[nodeLabel][0]]
                        else:
                            gMutationLocal[nodeLabel].append(str(counterGlobal) + ":-1:S:" + str(coor) + ":" + currentChar[nodeLabel][0])

                else:
                    if currentChar[node.get_parent().get_label()][0] in currentChar[nodeLabel]:
                        currentChar[nodeLabel] = [currentChar[node.get_parent().get_label()][0]]
                    else:
                        currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                        if nodeLabel not in gMutationLocal:
                            gMutationLocal[nodeLabel] = [str(counterGlobal) + ":-1:S:" + str(coor) + ":" + currentChar[nodeLabel][0]]
                        else:
                            gMutationLocal[nodeLabel].append(str(counterGlobal) + ":-1:S:" + str(coor) + ":" + currentChar[nodeLabel][0])

        counterGlobal += 1        

    # if (i == 4):
    #     print (gMutationLocal)
    #     break






