from treeswift import *
# from ete3 import Tree
from sys import argv
from readJson import *
import json

INF = 10000
VER = 2
HOR = 1
DIAG = 0

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

def blockIDConv (localBlock, globalBlock):
    if localBlock == -1:
        localBlock = 0
    localBlock  = bin(localBlock)[2:]
    localBlock = "0"*(32 - len(localBlock)) + localBlock
    globalBlock  = bin(globalBlock)[2:]
    globalBlock = "0"*(32 - len(globalBlock)) + globalBlock
    # print (localBlock + globalBlock)
    return  globalBlock + localBlock

def nucBinConv (listNuc):
    nucBin = ""
    for nuc in listNuc:
        nucBin = nucBin + char2bin(nuc)
    
    nucBin = "0" * (24 - len(nucBin)) + nucBin
    return nucBin

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

    # print (len(globalBlock))
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
        # print (len (globalBlock), eachSeq)
    
    return (globalBlock, idxList) 

f = open("/home/AD.UCSD.EDU/swalia/data/HIV/pangraph/HIV_5000.json","r")
# f = open("/home/AD.UCSD.EDU/swalia/sars_4k.json","r")
data = json.load(f)
f.close()

## Import Tree
# newick_read = open("/home/AD.UCSD.EDU/swalia/sars_4k.nwk", "r")
newick_read = open("/home/AD.UCSD.EDU/swalia/data/HIV/pangraph/HIV_5000.nwk", "r")
tree_read = newick_read.readlines()[0].split("\n")[0]
# tree_read = newick_read.readlines()[0].split("\n")[0] + ";"
# tree_new = Tree(tree_read)
# print (tree_new)

# for node in tree_new.traverse("preorder"):
#   # Do some analysis on node
#   print (node.name)

tree_old = tree_read
newick_read.close()
tree = read_tree_newick(tree_old)
# print(tree)
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
print ("global blocks generated")
print ("Global block:", gBlocks)
# print ("GIdx:", gIdx)

### Generating Gaps #####
gapListGlobal = []
for i in range (len (gIdx) + 1):
    # print (i)
    gapListLocal = []
    if i < len (gIdx):
        # print ("len of seqName", len(seqName))
        for j in range (len (seqName)):
            currIdx = gIdx[i][j]
            if (i == 0):
                prevIdx = 0
            else:
                prevIdx = gIdx[i-1][j] + 1

            gapListLocal, X, Y = blockAlignment(gapListLocal, seqPaths[seqName[j]][prevIdx: currIdx])
            # print (j, len(gapListLocal))
    else:

        for j in range (len (seqName)):
            currIdx = len(seqPaths[ seqName [j]]) 
            if (i == 0):
                prevIdx = 0
            else:
                prevIdx = gIdx[i-1][j] + 1
            

            gapListLocal, X, Y = blockAlignment(gapListLocal, seqPaths[seqName[j]][prevIdx: currIdx])


    gapListGlobal.append (gapListLocal)

print ("Gaps Generated")
# print ("GapList:", gapListGlobal)

### Appending gaps in each Sequence ###
blockList = dict()
for j in range (len (seqName)):
    
    blockListLocal = []

    for i in range (len (gIdx) + 1):
        if (i < len (gIdx)):
            currIdx = gIdx[i][j]
            if (i == 0):
                prevIdx = 0
            else:
                prevIdx = gIdx[i-1][j] + 1
        else:
            currIdx = len(seqPaths[ seqName [j]]) 
            if (i == 0):
                prevIdx = 0
            else:
                prevIdx = gIdx[i-1][j] + 1
            
        gapListLocal, X, Y = blockAlignment(gapListGlobal[i], seqPaths[seqName[j]][prevIdx: currIdx])

        blockListLocal.append(Y)
        if (i <= (len(gIdx) - 1)):
            blockListLocal.append(gBlocks[i])

    blockList[seqName[j]] = blockListLocal

print ("blockList generated")
# print (blockList["ON650488.1"])
# print (blockList["ON650470.1"])
# print (seqPaths["ON650488.1"])
# print (seqPaths["ON650470.1"])
c = 0
for i in blockList[seqName[0]]:
    if type(i) == list:
        c += len (i)
    else:
        c += 1
# print (c)

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


# print (seqName [0])
# print (blockListIdx ["ON650346.1"])
'''
### Check if Seq length is correct
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

    if (seqLen != seqLenWithoutGap):
        print (seqLen, seqLenWithoutGap, seqLenWithGap)
'''

### Apply Fitch on Gaps ###
gMutation = dict()
for globalIdx in range (len (blockList[seqName[0]])):
    idx = int(globalIdx/2)
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
                            gMutationLocal[nodeLabel] = str(idx) + ":" + str(i) + ":I:" + str(currentChar[nodeLabel][0])
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
                            gMutationLocal[child0Label] = str(idx) + ":" + str(i) + ":D:" + str(currentChar[child0Label][0])
                        else:
                            gMutationLocal[child0Label] = str(idx) + ":" + str(i) + ":I:" + str(currentChar[child0Label][0])
                    else:
                        currentChar[child0Label] = currentChar[nodeLabel]

                    Y = set(currentChar[child1Label])
                    if (list (Y & Z) == []):
                        currentChar[child1Label] = [currentChar[child1Label][0]]
                        if (currentChar[child1Label][0] == "-"):
                            gMutationLocal[child1Label] = str(idx) + ":" + str(i) + ":D:" + str(currentChar[child1Label][0])
                        else:
                            gMutationLocal[child1Label] = str(idx) + ":" + str(i) + ":I:" + str(currentChar[child1Label][0])
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
                    gMutation[nodeLabel].append(str(idx) + ":I:" + blockList[seqName[0]][globalIdx])
                else:
                    gMutation[nodeLabel] = [str(idx) + ":I:" + blockList[seqName[0]][globalIdx]]

print("Fitch on gap")


'''
# Re-generating blocks
positions = list()
gapLengths = list()

print(seqName[0], blockList[seqName[0]])
print (gapListGlobal)
pos = 0
for i in gapListGlobal:
    positions.append (pos)
    gapLengths.append (len(i))
    pos += 1

def copyList(inFile):
    outFile = []
    for eachElem in inFile:
        if type(eachElem) == list:
            localList = []
            for elem in eachElem:
                localList.append(elem)
            outFile.append(localList)
        else:
            outFile.append(eachElem)
    return outFile



blockListNode = dict()
mismatch = []
for eachNode in tree.traverse_preorder():
    nodeLabel = eachNode.get_label()
    if (eachNode.is_root()):
        blockListNode [nodeLabel] = []
        for i in range (len (positions)):
            if (gapLengths[i] > 0):
                blockListNode[nodeLabel].append (["-"]*gapLengths[i])
            if (i < len (positions) - 1):
                blockListNode[nodeLabel].append ("-")
        if nodeLabel in gMutation:
            for eachMut in gMutation[nodeLabel]:
                # print (eachMut)
                eachMut = eachMut.split(":")
                if (len (eachMut) == 4):
                    # Gap mutatation
                    idx = int(int(eachMut[0])*2)
                    gIdx = int (eachMut[1])
                    if (eachMut[2]=="I"):
                        blockListNode[nodeLabel][idx][gIdx] = eachMut[3]
                    else:
                        blockListNode[nodeLabel][idx][gIdx] = "-"

                else:
                    # Non-Gap mutation
                    idx = int(int(eachMut[0])*2 + 1)
                    if (eachMut[1]=="I"):
                        blockListNode[nodeLabel][idx] = eachMut[2]
                    else:
                        blockListNode[nodeLabel][idx] = "-"

        # print (nodeLabel, blockListNode[nodeLabel])
        # break
    else:
        parentLabel = eachNode.get_parent().get_label()
        blockListNode [nodeLabel] = copyList (blockListNode [parentLabel])
        # if (nodeLabel == 787):
        #     print ("787 parent:",parentLabel, blockListNode [789])
        #     break
        if nodeLabel in gMutation:
            for eachMut in gMutation[nodeLabel]:
                # print (eachMut)
                eachMut = eachMut.split(":")
                if (len (eachMut) == 4):
                    # Gap mutatation
                    idx = int(int(eachMut[0])*2)
                    gIdx = int (eachMut[1])
                    if (eachMut[2]=="I"):
                        blockListNode[nodeLabel][idx][gIdx] = eachMut[3]
                    else:
                        blockListNode[nodeLabel][idx][gIdx] = "-"

                else:
                    # Non-Gap mutation
                    idx = int(int(eachMut[0])*2 + 1)
                    if (eachMut[1]=="I"):
                        blockListNode[nodeLabel][idx] = eachMut[2]
                    else:
                        blockListNode[nodeLabel][idx] = "-"
    # print (nodeLabel, blockListNode[nodeLabel])

    if (eachNode.is_leaf()):
        regenBlockList = []
        for eachElem in blockListNode[nodeLabel]:
            if (type(eachElem) == list):
                for elem in eachElem:
                    if (elem != "-"):
                        regenBlockList.append(elem)
            else:
                regenBlockList.append(eachElem)
        if (regenBlockList != seqPaths[nodeLabel]):
            mismatch.append(nodeLabel)
            if (nodeLabel == "ON655904.1"):
                print (nodeLabel, regenBlockList, seqPaths[nodeLabel])

print (len (mismatch))

# exit(1)


### For a node
nodeToCheck = "ON650453.1"
for eachNode in tree.traverse_preorder():
    nodeLabel = eachNode.get_label()
    if (nodeLabel == nodeToCheck):
        node = eachNode
        break

while (node.is_root()==False):
    nodeLabel = node.get_label()
    # if (nodeLabel) in gMutation:
    print (nodeLabel)
    node = node.get_parent()
print (node.get_label(), gMutation[node.get_label()])
'''
blockData = getBlocks (data)
checkId = "ON651176.1"
# checkId = "ON651157.1"
checkBlock = "WYWTUMLSDN"
# print (seqPaths[checkId])
# blockInfo = getBlockInfo (checkBlock, blockData)
# # consSeq   = getBlockSeq (blockInfo)
# blockSubs = getBlockSubs (blockInfo)
# print (blockSubs[checkId])
# blockDel = getBlockDel (blockInfo)
# print (blockDel[checkId])
# blockIns = getBlockIns (blockInfo)
# print (blockIns[checkId])

#####################


# exit(1)
### Apply fitch on Substitutions ###
blockData = getBlocks (data)
blockMap = dict()
gSubs = dict()
gDel  = dict()
gIns  = dict()
gGaps = dict()
SUBS = True
DEL  = True
INS  = True



for node in tree.traverse_preorder():
    nodeLabel = node.get_label()
    gSubs[nodeLabel] = []

counterGlobal = 0

if SUBS == True & DEL == True:
    for i in range (len (blockList [seqName[0]])):
        counterLocal = 0
        # elem - just to check if it's a list or not 
        elem = blockList[seqName [0]][i]

        if type(elem) == list:
            # print ("It's a list", i)
            for y in range (len (elem)):
                element = gapListGlobal[int(i/2)][y]
                # print (element)
                blockInfo = getBlockInfo (element, blockData)
                consSeq   = getBlockSeq (blockInfo)
                blockSubs = getBlockSubs (blockInfo)
                blockDel = getBlockDel (blockInfo)

                # print (blockSubs)
                blockOcc = dict()
                for z in blockList:
                    if blockListIdx [z][i][y] != "-":
                        blockOcc [z] = blockListIdx [z][i][y]
                
                # print (len(blockOcc))
                # break
                coorList = set()
                coorDelList = set()
                for z in blockOcc:
                    coorList.update (blockSubs[z][blockOcc [z]][0])
                    coorDelList.update (blockDel[z][blockOcc [z]])
                
                # print (coorList, len(consSeq))
                for coor in coorList:
                    ## Apply fitch for each coordinate
                    coorDelPresent = False
                    currentChar = dict()
                    for node in tree.traverse_postorder():
                        nodeLabel = node.get_label()
                        if (node.is_root()):
                            rootLabel = nodeLabel
                        
                        if (node.is_leaf()):

                            if nodeLabel not in blockOcc:
                                currentChar[nodeLabel] = [consSeq[coor - 1]]
                            else:
                                idxList  = blockSubs [nodeLabel][blockOcc[nodeLabel]][0]
                                charList = blockSubs [nodeLabel][blockOcc[nodeLabel]][1]
                                idxDelList  = blockDel [nodeLabel][blockOcc[nodeLabel]]
                                if coor in idxList:
                                    charIdx = idxList.index(coor)
                                    currentChar[nodeLabel] = [charList[charIdx]]
                                elif coor in idxDelList:
                                    currentChar[nodeLabel] = ["-"]
                                    coorDelPresent = True
                                else:
                                    # print (coor, len(consSeq))
                                    currentChar[nodeLabel] = [consSeq[coor - 1]]
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
                            if consSeq[coor - 1] in currentChar[nodeLabel]:
                                currentChar[nodeLabel] = [consSeq[coor - 1]]
                            else:
                                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                                if currentChar[nodeLabel][0] == "-":
                                    if nodeLabel not in gSubs:
                                        gSubs[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor)]                                
                                    else:
                                        gSubs[nodeLabel].append (str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor))
                                    
                                else:
                                    if nodeLabel not in gSubs:
                                        gSubs[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor) + ":" + currentChar[nodeLabel][0]]
                                    else:
                                        gSubs[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor) + ":" + currentChar[nodeLabel][0])

                        else:
                            parentChar = currentChar[node.get_parent().get_label()][0]
                            if parentChar in currentChar[nodeLabel]:
                                currentChar[nodeLabel] = [parentChar]
                            else:
                                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                                ## Insertion
                                if parentChar == "-":
                                    if nodeLabel not in gSubs:
                                        gSubs[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":I:" + str(coor) + ":" + currentChar[nodeLabel][0]]                                
                                    else:
                                        gSubs[nodeLabel].append (str (counterGlobal) + ":" + str (counterLocal) + ":I:" + str(coor) + ":" + currentChar[nodeLabel][0])

                                else:
                                    ## Deletion
                                    if currentChar[nodeLabel][0] == "-":
                                        if nodeLabel not in gSubs:
                                            gSubs[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor)]                                
                                        else:
                                            gSubs[nodeLabel].append (str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor))
                                    
                                    ## Substitution
                                    else:
                                        if nodeLabel not in gSubs:
                                            gSubs[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor) + ":" + currentChar[nodeLabel][0]]                                
                                        else:
                                            gSubs[nodeLabel].append (str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor) + ":" + currentChar[nodeLabel][0])

                    if coorDelPresent:
                        coorDelList.remove(coor)

                for coor in coorDelList:
                    ## Apply fitch for each coordinate
                    currentChar = dict()
                    for node in tree.traverse_postorder():
                        nodeLabel = node.get_label()
                        if (node.is_root()):
                            rootLabel = nodeLabel
                        
                        if (node.is_leaf()):
                            if nodeLabel not in blockOcc:
                                currentChar[nodeLabel] = [consSeq[coor - 1]]
                            else:
                                idxList  = blockDel [nodeLabel][blockOcc[nodeLabel]]
                                if coor in idxList:
                                    currentChar[nodeLabel] = ["-"]
                                else:
                                    currentChar[nodeLabel] = [consSeq[coor - 1]]

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
                            if consSeq[coor - 1] in currentChar[nodeLabel]:
                                currentChar[nodeLabel] = [consSeq[coor - 1]]
                            else:
                                currentChar[nodeLabel] = ["-"]
                                if nodeLabel not in gDel:
                                    gDel[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor)]
                                else:
                                    gDel[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor))

                            
                        else:
                            parentChar = currentChar[node.get_parent().get_label()][0]
                            if parentChar in currentChar[nodeLabel]:
                                currentChar[nodeLabel] = [parentChar]
                            else:
                                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                                ## Insertion
                                if parentChar == "-":
                                    if nodeLabel not in gDel:
                                        gDel[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":I:" + str(coor) + ":" + currentChar[nodeLabel][0]]
                                    else:
                                        gDel[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":I:" + str(coor) + ":" + currentChar[nodeLabel][0])

                                ## Deletion
                                else:
                                    if nodeLabel not in gDel:
                                        gDel[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor)]
                                    else:
                                        gDel[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor))
                
                counterLocal += 1

        else:
            # print ("It's not a list", i)

            blockInfo = getBlockInfo (elem, blockData)
            consSeq   = getBlockSeq (blockInfo)
            blockSubs = getBlockSubs (blockInfo)
            blockDel = getBlockDel (blockInfo)
            
            blockOcc = dict()
            for z in blockList:
                blockOcc [z] = blockListIdx [z][i]

            coorList = set()
            coorDelList = set()
            for z in blockOcc:
                coorList.update (blockSubs[z][blockOcc [z]][0])
                coorDelList.update (blockDel[z][blockOcc [z]])

            for coor in coorList:
                ## Apply fitch for each coordinate
                coorDelPresent = False
                currentChar = dict()
                
                for node in tree.traverse_postorder():
                    nodeLabel = node.get_label()
                    if (node.is_root()):
                        rootLabel = nodeLabel
                    
                    if (node.is_leaf()):
                        idxList = blockSubs [nodeLabel][blockOcc[nodeLabel]][0]
                        charList = blockSubs [nodeLabel][blockOcc[nodeLabel]][1]
                        idxDelList  = blockDel [nodeLabel][blockOcc[nodeLabel]]
                        if coor in idxList:
                            charIdx = idxList.index(coor)
                            currentChar[nodeLabel] = [charList[charIdx]]
                        elif coor in idxDelList:
                            currentChar[nodeLabel] = ["-"]
                            coorDelPresent = True
                        else:
                            currentChar[nodeLabel] = [consSeq[coor - 1]]
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
                        if consSeq[coor - 1] in currentChar[nodeLabel]:
                            currentChar[nodeLabel] = [consSeq[coor - 1]]
                        else:
                            currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                            if currentChar[nodeLabel][0] == "-":
                                if nodeLabel not in gSubs:
                                    gSubs[nodeLabel] = [str (counterGlobal) + ":-1:D:" + str(coor)]                                
                                else:
                                    gSubs[nodeLabel].append (str (counterGlobal) + ":-1:D:" + str(coor))
                                
                            else:
                                if nodeLabel not in gSubs:
                                    gSubs[nodeLabel] = [str (counterGlobal) + ":-1:S:" + str(coor) + ":" + currentChar[nodeLabel][0]]
                                else:
                                    gSubs[nodeLabel].append(str (counterGlobal) + ":-1:S:" + str(coor) + ":" + currentChar[nodeLabel][0])

                    else:
                        parentChar = currentChar[node.get_parent().get_label()][0]
                        if parentChar in currentChar[nodeLabel]:
                            currentChar[nodeLabel] = [parentChar]
                        else:
                            currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                            ## Insertion
                            if parentChar == "-":
                                if nodeLabel not in gSubs:
                                    gSubs[nodeLabel] = [str (counterGlobal) + ":-1:I:" + str(coor) + ":" + currentChar[nodeLabel][0]]                                
                                else:
                                    gSubs[nodeLabel].append (str (counterGlobal) + ":-1:I:" + str(coor) + ":" + currentChar[nodeLabel][0])

                            else:
                                ## Deletion
                                if currentChar[nodeLabel][0] == "-":
                                    if nodeLabel not in gSubs:
                                        gSubs[nodeLabel] = [str (counterGlobal) + ":-1:D:" + str(coor)]                                
                                    else:
                                        gSubs[nodeLabel].append (str (counterGlobal) + ":-1:D:" + str(coor))
                                
                                ## Substitution
                                else:
                                    if nodeLabel not in gSubs:
                                        gSubs[nodeLabel] = [str (counterGlobal) + ":-1:S:" + str(coor) + ":" + currentChar[nodeLabel][0]]                                
                                    else:
                                        gSubs[nodeLabel].append (str (counterGlobal) + ":-1:S:" + str(coor) + ":" + currentChar[nodeLabel][0])
                if coorDelPresent:
                    coorDelList.remove(coor)
            
            for coor in coorDelList:
                ## Apply fitch for each coordinate
                currentChar = dict()
                for node in tree.traverse_postorder():
                    nodeLabel = node.get_label()
                    if (node.is_root()):
                        rootLabel = nodeLabel
                    
                    if (node.is_leaf()):
                        idxList = blockDel [nodeLabel][blockOcc[nodeLabel]]
                        if coor in idxList:
                            currentChar[nodeLabel] = ["-"]
                        else:
                            currentChar[nodeLabel] = [consSeq[coor - 1]]
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
                        if consSeq[coor - 1] in currentChar[nodeLabel]:
                            currentChar[nodeLabel] = [consSeq[coor - 1]]
                        else:
                            currentChar[nodeLabel] = ["-"]
                            if nodeLabel not in gDel:
                                gDel[nodeLabel] = [str(counterGlobal) + ":-1:D:" + str(coor)]
                            else:
                                gDel[nodeLabel].append(str(counterGlobal) + ":-1:D:" + str(coor))

                    else:
                        parentChar = currentChar[node.get_parent().get_label()][0]
                        if parentChar in currentChar[nodeLabel]:
                            currentChar[nodeLabel] = [parentChar]
                        else:
                            currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                            ## Insertion
                            if parentChar == "-":
                                if nodeLabel not in gDel:
                                    gDel[nodeLabel] = [str(counterGlobal) + ":-1:I:" + str(coor) + ":" + currentChar[nodeLabel][0]]
                                else:
                                    gDel[nodeLabel].append(str(counterGlobal) + ":-1:I:" + str(coor) + ":" + currentChar[nodeLabel][0])

                            ## Deletion
                            else:
                                if nodeLabel not in gDel:
                                    gDel[nodeLabel] = [str(counterGlobal) + ":-1:D:" + str(coor)]
                                else:
                                    gDel[nodeLabel].append(str(counterGlobal) + ":-1:D:" + str(coor))
            
            
            counterGlobal += 1        


print("Fitch on substitution and Deletion")
counterGlobal = 0


if INS == True:
    for i in range (len (blockList [seqName[0]])):
        counterLocal = 0
        # elem - just to check if it's a list or not 
        elem = blockList[seqName [0]][i]

        if type(elem) == list:
            # print ("It's a list", i)
            for y in range (len (elem)):
                element = gapListGlobal[int(i/2)][y]

                blockInfo = getBlockInfo (element, blockData)
                consSeq   = getBlockSeq (blockInfo)
                blockIns = getBlockIns (blockInfo)
                blockGaps = getBlockGap (blockInfo)
                
                for eachGap in blockGaps:
                    if counterGlobal not in gGaps:
                        gGaps[counterGlobal] = dict()
                    
                    if counterLocal not in gGaps[counterGlobal]:
                        gGaps[counterGlobal][counterLocal] = []

                    gGaps[counterGlobal][counterLocal].append([int(eachGap) , blockGaps[eachGap]])



                blockOcc = dict()
                for z in blockList:
                    if blockListIdx [z][i][y] != "-":
                        blockOcc [z] = blockListIdx [z][i][y]

                gapCoorList = set()
                for z in blockOcc:
                    for l in range (len (blockIns[z][blockOcc [z]][0])):
                        gapCoorList.add ((blockIns[z][blockOcc [z]][0][l][0], blockIns[z][blockOcc [z]][0][l][1]))

            
                for coor in gapCoorList:
                    ## Apply fitch for each coordinate
                    currentChar = dict()
                    for node in tree.traverse_postorder():
                        nodeLabel = node.get_label()
                        if (node.is_root()):
                            rootLabel = nodeLabel
                        
                        if (node.is_leaf()):
                            if nodeLabel not in blockOcc:
                                currentChar[nodeLabel] = ["-"]
                            else:
                                idxList = blockIns [nodeLabel][blockOcc[nodeLabel]][0]
                                # print (nodeLabel ,coor, idxList)
                                if [coor[0], coor[1]] in idxList:
                                    charIdx = idxList.index([coor[0], coor[1]])
                                    currentChar[nodeLabel] = [blockIns [nodeLabel][blockOcc[nodeLabel]][1][charIdx]]
                                else:
                                    currentChar[nodeLabel] = ["-"]

                        else:
                            childrenList = node.child_nodes()
                            X = set(currentChar[childrenList[0].get_label()])
                            Y = set(currentChar[childrenList[1].get_label()])
                            if (list (X & Y) == []):
                                currentChar[nodeLabel] = list (X | Y)
                            else:
                                currentChar[nodeLabel] = list (X & Y)

            #         # print (currentChar)
                    for node in tree.traverse_preorder():
                        nodeLabel = node.get_label()
                        
                        if node.is_root():
                            if "-" in currentChar[nodeLabel]:
                                currentChar[nodeLabel] = ["-"]
                            else:
                                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                                if nodeLabel not in gIns:
                                    gIns[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":I:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0]]
                                else:
                                    gIns[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":I:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0])
                            
                        else:
                            parentChar = currentChar[node.get_parent().get_label()][0]
                            if parentChar in currentChar[nodeLabel]:
                                currentChar[nodeLabel] = [parentChar]
                            else:
                                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                                ## Insertion
                                if parentChar == "-":
                                    if nodeLabel not in gIns:
                                        gIns[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":I:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0]]
                                    else:
                                        gIns[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":I:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0])

                                else:

                                    ## Deletion
                                    if currentChar[nodeLabel][0] == "-":

                                        if nodeLabel not in gIns:
                                            gIns[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor[0]) + ":" + str(coor[1])]
                                        else:
                                            gIns[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":D:" + str(coor[0]) + ":" + str(coor[1]))

                                    ## Substitution
                                    else:
                                        if nodeLabel not in gIns:
                                            gIns[nodeLabel] = [str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0]]
                                        else:
                                            gIns[nodeLabel].append(str (counterGlobal) + ":" + str (counterLocal) + ":S:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0])

                counterLocal += 1

        else:
            # print ("It's not a list", i)

            blockInfo = getBlockInfo (elem, blockData)
            consSeq   = getBlockSeq (blockInfo)
            # print (consSeq)
            blockIns = getBlockIns (blockInfo)
            blockGaps = getBlockGap (blockInfo)

            for eachGap in blockGaps:
                if counterGlobal not in gGaps:
                    gGaps[counterGlobal] = dict()
                
                if "-1" not in gGaps[counterGlobal]:
                    gGaps[counterGlobal]["-1"] = []

                gGaps[counterGlobal]["-1"].append([int(eachGap) , blockGaps[eachGap]])


            blockOcc = dict()
            for z in blockList:
                blockOcc [z] = blockListIdx [z][i]

            gapCoorList = set()
            for z in blockOcc:
                for l in range (len (blockIns[z][blockOcc [z]][0])):
                    gapCoorList.add ((blockIns[z][blockOcc [z]][0][l][0], blockIns[z][blockOcc [z]][0][l][1]))


            for coor in gapCoorList:
                ## Apply fitch for each coordinate
                currentChar = dict()
                for node in tree.traverse_postorder():
                    nodeLabel = node.get_label()
                    if (node.is_root()):
                        rootLabel = nodeLabel
                    
                    if (node.is_leaf()):
                        idxList = blockIns [nodeLabel][blockOcc[nodeLabel]][0]
                        # print (nodeLabel ,coor, idxList)
                        if [coor[0], coor[1]] in idxList:
                            charIdx = idxList.index([coor[0], coor[1]])
                            currentChar[nodeLabel] = [blockIns [nodeLabel][blockOcc[nodeLabel]][1][charIdx]]
                        else:
                            currentChar[nodeLabel] = ["-"]

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
                        if "-" in currentChar[nodeLabel]:
                            currentChar[nodeLabel] = ["-"]
                        else:
                            currentChar[nodeLabel] = currentChar[nodeLabel][0]
                            if nodeLabel not in gIns:
                                gIns[nodeLabel] = [str(counterGlobal) + ":-1:I:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0]]
                            else:
                                gIns[nodeLabel].append(str(counterGlobal) + ":-1:I:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0])
            
                    else:
                        parentChar = currentChar[node.get_parent().get_label()][0]
                        if parentChar in currentChar[nodeLabel]:
                            currentChar[nodeLabel] = [parentChar]
                        else:
                            currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                            ## Insertion
                            if parentChar == "-":
                                if nodeLabel not in gIns:
                                    gIns[nodeLabel] = [str(counterGlobal) + ":-1:I:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0]]
                                else:
                                    gIns[nodeLabel].append(str(counterGlobal) + ":-1:I:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0])

                            else:

                                ## Deletion
                                if currentChar[nodeLabel][0] == "-":

                                    if nodeLabel not in gIns:
                                        gIns[nodeLabel] = [str(counterGlobal) + ":-1:D:" + str(coor[0]) + ":" + str(coor[1])]
                                    else:
                                        gIns[nodeLabel].append(str(counterGlobal) + ":-1:D:" + str(coor[0]) + ":" + str(coor[1]))

                                ## Substitution
                                else:
                                    if nodeLabel not in gIns:
                                        gIns[nodeLabel] = [str(counterGlobal) + ":-1:S:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0]]
                                    else:
                                        gIns[nodeLabel].append(str(counterGlobal) + ":-1:S:" + str(coor[0]) + ":" + str(coor[1]) + ":" + currentChar[nodeLabel][0])



            counterGlobal += 1

print("Fitch on insertion")
# for node in tree.traverse_preorder():
#     if (node.get_label() == "ON650464.1"):
#         nodeQuery = node

# while nodeQuery.is_root() == False:
#     if nodeQuery.get_label() in gIns:
#         print (gIns[nodeQuery.get_label()])
#     nodeQuery = nodeQuery.get_parent()    
# if nodeQuery.get_label() in gIns:
#     print (gIns[nodeQuery.get_label()])


'''
# Regenrating the sequence
seqLabel = "ON650466.1"
listBlocks = []

for node in tree.traverse_preorder():

    nodeLabel = node.get_label()

    if node.is_root():
        root = node
        rootLabel = nodeLabel

        for eachGap in gapListGlobal:
            listBlocks.append( ["-"] * len(eachGap) )
            listBlocks.append("-")

        listBlocks = listBlocks[:-1]

        for eachMut in gMutation [rootLabel]:
            mut = eachMut.split(":")
            if len (mut) == 3:
                listBlocks[int (mut[0])] = mut[2]
            else:
                listBlocks[int (mut[0])] [int (mut[1])] = mut[3]

    else:
        if nodeLabel == seqLabel:
            currNode = node
            break

listNode = []
while currNode.get_label() != rootLabel:
    listNode.append(currNode.get_label())
    currNode = currNode.get_parent()
listNode.reverse()

for currNodeLable in listNode:
    if currNodeLable in gMutation:
        for eachMut in gMutation [currNodeLable]:
            mut = eachMut.split(":")
            if (mut[2] == "I"):
                listBlocks[int (mut[0])] [int (mut[1])] = mut[3]
            else:
                listBlocks[int (mut[0])] [int (mut[1])] = "-"

# print (listBlocks)
print (gGaps)

for i in range (len(listBlocks)):
    if (type(listBlocks[i]) == list):
        continue
    else:
        # list all the mutations
        element = listBlocks[i]
            # print (element)
        blockInfo = getBlockInfo (element, blockData)
        consSeq   = getBlockSeq (blockInfo)
        
        for eachSubs in gSubs[seqLabel]:
            eachSubs = eachSubs.split(":")
            if ((eachSubs[0] == str(int(i/2))) and (eachSubs[1] == '-1')):
                sliceIdx = int(eachSubs[3])
                if (sliceIdx == len(consSeq) - 1):
                    consSeq = consSeq[:int(eachSubs[3])] + eachSubs[4] 
                else:
                    consSeq = consSeq[:int(eachSubs[3])] + eachSubs[4] + consSeq[int(eachSubs[3])+1:]

        print (gDel)

'''

import mutation_annotation_pb2
from copy import deepcopy
import time
import sys


MAT = mutation_annotation_pb2.tree()

# Read the address book.
try:
  f = open(sys.argv[1], "rb")
  MAT.ParseFromString(f.read())
  f.close()
except IOError:
  print (sys.argv[1] + ": Could not open MAT.  Creating a new one.")


MAT.newick = tree_old

gapList = MAT.blockGaps

positions = list()
gapLengths = list()

pos = 0
for i in gapListGlobal:
    positions.append (pos)
    gapLengths.append (len(i))
    pos += 1

gapList.blockPosition.extend (positions)
gapList.blockGapLength.extend (gapLengths)

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
    # print ((int_list[-1]))
    return (int_list)

counterGlobal = 0

for i in range (len (blockList [seqName[0]])):
    counterLocal = 0
    # elem - just to check if it's a list or not 
    elem = blockList[seqName [0]][i]

    if type(elem) == list:
        for y in range (len (elem)):
            element = gapListGlobal[int(i/2)][y]
            blockInfo = getBlockInfo (element, blockData)
            consSeq   = getBlockSeq (blockInfo)

            newBlock = MAT.blocks.add()
            blockId  = blockIDConv(int(counterLocal), int(counterGlobal))
            newBlock.blockId = int(blockId, 2)
            newBlock.blockGapExist = True
            newBlock.consensusSeq.extend(seq2int(consSeq))
            # Check Open
            list_ = seq2int(consSeq)
            for num in list_:
                num = bin(num)[2:]
                num = "0" * (32 - len(num)) + num
                for zz in range (8):
                    if (num[4 * zz: 4 * (zz + 1)] == "1010"):
                        print ("found:", list_)
            # Check Close
            counterLocal += 1

    else:
        # print ("It's not a list", i)

        blockInfo = getBlockInfo (elem, blockData)
        consSeq   = getBlockSeq (blockInfo)
        
        newBlock = MAT.blocks.add()
        blockId  = blockIDConv(int(counterLocal), int(counterGlobal))
        newBlock.blockId = int(blockId, 2)
        newBlock.blockGapExist = False
        newBlock.consensusSeq.extend(seq2int(consSeq))

        list_ = seq2int(consSeq)
        # Check Open
        list_ = seq2int(consSeq)
        for num in list_:
            num = bin(num)[2:]
            num = "0" * (32 - len(num)) + num
            for zz in range (8):
                if (num[4 * zz: 4 * (zz + 1)] == "1010"):
                    print ("found:", list_)
            # Check Close

        counterGlobal += 1

# blockMap = dict()
# for i in range (len (blockList [seqName[0]])):
#     counterLocal = 0
#     # elem - just to check if it's a list or not 
#     elem = blockList[seqName [0]][i]

#     if type(elem) == list:
#         # print ("It's a list", i)
#         for y in range (len (elem)):
#             element = gapListGlobal[int(i/2)][y]
#             if element not in blockMap:
#                 blockMap[element] = [[],[]]

#             blockId  = blockIDConv(int(counterLocal), int(counterGlobal))
#             blockMap[element][0].append (blockId)
#             blockMap[element][1].append (True)
#             counterLocal += 1

#     else:
#         # print ("It's not a list", i)
#         if elem not in blockMap:
#             blockMap[elem] = [[],[]]
#         blockId  = blockIDConv(int(counterLocal), int(counterGlobal))
#         blockMap[elem][0].append(blockId)
#         blockMap[elem][1].append(False)
#         counterGlobal += 1 
# # print(blockMap)
# for elem in blockMap:
#     blockInfo = getBlockInfo (elem, blockData)
#     consSeq   = getBlockSeq (blockInfo)

#     newBlock = MAT.blocks.add()
#     newBlock.consensusSeq.extend(seq2int(consSeq))
#     for i in range (len (blockMap[elem][0])):
#         newBlockElem = newBlock.blocks.add()
#         newBlockElem.blockId = int(blockMap[elem][0][i], 2)
#         newBlockElem.blockGapExist = blockMap[elem][1][i]


# for i in range (len (blockList [seqName[0]])):
#     counterLocal = 0
#     # elem - just to check if it's a list or not 
#     elem = blockList[seqName [0]][i]

#     if type(elem) == list:
#         # print ("It's a list", i)
#         for y in range (len (elem)):
#             element = gapListGlobal[int(i/2)][y]
#             # print (element)
#             blockInfo = getBlockInfo (element, blockData)
#             consSeq   = getBlockSeq (blockInfo)

#             newBlock = MAT.blocks.add()
#             blockId  = blockIDConv(int(counterLocal), int(counterGlobal))
#             # print (int(counterLocal), int(counterGlobal), blockId)
#             newBlock.blockId = int(blockId, 2)
#             newBlock.blockGapExist = True
#             newBlock.consensusSeq.extend(seq2int(consSeq))

#             counterLocal += 1

#     else:
#         # print ("It's not a list", i)

#         blockInfo = getBlockInfo (elem, blockData)
#         consSeq   = getBlockSeq (blockInfo)
        
#         newBlock = MAT.blocks.add()
#         blockId  = blockIDConv(int(counterLocal), int(counterGlobal))
#         newBlock.blockId = int(blockId, 2)
#         newBlock.blockGapExist = False
#         newBlock.consensusSeq.extend(seq2int(consSeq))

#         counterGlobal += 1        

# print("MAT: blocks generated")

SsubsDict = {}
SinsDict = {}
SdelDict = {}
delDict = {}
insDict = {}
IsubsDict = {}
IdelDict = {}
IinsDict = {}

# for eachNode in gSubs:
for node in tree.traverse_preorder():
    eachNode = node.get_label()
    # print (eachNode)
    SsubsDict[eachNode] = {}
    SinsDict[eachNode] = {}
    SdelDict[eachNode] = {}
    delDict[eachNode] = {}
    insDict[eachNode] = {}
    IsubsDict[eachNode] = {}
    IdelDict[eachNode] = {}
    IinsDict[eachNode] = {}
    node = MAT.nodes.add()

    # block mutations
    if eachNode in gMutation:
        for eachBlockMut in gMutation[eachNode]:
            blockMut = node.blockMutation.add()
            eachBlockMut = eachBlockMut.split(":")
            
            if (len(eachBlockMut) == 3):
                counterLocal = "-1"
                counterGlobal = eachBlockMut[0]
            else:
                counterLocal = eachBlockMut[1]
                counterGlobal = eachBlockMut[0]
            blockId  = blockIDConv(int(counterLocal), int(counterGlobal))
            blockMut.blockId = int(blockId, 2)
            if (counterLocal == "-1"):
                blockMut.blockGapExist = False
            else:
                blockMut.blockGapExist = True

            if (eachBlockMut[-2] == "I"):
                blockMut.blockMutInfo = True
            else:
                blockMut.blockMutInfo = False
            
            # print (eachNode, eachBlockMut, int(counterLocal), int(counterGlobal), blockMut.blockId, blockMut.blockGapExist, blockMut.blockMutInfo)
            

    # print (eachNode)
    if eachNode in gSubs:
        for eachSubs in gSubs[eachNode]:
            subsInfo = eachSubs.split(":")

            if (subsInfo[2] == "S"):
                if subsInfo[0] not in SsubsDict[eachNode]:
                    SsubsDict[eachNode][subsInfo[0]] = {}
                
                if subsInfo[1] not in SsubsDict[eachNode][subsInfo[0]]:
                    SsubsDict[eachNode][subsInfo[0]][subsInfo[1]] = [[],[]]
                SsubsDict[eachNode][subsInfo[0]][subsInfo[1]][0].append(int(subsInfo[3]))
                SsubsDict[eachNode][subsInfo[0]][subsInfo[1]][1].append(subsInfo[4])

            elif (subsInfo[2] == "D"):
                if subsInfo[0] not in SdelDict[eachNode]:
                    SdelDict[eachNode][subsInfo[0]] = {}
                
                if subsInfo[1] not in SdelDict[eachNode][subsInfo[0]]:
                    SdelDict[eachNode][subsInfo[0]][subsInfo[1]] = []
                
                SdelDict[eachNode][subsInfo[0]][subsInfo[1]].append(int(subsInfo[3]))

            else:
                if subsInfo[0] not in SinsDict[eachNode]:
                    SinsDict[eachNode][subsInfo[0]] = {}
                
                if subsInfo[1] not in SinsDict[eachNode][subsInfo[0]]:
                    SinsDict[eachNode][subsInfo[0]][subsInfo[1]] = [[],[]]
                SinsDict[eachNode][subsInfo[0]][subsInfo[1]][0].append(int(subsInfo[3]))
                SinsDict[eachNode][subsInfo[0]][subsInfo[1]][1].append(subsInfo[4])
    
    if eachNode in gDel:
        for eachDel in gDel[eachNode]:
            delInfo = eachDel.split(":")

            if delInfo[2] == "D":
                if delInfo[0] not in delDict[eachNode]:
                    delDict[eachNode][delInfo[0]] = {}
                
                if delInfo[1] not in delDict[eachNode][delInfo[0]]:
                    delDict[eachNode][delInfo[0]][delInfo[1]] = []

                delDict[eachNode][delInfo[0]][delInfo[1]].append(int(delInfo[3]))

            else:
                if delInfo[0] not in insDict[eachNode]:
                    insDict[eachNode][delInfo[0]] = {}
                
                if delInfo[1] not in insDict[eachNode][delInfo[0]]:
                    insDict[eachNode][delInfo[0]][delInfo[1]] = [[],[]]

                insDict[eachNode][delInfo[0]][delInfo[1]][0].append(int(delInfo[3]))
                insDict[eachNode][delInfo[0]][delInfo[1]][1].append(delInfo[4])
        # print ("Deletion dict generated")

    if eachNode in gIns:
        for eachIns in gIns[eachNode]:
            insInfo = eachIns.split(":")
            
            if insInfo[2] == "D":
                if insInfo[0] not in IdelDict[eachNode]:
                    IdelDict[eachNode][insInfo[0]] = {}
                
                if insInfo[1] not in IdelDict[eachNode][insInfo[0]]:
                    IdelDict[eachNode][insInfo[0]][insInfo[1]] = {}

                if insInfo[3] not in IdelDict[eachNode][insInfo[0]][insInfo[1]]:
                    IdelDict[eachNode][insInfo[0]][insInfo[1]][insInfo[3]] = []
                
                
                IdelDict[eachNode][insInfo[0]][insInfo[1]][insInfo[3]].append(int(insInfo[4]))

            elif insInfo[2] == "I":
                if insInfo[0] not in IinsDict[eachNode]:
                    IinsDict[eachNode][insInfo[0]] = {}
                
                if insInfo[1] not in IinsDict[eachNode][insInfo[0]]:
                    IinsDict[eachNode][insInfo[0]][insInfo[1]] = {}

                if insInfo[3] not in IinsDict[eachNode][insInfo[0]][insInfo[1]]:
                    IinsDict[eachNode][insInfo[0]][insInfo[1]][insInfo[3]] = [[],[]]
                
                # print (insInfo)
                IinsDict[eachNode][insInfo[0]][insInfo[1]][insInfo[3]][0].append(int(insInfo[4]))
                IinsDict[eachNode][insInfo[0]][insInfo[1]][insInfo[3]][1].append(insInfo[5])

            else:
                if insInfo[0] not in IsubsDict[eachNode]:
                    IsubsDict[eachNode][insInfo[0]] = {}
                
                if insInfo[1] not in IsubsDict[eachNode][insInfo[0]]:
                    IsubsDict[eachNode][insInfo[0]][insInfo[1]] = {}

                if insInfo[3] not in IsubsDict[eachNode][insInfo[0]][insInfo[1]]:
                    IsubsDict[eachNode][insInfo[0]][insInfo[1]][insInfo[3]] = [[],[]]
                
                
                IsubsDict[eachNode][insInfo[0]][insInfo[1]][insInfo[3]][0].append(int(insInfo[4]))
                IsubsDict[eachNode][insInfo[0]][insInfo[1]][insInfo[3]][1].append(insInfo[5])

        # print ("Insertion dict generated")

    
    ## Substitutions -> Substitutions, Insertions and Deletions 
    for eachPos in SsubsDict[eachNode]:
        for eachGapPos in SsubsDict[eachNode][eachPos]:
            
            loc = deepcopy(SsubsDict[eachNode][eachPos][eachGapPos][0])
            loc.sort()
            # print (loc)
            i = 0
            while i < len(loc):
                nodeMut = node.nucMutation.add()
                currPos = loc[i]
                currIdx = loc[i]
                length = 1
                idxUnsort = SsubsDict[eachNode][eachPos][eachGapPos][0].index(currIdx)
                charList = [SsubsDict[eachNode][eachPos][eachGapPos][1][idxUnsort]]
                
                while True:
                    if (i+1) == len(loc):
                        # print (currPos, length, charList)
                        break
                
                    if currIdx == (loc[i+1] - 1):
                        currIdx = loc[i+1]
                        length += 1
                        i += 1
                        idxUnsort = SsubsDict[eachNode][eachPos][eachGapPos][0].index(currIdx)
                        charList.append(SsubsDict[eachNode][eachPos][eachGapPos][1][idxUnsort])
                        if (length >= 6):
                            break
                    else:
                        # print (currPos, length, charList)
                        break
                nodeMut.nucPosition = currPos - 1
                nodeMut.nucGapExist = False
                blockId  = blockIDConv(int(eachGapPos), int(eachPos))
                nodeMut.blockId = int(blockId, 2)
                if (eachGapPos == "-1"):
                    nodeMut.blockGapExist = False
                else:
                    nodeMut.blockGapExist = True
                lenBin = bin(length)[2:]
                lenMut = "0" * (4 - len (lenBin)) + lenBin 
                nucBin = nucBinConv (charList)
                # print (lenMut, nucBin)
                nodeMut.mutInfo = int( nucBin + lenMut  + "0000" ,2)

                if (nodeMut.mutInfo == 2576):
                    print ("Here:", charList, nucBin, lenMut, "0000")

                i += 1
                # if (nodeMut.nucPosition == 2):
                #     print (eachNode, nodeMut.nucPosition, nodeMut.nucGapExist, nodeMut.blockId, int(eachGapPos), int(eachPos), nodeMut.mutInfo, nucBin, lenMut, "0000")

    for eachPos in SdelDict[eachNode]:
        for eachGapPos in SdelDict[eachNode][eachPos]:
            
            loc = deepcopy(SdelDict[eachNode][eachPos][eachGapPos])
            # print (loc)
            loc.sort()
            # print (loc)
            i = 0
            while i < len(loc):
                nodeMut = node.nucMutation.add()
                currPos = loc[i]
                currIdx = loc[i]
                length = 1
                
                while True:
                    if (i+1) == len(loc):
                        # print (currPos, length)
                        break
                
                    if currIdx == (loc[i+1] - 1):
                        currIdx = loc[i+1]
                        length += 1
                        i += 1
                        if (length >= 6):
                            # print (currPos, length)
                            break
                    else:
                        # print (currPos, length)
                        break
                
                nodeMut.nucPosition = currPos - 1
                nodeMut.nucGapExist = False
                blockId  = blockIDConv(int(eachGapPos), int(eachPos))
                nodeMut.blockId = int(blockId, 2)
                if (eachGapPos == "-1"):
                    nodeMut.blockGapExist = False
                else:
                    nodeMut.blockGapExist = True
                lenBin = bin(length)[2:]
                lenMut = "0" * (4 - len (lenBin)) + lenBin 
                nucBin = "0" * 24
                # print (lenMut, nucBin)
                nodeMut.mutInfo = int( nucBin + lenMut + "0001" ,2)
                if (nodeMut.mutInfo == 2576):
                    print ("Here:", nucBin, lenMut, "0001")
                i += 1
                # print (eachNode, nodeMut.nucPosition, nodeMut.nucGapExist, nodeMut.blockId, int(eachGapPos), int(eachPos), nodeMut.mutInfo, nucBin, lenMut, "0000")

    for eachPos in SinsDict[eachNode]:
        for eachGapPos in SinsDict[eachNode][eachPos]:
            
            loc = deepcopy(SinsDict[eachNode][eachPos][eachGapPos][0])
            loc.sort()
            # print (loc)
            i = 0
            while i < len(loc):
                nodeMut = node.nucMutation.add()
                currPos = loc[i]
                currIdx = loc[i]
                length = 1
                idxUnsort = SinsDict[eachNode][eachPos][eachGapPos][0].index(currIdx)
                charList = [SinsDict[eachNode][eachPos][eachGapPos][1][idxUnsort]]
                
                while True:
                    if (i+1) == len(loc):
                        # print (currPos, length, charList)
                        break
                
                    if currIdx == (loc[i+1] - 1):
                        currIdx = loc[i+1]
                        length += 1
                        i += 1
                        idxUnsort = SinsDict[eachNode][eachPos][eachGapPos][0].index(currIdx)
                        charList.append(SinsDict[eachNode][eachPos][eachGapPos][1][idxUnsort])
                        if (length >= 6):
                            break
                    else:
                        # print (currPos, length, charList)
                        break
                nodeMut.nucPosition = currPos - 1
                nodeMut.nucGapExist = False
                blockId  = blockIDConv(int(eachGapPos), int(eachPos))
                nodeMut.blockId = int(blockId, 2)
                if (eachGapPos == "-1"):
                    nodeMut.blockGapExist = False
                else:
                    nodeMut.blockGapExist = True
                lenBin = bin(length)[2:]
                lenMut = "0" * (4 - len (lenBin)) + lenBin 
                nucBin = nucBinConv (charList)
                # print ("Insertion")
                nodeMut.mutInfo = int( nucBin + lenMut  + "0010" ,2)
                if (nodeMut.mutInfo == 2576):
                    print ("Here:", nucBin, lenMut, "0010")
                i += 1
    
    ## Deletions -> Insertions and Deletions 
    for eachPos in delDict[eachNode]:
        for eachGapPos in delDict[eachNode][eachPos]:
            
            loc = deepcopy(delDict[eachNode][eachPos][eachGapPos])
            loc.sort()
            # print (loc)
            i = 0
            while i < len(loc):
                nodeMut = node.nucMutation.add()
                currPos = loc[i]
                currIdx = loc[i]
                length = 1
                
                while True:
                    if (i+1) == len(loc):
                        # print (currPos, length)
                        break
                
                    if currIdx == (loc[i+1] - 1):
                        currIdx = loc[i+1]
                        length += 1
                        i += 1
                        if (length >= 6):
                            # print (currPos, length)
                            break
                    else:
                        # print (currPos, length)
                        break
                
                nodeMut.nucPosition = currPos - 1
                nodeMut.nucGapExist = False
                blockId  = blockIDConv(int(eachGapPos), int(eachPos))
                nodeMut.blockId = int(blockId, 2)
                if (eachGapPos == "-1"):
                    nodeMut.blockGapExist = False
                else:
                    nodeMut.blockGapExist = True
                lenBin = bin(length)[2:]
                lenMut = "0" * (4 - len (lenBin)) + lenBin 
                nucBin = "0" * 24
                # print (lenMut, nucBin)
                nodeMut.mutInfo = int( nucBin + lenMut + "0001" ,2)
                if (nodeMut.mutInfo == 2576):
                    print ("Here:", nucBin, lenMut, "0001")
                i += 1

    for eachPos in insDict[eachNode]:
        for eachGapPos in insDict[eachNode][eachPos]:
            
            loc = deepcopy(insDict[eachNode][eachPos][eachGapPos][0])
            loc.sort()
            # print (loc)
            i = 0
            while i < len(loc):
                nodeMut = node.nucMutation.add()
                currPos = loc[i]
                currIdx = loc[i]
                length = 1
                idxUnsort = insDict[eachNode][eachPos][eachGapPos][0].index(currIdx)
                charList = [insDict[eachNode][eachPos][eachGapPos][1][idxUnsort]]
                
                while True:
                    if (i+1) == len(loc):
                        # print (currPos, length, charList)
                        break
                
                    if currIdx == (loc[i+1] - 1):
                        currIdx = loc[i+1]
                        length += 1
                        i += 1
                        idxUnsort = insDict[eachNode][eachPos][eachGapPos][0].index(currIdx)
                        charList.append(insDict[eachNode][eachPos][eachGapPos][1][idxUnsort])
                        if (length >= 6):
                            break
                    else:
                        # print (currPos, length, charList)
                        break
                nodeMut.nucPosition = currPos - 1
                nodeMut.nucGapExist = False
                blockId  = blockIDConv(int(eachGapPos), int(eachPos))
                nodeMut.blockId = int(blockId, 2)
                if (eachGapPos == "-1"):
                    nodeMut.blockGapExist = False
                else:
                    nodeMut.blockGapExist = True
                lenBin = bin(length)[2:]
                lenMut = "0" * (4 - len (lenBin)) + lenBin 
                nucBin = nucBinConv (charList)
                # print ("Insertion")
                nodeMut.mutInfo = int( nucBin + lenMut  + "0010" ,2)
                if (nodeMut.mutInfo == 2576):
                    print ("Here:", nucBin, lenMut, "0010")
                i += 1

    ## Insertions -> Insertions, Deletions, and Substitution
    for eachPos in IinsDict[eachNode]:
        for eachGapPos in IinsDict[eachNode][eachPos]:
            for eachNucGap in IinsDict[eachNode][eachPos][eachGapPos]:
                
                loc = deepcopy(IinsDict[eachNode][eachPos][eachGapPos][eachNucGap][0])
                loc.sort()
                # print (loc)
                i = 0
                while i < len(loc):
                    nodeMut = node.nucMutation.add()
                    currPos = loc[i]
                    currIdx = loc[i]
                    length = 1
                    idxUnsort = IinsDict[eachNode][eachPos][eachGapPos][eachNucGap][0].index(currIdx)
                    charList = [IinsDict[eachNode][eachPos][eachGapPos][eachNucGap][1][idxUnsort]]
                    
                    while True:
                        if (i+1) == len(loc):
                            # print (currPos, length, charList)
                            break
                    
                        if currIdx == (loc[i+1] - 1):
                            currIdx = loc[i+1]
                            length += 1
                            i += 1
                            idxUnsort = IinsDict[eachNode][eachPos][eachGapPos][eachNucGap][0].index(currIdx)
                            charList.append(IinsDict[eachNode][eachPos][eachGapPos][eachNucGap][1][idxUnsort])
                            if (length >= 6):
                                break
                        else:
                            # print (currPos, length, charList)
                            break
                    nodeMut.nucGapPosition = currPos
                    nodeMut.nucPosition = int(eachNucGap)
                    nodeMut.nucGapExist = True
                    blockId  = blockIDConv(int(eachGapPos), int(eachPos))
                    nodeMut.blockId = int(blockId, 2)
                    if (eachGapPos == "-1"):
                        nodeMut.blockGapExist = False
                    else:
                        nodeMut.blockGapExist = True
                    lenBin = bin(length)[2:]
                    lenMut = "0" * (4 - len (lenBin)) + lenBin 
                    nucBin = nucBinConv (charList)
                    # print ("Insertion")
                    nodeMut.mutInfo = int( nucBin + lenMut  + "0010" ,2)
                    if (nodeMut.mutInfo == 2576):
                        print ("Here:", nucBin, lenMut, "0010")
                    i += 1

    for eachPos in IdelDict[eachNode]:
        for eachGapPos in IdelDict[eachNode][eachPos]:
            for eachNucGap in IdelDict[eachNode][eachPos][eachGapPos]:
                
                loc = deepcopy(IdelDict[eachNode][eachPos][eachGapPos][eachNucGap])
                loc.sort()
                # print (loc)
                i = 0
                while i < len(loc):
                    nodeMut = node.nucMutation.add()
                    currPos = loc[i]
                    currIdx = loc[i]
                    length = 1
                    
                    while True:
                        if (i+1) == len(loc):
                            # print (currPos, length)
                            break
                    
                        if currIdx == (loc[i+1] - 1):
                            currIdx = loc[i+1]
                            length += 1
                            i += 1
                            if (length >= 6):
                                # print (currPos, length)
                                break
                        else:
                            # print (currPos, length)
                            break
                    
                    nodeMut.nucGapPosition = currPos
                    nodeMut.nucPosition = int(eachNucGap)
                    nodeMut.nucGapExist = True
                    blockId  = blockIDConv(int(eachGapPos), int(eachPos))
                    nodeMut.blockId = int(blockId, 2)
                    if (eachGapPos == "-1"):
                        nodeMut.blockGapExist = False
                    else:
                        nodeMut.blockGapExist = True
                    lenBin = bin(length)[2:]
                    lenMut = "0" * (4 - len (lenBin)) + lenBin 
                    nucBin = "0" * 24
                    # print (lenMut, nucBin)
                    nodeMut.mutInfo = int( nucBin + lenMut + "0001" ,2)
                    if (nodeMut.mutInfo == 2576):
                        print ("Here:", nucBin, lenMut, "0001")
                    i += 1

    for eachPos in IsubsDict [eachNode]:
        for eachGapPos in IsubsDict[eachNode][eachPos]:
            for eachNucGap in IsubsDict[eachNode][eachPos][eachGapPos]:
                
                loc = deepcopy(IsubsDict[eachNode][eachPos][eachGapPos][eachNucGap][0])
                loc.sort()
                # print (loc)
                i = 0
                while i < len(loc):
                    nodeMut = node.nucMutation.add()
                    currPos = loc[i]
                    currIdx = loc[i]
                    length = 1
                    idxUnsort = IsubsDict[eachNode][eachPos][eachGapPos][eachNucGap][0].index(currIdx)
                    charList = [IsubsDict[eachNode][eachPos][eachGapPos][eachNucGap][1][idxUnsort]]
                    
                    while True:
                        if (i+1) == len(loc):
                            # print (currPos, length, charList)
                            break
                    
                        if currIdx == (loc[i+1] - 1):
                            currIdx = loc[i+1]
                            length += 1
                            i += 1
                            idxUnsort = IsubsDict[eachNode][eachPos][eachGapPos][eachNucGap][0].index(currIdx)
                            charList.append(IsubsDict[eachNode][eachPos][eachGapPos][eachNucGap][1][idxUnsort])
                            if (length >= 6):
                                break
                        else:
                            # print (currPos, length, charList)
                            break
                    nodeMut.nucGapPosition = currPos
                    nodeMut.nucPosition = int(eachNucGap)
                    nodeMut.nucGapExist = True
                    blockId  = blockIDConv(int(eachGapPos), int(eachPos))
                    nodeMut.blockId = int(blockId, 2)
                    if (eachGapPos == "-1"):
                        nodeMut.blockGapExist = False
                    else:
                        nodeMut.blockGapExist = True
                    lenBin = bin(length)[2:]
                    lenMut = "0" * (4 - len (lenBin)) + lenBin 
                    nucBin = nucBinConv (charList)
                    # print (lenMut, nucBin)
                    nodeMut.mutInfo = int( nucBin + lenMut  + "0000" ,2)
                    if (nodeMut.mutInfo == 2576):
                        print ("Here:", charList, nucBin, lenMut, "0000")
                    i += 1

    for counterGlobal in gGaps:
        for counterLocal in gGaps[counterGlobal]:
            newGap = MAT.gaps.add()
            blockId  = blockIDConv(int(counterLocal), int(counterGlobal))
            # print (counterLocal, counterGlobal)
            newGap.blockId = int(blockId, 2)
            if (counterLocal == "-1"):
                newGap.blockGapExist = False
            else:
                newGap.blockGapExist = True
            nucPos = []
            nucGapLen = []
            for eachGap in gGaps[counterGlobal][counterLocal]:
                nucPos.append(eachGap[0])
                nucGapLen.append(eachGap[1])
            newGap.nucPosition.extend(nucPos)
            newGap.nucGapLength.extend(nucGapLen)


# for node in tree.traverse_preorder():
#     if (node.get_label() == checkId):
#         nodeQuery = node

# typeCheck = insDict
# while nodeQuery.is_root() == False:
#     nodeLabel = nodeQuery.get_label()
#     if nodeLabel in typeCheck:
#         if "0" in typeCheck[nodeLabel]:
#             if "4" in typeCheck[nodeLabel]["0"]:   
#                 print ("Subs:", nodeLabel, typeCheck[nodeLabel]["0"]["4"])
#     nodeQuery = nodeQuery.get_parent() 

# nodeLabel = nodeQuery.get_label()
# print (nodeLabel)
# if nodeLabel in typeCheck:
#     if "0" in typeCheck[nodeLabel]:
#         if "4" in typeCheck[nodeLabel]["0"]:   
#             print ("Subs:", nodeLabel, typeCheck[nodeLabel]["0"]["4"])



# print (insDict[checkId]["0"])
# print (delDict[checkId]["0"])
# print (IinsDict[checkId]["0"])
# print (IdelDict[checkId]["0"])
f = open(sys.argv[1], "wb")
f.write(MAT.SerializeToString())
f.close()
