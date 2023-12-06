from treeswift import *
# from ete3 import Tree
from sys import argv
from readJson import *
import json
import time

MAX_CHAR_LEN = 6
SUB = "0000"
DEL = "0001"
INS = "0010"

start = time.time()
newick_read = open("/home/AD.UCSD.EDU/swalia/data/HIV/pangraph/HIV_5000.nwk", "r")
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

f = open("/home/AD.UCSD.EDU/swalia/data/HIV/HIV_msa.fasta", "r")
data = f.readlines()
f.close()

msa = dict()
sequence = ""
name = ""
for eachLine in data:
    if (eachLine[0] == ">"):
        # if name in msa:
            # print (len(msa[name]))
        name = eachLine.split("\n")[0][1:]
        msa[name] = ""
    else:
        msa[name] += eachLine.split("\n")[0]

mutations = dict()
mIns = dict()
mDel = dict()
mSub = dict()
mGap = dict()
gCoor = 0
lCoor = -1
consSeq = ""

for coor in range(len(msa[name])):
    currentChar = dict()
    for node in tree.traverse_postorder():
        nodeLabel = node.get_label()
        if (node.is_root()):
            rootLabel = nodeLabel
        
        if (node.is_leaf()):
            currentChar[nodeLabel] = [msa[nodeLabel][coor]]

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
                currentChar[nodeLabel] = ["-"]
            else:
                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]

            if currentChar[nodeLabel][0] == "-":
                lCoor += 1 #incrementing local coordinate within gap
                
            else:
                # ToDo: Store length not lCoor
                if (lCoor != -1):
                    if rootLabel not in mGap:
                        mGap[rootLabel] = [str(gCoor) + ":" + str(lCoor + 1)]
                    else:
                        mGap[rootLabel].append(str(gCoor) + ":" + str(lCoor + 1))

                gCoor += 1 #incrementing global coordinate
                lCoor = -1 #setting lCoor to -1
                consSeq += currentChar[nodeLabel][0]


        else:
            if currentChar[rootLabel][0] == "-":
                gCoorCurrent = gCoor
            else:
                gCoorCurrent = gCoor - 1 #use this gCoor in nodes other than root as we have already incremented it while reading root
            parentChar = currentChar[node.get_parent().get_label()][0]
            if parentChar in currentChar[nodeLabel]:
                currentChar[nodeLabel] = [parentChar]
            else:
                currentChar[nodeLabel] = [currentChar[nodeLabel][0]]
                ## Insertion
                if parentChar == "-":
                    if nodeLabel not in mIns:
                        mIns[nodeLabel] = ["I:" + str(gCoorCurrent) + ":" + str(lCoor) + ":" + currentChar[nodeLabel][0]]                                
                    else:
                        mIns[nodeLabel].append("I:" + str(gCoorCurrent) + ":" + str(lCoor) + ":" + currentChar[nodeLabel][0])                               


                else:
                    ## Deletion
                    if currentChar[nodeLabel][0] == "-":
                        if nodeLabel not in mDel:
                            mDel[nodeLabel] = ["D:" + str(gCoorCurrent) + ":" + str(lCoor)]                                
                        else:
                            mDel[nodeLabel].append("D:" + str(gCoorCurrent) + ":" + str(lCoor))
                    
                    ## Substitution
                    else:
                        if nodeLabel not in mSub:
                            mSub[nodeLabel] = ["S:" + str(gCoorCurrent) + ":" + str(lCoor) + ":" + currentChar[nodeLabel][0]]                               
                        else:
                            mSub[nodeLabel].append("S:" + str(gCoorCurrent) + ":" + str(lCoor) + ":" + currentChar[nodeLabel][0])
                        # if nodeLabel not in mSub:
                        #     mSub[nodeLabel] = dict()
                        # if gCoorCurrent not in mSub[nodeLabel]:
                        #     mSub[nodeLabel][gCoorCurrent] = dict()
                        # if lCoor not in mSub[nodeLabel][gCoorCurrent]:
                        #     mSub[nodeLabel][gCoorCurrent][lCoor] = list()
                                                    
                        # mSub[nodeLabel][gCoorCurrent][lCoor].append(currentChar[nodeLabel][0])

if (lCoor != -1):
    if rootLabel not in mGap:
        mGap[rootLabel] = [str(gCoor) + ":" + str(lCoor + 1)]
    else:
        mGap[rootLabel].append(str(gCoor) + ":" + str(lCoor + 1))


# exit(1)
import mutation_annotation_pb2
from copy import deepcopy
import time
import sys

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

def nucBinConv (listNuc):
    nucBin = ""
    for nuc in listNuc:
        nucBin = nucBin + char2bin(nuc)
    
    nucBin = "0" * (24 - len(nucBin)) + nucBin
    return nucBin

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

MAT = mutation_annotation_pb2.tree()

# Read the address book.
try:
  f = open(sys.argv[1], "rb")
  MAT.ParseFromString(f.read())
  f.close()
except IOError:
  print (sys.argv[1] + ": Could not open MAT.  Creating a new one.")



MAT.newick = tree_old

blocks = MAT.blocks.add()
blocks.blockId = 0
blocks.blockGapExist = False
blocks.consensusSeq.extend(seq2int(consSeq))

blockGap = MAT.blockGaps
blockPosition = [0]
blockGapLength = [0]
blockGap.blockPosition.extend(blockPosition)
blockGap.blockGapLength.extend(blockGapLength)

end = time.time()
print (end - start)

start = time.time()
for node in tree.traverse_preorder():
    newNode = MAT.nodes.add()
    nodeLabel = node.get_label()

    if (node.is_root()):
        blockMut = newNode.blockMutation.add()
        blockMut.blockId =  0
        blockMut.blockGapExist = False
        blockMut.blockMutInfo = True

        gaps = MAT.gaps.add()
        gaplist = []
        gaplen = []
        for eachGap in mGap[rootLabel]:
            eachGap = eachGap.split(":")
            if (eachGap[1] == "-1"):
                # gaplist.append(int(eachGap[0]))
                # gaplen.append(int(eachGap[1]))
                print (eachGap[0], eachGap[1])
            else:
                gaplist.append(int(eachGap[0]))
                gaplen.append(int(eachGap[1]))
            
            if (int(eachGap[0]) == 8741):
                print (int(eachGap[0]), int(eachGap[1]))
            if (int(eachGap[0]) == 8740):
                print (int(eachGap[0]), int(eachGap[1]))
            if (int(eachGap[0]) == 8742):
                print (int(eachGap[0]), int(eachGap[1]))
        print(len(gaplist), len(gaplen), len(mGap[rootLabel]))
        gaps.nucPosition.extend(gaplist)
        gaps.nucGapLength.extend(gaplen)
        gaps.blockId = 0
        gaps.blockGapExist = False




    # handle substitution
    if nodeLabel in mSub:
        globalSubsCoor = [] 
        globalSubsChar = []
        localSubsCoor = dict() 
        localSubsChar = dict()
        for eachSubs in mSub[nodeLabel]:
            eachSubs = eachSubs.split(":")
            gCoordinate = int(eachSubs[1])
            lCoordinate = int(eachSubs[2])
            char        = eachSubs[3]


            if (eachSubs[2] == "-1"):
                globalSubsCoor.append(gCoordinate)
                globalSubsChar.append(char)
            else:
                if gCoordinate not in localSubsCoor:
                    localSubsCoor[gCoordinate] = []
                    localSubsChar[gCoordinate] = []
                localSubsCoor[gCoordinate].append(lCoordinate)
                localSubsChar[gCoordinate].append(char)

        #  Global Char Concat
        loc = deepcopy(globalSubsCoor)
        loc.sort()
        
        i = 0
        while i < len(loc):
            nucMut = newNode.nucMutation.add()
            currPos = loc[i]
            currIdx = loc[i]
            length = 1
            idxUnsort = globalSubsCoor.index(currIdx)
            charList = [globalSubsChar[idxUnsort]]

            while True:
                if (i + 1) == len(loc):
                    break
                if currIdx == (loc[i+1] - 1):
                    currIdx = loc[i+1]
                    length += 1
                    i += 1
                    idxUnsort = globalSubsCoor.index(currIdx)
                    charList.append(globalSubsChar[idxUnsort])
                    if (length >= MAX_CHAR_LEN):
                        break
                else:
                    break
                
            nucMut.nucPosition = currPos
            nucMut.nucGapPosition = 0
            nucMut.nucGapExist = False
            nucMut.blockId = 0
            nucMut.blockGapExist = False

            lenBin = bin(length)[2:]
            lenMut = "0" * (4 - len (lenBin)) + lenBin
            nucBin = nucBinConv (charList)
            nucMut.mutInfo = int( nucBin + lenMut  + SUB ,2)

            # print (currPos, nucMut.nucGapExist, lenMut, length)
            i += 1

        #  Global Char Concat
        for globalSubsCoordinate in localSubsCoor:
            loc = deepcopy(localSubsCoor[globalSubsCoordinate])
            loc.sort()

            i = 0
            while i < len(loc):
                nucMut = newNode.nucMutation.add()
                currPos = loc[i]
                currIdx = loc[i]
                length = 1
                idxUnsort = localSubsCoor[globalSubsCoordinate].index(currIdx)
                charList = [localSubsChar[globalSubsCoordinate][idxUnsort]]

                while True:
                    if (i + 1) == len(loc):
                        break
                    if currIdx == (loc[i+1] - 1):
                        currIdx = loc[i+1]
                        length += 1
                        i += 1
                        idxUnsort = localSubsCoor[globalSubsCoordinate].index(currIdx)
                        charList.append(localSubsChar[globalSubsCoordinate][idxUnsort])
                        if (length >= MAX_CHAR_LEN):
                            break
                    else:
                        break
                    
                nucMut.nucPosition = globalSubsCoordinate
                nucMut.nucGapPosition = currPos
                nucMut.nucGapExist = True
                nucMut.blockId = 0
                nucMut.blockGapExist = False

                lenBin = bin(length)[2:]
                lenMut = "0" * (4 - len (lenBin)) + lenBin
                nucBin = nucBinConv (charList)
                nucMut.mutInfo = int( nucBin + lenMut  + SUB ,2)
                # print (currPos, nucMut.nucGapExist, lenMut)
                i += 1

            
            # nucMut = newNode.nucMutation.add()
            # nucMut.nucPosition = int(eachSubs[1])
            # nucMut.nucGapPosition = int(eachSubs[2])
            # if (eachSubs[2] == "-1"):
            #     nucMut.nucGapExist = False
            # else:
            #     nucMut.nucGapExist = True
            # nucBin = nucBinConv([eachSubs[3]])
            # nucMut.mutInfo = int( nucBin + "0001"  + "0000" ,2)
            # nucMut.blockId = 0
            # nucMut.blockGapExist = False
            
            # if (nucMut.nucPosition == 1478):
            #     print(nucMut.nucPosition, nucMut.nucGapPosition, nucMut.nucGapExist, nucMut.mutInfo)

    if nodeLabel in mDel:
        globalDelCoor = [] 
        localDelCoor = dict() 
        
        for eachDel in mDel[nodeLabel]:
            eachDel = eachDel.split(":")
            gCoordinate = int(eachDel[1])
            lCoordinate = int(eachDel[2])

            if (eachDel[2] == "-1"):
                globalDelCoor.append(gCoordinate)
            else:
                if gCoordinate not in localDelCoor:
                    localDelCoor[gCoordinate] = []
                localDelCoor[gCoordinate].append(lCoordinate)

        #  Global Char Concat
        loc = deepcopy(globalDelCoor)
        loc.sort()
        
        i = 0
        while i < len(loc):
            nucMut = newNode.nucMutation.add()
            currPos = loc[i]
            currIdx = loc[i]
            length = 1

            while True:
                if (i + 1) == len(loc):
                    break
                if currIdx == (loc[i+1] - 1):
                    currIdx = loc[i+1]
                    length += 1
                    i += 1
                    if (length >= MAX_CHAR_LEN):
                        break
                else:
                    break
                
            nucMut.nucPosition = currPos
            nucMut.nucGapPosition = 0
            nucMut.nucGapExist = False
            nucMut.blockId = 0
            nucMut.blockGapExist = False

            lenBin = bin(length)[2:]
            lenMut = "0" * (4 - len (lenBin)) + lenBin
            nucBin = "0" * 24
            nucMut.mutInfo = int( nucBin + lenMut  + DEL ,2)

            i += 1

        #  Global Char Concat
        for globalDelCoordinate in localDelCoor:
            loc = deepcopy(localDelCoor[globalDelCoordinate])
            loc.sort()

            i = 0
            while i < len(loc):
                nucMut = newNode.nucMutation.add()
                currPos = loc[i]
                currIdx = loc[i]
                length = 1

                while True:
                    if (i + 1) == len(loc):
                        break
                    if currIdx == (loc[i+1] - 1):
                        currIdx = loc[i+1]
                        length += 1
                        i += 1
                        if (length >= MAX_CHAR_LEN):
                            break
                    else:
                        break
                    
                nucMut.nucPosition = globalDelCoordinate
                nucMut.nucGapPosition = currPos
                nucMut.nucGapExist = True
                nucMut.blockId = 0
                nucMut.blockGapExist = False

                lenBin = bin(length)[2:]
                lenMut = "0" * (4 - len (lenBin)) + lenBin
                nucBin = "0" * 24
                nucMut.mutInfo = int( nucBin + lenMut  + DEL ,2)
                i += 1
            
            # nucMut = newNode.nucMutation.add()

            # eachDel = eachDel.split(":")
            # nucMut.nucPosition = int(eachDel[1])
            # nucMut.nucGapPosition = int(eachDel[2])
            # if (eachDel[2] == "-1"):
            #     nucMut.nucGapExist = False
            # else:
            #     nucMut.nucGapExist = True
            # nucBin = 24*"0"
            # nucMut.mutInfo = int( nucBin + "0001"  + "0001" ,2)
            # nucMut.blockId = 0
            # nucMut.blockGapExist = False

    if nodeLabel in mIns:
        globalInsCoor = [] 
        globalInsChar = []
        localInsCoor = dict() 
        localInsChar = dict()
        for eachIns in mIns[nodeLabel]:
            eachIns = eachIns.split(":")
            gCoordinate = int(eachIns[1])
            lCoordinate = int(eachIns[2])
            char        = eachIns[3]

            if (eachIns[2] == "-1"):
                globalInsCoor.append(gCoordinate)
                globalInsChar.append(char)
            else:
                if gCoordinate not in localInsCoor:
                    localInsCoor[gCoordinate] = []
                    localInsChar[gCoordinate] = []
                localInsCoor[gCoordinate].append(lCoordinate)
                localInsChar[gCoordinate].append(char)

        #  Global Char Concat
        loc = deepcopy(globalInsCoor)
        loc.sort()
        
        i = 0
        while i < len(loc):
            nucMut = newNode.nucMutation.add()
            currPos = loc[i]
            currIdx = loc[i]
            length = 1
            idxUnsort = globalInsCoor.index(currIdx)
            charList = [globalInsChar[idxUnsort]]

            while True:
                if (i + 1) == len(loc):
                    break
                if currIdx == (loc[i+1] - 1):
                    currIdx = loc[i+1]
                    length += 1
                    i += 1
                    idxUnsort = globalInsCoor.index(currIdx)
                    charList.append(globalInsChar[idxUnsort])
                    if (length >= MAX_CHAR_LEN):
                        break
                else:
                    break
                
            nucMut.nucPosition = currPos
            nucMut.nucGapPosition = 0
            nucMut.nucGapExist = False
            nucMut.blockId = 0
            nucMut.blockGapExist = False

            lenBin = bin(length)[2:]
            lenMut = "0" * (4 - len (lenBin)) + lenBin
            nucBin = nucBinConv (charList)
            nucMut.mutInfo = int( nucBin + lenMut  + INS ,2)

            # print (currPos, nucMut.nucGapExist, lenMut, length)
            i += 1

        #  Global Char Concat
        for globalInsCoordinate in localInsCoor:
            loc = deepcopy(localInsCoor[globalInsCoordinate])
            loc.sort()

            i = 0
            while i < len(loc):
                nucMut = newNode.nucMutation.add()
                currPos = loc[i]
                currIdx = loc[i]
                length = 1
                idxUnsort = localInsCoor[globalInsCoordinate].index(currIdx)
                charList = [localInsChar[globalInsCoordinate][idxUnsort]]

                while True:
                    if (i + 1) == len(loc):
                        break
                    if currIdx == (loc[i+1] - 1):
                        currIdx = loc[i+1]
                        length += 1
                        i += 1
                        idxUnsort = localInsCoor[globalInsCoordinate].index(currIdx)
                        charList.append(localInsChar[globalInsCoordinate][idxUnsort])
                        if (length >= MAX_CHAR_LEN):
                            break
                    else:
                        break
                    
                nucMut.nucPosition = globalInsCoordinate
                nucMut.nucGapPosition = currPos
                nucMut.nucGapExist = True
                nucMut.blockId = 0
                nucMut.blockGapExist = False

                lenBin = bin(length)[2:]
                lenMut = "0" * (4 - len (lenBin)) + lenBin
                nucBin = nucBinConv (charList)
                nucMut.mutInfo = int( nucBin + lenMut  + INS ,2)
                i += 1
            
            # nucMut = newNode.nucMutation.add()

            # eachIns = eachIns.split(":")
            # nucMut.nucPosition = int(eachIns[1])
            # nucMut.nucGapPosition = int(eachIns[2])
            # if (eachIns[2] == "-1"):
            #     nucMut.nucGapExist = False
            # else:
            #     nucMut.nucGapExist = True
            # nucBin = nucBinConv([eachIns[3]])
            # nucMut.mutInfo = int( nucBin + "0001"  + "0010" ,2)
            # nucMut.blockId = 0
            # nucMut.blockGapExist = False
            # if (nucMut.nucPosition == 471):
            #     print(nucMut.nucPosition, nucMut.nucGapPosition, nucMut.nucGapExist, nucMut.mutInfo)


end = time.time()
print (end - start)
f = open(sys.argv[1], "wb")
f.write(MAT.SerializeToString())
f.close()
