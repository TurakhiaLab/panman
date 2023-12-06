import json


# getPaths returns a map with leaf name as key and list of blocks corresponding to key as value
def getPaths (data):
    
    paths = data ['paths']
    seqPaths = dict()
    seqStrands = dict()

    for path in paths:
        seq = path['name']
        seqPaths[seq] = list()
        seqStrands[seq] = list()
        for block in path ["blocks"]:
            blockId = block ["id"]
            strand = block["strand"]
            seqPaths [seq].append (blockId)
            seqStrands[seq].append(strand)
    
    return seqPaths, seqStrands

def setOffset(data):
    paths = data ['paths']
    seqOffset = dict()

    for path in paths:
        c = path['offset']
        seq = path['name']
        seqOffset[seq] = c
    
    return seqOffset

def getBlocksLen (data):
    blocks = data['blocks']
    blockLen = dict()
    for block in blocks:
        blockLen[block['id']] = len(block['sequence'])
    return blockLen

def getBlocks (data):

    blocks = data ['blocks']
    return blocks

def getBlockInfo (block, blockData):

    for elem in blockData:
        if elem['id'] == block:
            return elem
        else:
            continue

def getBlockSeq (blockInfo):
    return blockInfo ["sequence"]

def getBlockSubs_helper (blockInfo, seq):
    blockSubs = blockInfo ['mutate']
    
    blockSubsData = dict()
    for elem in blockSubs:    
        if elem [0]['name'] != seq:
            continue
           
        blockSubsData [elem [0]['number']] = [elem [0]['strand']]

        coorList = []
        charList = []
        for mut in elem[1]:
            coorList.append (mut[0])
            charList.append (mut[1])
        blockSubsData [elem [0]['number']].append (coorList)
        blockSubsData [elem [0]['number']].append (charList)
    return blockSubsData

def getBlockSubs (blockInfo):
    blockSubs = blockInfo ['mutate']
    
    blockSubsData = dict()
    for elem in blockSubs:    
        if elem [0]['name'] not in blockSubsData:   
            blockSubsData [elem [0]['name']] = dict()
        blockSubsData [elem [0]['name']][elem [0]['number']] = []

        coorList = []
        charList = []
        for mut in elem[1]:
            coorList.append (mut[0])
            charList.append (mut[1])
        blockSubsData [elem [0]['name']][elem [0]['number']].append (coorList)
        blockSubsData [elem [0]['name']][elem [0]['number']].append (charList)
    return blockSubsData

def getBlockDel_helper (blockInfo, seq):
    blockDel = blockInfo ['delete']
    blockDelData = dict()
    for elem in blockDel:    
        if elem [0]['name'] != seq:   
            continue
        
        blockDelData [elem [0]['number']] = [elem [0]['strand']]

        coorList = []
        for mut in elem[1]:
            if mut != []:
                coorList += list(range(mut[0],mut[0] + mut[1]))
        blockDelData [elem [0]['number']].append(coorList)
    return blockDelData

def getBlockDel (blockInfo):
    blockDel = blockInfo ['delete']
    blockDelData = dict()
    for elem in blockDel:    
        if elem [0]['name'] not in blockDelData:   
            blockDelData [elem [0]['name']] = dict()
        blockDelData [elem [0]['name']][elem [0]['number']] = []

        coorList = []
        for mut in elem[1]:
            if mut != []:
                coorList += list(range(mut[0],mut[0] + mut[1]))
        blockDelData [elem [0]['name']][elem [0]['number']] += coorList
    return blockDelData

def getBlockIns_helper (blockInfo, seq):
    blockIns = blockInfo ['insert']

    blockInsData = dict()
    for elem in blockIns:
        if elem [0]['name'] != seq:
            continue   

        blockInsData [elem [0]['number']] = [elem [0]['strand']]

        coorList = []
        charList = []
        for mut in elem[1]:
            if mut != []:
                for z in range (len (mut [1])):
                    coorList    += [[mut[0][0], mut[0][1] + z]]
                    charList    += [mut[1][z]]
        
        blockInsData [elem [0]['number']] += [coorList, charList]
    return blockInsData

def getBlockIns (blockInfo):
    blockIns = blockInfo ['insert']

    blockInsData = dict()

    for elem in blockIns:
        if elem [0]['name'] not in blockInsData:   
            blockInsData [elem [0]['name']] = dict()
        blockInsData [elem [0]['name']][elem [0]['number']] = []

        coorList = []
        charList = []
        for mut in elem[1]:
            if mut != []:
                for z in range (len (mut [1])):
                    coorList    += [[mut[0][0], mut[0][1] + z]]
                    charList    += [mut[1][z]]
        
        blockInsData [elem [0]['name']][elem [0]['number']] += [coorList, charList]
    return blockInsData

def getBlockGap (blockInfo):
    return blockInfo ['gaps']


'''
f = open("/home/AD.UCSD.EDU/swalia/data/HIV/pangraph/HIV_5000.json","r")
data = json.load(f)
f.close()
# blocks = getBlocks (data)
getPath = getPaths (data)
count = 0
uniqueBlock = []
for eachSeq in getPath:
    for i in range(len(getPath[eachSeq])):
        if getPath[eachSeq][i] not in uniqueBlock:
            uniqueBlock.append(getPath[eachSeq][i])
         
print (len(uniqueBlock))
# for eachSeq in getPath:
#     if "RKCHLEBAON" in getPath[eachSeq]:
#         print (eachSeq, getPath[eachSeq])
# min_len = 1000
# max_len = -1
# for eachSeq in getPath:
#     curr_len = len(getPath[eachSeq])
#     if (min_len > curr_len):
#         min_len = curr_len
#         min_len_seq = eachSeq
#     if (max_len < curr_len):
#         max_len = curr_len
#         max_len_seq = eachSeq
# print (min_len, max_len)
# print (getPath[min_len_seq], getPath[max_len_seq])
'''