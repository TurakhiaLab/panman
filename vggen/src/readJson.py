import json


# getPaths returns a map with leaf name as key and list of blocks corresponding to key as value
def getPaths (data):
    
    paths = data ['paths']
    seqPaths = dict()

    for path in paths:
        seq = path['name']
        seqPaths[seq] = list()
        for block in path ["blocks"]:
            blockId = block ["id"]
            seqPaths [seq].append (blockId)
    
    return seqPaths

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
f = open("/home/AD.UCSD.EDU/swalia/data/sars/pangraph/sars_100.json","r")
data = json.load(f)
f.close()
blocks = getBlocks (data)
# getPath = getPaths (data)
# print (blocks[0])
# VVQWGYGSJE'], 'JOYTOADXGX
blockInfo = getBlockInfo ("VVQWGYGSJE", blocks)
# print (blocks)
print ( getBlockSeq (blockInfo))
# print (getBlockGap (blockInfo))

'''