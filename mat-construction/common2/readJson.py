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

'''
f = open("/home/AD.UCSD.EDU/swalia/data/ecoli/pangraph/ecoli_10.json","r")
data = json.load(f)
f.close()

blocks = getBlocks (data)
blockInfo = getBlockInfo ("QDKTNRNHWP", blocks)
print (getBlockSubs (blockInfo))
'''