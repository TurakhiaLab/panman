import json



# Returns paths associated with each sequence 
def getPaths( data ):
    # seqPaths: Path associated with each sequence is stored
    seqPaths = dict()
    paths = data['paths']
    noOfPaths = len( paths ) 

    # Iterating over each sequence path
    for path in paths:
        seq = path["name"]
        seqPaths[seq] = list()
        for block in path["blocks"]:
            blockId = block["id"]
            seqPaths[seq].append(blockId)

    return ( seqPaths )


def getMutations( data ):
    blocks   = data["blocks"]
    blockSub = dict()
    blockSubIdx = dict()
    blockDel = dict()
    blockGap = dict()
    blockIns = dict()
    blockSeq = dict()
    for block in blocks:
        blockId = block["id"] # block Identifier
        blockSeq[blockId] = block["sequence"].upper() # block's consensus sequence 
        
        # gaps in block
        blockIns[blockId] = dict()
        blockGap[blockId] = dict()
        for gapPos in block["gaps"]:
            blockGap[blockId][int(gapPos)] = block["gaps"][gapPos]
            blockIns[blockId][int(gapPos)] = dict()
            # blockGap.append([gapPos,block["gaps"][gapPos]])
        # blockGap = block["gaps"] 
        # print( "blockGap= ",blockGap )
        
        # substitutions in block
        blockSub[blockId] = dict()
        blockSubIdx[blockId] = set()
        for seqSubs in block["mutate"]:
            seqName = seqSubs[0]["name"]
            subsList = seqSubs[1]
            blockSub[blockId][seqName] = dict()
            for everySubs in subsList:
                blockSub[blockId][seqName][everySubs[0]-1] = everySubs[1].upper()
                blockSubIdx[blockId].add(everySubs[0]-1)


        # deletions in block
        blockDel[blockId] = dict()
        for seqSubs in block["delete"]:
            # ToDo: Does it handle multiple deletions
            seqName = seqSubs[0]["name"]
            subsList = seqSubs[1]
            blockDel[blockId][seqName] = set()
            for everySubs in subsList:
                for each_del in range(everySubs[0]-1, everySubs[0]-1 + everySubs[1]):
                    blockDel[blockId][seqName].add(each_del)

        # Insertions in block
        for each_gap in blockGap[blockId]:
            blockIns[blockId][each_gap] = dict()
            
        for seqSubs in block["insert"]:
            # ToDo: Does it handle multiple insertions
            seqName = seqSubs[0]["name"]
            subsList = seqSubs[1]
            
            for everySubs in subsList:
                if( len(everySubs) > 1 ):
                    blockIns[blockId][everySubs[0][0]][seqName] = dict()
                    gap_length = len( everySubs[1])
                    # print( "gap length", gap_length)
                    for each_insert in range(gap_length):
                        # print( blockIns[blockId][everySubs[0][0]-1][seqName], len(everySubs[1]), (each_insert))
                        blockIns[blockId][everySubs[0][0]][seqName][ everySubs[0][1]+each_insert ] = everySubs[1][each_insert].upper()
                    # print(everySubs[0],everySubs[1])
    

    return ( blockSeq, blockSub, blockSubIdx, blockDel, blockGap, blockIns )


"""
f = open("/home/AD.UCSD.EDU/swalia/data/sars/pangraph/sars_4k.json","r")
data = json.load(f)
f.close()
blockSeq, blockSub, blockSubIdx, blockDel, blockGap, blockIns = getMutations( data )

each_block = "SGBRFSJOXY"
each_seq = "ON655904.1"
print( blockSub[each_block][each_seq] )
print( blockDel[each_block][each_seq] )
print( blockGap[each_block] )
for each_gap in  blockGap[each_block]:
    print( each_gap, blockIns[each_block][each_gap][each_seq])
"""


