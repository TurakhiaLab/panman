
from treeswift import *
from sys import argv


# Mutations associated with subsitions (Issue: pangraph does noe list all the mutations check 975 index in all seqs and consensus seq)
def mutation_subs(blockSeq, blockSub, blockSubIdx, tree):
    
    mutate = dict()
    current_char = dict() 
    mutation_current_char = dict()
    k=0
    for block in blockSeq:
        mutate[ block ] = dict()
        for node in tree.traverse_postorder():
            mutate[block][node.get_label()] = dict()
        k+=1
    # print(k)

    k=0
    for block in blockSeq:
        for char_index in blockSubIdx[block]:
            blockChar = blockSeq[block][char_index].upper()
            ## Handling Substitutions
            for node in tree.traverse_postorder():

                if( node.is_root() == True ):
                    rootLabel = node.get_label()

                if ( node.is_leaf() == True ):
                    mutation_current_char[node.get_label()] = set()
                    ## Note: We do not need to store current char: crosscheck
                    current_char[node.get_label()] = set()
                    current_char[node.get_label()].add(blockChar)
                    if (node.get_label()) in blockSub[block]: 
                        if char_index in blockSub[block][node.get_label()]:
                            current_char[node.get_label()] = set()
                            current_char[node.get_label()].add(blockSub[block][node.get_label()][char_index])

                
                else:
                    mutation_current_char[node.get_label()] = set()
                    children = node.child_nodes()
                    current_char[node.get_label()] = set()
                    for child in children:
                        # Edit: change to the case when two nodes are joined choose if a common char_index exist, not always union
                        if ( (current_char[node.get_label()] & current_char[child.get_label()]) == set()):
                            current_char[node.get_label()] = set.union(current_char[node.get_label()], current_char[child.get_label()])
                        else:
                            current_char[node.get_label()] = current_char[node.get_label()] & current_char[child.get_label()]        
            
            
            # print("blockChar", blockChar, "current_char", current_char[rootLabel])
            if blockChar in current_char[rootLabel]:
                current_char[rootLabel] = set(blockChar)
            else:
                if (current_char[rootLabel] == set("N")):
                    current_char[rootLabel] = set("N")
                else:
                    current_char_dummy = current_char[rootLabel].pop()
                    if (current_char_dummy == "N"):
                        current_char[rootLabel] = set(current_char[rootLabel].pop())
                    else:
                        current_char[rootLabel] = set(current_char_dummy)

            mutation_current_char[node.get_label()] = set(blockChar)
            # print(current_char)

            for node in tree.traverse_preorder():
                if(node.is_root() == False):
                    parent_char = next(iter(current_char[node.get_parent().get_label()]))
                    if (parent_char) in current_char[node.get_label()]:
                        current_char[node.get_label()] = set(parent_char)
                        mutation_current_char[node.get_label()] = set(parent_char)
                    else:
                        mutation_current_char[node.get_label()] = set(parent_char)


            for each_node in current_char:
                # print(current_char[each_node])
                if(current_char[each_node] != mutation_current_char[each_node]):
                    mutate[block][each_node][char_index] = str(mutation_current_char[each_node].pop()) + ":" +  str(current_char[each_node].pop())
                    # print( each_node, mutate[block][each_node][char_index] )
        k += 1
    return ( mutate )

def mutation_del(blockSeq, blockDel, tree):
    # Mutations associated with deletions  
    deletions_index_list = dict()
    for block in blockDel:
        deletions_index_list[block] = set()
        current_block = blockDel[block]
        for each_seq in current_block:
            for each_idx in current_block[each_seq]:
                    deletions_index_list[block].add(each_idx)

    # print( (deletions_index_list) )

    # Fitch algorithm on deletion columns

    current_char = dict() 
    mutate = dict()
    mutation_current_char = dict()
    for block in blockSeq:
        mutate[ block ] = dict()
        for node in tree.traverse_postorder():
            mutate[block][node.get_label()] = dict()

    for block in deletions_index_list:
        for each_col in deletions_index_list[block]:
            current_col = each_col 
            blockChar = blockSeq[block][each_col].upper()
            # print( each_seq, each_col)
            # print(seqs[each_seq][current_col])
            for node in tree.traverse_postorder():
                if( node.is_root() == True ):
                    rootLabel = node.get_label()

                if ( node.is_leaf() == True ):
                    mutation_current_char[node.get_label()] = set()
                    current_char[node.get_label()] = set()
                    current_char[node.get_label()].add(blockChar)
                    if node.get_label() in blockDel[block]:
                        if each_col in blockDel[block][node.get_label()]:
                            current_char[node.get_label()] = set()
                            current_char[node.get_label()].add("N")
                
                else:
                    mutation_current_char[node.get_label()] = set()
                    children = node.child_nodes()
                    current_char[node.get_label()] = set()
                    for child in children:
                        # Edit: change to the case when two nodes are joined choose if a common char_index exist, not always union
                        if( (current_char[node.get_label()] & current_char[child.get_label()]) == set() ):
                            current_char[node.get_label()] = set.union(current_char[node.get_label()], current_char[child.get_label()])
                        else:
                            current_char[node.get_label()] = (current_char[node.get_label()] & current_char[child.get_label()])                           
            
            
            if blockChar in current_char[rootLabel]:
                current_char[rootLabel] = set(blockChar)
            else:
                current_char[rootLabel] = set("N")
                
            mutation_current_char[node.get_label()] = set(blockChar)

            for node in tree.traverse_preorder():
                if(node.is_root() == False):
                    parent_char = next(iter(current_char[node.get_parent().get_label()]))
                    
                    if (parent_char) in current_char[node.get_label()]:
                        current_char[node.get_label()] = set(parent_char)
                        mutation_current_char[node.get_label()] = set(parent_char)
                    else:
                        mutation_current_char[node.get_label()] = set(parent_char)

            for each_node in current_char:
                if(current_char[each_node] != mutation_current_char[each_node]):
                    ## Remember if mutation_current_char_value is "N", then the operation becomes insertion
                    ## Must take care out it in MAT.py
                    mutation_current_char_value = mutation_current_char[each_node].pop()
                    mutate[block][each_node][each_col] = str(mutation_current_char_value) + ":" + str(current_char[each_node].pop())
                    # print( each_node, mutate[block][each_node][each_col] )

    return (mutate)

def mutation_ins(blockSeq, blockIns, blockGap, tree):
    current_char = dict() 
    mutate = dict()
    mutation_current_char = dict()
    gap_list = list()
    # for block in blockSeq:
    #     mutate[ block ] = dict()


    # Not handling multiple blocks as blockgap gives an global reference, not particular to any block
    
    for block in blockIns:
        mutate[ block ] = dict()
        
        for each_gap in blockIns[block]:
            mutate[ block ][each_gap] = dict()
            # gap_list[ block ][each_gap] = dict()
            for node in tree.traverse_postorder():
                mutate[ block ][each_gap][node.get_label()] = dict()
                # gap_list[ block ][each_gap][node.get_label()] = dict()
            
            for insert_index in range(0, blockGap[block][each_gap]):
                insertions_index_list = dict()
                insertions_index_list[block] = dict()
                current_block = blockIns[block]
                blockChar = "N"
                for node in tree.traverse_postorder():
                    if( node.is_root() == True ):
                        rootLabel = node.get_label()

                    if ( node.is_leaf() == True ):
                        mutation_current_char[node.get_label()] = set()
                        current_char[node.get_label()] = set()
                        current_char[node.get_label()].add(blockChar)
                        
                        if node.get_label() in blockIns[block][each_gap]: 
                            # print( node.get_label() ,  blockIns[block][each_gap], insert_index, blockIns[block][each_gap][node.get_label()])
                            if insert_index in blockIns[block][each_gap][node.get_label()]:
                                current_char[node.get_label()] = set()
                                current_char[node.get_label()].add( blockIns[block][each_gap][node.get_label()][insert_index] )

                    else:
                        mutation_current_char[node.get_label()] = set()
                        children = node.child_nodes()
                        current_char[node.get_label()] = set()
                        for child in children:
                            # Edit: change to the case when two nodes are joined choose if a common char_index exist, not always union
                            if( (current_char[node.get_label()] & current_char[child.get_label()]) == set() ):
                                current_char[node.get_label()] = set.union(current_char[node.get_label()], current_char[child.get_label()])
                            else:
                                current_char[node.get_label()] = (current_char[node.get_label()] & current_char[child.get_label()])                           
                
                            
                            # current_char[node.get_label()] = set.union(current_char[node.get_label()], current_char[child.get_label()])


                ## blockChar is "N" so do not need to check for multiple "N" conditions
                if blockChar in current_char[rootLabel]:
                    current_char[rootLabel] = set(blockChar)
                else:
                    current_char[rootLabel] = set(current_char[rootLabel].pop())
                mutation_current_char[node.get_label()] = set(blockChar)


                for node in tree.traverse_preorder():
                    if(node.is_root() == False):
                        parent_char = next(iter(current_char[node.get_parent().get_label()]))
                        if (parent_char) in current_char[node.get_label()]:
                            current_char[node.get_label()] = set(parent_char)
                            mutation_current_char[node.get_label()] = set(parent_char)
                        else:
                            mutation_current_char[node.get_label()] = set(parent_char)

                for each_node in current_char:
                    if(current_char[each_node]!=mutation_current_char[each_node]):
                        mutate[block][each_gap][each_node][insert_index] = str(mutation_current_char[each_node].pop()) + ":" + str(current_char[each_node].pop())

            gap_list.append([block, each_gap, blockGap[block][each_gap]])
                # gap_list[block][each_gap][rootLabel]["G"] =  str(each_gap) + ":" + str(blockGap[block][each_gap])

    return (mutate, gap_list)

"""
from read_json import *
from read_seq import *
tree = read_tree_newick("/home/AD.UCSD.EDU/swalia/data/sars/pangraph/sars_4k.nwk")

label=0
for node in tree.traverse_postorder():
    if(node.get_label()==None):
        node.set_label(str(label))
        label += 1 

f = open("/home/AD.UCSD.EDU/swalia/data/sars/pangraph/sars_4k.json","r")
data = json.load(f)
f.close()
blockSeq, blockSub, blockSubIdx, blockDel, blockGap, blockIns = getMutations( data )
get_paths = getPaths( data )
print("getMutation end here")

# mutate = mutation_subs(blockSeq, blockSub, blockSubIdx, tree)
# print( mutate )
# mutate = mutation_del(blockSeq, blockDel, tree)
# print( mutate )
mutate, gap_list = mutation_ins(blockSeq, blockIns, blockGap, tree)
print(mutate)
"""

"""
each_block = "SGBRFSJOXY"
each_seq = "ON654166.1"
gap_index=112
blockSeq, blockSub, blockSubIdx, blockDel, blockGap, blockIns = getMutations( data )
get_paths = getPaths( data )
print("getMutation end here")

mutate, gap_list = mutation_ins(blockSeq, blockIns, blockGap, tree)
print("mutate end here")

for node in tree.traverse_postorder():
    if(node.get_label()==None):
        node.set_label(str(label))
        label += 1 
    if (node.get_label() == each_seq):
        save_node = node
    if (node.is_root()) == True:
        root_node = node

while (True):
    print(save_node.get_label(), mutate[each_block][gap_index][save_node.get_label()])        
    if(save_node == root_node):
        break
    save_node = save_node.get_parent()
"""    





