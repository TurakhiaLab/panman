
# from read_json import *
# from treeswift import *
# import sys


def map_blockid_fitch(map, count, blockid):
    if ( blockid ) not in map:
        count += 1
        
        map[blockid] = count
    return (map[blockid], count)


def fitch_block(tree, data, blockSeq, paths):

    label = 0
    for node in tree.traverse_postorder():
        if(node.get_label()==None):
            node.set_label(str(label))
            label += 1 


    ## Mapping blockid and Storing consensus sequence
    count_blockid = 0
    map = dict()
    for each_block in blockSeq:
        [block_id_compressed, count_blockid] = map_blockid_fitch(map, count_blockid, each_block)

    for each_seq in paths:
        new_block_list = list()
        for each_block in paths[each_seq]:
            new_block_list.append( map[each_block] )

        paths[each_seq] = new_block_list

    # print( paths )
    fitch_blocks = dict()
    current_char = dict()
    mutation_current_char = dict()
    for each_col in range(1, count_blockid + 1 ):
        for node in tree.traverse_postorder():
            each_node = node.get_label()
            current_char[each_node] = set()

            if (each_node) not in fitch_blocks:
                fitch_blocks[each_node] = list()


            if( node.is_root() == True ):
                rootLabel = each_node

            if ( node.is_leaf() == True ):
                if (each_col) in paths[each_node]:
                    # print( each_node )
                    current_char[each_node].add(each_col)
                else:
                    current_char[each_node].add(0)

            else:
                children = node.child_nodes()
                for each_child in children:
                    if ( (current_char[each_node] & current_char[each_child.get_label()]) == set()):
                        current_char[each_node] = set.union(current_char[each_node], current_char[each_child.get_label()])
                    else:
                        current_char[each_node] = current_char[each_node] & current_char[each_child.get_label()]

        # print( current_char[rootLabel], rootLabel )
        current_char_dummy = current_char[rootLabel].pop()
        current_char[rootLabel] = set()
        current_char[rootLabel].add(current_char_dummy)
        mutation_current_char[rootLabel] = set()
        mutation_current_char[rootLabel].add(0)

        for node in tree.traverse_preorder():
            if(node.is_root() == False):
                parent_char = next(iter(current_char[node.get_parent().get_label()]))
                mutation_current_char[node.get_label()] = set()
                mutation_current_char[node.get_label()].add(parent_char)
                if (parent_char) in current_char[node.get_label()]:
                    current_char[node.get_label()] = set()
                    current_char[node.get_label()].add(parent_char)

        for node in current_char:
            if(current_char[node] != mutation_current_char[node]):
                fitch_blocks[node].append( str(mutation_current_char[node].pop()) + ":" +  str(current_char[node].pop()) )

    return fitch_blocks
'''
newick_file = "/home/AD.UCSD.EDU/swalia/data/sars/pangraph/sars_4k.nwk"
json_file   = "/home/AD.UCSD.EDU/swalia/data/sars/pangraph/sars_4k.json"
fitch_blocks  = fitch_block(newick_file, json_file)
print (fitch_blocks)
# for node in fitch_blocks:
#     if( fitch_blocks[node] != []):
#         print(node, fitch_blocks[node])
newick_read = open(newick_file,"r") 
tree_old = newick_read.readlines()[0].split("\n")[0]
newick_read.close()
tree = read_tree_newick(tree_old)
label=0
for node in tree.traverse_postorder():
    if(node.get_label()==None):
        node.set_label(str(label))
        label += 1 
node = "ON655904.1"
for curr_node in tree.traverse_postorder():
    if (curr_node.get_label() == node):
        print (node, fitch_blocks[node])
        if(curr_node.is_root()):
            break
        node = curr_node.get_parent().get_label()
'''