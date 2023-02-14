import mutation_annotation_pb2
import sys
import time
from treeswift import *

def B2KB(B):
    return float(B/1024)

MAT = mutation_annotation_pb2.tree()

# Read the address book.
# f = open("sars_1k.mat", "rb")
# MAT.ParseFromString(f.read())
# f.close()
# # print ("newick:", MAT.newick.ByteSize())
# nodeSize=0
# for each_node in MAT.nodes:
#     nodeSize += each_node.ByteSize()


# blcokSize=0
# for each_node in MAT.blocks:
#     blcokSize += each_node.ByteSize()

# gapSize=0
# for each_node in MAT.gaps:
#     gapSize += each_node.ByteSize()

# print ("nodes Size:", B2KB(nodeSize), "KB")
# print ("blocks Size:", B2KB(blcokSize), "KB")
# print ("gaps Size:", B2KB(gapSize), "KB")
# print ("blockGaps Size:", MAT.blockGaps.ByteSize(), "KB")
# f.close()


def distance_from_root(node):
    curr_node = node
    distance = 0
    while (curr_node.is_root() == False):
        distance += 1
        curr_node = curr_node.get_parent()
    return distance

def height_from_base(tree):

    height_each_node = {}

    innerNode = 0
    for node in tree.traverse_preorder():
        if (node.get_label() == None):
            node.set_label(innerNode)
            innerNode += 1
        else:
            continue

    for node in tree.traverse_postorder():
        if (node.is_leaf()):
            height_each_node[node] = 0
        
        else:
            children = node.child_nodes()
            height_each_node[node] = max(height_each_node[children[0]], height_each_node[children[1]]) + 1

    return height_each_node

f = open("HIV_5k.mat", "rb")
MAT.ParseFromString(f.read())
f.close()

newick_read = open("/home/AD.UCSD.EDU/swalia/data/HIV/pangraph/HIV_5000.nwk", "r")
tree_read = newick_read.readlines()[0].split("\n")[0]

tree_old = tree_read
newick_read.close()
tree = read_tree_newick(tree_old)

nodes = MAT.nodes


height_each_node = height_from_base(tree)


heightwise_mut = dict()

height_list = []
for node in tree.traverse_preorder():
    height_list.append(height_each_node[node])

total_node_size = 0
for i in range (len(nodes)):
    curr_height = height_list[i]

    if curr_height not in heightwise_mut:
        heightwise_mut[curr_height] = 0
    
    heightwise_mut[curr_height] += nodes[i].ByteSize()
    total_node_size += nodes[i].ByteSize()

for h in heightwise_mut:
    heightwise_mut[h] = B2KB(heightwise_mut[h])


heightwise_mut["All Nodes"] =  B2KB(total_node_size)

blockSize=0
for each_node in MAT.blocks:
    blockSize += each_node.ByteSize()
heightwise_mut["Blocks"] =  B2KB(blockSize)

gapSize=0
for each_node in MAT.gaps:
    gapSize += each_node.ByteSize()
heightwise_mut["Gaps"] =  B2KB(gapSize)

blockGapSize = B2KB(MAT.blockGaps.ByteSize())
heightwise_mut["Block Gaps"] =  B2KB(blockGapSize)

print (heightwise_mut)
print ("nodes Size:", B2KB(total_node_size), "KB")
print ("blocks Size:", B2KB(blockSize), "KB")
print ("gaps Size:", B2KB(gapSize), "KB")
print ("blockGaps Size:", B2KB(MAT.blockGaps.ByteSize()), "KB")

# Plotting
import matplotlib.pyplot as plt
import numpy as np

x = np.array([])
y = np.array([])
for each_entry in heightwise_mut:
    print(heightwise_mut[each_entry])
    x = np.append(x, each_entry)
    y = np.append(y, heightwise_mut[each_entry])

plt.plot(x,y)
plt.savefig("plot")