from common.read_json import *
from common.read_seq import *
from common.fitch import *
from common.fitch_blocks import *
from treeswift import *

import v5.mutation_annotation_pb2
import sys
import time

MAT = v5.mutation_annotation_pb2.tree()


# Read the address book.
try:
  f = open(sys.argv[1], "rb")
  MAT.ParseFromString(f.read())
  f.close()
except IOError:
  print (sys.argv[1] + ": Could not open MAT.  Creating a new one.")


# Read newick tree
newick_read = open(sys.argv[2],"r") 
tree_old = newick_read.readlines()[0].split("\n")[0]
print (tree_old)
newick_read.close()
tree = read_tree_newick(tree_old)

label = 0
for node in tree.traverse_postorder():
    if(node.get_label()==None):
        node.set_label(str(label))
        label += 1 


# Read json file
f = open(sys.argv[3],"r")
data = json.load(f)
f.close()

## Reading Pangraph output json file
print("Reading pangraph's output json file ")
start_time = time.time()
blockSeq, blockSub, blockSubIdx, blockDel, blockGap, blockIns = getMutations( data )
paths = getPaths( data )
end_time = time.time()
print(str(end_time-start_time)+"\n")


# Applyting Fitch algorithm on substitutions, insertions and deletions
print("Applyting Fitch algorithm on substitutions ")
start_time = time.time()
mutate_sub = mutation_subs(blockSeq, blockSub, blockSubIdx,  tree) 
end_time = time.time()
print(str(end_time-start_time)+"\n")

print("Applyting Fitch algorithm on insertions ")
start_time = time.time()
mutate_ins, mutate_gap = mutation_ins(blockSeq, blockIns, blockGap,  tree)
end_time = time.time()
print(str(end_time-start_time)+"\n")

print("Applyting Fitch algorithm on deletions ")
start_time = time.time()
mutate_del = mutation_del(blockSeq, blockDel, tree)
end_time = time.time()
print(str(end_time-start_time)+"\n")

print("Applyting Fitch algorithm on blocks ")
start_time = time.time()
fitch_blocks = fitch_block(tree, data, blockSeq, paths)
end_time = time.time()
print(str(end_time-start_time)+"\n")


def map_blockid(map, count, blockid):
    if ( blockid ) not in map:
        count += 1
        value = bin(count)[2:]
        # Block_id -> 3B, therefore 24
        value = "0"*(24-len(value)) + value
        map[blockid] = value
    return (map[blockid], count)

def char2bin(char):
    bin = "1111"
    if ( char == "A"):
        bin = "0001"
    elif (char == "C"):
        bin = "0010"
    elif (char == "G"):
        bin = "0100"
    elif (char == "T"):
        bin = "1000"
    elif (char == "R"):
        bin = "0101"
    elif (char == "Y"):
        bin = "1010"
    elif (char == "S"):
        bin = "0110"
    elif (char == "W"):
        bin = "1001"
    elif (char == "K"):
        bin = "1100"
    elif (char == "M"):
        bin = "0011"
    elif (char == "B"):
        bin = "1110"
    elif (char == "D"):
        bin = "1101"
    elif (char == "H"):
        bin = "1011"
    elif (char == "V"):
        bin = "0111"
    else:
        bin = "1111"
    return (bin)

## Mapping blockid and Storing consensus sequence
count_blockid = 0
map = dict()

for each_block in blockSeq:
    [block_id_compressed, count_blockid] = map_blockid(map, count_blockid, each_block)
node_list = dict()

print(count_blockid)
new_count=0
max_block = 0
for each_block in map:
    value = int(map[each_block],2)
    if(value > max_block):
        max_block = value
    new_count += 1
print(max_block)
print(new_count)


## processing subsitution for protobuf format
print("processing subsitution for protobuf format ")
start_time = time.time()

count_subs = 0
count_subs_cont = 0
count_del = 0
count_del_cont = 0
count_ins = 0
count_ins_cont = 0
contiguous_block = 15

for each_block in mutate_sub:   
    for each_node in mutate_sub[each_block]:

        if (each_node) not in node_list:
            node_list[each_node] = list()
        
        sorted_index_list = list()
        for each_index in mutate_sub[each_block][each_node]:
            sorted_index_list.append(each_index)
        
        sorted_index_list.sort()
        if (sorted_index_list != []):
            index = 0

        else:
            continue

        prev_idx = 0
       
        while (index < len(sorted_index_list)):
            mut_list = list()
            idx = sorted_index_list[index]
            prev_idx = idx
            while (idx+1) in sorted_index_list:
                index += 1
                idx += 1
                count_subs += 1
                count_subs_cont += 1
                if ((idx - prev_idx) == contiguous_block):
                    break

            mut_list.append("S")
            mut_list.append(prev_idx)
            block_id_compressed = map[each_block]
           
            mut_list.append(block_id_compressed)
            # print("subs:", idx + 1 - prev_idx)
            mut_list.append(idx + 1 - prev_idx)

            mutation_chars = ""

            for i in range(prev_idx, idx + 1):
                mutation_string = mutate_sub[each_block][each_node][i]
                mutation_chars += char2bin(mutation_string[2])

            mut_list.append(mutation_chars)
            
            node_list[each_node].append(mut_list)

            index += 1
            count_subs += 1
            
            # if (block_id_compressed == "000000000000001010010101"):
            #     print("S:", prev_idx)    

end_time = time.time()
print(str(end_time-start_time)+"\n")


# # processing deletions for protobuf format
print("processing deletions for protobuf format ")
start_time = time.time()
for each_block in mutate_del:
    
    # print ("each block:", each_block)
    for each_node in mutate_del[each_block]:
        
        if (each_node) not in node_list:
            node_list[each_node] = list()         

        sorted_index_list = list()
        for each_index in mutate_del[each_block][each_node]:
            sorted_index_list.append(each_index)

        sorted_index_list.sort()
        if (sorted_index_list != []):
            index = 0

        else:
            continue
        # print (sorted_index_list)
        prev_idx = 0

        while (index < len(sorted_index_list)):
            mut_list = list()
            idx = sorted_index_list[index]
            prev_idx = idx
            mutation_string = mutate_del[each_block][each_node][idx]
            mutation_char = mutation_string[0] 
            if (mutation_char == "N"):
                mutation_type = "I"
            else:
                mutation_type = "D"

            while (idx+1) in sorted_index_list:
                mutation_string = mutate_del[each_block][each_node][idx + 1]
                if ((mutation_char == "N") & (mutation_string[0] != "N")):
                    break
                elif ((mutation_char != "N") & (mutation_string[0] == "N")):
                    break
                else:
                    index += 1
                    idx += 1
                    count_del += 1
                    count_del_cont += 1
                    if ((idx - prev_idx) == contiguous_block):
                        break

            mut_list.append(mutation_type)
            mut_list.append(prev_idx)
            
            block_id_compressed = map[each_block]
            mut_list.append(block_id_compressed)
            # print("del:", idx + 1 - prev_idx)
            
            mut_list.append(idx + 1 - prev_idx)

            mutation_chars = ""
            for i in range(prev_idx, idx + 1):
                mutation_string = mutate_del[each_block][each_node][i]
                if (mutation_type == "I"):
                    mutation_chars += char2bin(mutation_string[2])
            
            if (mutation_type == "I"):
                mut_list.append(mutation_chars)

            node_list[each_node].append(mut_list)

            index += 1
            count_del += 1
            # if (block_id_compressed == "000000000000001010010101"):
            #     print("D:", prev_idx)

end_time = time.time()
print(str(end_time-start_time)+"\n")


## processing insertions for protobuf format
print("processing insertions for protobuf format ")
start_time = time.time()
for each_block in mutate_ins:
    for each_gap in mutate_ins[each_block]:
        for each_node in mutate_ins[each_block][each_gap]:
            if (each_node) not in node_list:
                node_list[each_node] = list()
            
            sorted_index_list = list()
            for each_index in mutate_ins[each_block][each_gap][each_node]:
                sorted_index_list.append(each_index)

            sorted_index_list.sort()
            if (sorted_index_list != []):
                index = 0
            else:
                continue

            prev_idx = 0
            while (index < len(sorted_index_list)):
                mut_list = list()
                idx = sorted_index_list[index]
                prev_idx = idx
                mutation_string = mutate_ins[each_block][each_gap][each_node][idx]
                mutation_char = mutation_string[0] 
                if (mutation_char == "N"):
                    mutation_type = "I"
                else:
                    mutation_type = "D"

                while (idx+1) in sorted_index_list:
                    mutation_string = mutate_ins[each_block][each_gap][each_node][idx + 1]
                                        
                    if ((mutation_char == "N") & (mutation_string[0] != "N")):
                        break
                    elif ((mutation_char != "N") & (mutation_string[0] == "N")):
                        break
                    else:
                        index += 1
                        idx += 1
                        count_ins += 1
                        count_ins_cont += 1
                        if ((idx - prev_idx) == contiguous_block):
                            break
                
                mut_list.append(mutation_type)
                mut_list.append(each_gap)
                mut_list.append(prev_idx)

                block_id_compressed = map[each_block]
                mut_list.append(block_id_compressed)

                # print("ins: ",idx + 1 - prev_idx)
                mut_list.append(idx + 1 - prev_idx)                

                mutation_chars = ""

                for i in range(prev_idx, idx + 1):
                    mutation_string = mutate_ins[each_block][each_gap][each_node][i]
                    if (mutation_type == "I"):
                        mutation_chars += char2bin(mutation_string[2])
                
                if (mutation_type == "I"):
                    mut_list.append(mutation_chars)
                
                node_list[each_node].append(mut_list)
                    
                index += 1
                count_ins += 1
                # if (block_id_compressed == "000000000000001010010101"):
                #     print("I:", prev_idx, each_gap, prev_idx)


end_time = time.time()
print(str(end_time-start_time)+"\n")


# print("No of substitutions", count_subs, count_subs_cont, "No of deletions", count_del, count_del_cont, "No of insertions", count_ins, count_ins_cont )


def condensed_mut(blockid, length, type):
    length_bin = bin(length)[2:]
    length = "0"*(5-len(length_bin)) + length_bin
    # print("Non SNP: ", length_bin)
    return (blockid + length + type)

def condensed_mut_snp(blockid, length, char, type):
    length = "1"

    return (blockid + length + char +  type)


def int2blockid(num):
    num = bin(num)[2:]
    num = "0" * (24 - len(num)) + num
    return num

def int2gaplen(num):
    num = bin(num)[2:]
    # num = "0" * (8 - len(num)) + num
    return num

# print(map)

gap_list_position = list()
gap_list_block_id = list()
gap_list_gap_length = list()
for each_gap_mut in mutate_gap:
    gap_list_position.append(each_gap_mut[1])
    gap_list_block_id.append (int (map[each_gap_mut[0]], 2))
    gap_list_gap_length.append( int( int2gaplen(each_gap_mut[2]), 2))
    
    value = int (map[each_gap_mut[0]], 2)
    if (value > 11720):    
        print ("G: ", value, each_gap_mut[0])
# print( gap_list_position )
# print( gap_list_condensed )

## Storing mutations to MAT
# position: 32 bits
# gap_position: 32 bits
# condensed: 
    # subs/del/ins: block_id (24 bits) + 5 bits (len of nuc) + 3 bits (type of mutation)
    # subs(1)/ins(1): block_id (24 bits) + 1 bits (len of nuc) + 4 bits (mutated char) + 3 bits (type of mutation)
    # subs:     000
    # del:      001
    # ins:      010
    # subs(1):  011
    # ins(1):   100
# nuc: 64 bits (16 nucs)
count_block_mut = 0
SUBS    =   "000"
DEL     =   "001"
INS     =   "010"
SUBS1   =   "011"
INS1    =   "100"
DEL1    =   "101"
MAT.newick = tree_old 

count_subs = 0
count_subs1 = 0
count_del = 0
count_del1 = 0
count_ins = 0
count_ins1 = 0

for node in tree.traverse_postorder():
    each_node = node.get_label()

    new_node = MAT.nodes.add()

    # Listing gap mutations
    if (node.is_root() == True):
        MAT.gaps.position.extend(gap_list_position)
        MAT.gaps.block_id.extend(gap_list_block_id)
        MAT.gaps.gap_length.extend(gap_list_gap_length)

    # Listing nuc_mutations
    for each_mutation in node_list[each_node]:
        nuc_mut = new_node.nuc_mutation.add()
        

        ## Substitions
        if (each_mutation[0] == "S"):
            nuc_mut.position = each_mutation[1]
            if(each_mutation[3] == 1):
                count_subs1 += 1
                nuc_mut.condensed = int(condensed_mut_snp(each_mutation[2], each_mutation[3], each_mutation[4], SUBS1), 2 )
            else:
                count_subs += each_mutation[3]
                nuc_mut.condensed = int(condensed_mut(each_mutation[2], each_mutation[3], SUBS), 2 )
                # print(each_mutation[2], each_mutation[3], bin(nuc_mut.condensed))
                nuc_mut.nucs = int(each_mutation[4], 2)

        elif (each_mutation[0] == "D"):
            
            if ( len(each_mutation) == 4 ):
                if (each_mutation[3] == 1):
                    count_del1 += 1
                    nuc_mut.position = each_mutation[1]
                    nuc_mut.condensed = int( condensed_mut_snp( each_mutation[2], each_mutation[3], "0000" , DEL1 ), 2 )
                else:
                    count_del += each_mutation[3]
                    nuc_mut.position = each_mutation[1]
                    nuc_mut.condensed = int( condensed_mut( each_mutation[2], each_mutation[3] , DEL ), 2 )
                
            else:
                if (each_mutation[3] == 1):
                    count_del1 += 1
                    nuc_mut.position = each_mutation[1]
                    nuc_mut.gap_position = each_mutation[2]
                    nuc_mut.condensed = int( condensed_mut_snp( each_mutation[3], each_mutation[4], "0000" ,DEL1 ), 2 )
                else:
                    count_del += each_mutation[4]
                    nuc_mut.position = each_mutation[1]
                    nuc_mut.gap_position = each_mutation[2]
                    nuc_mut.condensed = int( condensed_mut( each_mutation[3], each_mutation[4], DEL ), 2 )
        else:
            
            if ( len(each_mutation) == 5 ):
                block__id = int(each_mutation[2],2)
                #if (block__id > new_count):
                    #print (block__id, each_mutation)
                if (each_mutation[3] == 1):
                    count_ins1 += 1
                    nuc_mut.position = each_mutation[1]
                    nuc_mut.condensed = int( condensed_mut_snp( each_mutation[2], each_mutation[3], each_mutation[4], INS1 ), 2 )
                else:
                    count_ins += each_mutation[3]
                    nuc_mut.position = each_mutation[1]
                    nuc_mut.condensed = int( condensed_mut( each_mutation[2], each_mutation[3], INS ), 2 )
                    nuc_mut.nucs = int( each_mutation[4], 2)
            else:
                # print ("I: ", each_mutation[4])
                block__id = int(each_mutation[3],2)
                #if (block__id > new_count):
                    #print (block__id, each_mutation)
                if (each_mutation[4] == 1):
                    count_ins1 += 1
                    nuc_mut.position = each_mutation[1]
                    nuc_mut.gap_position = each_mutation[2]
                    nuc_mut.condensed = int( condensed_mut_snp( each_mutation[3], each_mutation[4], each_mutation[5], INS1 ), 2 )

                else:
                    count_ins += each_mutation[4]
                    nuc_mut.position = each_mutation[1]
                    nuc_mut.gap_position = each_mutation[2]
                    nuc_mut.condensed = int( condensed_mut( each_mutation[3], each_mutation[4], INS ), 2 )
                    nuc_mut.nucs = int( each_mutation[5], 2)

    
    #  Listing block mutations
    block_mut_list = list()
    # print(each_node, fitch_blocks[each_node])
    for each_mut in fitch_blocks[each_node]:
        count_block_mut += 1
        mut_values = each_mut.split(":")
        # print(mut_values)
        # Insertion: blockid + 00000 + 010
        # Deletion: blockid + 00000 + 001
        
        if (mut_values[0] == "0"):
            # print("Insertion", mut_values[1] )
            blockid = int2blockid(int(mut_values[1]))
            block_mut_list.append( int((blockid + "00000010"), 2) )
        else:
            # print("Insertion", mut_values[0]) 
            blockid = int2blockid(int(mut_values[0]))
            block_mut_list.append( int((blockid + "00000001"), 2))
    
    if (( fitch_blocks[each_node] ) != []):
        # print (block_mut_list)
        new_node.block_mutation.condensed_block_mut.extend(block_mut_list)

print ("Total Substitutions: ", count_subs)
print ("Total Insertions: ", count_ins)
print ("Total Deletions: ", count_del)
print ("Total  SNP Substitutions: ", count_subs1)
print ("Total SNP Insertions: ", count_ins1)
print ("Total SNP Deletions: ", count_del1)


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
            dec = int( bin_compact, 2 )
            bin_compact = ""
            int_list.append(dec)
        
        else:
            continue
    return (int_list)

blockid_list = list()
consensus_list = list()

for each_block in blockSeq:
    new_block = MAT.blocks.add()
    new_block.block_id = int(("00000000" + map[each_block]), 2)
    consensus_seq_int_list = seq2int( blockSeq[each_block] )
    new_block.consensus_seq.extend(consensus_seq_int_list)


# print( count_block_mut )
# Write the new address book back to disk.
f = open(sys.argv[1], "wb")
f.write(MAT.SerializeToString())
f.close()







































