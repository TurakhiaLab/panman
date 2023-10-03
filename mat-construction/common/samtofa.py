import sys
import re

# Error Message
try:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
except:
    print ("Input and output files are not passed")
    exit()
    
# Read sam file
sam_file = open(input_file,"r")
sam_data = sam_file.readlines()
sam_file.close()

# msa dictionary
msa = dict()

# msa generation
count = 0
for each_line in sam_data:
    if ( each_line[0] == "@" ):
        continue
    elif ( each_line == "" or each_line == "\n"):
        break
    else:
        line_split = each_line.split("\t")
        seq_name = line_split[0] 
        tranformation = line_split[5]
        consensus_seq = line_split[9]
        if(consensus_seq == "*"):
            continue
        else:
            tranformation_array = re.split('(\d+)', tranformation)[1:]
            original_seq = ""
            flag = 0
            if ( seq_name ) in msa:
                original_seq = msa[ seq_name ]
                flag = 1
            
            original_seq_index = 0
            for (index) in range( 1, len(tranformation_array), 2 ):
                mutation = tranformation_array[index]
                mutation_len = int(tranformation_array[index-1])
                if( mutation == "H" ):
                    continue
                
                elif( mutation == "S"):
                    if( flag == 1 ):
                        # original_seq = original_seq[ : original_seq_index ] + consensus_seq[original_seq_index : original_seq_index + mutation_len ] + original_seq[ original_seq_index + mutation_len : ]
                        continue
                    else:
                        original_seq += consensus_seq[original_seq_index : original_seq_index + mutation_len ]

                elif (mutation == "M" or mutation == "I"):
                    if( flag == 1 ):
                        original_seq = original_seq[ : original_seq_index ] + consensus_seq[original_seq_index : original_seq_index + mutation_len ] + original_seq[ original_seq_index + mutation_len : ]
                    else:
                        original_seq += consensus_seq[original_seq_index : original_seq_index + mutation_len ]
                
                elif( mutation == "D"):
                    if( flag == 1 ):
                        original_seq = original_seq[ : original_seq_index ] + ("-" * mutation_len) + original_seq[ original_seq_index + mutation_len : ]

                    else:
                        original_seq += ("-" * mutation_len)
                

                original_seq_index += mutation_len

            msa[ seq_name ] = original_seq
            count += 1

f = open(output_file, "w")
for each_seq in msa:
    f.write(">")
    f.write(each_seq)
    f.write("\n")
    f.write(msa[each_seq])
    f.write("\n")

f.close()