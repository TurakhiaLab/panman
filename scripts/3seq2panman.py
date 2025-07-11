import re
from sys import argv

def read_msa(msa_file):
    f=open(msa_file, "r")
    d=f.readlines()
    f.close()

    seqs=dict()

    seq_name=""
    seq=""
    for l in d:
        if (l[0]=='>'):
            if (seq!=""):
                seqs[seq_name]=seq
            seq_name=l[1:].split("\n")[0]
            seq=""
        else:
            seq+=l.split("\n")[0]
    seqs[seq_name]=seq
    
    return seqs

def msatorefcoord(position, seq):
    ref_seq=0
    for i in range(len(seq)):
        if (i == position):
            return ref_seq
        if (seq[i]=='-'):
            continue
        else:
            ref_seq+=1
    return ref_seq



def process_file(input_file, output_file, seqs):
    with open(input_file, 'r') as infile:
        data = infile.readlines()

    output = []
    map_ = dict()
    for line in data[1:]:
        # Split the line by whitespace (tab-separated format)
        row = line.strip().split()
        
        if len(row) < 13:
            continue
        
        breakpoint1 = row[12]
        breakpoint2 = row[14]

        # split_breakpoints = re.split(r' & |-', breakpoints)
        split_breakpoints1 = re.split('-', breakpoint1)
        split_breakpoints2 = re.split('-', breakpoint2)

        seq = seqs[row[2]]
        if (row[0] in map_):
            continue
        if (row[1] in map_):
            continue
        if (row[2] in map_):
            continue
        map_[row[0]]=1
        map_[row[1]]=1
        map_[row[2]]=1

        # Format the new row based on the transformation rules
        new_row = [
            "R",  # First column R
            0,
            row[0], 
            0, 
            row[1],  
            msatorefcoord(int(split_breakpoints1[0]), seq),  # 6th split value (start)
            msatorefcoord(int(split_breakpoints1[1]), seq),  # 7th split value (end)
            msatorefcoord(int(split_breakpoints2[0]), seq),  # 8th split value (start)
            msatorefcoord(int(split_breakpoints2[1]), seq),  # 9th split value (end)
            0,
            row[2],  
            'B'
        ]

        # Add the formatted row to the output list
        output.append("\t".join(map(str, new_row)))

    # If an output file is provided, write the results to that file
    if output_file:
        with open(output_file, 'w') as outfile:
            outfile.write("\n".join(output))
    else:
        # Otherwise, print the output to the console
        print("\n".join(output))

# Example usage:
if (len(argv)!=4):
    print("Usage python3 3seq2panman.py 3seq.rec 3seq_panman.rec in.msa")
input_file = argv[1]
output_file = argv[2]
msa_file = argv[3]

seqs = read_msa(msa_file)
process_file(input_file, output_file, seqs)
