#!/usr/bin/env python
from Bio import SeqIO
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument('--in_fasta', type=str, help='file path to the input fasta')
parser.add_argument('--out-dir', type=str, help='output dir path')
args = parser.parse_args()

fasta_file = args.in_fasta

for record in SeqIO.parse(fasta_file, "fasta"):
    name = record.id
    title = name.replace(' ', '_').replace('/', '_')
    output_filename = f"{args.out_dir}/{title}.fa"
    with open(output_filename, 'w', encoding='utf-8') as outfile:
        outfile.write(">"+name+"\n")
        outfile.write(str(record.seq))