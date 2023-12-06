import subprocess
import time
import random
import sys
import os

from Bio import SeqIO

input_fastas = ['build/fasta/2k.fasta']

for f in input_fastas:
    fasta_sequences = list(SeqIO.parse(open(f),'fasta'))

    # pick 10 random nodes
    idxs = [random.randint(0, len(fasta_sequences)) for _ in range(100)]

    for idx in idxs:
        fasta = fasta_sequences[idx]
        name, sequence = fasta.id, str(fasta.seq)
        with open(f'fastq_temp/{name}.fa', 'w+') as f:
            f.write(f'>{name}\n{sequence}\n')
        subprocess.run(f'iss generate --model NovaSeq --genomes fastq_temp/{name}.fa -n 1000 --output fastq/{name}.1kreads.fastq --cpus 8', shell=True)
      
        # subprocess.run(f'seqtk sample -s100 fastq_temp/{name}_R1.fastq {num_reads} > fastq/{name}_sampled_{num_reads}.fastq', shell=True)
