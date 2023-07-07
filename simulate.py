import subprocess
import time
import random
import sys
import os

from Bio import SeqIO

input_fasta = 'build/fasta/ten.fasta'
fasta_sequences = SeqIO.parse(open(input_fasta),'fasta')

for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    num_reads = 800
    with open(f'fastq_temp/{name}.fa', 'w') as f:
        f.write(f'>{name}\n{sequence}\n')
    subprocess.run(f'iss generate --genomes fastq_temp/{name}.fa --model miseq -n 5000 --output fastq_temp/{name} --cpus 100', shell=True)
    subprocess.run(f'seqtk sample -s100 fastq_temp/{name}_R1.fastq {num_reads} > fastq/{name}_sampled_{num_reads}.fastq', shell=True)


# while True:
#     node = bfs[random.randint(0, len(bfs)-1)]
#     ct += 1
#     if ct > 30:
#         break
#     name_raw = node.id
#     name = node.id.replace('|', '_').replace('/', '_')
#     genome = reconstruct_genome(node.id, t, ref)    
#     with open(f'{name}.fa', 'w') as f:
#         f.write(f'>{name}\n{genome}\n')
#     subprocess.run(f'iss generate --genomes {name}.fa --model miseq -n 5000 --output {name}_reads --cpus 100', shell=True)
#     subprocess.run(f'seqtk sample -s100 {name}_reads_R1.fastq {num_reads} > {name}_1_sampled_{num_reads}.fastq', shell=True)
#     subprocess.run(f'seqtk sample -s100 {name}_reads_R2.fastq {num_reads} > {name}_2_sampled_{num_reads}.fastq', shell=True)

