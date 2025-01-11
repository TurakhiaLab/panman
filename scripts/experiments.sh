# Constructing PanGraph alignment and Tree from raw genome sequences
./pangraph build -k mmseqs -a 200 -b 30 data.fa | ./pangraph polish

# Constructing GFA from raw genome sequences using PGGB
./pggb -i data.fa -t 32 -n num_sequences -o output_dir

# Constructing VG
./vg convert -p data.gfa -t 32 -v > data.vg 

# Constructing GBZ
./vg gbwt -G data.gfa --num-threads 32 --gbz-format -g data.gbz

# Constructing MSA from raw genome sequences
./mafft --auto --keeplength --addfragments data.fa > data.msa

# Constructing PanMANs from Pangraph alignment
./panmanUtils -P data.pangraph -N data.nwk -o data

# Constructing PanMANs from GFA
./panmanUtils -G data.gfa -N data.nwk -o data

# Constructing PanMANs from MSA
./panmanUtils -M data.msa -N data.nwk -o data

##################### panmanUtils experiements for a PanMAN of 4000 RSV sequences ############################################

# Download the PanMAN



# Extracting summary statistics from PanMANs
./panmanUtils -I rsv_4000.panman --summary

# Extracting raw sequences in FASTA format
./panmanUtils -I rsv_4000.panman --fasta -o rsv_4000

# Extracting MSA in FASTA format
./panmanUtils -I rsv_4000.panman --fasta-aligned -o rsv_4000

# Extracting variations in VCF format
./panmanUtils -I rsv_4000.panman --vcf -o rsv_4000

# Converting PanMAN into a  variations in VCF format
./panmanUtils -I rsv_4000.panman --gfa -o rsv_4000






