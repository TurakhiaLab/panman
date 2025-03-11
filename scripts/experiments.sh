
# Uncomment one of the following lines 
DEST_FILE="tb_400.fa"
# DEST_FILE="rsv_4000.fa"
# DEST_FILE="sars_20000.fa"
# DEST_FILE="HIV_20000.fa"
# DEST_FILE="ecoli_1000.fa"
# DEST_FILE="klebs_1000.fa"


gdown --folder https://drive.google.com/drive/folders/1bJ1GWOACNswomgK001WJCeQEkcoZhfGZ?usp=sharing
cd RawData
gunzip *.gz

# Constructing PanGraph from the Raw Sequences
pangraph build $DEST_FILE | pangraph polish >out.json 2>out.nwk
awk '/tree/ {{split($0,a,"tree:  "); print a[2]}}' out.nwk > temp.newick && mv temp.newick out.nwk

# Constructing GFA from raw genome sequences using PGGB
pggb -i $DEST_FILE -t 32 -n num_sequences -o output_dir

# Constructing VG
vg convert -p out.gfa -t 32 -v > out.vg 

# Constructing GBZ
vg gbwt -G out.gfa --num-threads 32 --gbz-format -g out.gbz

# Constructing MSA from raw genome sequences
mafft --auto --keeplength --addfragments $DEST_FILE > out.msa

# Constructing PanMANs from Pangraph alignment
panmanUtils -P out.json -N out.nwk -o data

# Constructing PanMANs from GFA
panmanUtils -G out.gfa -N out.nwk -o data

# Constructing PanMANs from MSA
panmanUtils -M out.msa -N out.nwk -o data

# Extracting summary statistics from PanMANs
panmanUtils -I out.panman --summary

# Extracting raw sequences in FASTA format
panmanUtils -I out.panman --fasta -o out

# Extracting MSA in FASTA format
panmanUtils -I out.panman --fasta-aligned -o out

# Extracting variations in VCF format
panmanUtils -I out.panman --vcf -o out

# Converting PanMAN into a  variations in VCF format
panmanUtils -I out.panman --gfa -o out

################ Experiments with 8M SARS COV-2 genomes #################

wget https://zenodo.org/records/14633185/files/sars_8M.aln.gz.xz?download=1
xz -d sars_8M.aln.gz.xz
gunzip sars_8M.aln.gz

TEMP_FILE="sars_8M_updated.fa"
while read -r line; do
    if [[ $line == ">"* ]]; then
        echo $line | awk '{print $1}' | sed 's/[|].*//g' >> "$TEMP_FILE"
        # replace_names "$line" >> "$TEMP_FILE"
    else
        echo "$line" >> "$TEMP_FILE"
    fi
done < "sars_8M.aln"

gdown --folder https://drive.google.com/drive/folders/10Vlycv5KBjwMB1RQfLR_CMBaqmniNpbf
panmanUtils --input-msa $TEMP_FILE --input-newick SARS_8M/sars_8M.nwk -o sars_8M --low-mem-mode --refFile SARS_8M/wuhan.aln
panmanUtils -I panman/sars_8M.panman --summary

panmanUtils --create-network panman/sars_8M.panman --input-file SARS_8M/recombinations.txt -o sars_8M_network
panmanUtils -I panman/sars_8M_network.panman --summary


