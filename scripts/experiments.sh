
ARCH=$(uname -m)

# Uncomment one of the following lines 
DEST_FILE="rsv_4000.fa"
# DEST_FILE="tb_400.fa"
# DEST_FILE="sars_20000.fa"
# DEST_FILE="HIV_20000.fa"
# DEST_FILE="ecoli_1000.fa"
# DEST_FILE="klebs_1000.fa"


gdown --folder https://drive.google.com/drive/folders/1bJ1GWOACNswomgK001WJCeQEkcoZhfGZ?usp=sharing
cd RawData
gunzip *.gz

# Constructing PanGraph from the Raw Sequences
pangraph build $DEST_FILE >out.json 2>out.nwk
awk '/tree/ {{split($0,a,"tree:  "); print a[2]}}' out.nwk > temp.newick && mv temp.newick out.nwk

if [ "$ARCH" = "x86_64" ]; then
    # Constructing GFA from raw genome sequences using PGGB
    samtools faidx $DEST_FILE
    pggb -i $DEST_FILE -t 32 -n 4000 -o output_dir # change 4000 to number of genomes for different datasets
    mv output_dir/*.gfa out.gfa

    # Constructing VG
    vg convert -g out.gfa -t 32 -v > out.vg 

    # Constructing GBZ
    vg gbwt -G out.gfa --num-threads 32 --gbz-format -g out.gbz
fi

# Constructing MSA from raw genome sequences
mafft --auto $DEST_FILE > out.msa

# Constructing PanMANs from Pangraph alignment
panmanUtils -P out.json -N out.nwk -o data

if [ "$ARCH" = "x86_64" ]; then
    # Constructing PanMANs from GFA
    panmanUtils -G out.gfa -N out.nwk -o data

    # Constructing PanMANs from MSA
    panmanUtils -M out.msa -N out.nwk -o data
fi

# Extracting summary statistics from PanMANs
panmanUtils -I panman/out.panman --summary

# Extracting raw sequences in FASTA format
panmanUtils -I panman/out.panman --fasta -o out

# Extracting MSA in FASTA format
panmanUtils -I panman/out.panman --fasta-aligned -o out

# Extracting variations in VCF format
panmanUtils -I panman/out.panman --vcf -o out

# Converting PanMAN into a  variations in VCF format
panmanUtils -I panman/out.panman --gfa -o out

# Recombinations in tb, klebs, ecoli
panmanUtils -I panman/out.panman --fasta-aligned -o out
3seq -gen-p ptable 40 # change 40 to number of genomes for different datasets
3seq -f info/out_0.msa -ptable ptable -id out
python3 3seq2panman.py out.3s.rec out.3s.panman.rec info/out_0.msa

################ Experiments with 8M SARS COV-2 genomes #################
T_OUT=twilight_out
gdown --folder https://drive.google.com/drive/folders/10Vlycv5KBjwMB1RQfLR_CMBaqmniNpbf
mkdir $T_OUT
wget -O sars_8M.aln.gz.xz https://zenodo.org/records/15059329/files/sars_8M.aln.gz.xz?download=1
xz -d sars_8M.aln.gz.xz
gunzip sars_8M.aln.gz
sed '/^>/!s/-//g' sars_8M.aln > sars_8M.fa
twilight -t SARS_8M/sars_8M.nwk -i sars_8M.fa -p y -v -r 1 --trans -8 -o $T_OUT/sars_8M_without_wuhan.msa --gap-open -60
cp SARS_8M/wuhan.fa $T_OUT
twilight -f $T_OUT -a l -p y -v -r 1 --trans -8 -o sars_8M.msa --gap-open -60
grep -a -A 1 Wuhan-Hu-1 sars_8M.msa > wuhan.aln
panmanUtils --input-msa sars_8M.msa --input-newick SARS_8M/sars_8M.nwk -o sars_8M --low-mem-mode --refFile wuhan.aln

panmanUtils -I panman/sars_8M.panman --summary
panmanUtils -I panman/sars_8M.panman --fasta
panmanUtils -I panman/sars_8M.panman --vcf

panmanUtils --create-network panman/sars_8M.panman --input-file SARS_8M/sars_8M_cmut -o sars_8M_network
panmanUtils -I panman/sars_8M_network.panman --summary


