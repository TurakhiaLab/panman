##################### Experiements for 4000 RSV sequences ############################################

# Download the PanMAN
DEST_FILE="rsv_4000.fa"
FILE_ID="1DAvh5mNOPTl5KL06QONy47oFgHKTo7uT"
CONFIRM=$(wget --quiet --save-cookies cookies.txt --keep-session-cookies --no-check-certificate "https://drive.google.com/uc?export=download&id=$FILE_ID" -O - | grep -o 'confirm=[^&]*' | sed 's/confirm=//')
if [ "$CONFIRM" ]; then
  wget --load-cookies cookies.txt -O "$DEST_FILE" "https://drive.google.com/uc?export=download&id=$FILE_ID&confirm=$CONFIRM"
else
  wget --load-cookies cookies.txt -O "$DEST_FILE" "https://drive.google.com/uc?export=download&id=$FILE_ID"
fi
rm cookies.txt

# Constructing PanGraph from the Raw Sequences
./pangraph build -k mmseqs $DEST_FILE | ./pangraph polish > rsv_4000.json 2> rsv_4000.nwk
awk '/tree/ {{split($0,a,"tree:  "); print a[2]}}' rsv_4000.nwk > temp.newick && mv temp.newick rsv_4000.nwk

# Constructing GFA from raw genome sequences using PGGB
./pggb -i $DEST_FILE -t 32 -n num_sequences -o output_dir

# Constructing VG
./vg convert -p rsv_4000.gfa -t 32 -v > rsv_4000.vg 

# Constructing GBZ
./vg gbwt -G rsv_4000.gfa --num-threads 32 --gbz-format -g rsv_4000.gbz

# Constructing MSA from raw genome sequences
./mafft --auto --keeplength --addfragments $DEST_FILE > rsv_4000.msa

# Constructing PanMANs from Pangraph alignment
./panmanUtils -P rsv_4000.json -N rsv_4000.nwk -o data

# Constructing PanMANs from GFA
./panmanUtils -G rsv_4000.gfa -N rsv_4000.nwk -o data

# Constructing PanMANs from MSA
./panmanUtils -M rsv_4000.msa -N rsv_4000.nwk -o data

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






