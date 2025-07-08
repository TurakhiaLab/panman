# Script to generate Smaller PanMANs

set -x
ARCH=$(uname -m)
if [ "$ARCH" = "linux/amd64" ] || [ "$ARCH" = "x86_64" ]; then
    echo "Running on "$ARCH" architecture"
else
    echo "ERROR: Running on "$ARCH" architecture. Cannot proceed as the baseline tools (Pangraph v0.7.3, PGGB v0.6.0) are not supported on this architecture."
    echo "Please use this script with an x86_64  processor."
    return;
fi

gdown --folder https://drive.google.com/drive/folders/1pmSTV4akkB659P2RgNaMn5Tk9QetvoV-?usp=sharing
cd RawData_small

DEST_FILES="rsv_400.fa tb_40.fa sars_2000.fa HIV_2000.fa ecoli_100.fa klebs_100.fa"

for DEST_FILE in $DEST_FILES; do
    num=$(echo $DEST_FILE | cut -d '_' -f2 | cut -d '.' -f1)
    type=$(echo $DEST_FILE | cut -d '_' -f1)
    echo "Processing $DEST_FILE with $num genomes"
    
    
    # Constructing PanGraph from the Raw Sequences
    pangraph build $DEST_FILE > ${type}.json 2> ${type}_pangraph.out
    awk '/tree/ {{split($0,a,"tree:  "); print a[2]}}' ${type}_pangraph.out > temp.newick && mv temp.newick ${type}.nwk

    
    # Constructing GFA from raw genome sequences using PGGB
    samtools faidx $DEST_FILE
    if [ "$type" = "tb" ] || [ "$type" = "ecoli" ] || [ "$type" = "klebs" ]; then
        pggb -i $DEST_FILE -t 32 -n $num -o output_dir -X -b -v -T 16 -B 10K
    else
        pggb -i $DEST_FILE -t 32 -n $num -o output_dir -X -b -v -T 32 -B 10K
    fi
    mv output_dir/*.gfa ${type}.gfa

    # Constructing VG
    vg convert -g ${type}.gfa -t 32 -v > ${type}.vg

    # Constructing GBZ
    vg gbwt -G ${type}.gfa --num-threads 32 --gbz-format -g ${type}.gbz

    # Constructing PanMANs from Pangraph alignment
    panmanUtils -P ${type}.json -N ${type}.nwk -o ${type}

    # Extracting summary statistics from PanMANs
    panmanUtils -I panman/${type}.panman --summary

    # Extracting raw sequences in FASTA format
    # panmanUtils -I panman/${type}.panman --fasta -o ${type}

    # Extracting MSA in FASTA format
    # panmanUtils -I panman/${type}.panman --fasta-aligned -o ${type}

    # Extracting variations in VCF format
    # panmanUtils -I panman/${type}.panman --vcf -o ${type}

    # Converting PanMAN into a  variations in VCF format
    # panmanUtils -I panman/${type}.panman --gfa -o ${type}

done

