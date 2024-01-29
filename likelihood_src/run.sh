#!/bin/bash
# Execute the file with sudo privaleges


fp64_exec="$PWD/../build_phylo/core_likelihood_fp64"
fp32_exec="$PWD/../build_phylo/core_likelihood_fp32"
trees_path="/home/swalia@AD.UCSD.EDU/work/phylo-accel/dataset/singlegene_trees"
alignments_path="/home/swalia@AD.UCSD.EDU/work/phylo-accel/dataset/singlegene_alignments"
sars_path="/home/swalia@AD.UCSD.EDU/work/phylo-accel/dataset"
RESULT_FOLDER="$PWD/results/$(date +%Y-%m-%d)"

dataset=""
dataset+="SongD1 "
dataset+="JarvD5a "
dataset+="JarvD5b "
dataset+="MisoD2a "
dataset+="MisoD2b "
dataset+="PrumD6 "
# dataset+="sars_10k "

valid_list=""
valid_list+="aln "
valid_list+="fasta "


for ds in $dataset; do
    echo "$ds"
    mkdir -p $RESULT_FOLDER/$ds
    if [[ $ds == "sars_10k" ]]; then
        fasta=$sars_path/${ds}.fasta
        tree=$sars_path/${ds}.nwk
        $fp64_exec $fasta $tree  >> $RESULT_FOLDER/$ds/fp_64
        $fp32_exec $fasta $tree  >> $RESULT_FOLDER/$ds/fp_32
    else 
        for i in "$alignments_path/$ds"/*; do
            bname=$(basename $i)
            valid=(${bname//./ })
            if [[ "$valid_list" != *"${valid[-1]}"* ]]; then
                continue;
            fi
            fasta=$alignments_path/$ds/$bname
            bname="${bname%.*}"
            tree=$trees_path/$ds/Best_observed/$bname.best.tre
            $fp64_exec $fasta $tree >> $RESULT_FOLDER/$ds/fp_64
            $fp32_exec $fasta $tree >> $RESULT_FOLDER/$ds/fp_32
        done
    fi
done



