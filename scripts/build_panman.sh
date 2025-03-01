#!/bin/bash
  
## Defines
PANMAN_HOME=/home/panman
PANMAN_BUILD=/home/panman/build
DATASET_PATH=/home/dataset
DATASET=sars_20
PANGRAPH_HOME=/home/pangraph/pangraph.sh
PANGRAPH_OUTPUT=$PANMAN_HOME/build/pangraph
panmanUtils=$PANMAN_BUILD/panmanUtils

# method
method=$1

cd $PANMAN_BUILD

if [[ $method=="pangraph" ]] then
    ##### Commands generate PanGraph (JSON) and Tree Topology (Newick) from raw sequences in FASTA format ####
    mkdir -p pangraph
    echo "Building PanGraph..."
    $PANGRAPH_HOME "$DATASET_PATH/$DATASET.fa" "$PANGRAPH_OUTPUT/$DATASET.json" 2> "$PANGRAPH_OUTPUT/$DATASET.nwk"
    echo $(cat "$PANGRAPH_OUTPUT/$DATASET.nwk" | grep "tree" | awk '{split($0,a,"tree:  "); print a[2]}') > $PANGRAPH_OUTPUT/$DATASET.nwk

    #### Run panmanUtils to construct PanMAN using PanGraph ####
    echo "Building PanMAN from PanGraph..."
    $panmanUtils -P $PANGRAPH_OUTPUT/$DATASET.json -N $PANGRAPH_OUTPUT/$DATASET.nwk -o $DATASET
elif [[ $method=="msa" ]] then
    ##### Commands generate GFA (JSON) and Tree Topology (Newick) from raw sequences in FASTA format ####
    

    #### Run panmanUtils to construct PanMAN using PanGraph ####
    echo "Building PanMAN from MSA..."
    $panmanUtils -M $PANGRAPH_OUTPUT/$DATASET.msa -N $PANGRAPH_OUTPUT/$DATASET.nwk -o $DATASET
elif [[ $method=="gfa" ]] then
    ##### Commands generate GFA (JSON) and Tree Topology (Newick) from raw sequences in FASTA format ####
    

    #### Run panmanUtils to construct PanMAN using PanGraph ####
    echo "Building PanMAN from GFA..."
    $panmanUtils -G $PANGRAPH_OUTPUT/$DATASET.gfa -N $PANGRAPH_OUTPUT/$DATASET.nwk -o $DATASET
else
    echo "Error: Either pass pangraph or gfa or msa as an argument";
fi



