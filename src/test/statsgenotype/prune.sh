for file in pruned/* ; do
    echo $file
    base_name=$(basename "$file")
    sample="${base_name%.*}"
    echo $sample

    echo "subtree --input-file $file --output-file $sample --node-ids \"\"\nexit\n" | panmat-utils /home/azhang/rotations/rotation_2/pangenome-mat/sars2k.pmat
    #echo "place $file\nexit\n" | build/panmat-utils results/sars_16k.mat > ${file}.placement.stdout 2> ${file}.placement.stderr
done