for file in fastq/*_R1.fastq ; do
    echo "$file"
    echo "place $file\nexit\n" | build/panmat-utils results/sars_16k.mat > ${file}.placement.stdout 2> ${file}.placement.stderr
done