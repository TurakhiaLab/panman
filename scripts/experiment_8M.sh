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

twilight -f $T_OUT --gap-ends 0 -p y -v -r 1 --trans -8 -o sars_8M.msa --gap-open -60
grep -a -A 1 Wuhan-Hu-1 sars_8M.msa > wuhan.aln
panmanUtils --input-msa sars_8M.msa --input-newick SARS_8M/sars_8M.nwk -o sars_8M --low-mem-mode --refFile wuhan.aln

panmanUtils -I panman/sars_8M.panman --summary
panmanUtils -I panman/sars_8M.panman --fasta
panmanUtils -I panman/sars_8M.panman --vcf

panmanUtils --create-network panman/sars_8M.panman --input-file SARS_8M/sars_8M_cmut -o sars_8M_network
panmanUtils -I panman/sars_8M_network.panman --summary