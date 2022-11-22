# pangenome-mat

Updated version of MAT from [Usher]

## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake -D Protobuf_PROTOC_EXECUTABLE=/usr/bin/protoc ..
```

## Run Instructions
```
cd build
./panmat-utils ../results/<mat file name>
```

## MAT Summary (Example)

```
> summary
Total Nodes in Tree: 187
Total Samples in Tree: 94
Total Substitutions: 146453
Total Insertions: 93396
Total Deletions: 101687
Total SNP Substitutions: 1815773
Total SNP Insertions: 24660
Total SNP Deletions: 10597
Max Tree Depth: 7
Mean Tree Depth: 6.6383

Summary creation time: 57906170
```

## MAT Writer
Writes MAT to protobuf format in the `build/pmat` directory
```
> write <filename>
```

## Newick
Prints newick string of tree
```
> newick
```

## Subtree Extract (Protobuf)
Extracts consolidated version of subtree containing given nodes and writes it to `build/pmat` directory in protobuf format. Supports reading identifier list from file.
```
> subtree --input=<input_file_name> <filename>
> subtree <filename> <node identifier 1> <node identifier 2> <node identifier 3> ...
```

## Subtree Extract (Newick)
Extracts consolidated version of subtree containing given nodes and writes its newick string to `build/pmat` directory. Supports reading identifier list from file.
```
> subtree-newick --input=<input_file_name> <filename>
> subtree-newick <filename> <node identifier 1> <node identifier 2> <node identifier 3> ...
```

## FASTA Writer (In development)
Extracts sequences from the tree and writes them in FASTA format in the `build/fasta` directory. Currently doesn't incorporate non-SNP indels.
```
> fasta [--aligned] <filename>
```


   [Usher]: <https://github.com/yatisht/usher>
