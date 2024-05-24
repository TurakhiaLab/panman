# Welcome to PanMAN Wiki

## What is a PanMAN?
PanMAN or Pangenome Mutation-Annotated Network is a novel data representation for pangenomes that provides massive leaps in both representative power and storage efficiency. Specifically, PanMANs are composed of mutation-annotated trees, called PanMATs, which, in addition to substitutions, also annotate inferred insertions, deletions, and even structural mutations on the different branches. Multiple PanMATs are connected in the form of a network using edges to generate a PanMAN. PanMANs are the most compressible pangenomic format for the different microbial datasets (SARS-CoV-2, RSV, HIV, Mycobacterium. Tuberculosis, E. Coli, and Klebsiella pneumoniae), providing 2.9 to 559-fold compression over standard pangenomic formats. 

## panmanUtils
panmanUtils includes multiple algorithms to construct PanMANs and to support various functionalities to modify and extract useful information from PanMANs. 

### Installation
panmanUtils can be installed using two different options (i) Installation script and (ii) Docker, as described below.

#### Installation Scripts
```
git clone https://github.com/TurakhiaLab/panman.git
cd panman/install
./installUbuntu.sh
```

#### Docker
```
ToDo
```

### Construction of PanMANs using panmanUtils
Since PanMAN can be composed of a single or multiple PanMATs, to construct a PanMAN, we start with a single tree PanMAN (or PanMAT) and then split it up into a network of multiple PanMATs using the inferred complex mutations provided as input. To construct the starting single-tree PanMAN representing a collection of sequences, panmanUtils require two inputs:
1. Tree topology representing the phylogenetic relationship of the input sequences
2. A pangenomic data structure (PanGraph or GFA or FASTA) representing the multiple-sequence alignment (MSA) corresponding to the sequence collection.

Construct a single tree PanMAN from PanGraph

```
$ ./panmanUtils --pangraph-in=<path to PanGraph JSON file> --newick-in=<path to newick file>
```
Similarly, if you'd like to construct a PanMAN using a GFA or an MSA as an input, use the `--gfa-in` or the `--msa-in` options instead of `--pangraph-in`.


> **NOTE:** Currently, we only support GFAv1.1 consisting of Segments, un-overlapping Links and Paths.

### Functionalities in panmanUtils
All panmanUtils functionality commands manipulates input PanMAN file (.panman).
```
$ ./panmanUtils -I <path to PanMAN file> {opt}
```

Specific options:
```
--summary: Prints summary of a PanMAN.
--fasta: Prints tip/internal sequences (FASTA format)
--fasta-aligned: Prints MSA of sequences for each PanMAT in a PanMAN (FASTA format)
--maf: Prints m-WGA for each PanMAT in a PanMAN (MAF format)
--vcf: Prints variations of all sequences from any PanMAT in a PanMAN (VCF format)
--gfa: Converts any PanMAT in a PanMAN to a GFA file
-o: path to output file
--threads: Number of threads to use when possible (Default: use all available cores).
--help: Print help messages.
```

#### Summary extract
The summary feature extract node and tree level statistics of a PanMAN, which contains a summary of its geometric and parsimony information.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --summary
```
```
$ ./panmanUtils -I ecoli_10.panman --summary
```

#### Tip/internal node sequences extract
Extract tip and internal node sequences from a PanMAN in a FASTA format.
* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --fasta -o <path to output file>
```
```
$ ./panmanUtils -I ecoli_10.panman --fasta -o ecoli_10.fasta
```

#### Multiple Sequence Alignment (MSA) extract
Extract MSA of sequences for each PanMAT in a PanMAN in a FASTA format.
* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --fasta-aligned -o <path to output file>
```
```
$ ./panmanUtils -I ecoli_10.panman --fasta-aligned -o ecoli_10.msa
```

#### Multiple Whole Genome Alignment (m-WGA) extract
Extract m-WGA for each PanMAT in a PanMAN in the form of a UCSC multiple alignment format (MAF).
* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --maf -o <path to output file>
```
```
$ ./panmanUtils -I ecoli_10.panman --maf -o ecoli_10.maf
```

#### Variant Call Format (VCF) extract
Extract variations of all sequences from any PanMAT in a PanMAN in the form of a VCF file with respect to <i>any</i> reference sequence (ref) in the PanMAT.
* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --vcf -o <path to output file> -r ref
```
```
$ ./panmanUtils -I ecoli_10.panman --vcf -o ecoli_10.vcf -r S5
```

#### Graphical fragment assembly (GFA) extract
Convert any PanMAT in a PanMAN to a Graphical fragment assembly (GFA) file representing the pangenome.
* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --gfa -o <path to output file>
```
```
$ ./panmanUtils -I ecoli_10.panman --gfa -o ecoli_10.gfa
```

#### Annotate
Annotate any node in a PanMAN with a custom string (nodes can be later searched by these annotations)
* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --annotate <path to file containing list of annotations>
```
```
$ ./panmanUtils -I ecoli_10.panman --annotate info.txt
```



