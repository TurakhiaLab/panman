# Welcome to PanMAN Wiki

## <b>What is a PanMAN?</b>
PanMAN or Pangenome Mutation-Annotated Network is a novel data representation for pangenomes that provides massive leaps in both representative power and storage efficiency. Specifically, PanMANs are composed of mutation-annotated trees, called PanMATs, which, in addition to substitutions, also annotate inferred indels (Fig. 1b), and even structural mutations (Fig. 1a) on the different branches. Multiple PanMATs are connected in the form of a network using edges to generate a PanMAN (Fig. 1c). PanMAN's representative power is compared against existing pangenomic formats in Fig. 1d. PanMANs are the most compressible pangenomic format for the different microbial datasets (SARS-CoV-2, RSV, HIV, Mycobacterium. Tuberculosis, E. Coli, and Klebsiella pneumoniae), providing 2.9 to 559-fold compression over standard pangenomic formats. 

<b>Figure 1: Overview of the PanMAN data structure</b>
<img src="images/panman.svg" width="1200" height="1200"/>

## <i><b>panmanUtils</b></i>
<i>panmanUtils</i> includes multiple algorithms to construct PanMANs and to support various functionalities to modify and extract useful information from PanMANs (Fig. 2). 

<b>Figure 2:  Overview of panmanUtils' functionalities</b>
<img src="images/utility.svg" width="1200" height="1200"/>

### <b><i>panmanUtils</i> Video Tutorial</b>
TBA


### <b>Installation</b>
panmanUtils can be installed using two different options, as described below: <br>
1. Installation script <br>
2. Docker

#### Installation Scripts
```
git clone https://github.com/TurakhiaLab/panman.git
cd panman/install
./installUbuntu.sh
```

#### Docker
```
TBA
```

### <b>Construction of PanMANs using <i>panmanUtils</i></b>
Since PanMAN can be composed of a single or multiple PanMATs, to construct a PanMAN, we start with a single tree PanMAN (or PanMAT) and then split it up into a network of multiple PanMATs using the inferred complex mutations provided as input. To construct the starting single-tree PanMAN representing a collection of sequences, panmanUtils require two inputs:<br>
1. Tree topology representing the phylogenetic relationship of the input sequences <br>
2. A pangenomic data structure (PanGraph or GFA or FASTA) representing the multiple-sequence alignment (MSA) corresponding to the sequence collection.

* Example syntax and Usage 
```
$ ./panmanUtils --pangraph-in=<path to PanGraph JSON file> --newick-in=<path to newick file> --output-file=<prefix of panman's file name>
```
```
$ ./panmanUtils --pangraph-in=ecoli_10.json --newick-in=ecoli_10.nwk --output-file=ecoli_10
```
Similarly, if you'd like to construct a PanMAN using a GFA or an MSA as an input, use the `--gfa-in` or the `--msa-in` options instead of `--pangraph-in`.
> **NOTE:** Currently, we only support GFAv1.1 consisting of Segments, un-overlapping Links and Paths.

### <b>Functionalities in <i>panmanUtils</i></b>
All panmanUtils functionality commands manipulate the input PanMAN file.
```
$ ./panmanUtils -I <path to PanMAN file> {opt}
```

Specific options:
```
-s [ --summary ]            Print PanMAN summary
-t [ --newick ]             Print newick string of all trees in a PanMAN
-f [ --fasta ]              Print tip/internal sequences (FASTA format)
-m [ --fasta-aligned ]      Print MSA of sequences for each PanMAT in a PanMAN (FASTA format)
-b [ --subnet ]             Extract subnet of given PanMAN to a new PanMAN file based on the list of nodes provided in the input-file
-v [ --vcf ]                Print variations of all sequences from any PanMAT in a PanMAN (VCF format)
-g [ --gfa ]                Convert any PanMAT in a PanMAN to a GFA file
-w [ --maf ]                Print m-WGA for each PanMAT in a PanMAN (MAF format)
-a [ --annotate ]           Annotate nodes of the input PanMAN based on the list provided in the input-file
-r [ --reroot ]             Reroot a PanMAT in a PanMAN based on the input sequence id (--reference)
-v [ --aa-translation ]     Extract amino acid translations in tsv file
-n [ --reference ] arg      Identifier of reference sequence for PanMAN construction (optional), VCF extract (required), or reroot (required)
-s [ --start ] arg          Start coordinate of protein translation
-e [ --end ] arg            End coordinate of protein translation
-i [ --input-file ] arg     Path to the input file, required for --subnet and --annotate
-o [ --output-file ] arg    Prefix of the output file name
```

> **NOTE:** When output-file argument is optional and is not provided to <i>panmanUtils</i>, the output will be printed in the terminal.

#### Summary extract
The summary feature extracts node and tree level statistics of a PanMAN, that contains a summary of its geometric and parsimony information.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --summary --output-file=<prefix of output file> (optional)
```
```
$ ./panmanUtils -I ecoli_10.panman --summary --output-file=ecoli_10
```

#### Newick extract
Extract Newick string of all trees in a PanMAN.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --newick --output-file=<prefix of output file> (optional)
```
```
$ ./panmanUtils -I ecoli_10.panman --newick --output-file=ecoli_10
```

#### Tip/internal node sequences extract
Extract tip and internal node sequences from a PanMAN in a FASTA format.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --fasta --output-file=<prefix of output file> (optional)
```
```
$ ./panmanUtils -I ecoli_10.panman --fasta --output-file=ecoli_10
```

#### Multiple Sequence Alignment (MSA) extract
Extract MSA of sequences for each PanMAT in a PanMAN in a FASTA format.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --fasta-aligned --output-file=<prefix of output file> (optional)
```
```
$ ./panmanUtils -I ecoli_10.panman --fasta-aligned --output-file=ecoli_10
```

#### Multiple Whole Genome Alignment (m-WGA) extract
Extract m-WGA for each PanMAT in a PanMAN in the form of a UCSC multiple alignment format (MAF).

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --maf --output-file=<prefix of output file> (optional)
```
```
$ ./panmanUtils -I ecoli_10.panman --maf --output-file=ecoli_10
```

#### Variant Call Format (VCF) extract
Extract variations of all sequences from any PanMAT in a PanMAN in the form of a VCF file with respect to <i>any</i> reference sequence (ref) in the PanMAT.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --vcf -reference=ref --output-file=<prefix of output file> (optional) 
```
```
$ ./panmanUtils -I ecoli_10.panman --vcf -reference=NC_000913.3 --output-file=ecoli_10 
```

#### Graphical fragment assembly (GFA) extract
Convert any PanMAT in a PanMAN to a Graphical fragment assembly (GFA) file representing the pangenome.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --gfa --output-file=<prefix of output file> (optional)
```
```
$ ./panmanUtils -I ecoli_10.panman --gfa --output-file=ecoli_10 
```

#### Subnetwork extract
Extract a subnetwork from given PanMAN and write it to a new PanMAN file based on the list of nodes provided in the input-file.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --subnet --input-file=<path to a file containing list of nodes> --output-file=<prefix of output file>
```
```
$ ./panmanUtils -I ecoli_10.panman --subnet --input-file=nodes.txt --output-file=ecoli_10_subnet
```

#### Annotate
Annotate nodes in a PanMAN with a custom string, later searched by these annotations, using a input tsv file containting list of nodes and their corresponding custom annotations. 

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --annotate <path to file containing list of annotations> --output-file=ecoli_10_annotate
```
```
$ ./panmanUtils -I ecoli_10.panman --annotate --input-file=annotations.tsv --output-file=ecoli_10_annotate
```
> **NOTE:** If output-file is not provided to <i>panmanUtils</i>, the annotated PanMAN will be written to the same file.

#### Amino Acid Translation
Extract amino acid translations from a PanMAN in tsv file.

* Example syntax and Usage
```
$ ./panmanUtils -I <path to PanMAN file> --aa-translations --output-file=<prefix of output file> (optional)
```
```
$ ./panmanUtils -I ecoli_10.panman --aa-translations --output_file=ecoli_10
```

## <b>Contributions</b>
We welcome contributions from the community to enhance the capabilities of PanMAN and panmanUtils. If you encounter any issues or have suggestions for improvement, please open an issue on [PanMAN GitHub page](https://github.com/TurakhiaLab/panman). For general inquiries and support, reach out to our team.

## <b>Citing PanMAN</b>
If you use the PanMANs or panmanUtils in your research or publications, we kindly request that you cite the following paper: XXX
