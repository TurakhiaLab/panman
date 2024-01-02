# PanMAT Utilities (panmat-utils)
panmat-utils is a suite of tools used to analyze and manipulate pangenome mutation annotated tree (.pmat) files. Use the following command to load the PanMAT file into the utility:
```
$ ./panmat-utils --panmat-in=<path to PanMAT file>
```
## 1. Extract Summary of a PanMAT
panmat-utils _summary_ can be used to get geometric and parsimonious details of the tree. This summary can be obtained by simply using the _summary_ command.
```
> summary
Total Nodes in Tree: 3
Total Samples in Tree: 2
Total Substitutions: 6910
Total Insertions: 8773
Total Deletions: 5120
Total Inversions: 1487
Total SNP Substitutions: 60198
Total SNP Insertions: 1400
Total SNP Deletions: 1014
Max Tree Depth: 1
Mean Tree Depth: 1

Summary creation time: 2952753 nanoseconds
```

## 2. Extracting Sequences from a PanMAT

Sequences represented by the PanMAT can be extracted from it using the PanMAT utility.

### 2.1 Extract all Raw Sequences

If you wish to extract all the raw sequences represented by the PanMAT, you can do so using the _fasta_ feature of the PanMAT.

```
> fasta --output-file=<output file name>
```

This will create a FASTA file called `<output file name>.fasta` in the `fasta` subdirectory consisting of all the raw sequences.

### 2.2 Extracting Selected sequences
You can also extract selected sequences by listing their identifiers using _sequences_ command:
```
> sequences --list=<space separated list of sequences> --output-file=<output file name>
```
This will also create a FASTA file called `<output file name>.fasta` in the `fasta` directory consisting the raw sequences for the given identifiers.

### 2.3 Extracting a Multiple Sequence Alignment (MSA) of all Sequences
A PanMAT inherently stores a MSA of all the sequences that it represents. This MSA can be extracted using the _fasta_ feature using the `--aligned` option.
```
fasta --aligned --output-file=<output file name>
```

## 3. Extract a Subtree from a PanMAT
You can extract the most compact subtree containing a given set of nodes from a PanMAT into another PanMAT file using the PanMAT utility. To do this, you can use the _subtree_	command of the PanMAT utility as follows:
```
> subtree --node-ids=<space_separated_list_of_node_ids> --output-file=<output_file_name>
```
This command creates a PanMAT file called `<output_file_name>.pmat` in the `pmat` subdirectory representing the nodes specified in the command. This PanMAT file can be analyzed separately using the utility.

### 3.1 Extract Subtree newick string 
If only the newick string of the minimal subtree is to be extracted, you can specify the option `--newick=true` to the above command.


## 4. Annotate and Search
You can use the PanMAT utility to annotate nodes of a PanMAT. For this, you need to create a comma separated annotations file in which the first column stores the node identifiers of the nodes that need to be annotated and the subsequent columns represent the annotations. A sample file would look like this:
```
seq_1,Winter-2020,Rome
seq_2,Fall-2020,NYC
seq_3,Winter-2020,NYC
```
This file would annotate the node representing `seq_1` with `Winter-2020` and `Rome`, `seq_2` with `Fall-2020` and `NYC` and `seq_3` with `Winter-2020` and `NYC`.

You can _annotate_ command to annotate the PanMAT as follows:
```
> annotate --input-file=<path_to_annotation_file>
```

You can also search nodes of a PanMAT by annotations using the _search_ utility as follows:
```
> search --keywords=Winter-2020 NYC
Winter-2020: seq_1; seq_3;
NYC: seq_2; seq_3;
```

## 5. Newick String
You can get the newick string representing the PanMAT by using the _newick_ command. The newick string will be printed to the console.
```
> newick
(NZ_CP007391.1,NZ_CP007265.1)node_1;
```

## 6. Extracting Other File Formats from a PanMAT
Other file formats can be extracted from a PanMAT using the PanMAT utility. These include MAF, VCF, GFA, and a custom Amino Acid translation file format.

### 6.1 Multiple Alignment Format (MAF)
The [MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) file format is a sequence alignment file format in which common sections of sequences (called *blocks*) are aligned together. The PanMAT file format inherently stores sequences split into blocks and therefore, the utility supports MAF extraction. It can be done use _maf_ command as follows:
```
> maf --output-file=<output_file_name>
```
This will create a MAF file with the name `<output_file_name>.maf` in the `maf` subdirectory.

### 6.2 VCF
A VCF file representing all the sequences in the PanMAT can be extracted with respect to **any chosen reference sequence** can be extracted using the _vcf_ command as follows:

```
> vcf --reference=<reference_sequence_id> --output-file=<output_file_name>
```
This will create a VCF file with the name `<output_fiile_name.vcf>` in the `vcf` subdirectory.

### 6.3 GFA

The PanMAT utility can also be used to extract a GFA file consisting of all the sequences in the PanMAT using _genGFA_ command as follows:
```
genGFA --output-file=<output_file_name>
```
This creates a GFA file with the name `<output_file_name>.gfa` in the `gfa` subdirectory.

### 6.4 Amino Acid Translations

You can also extract the list of amino acid translations in the sequence represented by each node in the PanMAT with respect to the sequence represented by the root between a given pair of `start` and `end` coordinates. These are returned in a custom tab separated (tsv) file format. Here's an example of a sample amino acid translations file:
```
node_id	aa_mutations
node_2	S:10:Ala;S:23:Cys;S:32:Tyr;
ON817448.1	I:34:Trp;S:57:Met;
...
```
Here, the amino acid mutations are separated by semi-colons. Mutations starting with `S`, `I`, and `D` represent substitutions, insertions, and deletions, respectively.

This file can be extracted using _aaTranslations_ command in the PanMAT utility.
```
> aaTranslations --start=<start_coordinate> --end=<end_coordinate> --output-file=<output_file_name>
```
This will create a file called `<output_file_name>.tsv` in the `aa_translations` subdirectory.

## 7. Write PanMAT to a file
The PanMAT can be written to a new file using _write_ command of the PanMAT utility as follows:
```
> write <file_name>
```
This creates a file called `<filename>.pmat` in the `pman`subdirectory. This subdirectory is intended to contain all the PanMAT files.