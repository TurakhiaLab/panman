# Exploring utilities in <i>panmanUtils</i>

Here, we will learn to use exploit various functionalities provided in <i>panmanUtils</i> software for downstream applications in epidemiological, microbiological, metagenomic, ecological, and evolutionary studies.

**Step 0:** The Steps below require panmanUtils and a PanMAN. If not done so far, refer to [installation guide](install.md) to install panmanUtils and [construction](construction.md) instructions to build a PanMAN. Alternatively, users can download pre-built PanMANs using the following command
```bash
# Assuming $PANMAN directs to the panman repository directory
cd $PANMAN_HOME
mkdir -p build/panman && cd build/panman
ToDO
```

### Functionalities in <i>panmanUtils</i>
All panmanUtils functionality commands manipulate the input PanMAN file.
```bash
cd $PANMAN_HOME/build
./panmanUtils -I <path to PanMAN file> {opt}
```



| **Option**                       | **Description**                                                                                                   |
|----------------------------------|-------------------------------------------------------------------------------------------------------------------| 
|`-I`, `--input-panman`            | Input PanMAN file path                                                                                            |
| `-s`, `--summary`                | Print PanMAN summary                                                                                              |
| `-t`, `--newick`                 | Print Newick string of all trees in a PanMAN                                                                      |
| `-f`, `--fasta`                  | Print tip/internal sequences (FASTA format)                                                                       |
| `-m`, `--fasta-aligned`          | Print MSA of sequences for each PanMAT in a PanMAN (FASTA format)                                                 |
| `-b`, `--subnet`                 | Extract subnet of given PanMAN to a new PanMAN file based on the list of nodes provided in the input file         |
| `-v`, `--vcf`                    | Print variations of all sequences from any PanMAT in a PanMAN (VCF format)                                        |
| `-g`, `--gfa`                    | Convert any PanMAT in a PanMAN to a GFA file                                                                      |
| `-w`, `--maf`                    | Print m-WGA for each PanMAT in a PanMAN (MAF format)                                                              |
| `-a`, `--annotate`               | Annotate nodes of the input PanMAN based on the list provided in the input file                                   |
| `-r`, `--reroot`                 | Reroot a PanMAT in a PanMAN based on the input sequence id (`--reference`)                                        | 
| `-v`, `--aa-translation`         | Extract amino acid translations in tsv file                                                                       | 
| `-e`, `--extended-newick`        | Print PanMAN's network in extended-newick format                                                                  |
| `-k`, `--create-network`         | Create PanMAN with network of trees from single or multiple PanMAN files                                          |
| `-p`, `--printMutations`         | Create PanMAN with network of trees from single or multiple PanMAN files                                          |
| `-q`, `--acr`                    | ACR method `[fitch(default), mppa]`                                                                                 |
| `-n`, `--reference`              | Identifier of reference sequence for PanMAN construction (optional), VCF extract (required), or reroot (required) | 
| `-s`, `--start`                  | Start coordinate of protein translation                                                                           | 
| `-e`, `--end`                    | End coordinate of protein translation                                                                             |
| `-d`, `--treeID`                 | Tree ID, required for `--vcf`                                                                                     |
| `-i`, `--input-file`             | Path to the input file, required for `--subnet`, `--annotate`, and `--create-network`                             |
| `-o`, `--output-file`            | Prefix of the output file name                                                                                    |


> **NOTE:** When output-file argument is optional and is not provided to <i>panmanUtils</i>, the output will be printed in the terminal.

!!!Note
    For all the examples below, `sars_20.panman` will be used as input panman. Alternatively, users can provide custom build panman using the instructions provided [here](construction.md).

#### Summary extract
The summary feature extracts node and tree level statistics of a PanMAN, that contains a summary of its geometric and parsimony information.

* Usage Syntax
```bash
./panmanUtils -I <path to PanMAN file> --summary --output-file=<prefix of output file> (optional)
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman  --summary --output-file=sars_20
```

#### Newick extract
Extract Newick string of all trees in a PanMAN.

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --newick --output-file=<prefix of output file> (optional)
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --newick --output-file=sars_20
```

#### Extended Newick extract
Extract network in Extended Newick format.

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> ----extended-newick --output-file=<prefix of output file> (optional)
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman ----extended-newick --output-file=sars_20
```

#### Tip/internal node sequences extract
Extract tip and internal node sequences from a PanMAN in a FASTA format.

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --fasta --output-file=<prefix of output file> (optional)
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --fasta --output-file=sars_20
```

#### Multiple Sequence Alignment (MSA) extract
Extract MSA of sequences for each PanMAT (with pseduo-root coordinates) in a PanMAN in a FASTA format.

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --fasta-aligned --output-file=<prefix of output file> (optional)
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --fasta-aligned --output-file=sars_20
```

#### Multiple Whole Genome Alignment (m-WGA) extract
Extract m-WGA for each PanMAT in a PanMAN in the form of a UCSC multiple alignment format (MAF).

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --maf --output-file=<prefix of output file> (optional)
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --maf --output-file=sars_20
```

#### Variant Call Format (VCF) extract
Extract variations of all sequences from any PanMAT in a PanMAN in the form of a VCF file with respect to <i>any</i> reference sequence (ref) in the PanMAT.

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --vcf -reference=ref --output-file=<prefix of output file> (optional) 
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --vcf -reference="Switzerland/SO-ETHZ-500145/2020|OU000199.2|2020-11-12" --output-file=sars_20 
```

#### Graphical fragment assembly (GFA) extract
Convert any PanMAT in a PanMAN to a Graphical fragment assembly (GFA) file representing the pangenome.

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --gfa --output-file=<prefix of output file> (optional)
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --gfa --output-file=sars_20 
```

#### Subnetwork extract
Extract a subnetwork from a given PanMAN and write it to a new PanMAN file based on the list of nodes provided in the input-file.

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --subnet --input-file=<path to a file containing list of nodes> --output-file=<prefix of output file>
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --subnet --input-file=nodes.txt --output-file=ecoli_10_subnet
```

#### Annotate
Annotate nodes in a PanMAN with a custom string, later searched by these annotations, using an input TSV file containing a list of nodes and their corresponding custom annotations. 

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --annotate --input-file=<path to file containing list of annotations> --output-file=ecoli_10_annotate
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --annotate --input-file=annotations.tsv --output-file=ecoli_10_annotate
```
> **NOTE:** If output-file is not provided to <i>panmanUtils</i>, the annotated PanMAN will be written to the same file.

#### Amino Acid Translation
Extract amino acid translations from a PanMAN in TSV file.

* Usage syntax
```bash
./panmanUtils -I <path to PanMAN file> --aa-translations --output-file=<prefix of output file> (optional)
```
* Example
```bash
cd $PANMAN_HOME/build
./panmanUtils -I panman/sars_20.panman --aa-translations --output_file=sars_20
```
