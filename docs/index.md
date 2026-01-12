# <b>Welcome to PanMAN Wiki</b>
<div align="center">
    <img src="images/logo.png"/>
</div>

## <b>Introduction</b> 
### What are PanMANs?
<iframe width="1000" height="600" src="https://www.youtube.com/embed/VZFIg7x9tVw" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
PanMAN or Pangenome Mutation-Annotated Network is a novel data representation for pangenomes that provides massive leaps in both representative power and storage efficiency. Specifically, PanMANs are composed of mutation-annotated trees, called PanMATs, which, in addition to substitutions, also annotate inferred indels (Fig. 2b), and even structural mutations (Fig. 2a) on the different branches. Multiple PanMATs are connected in the form of a network using edges to generate a PanMAN (Fig. 2c). PanMAN's representative power is compared against existing pangenomic formats in Fig. 1. PanMANs are the most compressible pangenomic format for the different microbial datasets (SARS-CoV-2, RSV, HIV, Mycobacterium. Tuberculosis, E. Coli, and Klebsiella pneumoniae), providing 2.9 to 559-fold compression over standard pangenomic formats. 

<div align="center">
    <img src="images/figure1.svg" width="600" height="600"/><br>
    <b>Figure 1: Comparison of representative power of PanMAN against other pangenomic formats (yellow ticks indicate partial representative ability)</b><br>
</div>
<br>
<br>
<div align="center">
    <img src="images/panman.svg" width="600" height="600"/><br>
    <b>Figure 2: Overview of the PanMAN data structure</b><br>
</div>




### PanMAN's Protocol Buffer file format
PanMAN utilizes Google’s protocol buffer (protobuf, [https://protobuf.dev/](https://protobuf.dev/)), a binary serialization file format, to compactly store PanMAN's data structure in a file. Fig. 3 provides the .proto file defining the PanMAN’s structure. At the top level, the file format of PanMANs encodes a list (declared as a repeated identifier in the .protof file) of PanMATs. Each PanMAT object stores the following data elements: (a) a unique identifier, (b) a phylogenetic tree stored as a string in Newick format, (c) a list of mutations on each branch ordered according to the pre-order traversal of the tree topology, (d) a block mapping object to record homologous segments identified as duplications and rearrangements, which are mapped against their common consensus sequence; the block-mapping object is also used to derive the pseudo-root, e) a gap list to store the position and length of gaps corresponding to each block's consensus sequence. Each mutation object encodes the node's block and nucleotide mutations that are inferred on the branches leading to that node. If a block mutation exists at a position described by the Block-ID field (int32), the block mutation field (bool) is set to 1, otherwise set to 0, and its type is stored as a substitution to and from a gap in Block mutation type field (bool), encoded as 0 or 1, respectively. In PanMAN, each nucleotide mutation within a block inferred on a branch has four pieces of information, i.e., position (middle coordinate), gap position (last coordinate), mutation type, and mutated characters. To reduce redundancy in the file, consecutive mutations of the same type are packed together and stored as a mutation info (int32) field, where mutation type, mutation length, and mutated characters use 3, 5, and 24 bits, respectively. PanMAN stores each character using one-hot encoding, hence, one "Nucleotide Mutations" object can store up to 6 consecutive mutations of the same type. PanMAN's file also stores the complex mutation object to encode the type of complex mutation and its metadata such as PanMATs' and nodes' identifiers, breakpoint coordinates, etc. The entire file is then compressed using XZ ([https://github.com/tukaani-project/xz](https://github.com/tukaani-project/xz)) to enhance storage efficiency.

<div align="center">
    <img src="images/pb.svg" width="600" height="600"/><br>
    <b>Figure 3: PanMAN's file format</b>
</div>

### <i>panmanUtils</i>
<i>panmanUtils</i> includes multiple algorithms to construct PanMANs and to support various functionalities to modify and extract useful information from PanMANs (Fig. 4). 

<div align="center">
    <img src="images/utility.svg" width="600" height="600"/><br>
    <b>Figure 4:  Overview of panmanUtils' functionalities</b>
</div>

### Video Tutorial
<iframe width="1000" height="600" src="https://www.youtube.com/embed/eh9zQElrmLI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

<a name="install"></a>
## <b><i>panmanUtils</i> Installation Methods</b>

<b><i>panmanUtils</i></b> software can be installed using four different methods:
 
1. Conda (Recommended) 
2. Docker Image
3. Dockerfile
4. Installation scripts 

### 1. Using conda (recommended)
Users can install <i>panmanUtils</i> through installation of [panman conda package](https://bioconda.github.io/recipes/panman/README.html#package-package%20&#x27;panman&#x27;), compatible with `linux-64` and `osx-64`. For modern macs using Apple silicon (arm64), you need to install Rosetta 2.

#### i. Dependencies
1. [Conda](https://docs.conda.io/en/latest/)
#### ii. Install panman conda package 
```
# Create and activate a new environment for panman
conda create -n panman-env python=3.11 -y
conda activate panman-env

# Set up channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# On macOS ARM:
# conda config --env --set subdir osx-64

# Install the panman package
conda install panman -y
```
#### iii. Run <i>panmanUtils</i> 
```
panmanUtils --help
```
### 2. Using Docker Image

To use <i>panmanUtils</i> in a docker container, users can create a docker container from a docker image, by following these steps (compatible with `linux-64` and `osx-64`).

#### i. Dependencies
1. [Docker](https://docs.docker.com/engine/install/)
#### ii. Pull and build the PanMAN docker image from DockerHub
```bash
## Note: If the Docker image already exists locally, make sure to pull the latest version using 
## docker pull swalia14/panman:latest

## If the Docker image does not exist locally, the following command will pull and run the latest version
docker run -it swalia14/panman:latest
```
#### iii. Run <i>panmanUtils</i>
```bash
# Insider docker container
panmanUtils --help
```
!!!Note
    <b>The docker image comes with preinstalled <i>panmanUtils</i> and other tools such as PanGraph, PGGB, and RIVET.</b>

### 3. Using DockerFile
Docker container with preinstalled <i>panmanUtils</i> can also be built from DockerFile by following these steps (compatible with `linux-64` and `osx-64`).

#### i. Dependencies
1. [Docker](https://docs.docker.com/engine/install/)
2. [Git](https://git-scm.com/downloads)
#### ii. Clone the repository and build a docker image
```bash
git clone https://github.com/TurakhiaLab/panman.git
cd panman/docker
docker build -t panman .
```
#### iii. Build and run the docker container
```bash
docker run -it panman
```
#### iv. Run <i>panmanUtils</i>
```bash
# Insider docker container
panmanUtils --help
```

### 4. Using installation script (Least recommended)
We provide scripts to install panmanUtils from source code (requires `sudo` access, compatible with `Linux` only). `Mac` users can use [MacOS specific installation script](https://github.com/TurakhiaLab/panman/blob/main/install/installationMacOS.sh), that uses `conda` to install <i>panmanUtils</i>.
#### 1. Dependencies
1. [Git](https://git-scm.com/downloads)

#### 2. Clone the repository
```bash
git clone https://github.com/TurakhiaLab/panman.git
cd panman
```
#### 3. Run the installation script
```bash
chmod +x install/installationUbuntu.sh
./install/installationUbuntu.sh
```
#### 4. Run <i>panmanUtils</i>
```bash
cd build
./panmanUtils --help
```
!!!Note
    <b><i>panmanUtils</i> is built using CMake and depends upon libraries such as Boost, cap'n proto, etc, which are also installed in `installationUbuntu.sh`. If users face version issues, try using the conda or docker methods detailed above.</b>

<a name="construction"></a>
## <b>PanMAN Construction</b>

Here, we will learn to build PanMAN from various input formats.

**Step 0:** The Steps below require <i>panmanUtils</i>, if not done so far, refer to [installation guide](#install) to install <i>panmanUtils</i>. To check if <i>panmanUtils</i> is properly installed or not, run the following command, and it should execute without error
```bash
panmanUtils --help
```
### Building PanMAN from Alignments (PanGraph/GFA/MSA)
#### Building PanMAN from PanGraph
**Step 1:** Check if `sars_20.json` and `sars_20.nwk` files exist in `test` directory. Alternatively, users can provide custom PanGraph (JSON) and tree topology (Newick format) files to build a panman. 

**Step 2:** Run <i>panmanUtils</i> with the following command to build a panman from PanGraph:

```bash
panmanUtils -P $PANMAN_HOME/test/sars_20.json -N $PANMAN_HOME/test/sars_20.nwk -o sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

#### Building PanMAN from GFA

**Step 1:** Check if `sars_20.gfa` and `sars_20.nwk` files exist in `test` directory. Alternatively, users can provide custom GFA and tree topology (Newick format) files to build a panman. 

**Step 2:** Run <i>panmanUtils</i> with the following command to build a panman from GFA:

```bash
panmanUtils -G $PANMAN_HOME/test/sars_20.gfa -N $PANMAN_HOME/test/sars_20.nwk -o sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

#### Building PanMAN from MSA (FASTA format)

**Step 1:** Check if `sars_20.msa` and `sars_20.nwk` files exist in `test` directory. Alternatively, users can provide custom MSA (FASTA format) and tree topology (Newick format) files to build a panman. 

**Step 2:** Run <i>panmanUtils</i> to build a panman from GFA using the following commands:

```bash
panmanUtils -M $PANMAN_HOME/test/sars_20.msa -N $PANMAN_HOME/test/sars_20.nwk -o sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

### Building PanMAN from raw genome sequences or fragment assemblies using Snakemake Workflow
We provide a Snakemake workflow to construct PanMANs from raw sequences (FASTA format) or from fragment assemblies.

!!!Note
    <b>The Snakemake workflow uses various tools such as PanGraph tool, PGGB, MAFFT, and MashTree to build input PanGraph, GFA, MSA, and Tree topology files, respectively and it is particularly designed to be used in the docker container build from either the provided docker image or the DockerFile (instructions provided [here](#install)).</b>

#### Building PanMAN from raw genome sequences
**Step 1:** Run the following command to construct a panman from raw sequences.

```bash
cd $PANMAN_HOME/workflows
snakemake --use-conda --cores 8 --config RUNTYPE="pangraph/gfa/msa" FASTA="[user_input]" SEQ_COUNT="Number of sequences" ASSEM="NONE" REF="NONE" TARGET="NONE"
```

#### Building PanMAN from fragment assemblies
**Step 1:** Run the following command to construct a panman from fragment assemblies.

```bash 
cd $PANMAN_HOME/workflows
snakemake --use-conda --cores 8 --config RUNTYPE="pangraph/gfa/msa" FASTA="None" SEQ_COUNT="Number of sequences" ASSEM="frag" REF="reference_file" TARGET="target.txt"
```
Here, target.txt includes a list of files that contain the fragmented assemblies.

## <b>Exploring utilities in <i>panmanUtils</i></b>

Here, we will learn to use various functionalities provided in <i>panmanUtils</i> software for downstream applications in epidemiological, microbiological, metagenomic, ecological, and evolutionary studies.

**Step 0:** The Steps below require panmanUtils and a PanMAN. We provide a pre-built panman (`sars_20.panman`), otherwise, refer to [installation guide](#install) to install panmanUtils and [construction](#construction) instructions to build a PanMAN. 
<!-- ```bash
# Assuming $PANMAN directs to the panman repository directory
cd $PANMAN_HOME
mkdir -p build/panman && cd build/panman
ToDO
``` -->

### Functionalities in <i>panmanUtils</i>
All panmanUtils functionality commands manipulate the input PanMAN file.
```bash
cd $PANMAN_HOME/build
panmanUtils -I <path to PanMAN file> {opt}
```
<a name="table1"></a>
<div align="center"> <b>Table 1:</b> List of functionalities supported by <i>panmanUtils</i> </div>

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
| `-V`, `--version`            | panmanUtils Version                                     |



> **Important:** When output-file argument is optional and is not provided to <i>panmanUtils</i>, the output will be printed in the terminal.

!!!Note
    <b>For all the examples below, `sars_20.panman` will be used as input panman. Alternatively, users can provide custom build panman using the instructions provided [here](#construction).</b>
!!! Note
    <b>Users can reduce memory consumption by lowering the number of CPU threads (default set to 32) through the --threads option in panmanUtils, at a cost of higher latency.</b>

#### Summary extract
The summary feature extracts node and tree level statistics of a PanMAN, that contains a summary of its geometric and parsimony information.

* Usage Syntax
```bash
panmanUtils -I <path to PanMAN file> --summary --output-file=<prefix of output file> (optional)
```
* Example
```bash
panmanUtils -I panman/sars_20.panman  --summary --output-file=sars_20
```

#### Newick extract
Extract Newick string of all trees in a PanMAN.

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --newick --output-file=<prefix of output file> (optional)
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --newick --output-file=sars_20
```

#### Extended Newick extract
Extract network in Extended Newick format.

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --extended-newick --output-file=<prefix of output file> (optional)
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --extended-newick --output-file=sars_20
```

#### Tip/internal node sequences extract
Extract tip and internal node sequences from a PanMAN in a FASTA format.

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --fasta --output-file=<prefix of output file> (optional)
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --fasta --output-file=sars_20
```

#### Multiple Sequence Alignment (MSA) extract
Extract the MSA of sequences for each PanMAT (with pseudo-root  coordinates) in a PanMAN in a FASTA format.

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --fasta-aligned --output-file=<prefix of output file> (optional)
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --fasta-aligned --output-file=sars_20
```

#### Multiple Whole Genome Alignment (m-WGA) extract
Extract m-WGA for each PanMAT in a PanMAN in the form of a UCSC multiple alignment format (MAF).

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --maf --output-file=<prefix of output file> (optional)
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --maf --output-file=sars_20
```

#### Variant Call Format (VCF) extract
Extract variations of all sequences from any PanMAT in a PanMAN in the form of a VCF file with respect to <i>any</i> reference sequence (ref) in the PanMAT.

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --vcf -reference=ref --output-file=<prefix of output file> (optional) 
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --vcf -reference="Switzerland/SO-ETHZ-500145/2020|OU000199.2|2020-11-12" --output-file=sars_20 
```

#### Graphical fragment assembly (GFA) extract
Convert any PanMAT in a PanMAN to a Graphical fragment assembly (GFA) file representing the pangenome.

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --gfa --output-file=<prefix of output file> (optional)
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --gfa --output-file=sars_20 
```

#### Subnetwork extract
Extract a subnetwork from a given PanMAN and write it to a new PanMAN file based on the list of nodes provided in the input-file.

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --subnet --input-file=<path to a file containing list of nodes> --output-file=<prefix of output file>
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --subnet --input-file=nodes.txt --output-file=sars_20_subnet
```

#### Annotate
Annotate nodes in a PanMAN with a custom string, later searched by these annotations, using an input TSV file containing a list of nodes and their corresponding custom annotations. 

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --annotate --input-file=<path to file containing list of annotations> --output-file=sars_20_annotate
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --annotate --input-file=annotations.tsv --output-file=sars_20_annotate
```
> **NOTE:** If output-file is not provided to <i>panmanUtils</i>, the annotated PanMAN will be written to the same file.

#### Amino Acid Translation
Extract amino acid translations from a PanMAN in a TSV file.

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --aa-translations --output-file=<prefix of output file> (optional)
```
* Example
```bash
panmanUtils -I panman/sars_20.panman --aa-translations --output_file=sars_20
```

#### Range Query
<i>panmanUtils</i> allow extracting alignment of all the sequences of a single PanMAT in a PanMAN (FASTA format) with respect to a user-defined reference sequence between positions [start:end]

* Usage syntax
```bash
panmanUtils -I <path to PanMAN file> --index no -x start -y end --reference=<ref sequence name>
```
* Example
```bash
panmanUtils -I <path to PanMAN file> --index no -x 10 -y 100 --reference="Switzerland/SO-ETHZ-500145/2020|OU000199.2|2020-11-12"
```

### <i>panmanUtils</i> Interactive mode
**Step 1:** Users can enter <i>panmanUtils</i>'s interactive mode by passing input panman as input using the following command:

```bash
panmanUtils -I <path to PanMAN file>
## Example
panmanUtils -I panman/sars_20.panman
```

!!! Note
    <b>The interactive mode should look like the image attached below</b>

 ![Interactive Mode](images/interactiveMode.png)

**Step 2:** Use the commands listed in [Table 1](#table1) to perform desired operation

## <b>Contributions</b>
We welcome contributions from the community to enhance the capabilities of PanMAN and panmanUtils. If you encounter any issues or have suggestions for improvement, please open an issue on [PanMAN GitHub page](https://github.com/TurakhiaLab/panman/issues). For general inquiries and support, reach out to our team.

## <b>Citing PanMAN</b>
If you use the PanMANs or <i>panmanUtils</i> in your research or publications, we kindly request that you cite the following paper:
* Walia S., Motwani H., Tseng YH., Smith K., Corbett-Detig R., Turakhia Y., <i>Compressive pangenomics using mutation-annotated networks</i>. <b>Nat Genet</b> (2026). [https://doi.org/10.1038/s41588-025-02478-7](https://doi.org/10.1038/s41588-025-02478-7)
