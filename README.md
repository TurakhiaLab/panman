[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg
[license-link]: [https://github.com/TurakhiaLab/panman/LICENSE](https://github.com/TurakhiaLab/panman/blob/main/LICENSE)
[![License][license-badge]][license-link]
[![DOI](https://img.shields.io/badge/DOI-https://zenodo.org/records/17781629-beige)](https://zenodo.org/records/17781629)
[<img src="https://img.shields.io/badge/Install with-Conda-aquamarine.svg?logo=Anaconda">](https://bioconda.github.io/recipes/panman/README.html#package-package%20&#x27;panman&#x27;)
[<img src="https://img.shields.io/badge/Install with-Docker-informational.svg?logo=Docker">](https://hub.docker.com/r/swalia14/panman)
[<img src="https://img.shields.io/badge/Published in-Nature Genetics-red">](https://doi.org/10.1038/s41588-025-02478-7)
[<img src="https://img.shields.io/badge/Build with-CMake-green.svg?logo=CMake">](https://github.com/TurakhiaLab/panman/blob/main/install/installationUbuntu.sh)
[<img src="https://img.shields.io/badge/Watch it on-Youtube-FF0000.svg?logo=YouTube">](https://www.youtube.com/watch?v=eh9zQElrmLI)
[![Build Status](https://github.com/TurakhiaLab/panman/actions/workflows/ci.yml/badge.svg)](https://github.com/TurakhiaLab/panman/actions)

<!-- [<img src="https://img.shields.io/badge/Made with-Snakemake-aquamarine.svg?logo=snakemake">](https://snakemake.readthedocs.io/en/v7.19.1/index.html) -->

 
# Pangenome Mutation-Annotated Network (PanMAN)
<div align="center">
  <img src="docs/images/logo.svg"/>
</div>

## Table of Contents
- [Introduction](#intro) ([Wiki](https://turakhia.ucsd.edu/panman/))
  - [PanMANs](#panman)
  - [<i>panmanUtils</i>](#panmanUtils)
- [Installation](#install)
  - [Using Conda](#conda) (Recommended)
  - [Using Docker Image](#image)
  - [Using DockerFile](#file)
  - [Using Installation Script](#script)
- [PanMAN Construction](#construct)
  - [From PanGraph/GFA/MSA](#pangraph)
  - [From Raw Sequences](#raw)
  - [From Fragment Assemblies](#frag)
- [<i>panmanUtils</i> functionalities](#function)
- [Contribute](#contributions)
- [Citing PanMAN](#cite_panman)

## <a name="intro"></a> Introduction
Here we provide an overview of PanMAN, <i>panmanUtils</i>, and its installation methods and usage. For more information please see our [Wiki](https://turakhia.ucsd.edu/panman/).
### <a name="panman"></a> What is a PanMAN?
PanMAN or Pangenome Mutation-Annotated Network is a novel data representation for pangenomes that provides massive leaps in both representative power and storage efficiency. Specifically, PanMANs are composed of mutation-annotated trees, called PanMATs, which, in addition to substitutions, also annotate inferred indels (Fig. 1b), and even structural mutations (Fig. 1a) on the different branches. Multiple PanMATs are connected in the form of a network using edges to generate a PanMAN (Fig. 1c). PanMAN's representative power is compared against existing pangenomic formats in Fig. 1d. PanMANs are the most compressible pangenomic format for the different microbial datasets (SARS-CoV-2, RSV, HIV, Mycobacterium. Tuberculosis, E. Coli, and Klebsiella pneumoniae), providing 2.9 to 559-fold compression over standard pangenomic formats. 
<div align="center">
    <div><b>Figure 1: Overview of the PanMAN data structure</b></div>
    <img src="docs/images/panman.svg" width="500"/>
</div>

### <a name="panmanUtils"></a> <i>panmanUtils</i>
<i>panmanUtils</i> includes multiple algorithms to construct PanMANs and to support various functionalities to modify and extract useful information from PanMANs (Fig. 2).

<div align="center">
    <div><b>Figure 2: Overview of <i>panmanUtils</i>' functionalities</b></div>
    <img src="docs/images/utility.svg" width="500"/>
</div>


## <a name="install"></a> Installation
<b><i>panmanUtils</i></b> software can be installed using four different methods: 
1. Conda (Recommended) 
2. Docker Image
3. Dockerfile
4. Installation scripts 
### 1. Using conda (recommended) <a name="conda"></a>
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
### 2. Using Docker Image <a name="image"></a>

To use <i>panmanUtils</i> in a docker container, users can create a docker container from a docker image, by following these steps (compatible with `linux-64` and `osx-64`).

#### i. Dependencies
1. [Docker](https://docs.docker.com/engine/install/)
#### ii. Pull and build the PanMAN docker image from DockerHub
```bash
## Note: If the Docker image already exist locally, make sure to pull the latest version using 
## docker pull swalia14/panman:latest

## If the Docker image does not exist locally, the following command will pull and run the latest version
docker run -it swalia14/panman:latest
```
#### iii. Run <i>panmanUtils</i>
```bash
# Insider docker container
panmanUtils --help
```

### 3. Using DockerFile <a name="file"></a>
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
#### i. Dependencies
1. [Git](https://git-scm.com/downloads)

#### ii. Clone the repository
```bash
git clone https://github.com/TurakhiaLab/panman.git
cd panman
```
#### iii. Run the installation script <a name="script"></a>
```bash
chmod +x install/installationUbuntu.sh
./install/installationUbuntu.sh
```
#### iv. Run <i>panmanUtils</i>
```bash
cd build
./panmanUtils --help
```

## <a name="construct"></a> PanMAN Construction
Once the package is installed, PanMANs can be constructed from PanGraph [or GFA or MSA] and Tree topology (Newick format) using <i>panmanUtils</i>. Here we provide examples for constructing PanMANs from PanGraph (JSON) and custom dataset. Alternatively, users can follow the instructions provided in [wiki](https://turakhia.ucsd.edu/panman/) for other methods.
### <a name="pangraph">Building PanMAN from PanGraph

**Step 1:** Check if `sars_20.json` and `sars_20.nwk` files exist in `test` directory. 

**Step 2:** Run <i>panmanUtils</i> with the following command to build a panman from PanGraph:

```bash
panmanUtils -P $PANMAN_HOME/test/sars_20.json -N $PANMAN_HOME/test/sars_20.nwk -o sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

### <a name="raw"></a> Building PanMAN from raw sequences or fragment assemblies using Snakemake Workflow
We provide a Snakemake workflow to construct PanMANs from raw sequences (FASTA format) or from fragment assemblies.

!!!Note
    The Snakemake workflow uses various tools such as PanGraph tool, PGGB, MAFFT, and MashTree to build input PanGraph, GFA, MSA, and Tree topology files, respectively and it is particularly designed to be used in the docker container build from either the provided docker image or the DockerFile (instructions provided [here](#install)).

#### Building PanMAN from raw genome sequences
**Step 1:** Run the following command to construct a panman from raw sequences.
* Usage
```bash
cd $PANMAN_HOME/workflows
snakemake --use-conda --cores 8 --config RUNTYPE="pangraph/gfa/msa" FASTA="[user_input]" SEQ_COUNT="Number of sequences" ASSEM="NONE" REF="NONE" TARGET="NONE"
```
* Example
```bash
cd $PANMAN_HOME/workflows
snakemake --use-conda --cores 8 --config RUNTYPE="pangraph" FASTA="$PANMAN_HOME/test/sars_20.fa" SEQ_COUNT="20" ASSEM="NONE" REF="NONE" TARGET="NONE"
```


#### <a name="frag"></a> Building PanMAN from fragment assemblies
**Step 1:** Run the following command to construct a panman from fragment assemblies.

```bash
cd $PANMAN_HOME/workflows
snakemake --use-conda --cores 8 --config RUNTYPE="pangraph/gfa/msa" FASTA="None" SEQ_COUNT="Number of sequences" ASSEM="frag" REF="reference_file" TARGET="target.txt"
```
Here, target.txt includes a list of files that contain the fragmented assemblies.

## <a name="function"></a> <i>panmanUtils</i> functionalities
<i>panmanUtils</i> provide various functionalities such as summary, [Raw sequence, MSA, VCF, GFA] extract, sub-network pruning, and many more. Please refer to [wiki](https://turakhia.ucsd.edu/panman/) for detailed information. Here we provide usage syntax and examples for summary and VCF extract.

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


## <a name="contri"></a> Contribute <br>
We welcome contributions from the community to enhance the capabilities of PanMAN and <i>panmanUtils</i>. If you encounter any issues or have suggestions for improvement, please open an issue on [PanMAN GitHub page](https://github.com/TurakhiaLab/panman/issues). For general inquiries and support, reach out to our team.

## <a name="cite_panman"></a> Citing PanMAN <br>
If you use the PanMANs or <i>panmanUtils</i> in your research or publications, we kindly request that you cite the following paper:
* Walia S., Motwani H., Tseng YH., Smith K., Corbett-Detig R., Turakhia Y., <i>Compressive pangenomics using mutation-annotated networks</i>. <b>Nat Genet</b> (2026). [https://doi.org/10.1038/s41588-025-02478-7](https://doi.org/10.1038/s41588-025-02478-7)

