[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg
[license-link]: [https://github.com/TurakhiaLab/panman/LICENSE](https://github.com/TurakhiaLab/panman/blob/main/LICENSE)
[![License][license-badge]][license-link]
[![DOI](https://img.shields.io/badge/DOI-https://zenodo.org/records/12630607-beige)](https://doi.org/10.5281/zenodo.14633185)
[<img src="https://img.shields.io/badge/Install with-DockerHub-informational.svg?logo=Docker">](https://hub.docker.com/r/swalia14/panman)
[<img src="https://img.shields.io/badge/Submitted to-bioRxiv-critical.svg?logo=LOGO">](https://doi.org/10.1101/2024.07.02.601807)
[<img src="https://img.shields.io/badge/Build with-CMake-green.svg?logo=snakemake">](https://cmake.org)
[<img src="https://img.shields.io/badge/Made with-Snakemake-aquamarine.svg?logo=snakemake">](https://snakemake.readthedocs.io/en/v7.19.1/index.html)
[<img src="https://img.shields.io/badge/Watch it on-Youtube-FF0000.svg?logo=YouTube">](https://www.youtube.com/watch?v=eh9zQElrmLI)
[![Build Status](https://github.com/TurakhiaLab/panman/actions/workflows/ci.yml/badge.svg)](https://github.com/TurakhiaLab/panman/actions)

 
# Pangenome Mutation Annotated Network (PanMAN)
<div align="center">
  <img src="docs/images/logo.svg"/>
</div>

This fork is for implementing and testing imputation for missing/`N` sequences.

## Table of Contents
- [Introduction](#intro) ([Wiki](https://turakhia.ucsd.edu/panman/))
  - [PanMANs](#panman)
  - [<i>panmanUtils</i>](#panmanUtils)
- [Installation](#install)
  - [Using Installation Script](#script)
  - [Using Docker Image](#image)
  - [Using DockerFile](#file)
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
### <a name="script"></a> Using installation script (requires sudo access)

**Step 0:** Dependencies
```bash
Git
```

**Step 1:** Clone the repository
```bash
git https://github.com/TurakhiaLab/panman.git
cd panman
```
**Step 2:** Run the installation script
```bash
chmod +x install/installationUbuntu.sh
./install/installationUbuntu.sh
```
**Step 3:** Run panmanUtils
```bash
cd build
./panmanUtils --help
```
### <a name="image"></a> Using Docker Image

To use <i>panmanUtils</i> in a docker container, users can create a docker container from a docker image, by following these steps

**Step 0:** Dependencies
```bash
Docker
```
**Step 1:** Pull the PanMAN docker image from DockerHub
```bash
docker pull swalia14/panman:latest
```
**Step 2:** Build and run the docker container
```bash
docker run -it swalia14/panman:latest
```
**Step 3:** Run panmanUtils
```bash
# Insider docker container
cd /home/panman/build
./panmanUtils --help
```

###  <a name="file"></a> Using DockerFile
Docker container with preinstalled <i>panmanUtils</i> can also be built from DockerFile by following these steps

**Step 0:** Dependencies
```bash
Docker
Git
``` 
**Step 1:** Clone the repository
```bash
git https://github.com/TurakhiaLab/panman.git
cd panman
```
**Step 2:** Build a docker image
```bash
cd docker
docker build -t panman .
```
**Step 3:** Build and run docker container
```bash
docker run -it panman
```
**Step 4:** Run panmanUtils
```bash
# Insider docker container
cd /home/panman/build
./panmanUtils --help
```

## <a name="construct"></a> PanMAN Construction
Once the package is installed, PanMANs can be constructed from PanGraph [or GFA or MSA] and Tree topology (Newick format) using <i>panmanUtils</i>. Here we provide examples for constructing PanMANs from PanGraph (JSON) and custom dataset. Alternatively, users can follow the instructions provided in [wiki](https://turakhia.ucsd.edu/panman/) for other methods.
### Building PanMAN from PanGraph

**Step 1:** Check if `sars_20.json` and `sars_20.nwk` files exist in `test` directory. 

**Step 2:** Run <i>panmanUtils</i> with the following command to build a panman from PanGraph:

```bash
cd $PANMAN_HOME/build
./panmanUtils -P $PANMAN_HOME/test/sars_20.json -N $PANMAN_HOME/test/sars_20.nwk -O sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

### <a name="raw"></a> Building PanMAN from raw sequences or fragment assemblies using Snakemake Workflow
We provide a Snakemake workflow to construct PanMANs from raw sequences (FASTA format) or from fragment assemblies.

!!!Note
    The Snakemake workflow uses various tools such as PanGraph tool, PGGB, MAFFT, and MashTree to build input PanGraph, GFA, MSA, and Tree topology files, respectively and it is particularly designed to be used in the docker container build from either the provided docker image or the DockerFile (instructions provided [here](#install)).

#### Building PanMAN from raw genome sequences
**Step 1:** Run the following command to construct a panman from raw sequences.

```bash
cd $PANMAN_HOME/workflows
conda activate snakemake
snakemake --use-conda --cores 8 --config RUNTYPE="pangraph/gfa/msa" FASTA="[user_input]" SEQ_COUNT="Number of sequences" ASSEM="NONE" REF="NONE" TARGET="NONE"
```

#### <a name="frag"></a> Building PanMAN from fragment assemblies
**Step 1:** Run the following command to construct a panman from fragment assemblies.

```bash
cd $PANMAN_HOME/workflows
conda activate snakemake
snakemake --use-conda --cores 8 --config RUNTYPE="pangraph/gfa/msa" FASTA="None" SEQ_COUNT="Number of sequences" ASSEM="frag" REF="reference_file" TARGET="target.txt"
```
Here, target.txt includes a list of files that contain the fragmented assemblies.

## <a name="function"></a> <i>panmanUtils</i> functionalities
<i>panmanUtils</i> provide various functionalities such as summary, [Raw sequence, MSA, VCF, GFA] extract, sub-network pruning, and many more. Please refer to [wiki](https://turakhia.ucsd.edu/panman/) for detailed information. Here we provide usage syntax and examples for summary and VCF extract.

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


## <a name="contri"></a> Contribute <br>
We welcome contributions from the community to enhance the capabilities of PanMAN and <i>panmanUtils</i>. If you encounter any issues or have suggestions for improvement, please open an issue on [PanMAN GitHub page](https://github.com/TurakhiaLab/panman). For general inquiries and support, reach out to our team.

## <a name="cite_panman"></a> Citing PanMAN <br>
If you use the PanMANs or <i>panmanUtils</i> in your research or publications, we kindly request that you cite the following paper: 
* Sumit Walia, Harsh Motwani, Kyle Smith, Russell Corbett-Detig, Yatish Turakhia, "<i>Compressive Pangenomics Using Mutation-Annotated Networks</i>", bioRxiv 2024.07.02.601807; doi: [10.1101/2024.07.02.601807](https://doi.org/10.1101/2024.07.02.601807)

