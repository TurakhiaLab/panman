# Pangenome Mutation Annotated Newtork (PanMAN) and <i>panmanUtils</i>

## What is a PanMAN?
PanMAN or Pangenome Mutation-Annotated Network is a novel data representation for pangenomes that provides massive leaps in both representative power and storage efficiency. Specifically, PanMANs are composed of mutation-annotated trees, called PanMATs, which, in addition to substitutions, also annotate inferred indels (Fig. 1b), and even structural mutations (Fig. 1a) on the different branches. Multiple PanMATs are connected in the form of a network using edges to generate a PanMAN (Fig. 1c). PanMAN's representative power is compared against existing pangenomic formats in Fig. 1d. PanMANs are the most compressible pangenomic format for the different microbial datasets (SARS-CoV-2, RSV, HIV, Mycobacterium. Tuberculosis, E. Coli, and Klebsiella pneumoniae), providing 2.9 to 559-fold compression over standard pangenomic formats. 
<div align="center">
    <div><b>Figure 1: Overview of the PanMAN data structure</b></div>
    <img src="docs/images/panman.svg" width="300"/>
</div>

## <i><b>panmanUtils</b></i>
<i>panmanUtils</i> includes multiple algorithms to construct PanMANs and to support various functionalities to modify and extract useful information from PanMANs (Fig. 2).
<div align="center">
    <div><b>Figure 2: Overview of panmanUtils' functionalities</b></div>
    <img src="docs/images/utility.svg" width="300"/>
    ![](docs/images/utility.svg)
</div>

## Installation
panmanUtils can be installed using two different options (i) installation script and (ii) Docker, as described below.

### Installation Scripts
```
git clone https://github.com/TurakhiaLab/panman.git
cd panman/install
./installUbuntu.sh
```

### Docker
```
ToDo
```

