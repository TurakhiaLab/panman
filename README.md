[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/TurakhiaLab/panman/LICENSE
[![License][license-badge]][license-link]
<!-- [![DOI](https://zenodo.org/badge/1.svg)](https://zenodo.org/badge/latestdoi/1) -->

# Pangenome Mutation Annotated Network (PanMAN)
<div align="center">
  <img src="docs/images/logo.svg"/>
</div>

## Table of Contents
- [Overview of PanMANs and <i>panmanUtils</i>](#overview)
- [Installation and Usage](#install) ([Documentation](https://turakhia.ucsd.edu/panman/))

- [Contribute](#contributions)
- [Citing PanMAN](#cite_panman)

## <a name="overview"></a> Overview of PanMAN and <i>panmanUtils</i> <br>
### What is a PanMAN?
PanMAN or Pangenome Mutation-Annotated Network is a novel data representation for pangenomes that provides massive leaps in both representative power and storage efficiency. Specifically, PanMANs are composed of mutation-annotated trees, called PanMATs, which, in addition to substitutions, also annotate inferred indels (Fig. 1b), and even structural mutations (Fig. 1a) on the different branches. Multiple PanMATs are connected in the form of a network using edges to generate a PanMAN (Fig. 1c). PanMAN's representative power is compared against existing pangenomic formats in Fig. 1d. PanMANs are the most compressible pangenomic format for the different microbial datasets (SARS-CoV-2, RSV, HIV, Mycobacterium. Tuberculosis, E. Coli, and Klebsiella pneumoniae), providing 2.9 to 559-fold compression over standard pangenomic formats. 
<div align="center">
    <div><b>Figure 1: Overview of the PanMAN data structure</b></div>
    <img src="docs/images/panman.svg" width="500"/>
</div>

### <i><b>panmanUtils</b></i>
<i>panmanUtils</i> includes multiple algorithms to construct PanMANs and to support various functionalities to modify and extract useful information from PanMANs (Fig. 2).

<!-- #### PanMAN constrution

<div align="center">
 <div><b>Figure 2: PanMAN construction pipeline using panmanUtils</b></div>
 <img src="docs/images/construct.svg" width="500"/>
</div> -->

<div align="center">
    <div><b>Figure 2: Overview of <i>panmanUtils</i>' functionalities</b></div>
    <img src="docs/images/utility.svg" width="500"/>
</div>


## <a name="install"></a> Installation and Usage <br>
For information on pnamanUtils installation and usage, please see our documentation page available [here](https://turakhia.ucsd.edu/panman/)

## <a name="contri"></a> Contribute <br>
We welcome contributions from the community to enhance the capabilities of PanMAN and <i>panmanUtils</i>. If you encounter any issues or have suggestions for improvement, please open an issue on [PanMAN GitHub page](https://github.com/TurakhiaLab/panman). For general inquiries and support, reach out to our team.

## <a name="cite_panman"></a> Citing PanMAN <br>
If you use the PanMANs or <i>panmanUtils</i> in your research or publications, we kindly request that you cite the following paper: XXX