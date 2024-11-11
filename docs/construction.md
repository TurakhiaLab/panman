# PanMAN Construction

Here, we will learn to build PanMAN from various input formats.

**Step 0:** The Steps below require <i>panmanUtils</i>, if not done so far, refer to [installation guide](install.md) to install <i>panmanUtils</i>. To check if <i>panmanUtils</i> is properly installed or not, run the following command, and it should execute without error
```bash
# enter into the panman directory (assuming $PANMAN directs to the panman repository directory)
cd $PANMAN_HOME
```
```bash
cd $PANMAN_HOME/build
./panmanUtils --help
```
### Building PanMAN from PanGraph

**Step 1:** Check if `sars_20.json` and `sars_20.nwk` files exist in `test` directory. Alternatively, users can provide custom PanGraph (JSON) and tree topology (Newick format) files to build a panman. 

<!-- ```bash
cd $PANMAN_HOME/dataset
TODO
``` -->

**Step 2:** Run <i>panmanUtils</i> with the following command to build a panman from PanGraph:

```bash
cd $PANMAN_HOME/build
./panmanUtils -P $PANMAN_HOME/test/sars_20.json -N $PANMAN_HOME/test/sars_20.nwk -O sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

### Building PanMAN from GFA

**Step 1:** Check if `sars_20.gfa` and `sars_20.nwk` files exist in `test` directory. Alternatively, users can provide custom GFA and tree topology (Newick format) files to build a panman. 

<!-- ```bash
cd $PANMAN_HOME/dataset
TODO
``` -->

**Step 2:** Run <i>panmanUtils</i> with the following command to build a panman from GFA:

```bash
cd $PANMAN_HOME/build
./panmanUtils -G $PANMAN_HOME/test/sars_20.gfa -N $PANMAN_HOME/test/sars_20.nwk -O sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

### Building PanMAN from MSA (FASTA format)

**Step 1:** Check if `sars_20.msa` and `sars_20.nwk` files exist in `test` directory. Alternatively, users can provide custom MSA (FASTA format) and tree topology (Newick format) files to build a panman. 
<!-- 
```bash
cd $PANMAN_HOME/dataset
TODO
``` -->

**Step 2:** Run <i>panmanUtils</i> with the following command to build a panman from GFA:

```bash
cd $PANMAN_HOME/build
./panmanUtils -M $PANMAN_HOME/test/sars_20.msa -N $PANMAN_HOME/test/sars_20.nwk -O sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

### Building PanMAN from raw genome sequences
We provide scripts to construct <i>panmanUtils</i> inputs (PanGraph/GFA/MSA and Newick) from raw sequences (FASTA format), followed by building a panman.

!!!Note
    This script uses various tools such as PanGraph tool, PGGB, MAFFT, and MashTree to build input PanGraph, GFA, MSA, and Tree topology files, respectively. The script is particularly designed to be used in the docker container build from either the provided docker image or the DockerFile (instructions provided [here](install.md)). ß

**Step 1:** Check if the `sars_20.fa` file exists in `test` directory. Alternatively, users can provide custom raw sequences (FASTA format) to build a panman. 

<!-- ```bash
cd $PANMAN_HOME/dataset
TODO
``` -->

**Step 2:** Run the following command to construct a panman from raw sequences.

```bash
cd $PANMAN_HOME/scripts
chmod +x build_panman.sh
./build_panman.sh [pangraph or gfa or msa]
```
    