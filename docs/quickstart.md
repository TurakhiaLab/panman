# Quick start

Here, we will learn to build PanMAN from various input formats.

**Step 0:** The Steps below require panmanUtils, if not done so far, refer to [installation guide](install.md) to install panmanUtils. To check if panmanUtils is properly installed or not, run the following command, and it should execute without error
```bash
# enter into the panman directory (assuming $PANMAN directs to the panman repository directory)
cd $PANMAN_HOME
```
```bash
cd $PANMAN_HOME/build
./panmanUtils --help
```
### Building PanMAN from PanGraph

**Step 1:** Check if `sars_20.json` and `sars_20.nwk` files exist in `test` directory. Otherwise, follow the instructions to download the dataset.

```bash
cd $PANMAN_HOME/dataset
TODO
```

**Step 2:** Run panmanUtils with the following command to build a panman from PanGraph:

```bash
cd $PANMAN_HOME/build
./panmanUtils -P $PANMAN_HOME/test/sars_20.json -N $PANMAN_HOME/test/sars_20.nwk -O sars_20
```
The above command will run <i>panmanUtils</i> program and build `sars_20.panman` in `$PANMAN_HOME/build/panman` directory.

### Building PanMAN from raw genome sequences
We provide scripts to first construct PanGraph from raw sequences, followed by building a panman.
**Step 1:** Check if the `sars_20.fa` file exists in `test` directory. Otherwise, follow the instructions to download the dataset.

```bash
cd $PANMAN_HOME/dataset
TODO
```

**Step 2:** Run the following command to construct a panman from raw sequences.

```bash
cd $PANMAN_HOME/scripts
chmod +x build_panman.sh
./build_panman.sh
```
!!!Note
    The above script is particularly designed to be used in the docker container build either from the provided docker image or the DockerFile (instructions provided [here](install.md))

