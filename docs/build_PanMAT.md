# Building a PanMAT
To build a PanMAT, the PanMAT utility (panmat-utils) requires two inputs: 
1. The tree topology of the input sequences in the form a newick string 
2. A data structure consisting of an alignment of sequences to be represented by the PanMAT. This data-structure can be represented by a PanGraph file, a GFA file or an MSA file representing the PanGenome. 

Use the following command to build the PanMAT

```
$ ./panmat-utils --pangraph-in=<path to PanGraph JSON file> --newick-in=<path to newick file>
```
If you'd like to give a GFA or an MSA file as an input, you could use the `--gfa-in` or the `--msa-in` options.

---
> **_NOTE:_** Currently, we only support a simple version of the GFAv1.1 consisting of Segments, un-overlapping Links and Paths.
---

This will generate the PanMAT and give you access to the command line interface of the PanMAT utility where you can enter instructions to analyze and manipulate PanMAT files, (see [PanMAT utilities](utilities.md)).
```
Data load time: 40045364 nanoseconds
>
```
<!-- ## Writing the PanMAT to a file
The PanMAT that you build using the source files can be written to a new file for fast and easy access from the utility's command line interface as follows:
```
> write <file_name>
```
This creates a file called `<filename>.pmat` in the `pman`subdirectory. This subdirectory is intended to contain all the PanMAT files.

## Reading a PanMAT
You can load the PanMAT file using the PanMAT utility directly using the `--panmat-in` option.
```
~$ ./panmat-utils --panmat-in=<path to PanMAT file>
``` -->