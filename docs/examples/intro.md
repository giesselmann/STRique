# Usage

STRique works on raw nanopore read data in either single or bulk fast5 format. Batches of single reads in tar archives work as well. Before starting the repeat detection the raw data folder must be indexed to enable extraction of single reads. The index file contains relative paths to the reads and must be saved in the indexed directory.

	python3 scripts/STRique.py index [OPTIONS] input

	positional arguments:
	  input                 	Input batch or directory of batches

	optional arguments:
	  --recursive           	Recursively scan input
	  --out_prefix OUT_PREFIX	Prefix for file paths in output
	  --tmp_prefix TMP_PREFIX	Prefix for temporary data

The repeat detection requires an indexed raw data archive and the alignment of the reads:

	python3 scripts/STRique.py count [OPTIONS] f5 model repeat

	positional arguments:
	  f5Index          Fast5 index
	  model            pore model
	  repeat           repeat region config file

	optional arguments:
	  --out OUT        output file name, if not given print to stdout
	  --algn ALGN      alignment in sam format, if not given read from stdin
	  --config CONFIG  Config file with HMM transition probabilities
	  --t T            Number of processes to use in parallel


## Configuration

Targeted repeats are configured in a tsv file with columns

```
chr	begin	end	name	repeat	prefix	suffix
```

```
e.g. c9orf72
chr9	27573527	27573544	c9orf72	GGCCCC	...GCCCCGACCACGCCCC	TAGCGCGCGACTCCTG...
```
The longer the prefix/ suffix sequences, the more reliable the signal alignment. Repeat, prefix and suffix sequence are always on template strand
