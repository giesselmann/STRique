# STRique
[![Latest GitHub release](https://img.shields.io/github/release-pre/giesselmann/STRique.svg)](https://github.com/giesselmann/STRique/releases/latest) 
[![Build Status](https://travis-ci.org/giesselmann/STRique.svg?branch=master)](https://travis-ci.org/giesselmann/STRique) 
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/giesselmann/strique.svg)](https://hub.docker.com/r/giesselmann/strique/)

STRique is a nanopore raw signal repeat detection pipeline

## Dependencies
**Python 3.5 or higher**

	pomegranate
	numpy, scipy, scikit-image
	h5py

**C++**

	g++ > 5 (C++14 support)
	SeqAn2
	Pybind11

Dependencies get downloaded and build by the setup script.

## Installation
In order to download, build and install STRique, execute the following commands (Consider using a python virtual environment):

    git clone --recursive https://github.com/giesselmann/STRique
    cd STRique
	pip3 install -r requirements.txt
    python3 setup.py install

Alternatively we provide a Docker container containing the STRique pipeline:

	docker pull giesselmann/strique

Installation time < 5min.

## Usage

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

	chr	begin	end	name	repeat	prefix	suffix

	e.g. c9orf72
	chr9	27573527	27573544	c9orf72	GGCCCC	...GCCCCGACCACGCCCC	TAGCGCGCGACTCCTG...

The longer the prefix/ suffix sequences, the more reliable the signal alignment. Repeat, prefix and suffix sequence are always on template strand

## Test & Examples
Test the pipeline with the following commands in the cloned repository:

	cd STRique
	python3 scripts/STRique_test.py
	python3 scripts/STRique.py index data/ > data/reads.fofn
	cat data/c9orf72.sam | python3 scripts/STRique.py count ./data/reads.fofn ./models/template_median68pA6mer.model ./configs/repeat_config.tsv --config ./configs/STRique.json

You should see output similar to

	ID      target strand count score_prefix score_suffix log_p offset ticks
	ce47b364-ed6e-4409-808a-1041c0b5aac2 c9orf72 - 735 6.3155927807600545 6.031860427335506 -119860.52066647023 1633 40758

Run time <1 min on a typical desktop computer, multiprocessing supported
