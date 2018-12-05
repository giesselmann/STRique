# STRique [![Build Status](https://travis-ci.org/giesselmann/STRique.svg?branch=master)](https://travis-ci.org/giesselmann/STRique)
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
In order to download, build and install STRique , execute the following commands:

    git clone --recursive https://github.com/giesselmann/STRique
    cd STRique
    python3 setup.py install
    
Installation time < 5min

## Usage	  
	usage: python3 STRique.py [-h] [--out OUT] [--algn ALGN] [--config CONFIG] [--t T] f5 model repeat

	STR Detection in raw nanopore data

	positional arguments:
	  f5               fast5 file directory
	  model            pore model
	  repeat           repeat region config file

	optional arguments:
	  -h, --help       show this help message and exit
	  --out OUT        output file name, if not given print to stdout
	  --algn ALGN      alignment in sam format, if not given read from stdin
	  --config CONFIG  Config file with HMM transition probabilities
	  --t T            Number of processes to use in parallel

## Configuration
	targeted repeats are configured in a tsv file with columns
	
	chr	begin	end	name	repeat	prefix	suffix
	
	e.g. c9orf72
	chr9	27573527	27573544	c9orf72	GGCCCC	...GCCCCGACCACGCCCC	TAGCGCGCGACTCCTG...
	
	the longer the prefix/ suffix sequences, the more reliable the signal alignment. 
	Repeat, prefix and suffix sequence are always on template strand
	
## Test
Test the pipeline with the following commands in the cloned repository:

	python scripts/STRique_test.py	
	cat data/c9orf72.sam | python3 scripts/STRique.py ./data ./models/template_median68pA6mer.model ./configs/repeat_config.tsv --config ./configs/STRique.json

You should see output similar to 

	ID      target strand count score_prefix score_suffix log_p offset ticks
	ce47b364-ed6e-4409-808a-1041c0b5aac2 c9orf72 - 735 6.3155927807600545 6.031860427335506 -119860.52066647023 1633 40758

Run time <1 min on a typical desktop computer, multiprocessing supported
