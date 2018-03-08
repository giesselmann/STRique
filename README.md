# repetitION
repetitION is a nanopore raw signal repeat detection pipeline
## Dependencies
**Python**

	Python 3.5 or higher
	pomegranate
	biopython
	numpy, scipy
	h5py
		
**C++**
	
	g++ > 5 (C++14 support)
	SeqAn2 (downloaded by install script)
	Pybind11 (downloaded by install script)

## Installation
In order to download, build and install repetitION , execute the following commands:

    git clone --recursive https://github.com/giesselmann/repetitION
    cd repetitION
	pip install .

## Usage
	python repetition.py [--sam SAM] [--ID ID] [--t T] f5 model config out

	positional arguments:
	  f5          fast5 file directory
	  model       pore model
	  config      repeat config file
	  out         output directory

	optional arguments:
	  --sam SAM   alignment file in sam format
	  --ID ID     Read ID filter, one ID per line
	  --t T       Number of processes to use in parallel
	  
## Test
Test the pipeline with the following commands in the cloned repository:

	python scripts/repetition.py data models/template_median68pA6mer.model configs/c9orf72.json ./

A File *0.tsv* with similar content should have been created:

	ID	ref	flag	count	score_prefix	score_suffix	log_p	ticks	offset
	ce47b364-ed6e-4409-808a-1041c0b5aac2	chr9	16	147	4713.84716796875	4765.23779296875	-28627.333236486844	40978	1635
