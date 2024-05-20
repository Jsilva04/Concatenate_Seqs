# Sequence Concatenator

## Description
This script is written in Python and concatenates sequence data from files in Fasta, Phylip or Nexus formats. 
The program can handle missing sequences by either removing individuals not represented in all fragments or adding sequences filled with 'N'.

## Features
- Supports input files in Fasta, Phylip and Nexus.
- Concatenates sequences for each individual across multiple files.
- Handles missing sequences.

## Requirements
* Python 3.x
* Biopython

## Installation
Install Biopython using pip:
```bash
pip install biopython 
```
Install Biopython using conda:
```bash
conda install conda-forge::biopython
```

## Usage
```bash
python concatenate_sequences.py [-h] -f {format} -o [OUTPUT] [--fill-missing] [--sequence-length SEQUENCE_LENGTH] input_files
```

## License
MIT License
