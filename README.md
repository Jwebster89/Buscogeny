# Buscogeny: A Busco Phylogeny generator
Create a phylogenetic tree from BUSCO identified single copy genes

## Requirements
- Biopython
- pandas
- numpy
- alive bar
- mafft
- iqtree

## Quick Usage
Buscogeny can be run as follows with the hypocreales_odb10 database as an example

`buscogeny.py -i ./Folder_of_BUSCO_outputs/ -d ./hypocreales_odb10`

## Running BUSCO
Below is an example loop for running BUSCO prior to running Buscogeny. Note, the output folder will be named Genome1.fna_busco512_odb10 in the below example, this is important (for now). Please use this format to run BUSCO.
```
for file in ./input/*.fna ; do base=$(basename $file) ; busco -i $file -o ${base}_busco512_odb10 -m genome -l ./hypocreales_odb10 ; done
```

## Options and usage
```

usage: buscogeny.py -i INPUT -d DB [-h]

Create Phylogenies from BUSCO output

Required Arguments:
  -i INPUT, --input INPUT
                        input folder of Busco output directories
  -d DB, --db DB        Location of odb database

Optional Arguments:
  -h, --help            show this help message and exit

```
