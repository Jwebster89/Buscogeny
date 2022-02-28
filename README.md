# Buscogeny: A Busco Phylogeny generator
Create a phylogenetic tree from BUSCO identified single copy genes

## Requirements
- Biopython
- pandas
- numpy
- alive bar
- mafft
- iqtree
- clipkit

## Installation
```
git clone https://github.com/Jwebster89/Buscogeny.git
cd Buscogeny
conda create -n buscogeny biopython pandas numpy mafft iqtree busco=5.1.2 clipkit
conda activate buscogeny
pip install alive-progress
```


## Quick Usage
Buscogeny can be run as follows with the hypocreales_odb10 database as an example

`buscogeny.py -i ./Folder_of_BUSCO_outputs/ -d ./hypocreales_odb10`

	
## Options and usage
```

        Create Phylogenies from BUSCO output.

        Version: 0.1.1

Required Arguments:
  -i INPUT, --input INPUT
                        Input folder of Genomes
  -d DB, --db DB        Location of odb database

Optional Arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use. Default 8


```
