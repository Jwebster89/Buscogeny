# Buscogeny: A Busco Phylogeny generator
Buscogeny (Busco Phylogeny) is a tool designed to create phylogenetic trees from BUSCO (Benchmarking Universal Single-Copy Orthologs) identified single-copy genes. It simplifies the process of phylogenetic analysis by automating the identification and alignment of orthologs and generating a phylogenetic tree using these aligned sequences.

## Requirements
Buscogeny requires the following software and libraries:
- Biopython
- pandas
- matplotlib
- plotly
- numpy
- alive bar
- mafft
- iqtree
- clipkit
- ClonalFrameML
- maskrc-svg
- BUSCO (5.4.0+ for -s correct behaviour) - Note: For odb10 databases please use BUSCO > 5.4.0 and < 5.8.0. For odb12, please use BUSCO 5.8.0+

## Installation
To install Buscogeny and its dependencies, follow these steps:
```
git clone https://github.com/Jwebster89/Buscogeny.git
cd Buscogeny
conda create -n buscogeny biopython pandas numpy mafft iqtree busco clipkit matplotlib plotly maskrc-svg clonalframeml
conda activate buscogeny
pip install alive-progress
```


## Quick Usage
Buscogeny can be run as follows with the hypocreales_odb10 database and 60 threads as an example

`buscogeny.py -i ./Folder_of_genomes -d ./hypocreales_odb10 -t 60`

## Output
```
Buscogeny_out/
│
├── alignments/
│   ├── alignment1.aln
│   ├── alignment2.aln
│   └── ... (other alignment files)
│
├── BUSCO/
│   ├── genome1_busco_out/
│   ├── genome2_busco_out/
│   └── ... (other BUSCO output directories for each genome)
│
├── ClonalFrameML/ (--rc_filt true)
│   ├── BUSCO_alignment.cfml.em.em.txt
│   ├── ... (other ClonalFrameML output)
│   ├── BUSC_alignment.rc_masked.aln
│   └── ... (other maskrc_svg output)
│
├── Genome_ortholog_counts.png
├── Genome_ortholog_counts.html
│
├── iqtree/
│   ├── iqtree.treefile
│   ├── iqtree.log
│   ├── iqtree.iqtree
│   └── ... (other IQ-TREE output files)
│
├── supermatrix/
│   ├── superaln.clstl
│   ├── superaln.degapped.clstl
│   ├── concatenated_alignment.fasta (--rc_filt true)
│   └── core_aln.xmfa (--rc_filt true)
│
└── targets/
    ├── target1_all.fna
    ├── target2_all.fna
    └── ... (other target files)

```

## Options and usage
```
usage: buscogeny.py -i INPUT -d DB -o OUTPUT [-h] [-t THREADS] [-g GAPPY_THRESHOLD] [-e EXCLUDE_THRESHOLD] [-s {prot,nucl}] [-r]

        Create Phylogenies from BUSCO output. 

        Version: 2.0.0

Required Arguments:
  -i INPUT, --input INPUT
                        Input folder of Genomes
  -d DB, --db DB        Location of odb database
  -o OUTPUT, --output OUTPUT
                        Location of output [prefix]_Buscogeny_out

Optional Arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use. Default 8
  -g GAPPY_THRESHOLD, --gappy_threshold GAPPY_THRESHOLD
                        Specifies gaps threshold used by clipkit. Default 0.05
  -e EXCLUDE_THRESHOLD, --exclude_threshold EXCLUDE_THRESHOLD
                        Specifies proportion of alignments an isolate is allowed to be missing from.
  -s {prot,nucl}, --seq_type {prot,nucl}
                        Alignment using: 'prot' for protein, 'nucl' for nucleotide. Default is 'prot'.
  -r, --rc_filt         Enable recombination filtering
```
## Notes
- Please ensure that when specifying the odb database to use a relative or absolute path. E.g. `./bacteria_odb10`, and not `bacteria_odb10`
- Buscogeny will skip folders if it detects previous output. This can be useful for rerunning if the run is interrupted or you would like to reuse previous BUSCO calculations. However, care should be taken to remove the interupted processes outputs.
- Recombination filtering is only available on nucleotide sequences, run Buscogeny with `-s nucl` if you plan to also use `-r`.
- A current issue between BUSCO and python results in a "SyntaxWarning: invalid escape sequence '\w'", this is planned to be addressed in the next version of BUSCO.
- BUSCO odb databases can be downloaded from https://busco-data.ezlab.org/v5/data/lineages/