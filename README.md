# GENOMES TO MSA 

## Description

This tool converts a list of genome files into a Multiple Sequence Alignment (MSA) file aligned with the Sibeliaz tool, ready for phylogenetic analysis. It also includes an option to generate a RapidNJ tree.

## Usage Example

```
bash
python3 Genome_to_MSA.py -g Example_dataset/Yersinia_genome.tab -o Yersinia_Genome_to_MSA_dir -nsites 10000 -K 15 -t 16 --Run_NJ_phylogeny yes
```


## Available options :
```
usage: Genome_to_MSA.py [-h] [-g GENOME_PATH] [-o OUTPUT_DIR] [-nsites NSITES] [-K K] [-t THREADS] [-Taxon TAXON] [-NJ RUN_NJ_PHYLOGENY]

Allow to generate a MSA from a list of genomes using Sibeliaz.

options:
  -h, --help            show this help message and exit
  -g GENOME_PATH, --genome_path GENOME_PATH
                        A txt file with one column with the full path of the genomes to align
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The output directory where will be saved files including the alignement file
  -nsites NSITES, --nsites NSITES
                        The number of desired sites in the final MSA (taken randomly)
  -K K, --K K           For small datasets, like bacteria, Sibeliaz recommend k=15, and for mammalian-sized genomes k=25. The default is 25.
  -t THREADS, --threads THREADS
                        Number of threads to use.
  -Taxon TAXON, --Taxon TAXON
                        comma separated list of taxa for wich we want to force to include all sites in the MSA, in that case number of final sites >=
                        nsites. Name of Taxa should be the name of the Genome file without the full path and remove everythin after the '.'
  -NJ RUN_NJ_PHYLOGENY, --Run_NJ_phylogeny RUN_NJ_PHYLOGENY
                        Set to yes if you want a quick phylogeny made with rapidnj.
```

Note : you will need Sibeliaz and optionnaly rapidnj if asked to be installed as well as few python packages (see header of the python script to install them).

______________

## Expected input table file format 
```
cat Yersinia_genome.tab

/path1/path2/Genome1.fna
/path1/path2/Genome2.fna
/path1/path2/Genome3.fna
/path1/path2/Genome4.fna
```
______________

## Outputfiles :

- ***blocks_coords.gff*** :  GFF file output from Sibeliaz.
- ***alignment.maf*** : Maf file output from Sibeliaz.
- ***output_sequences.fasta*** : MSA file containing all the homologous sites found among species.
- ***output_sequences_trimmed.fasta*** : MSA file containing a set o X sites selected (n=nsites).
- ***output_sequences_trimmed.nwk*** : Newik file if option --Run_NJ_phylogeny was set to yes (default is 100 bootstraps done).

The contig names in the final fasta files will be the name of your Genomes files without the extension, so in this example the expected output would be : 

```
>Genome1
AAAXXXXX
>Genome2
AAAXXXXX
>Genome3
AAAXXXXX
>Genome4
AAAXXXXX
```


A full run including NJ tree on 14 Yersinia genomes with 100,000 selected sites took around 30min of run. 





