
██████╗ ███████╗███╗   ██╗ ██████╗ ███╗   ███╗███████╗███████╗    ████████╗ ██████╗     ███╗   ███╗███████╗ █████╗ 
██╔════╝ ██╔════╝████╗  ██║██╔═══██╗████╗ ████║██╔════╝██╔════╝    ╚══██╔══╝██╔═══██╗    ████╗ ████║██╔════╝██╔══██╗
██║  ███╗█████╗  ██╔██╗ ██║██║   ██║██╔████╔██║█████╗  ███████╗       ██║   ██║   ██║    ██╔████╔██║███████╗███████║
██║   ██║██╔══╝  ██║╚██╗██║██║   ██║██║╚██╔╝██║██╔══╝  ╚════██║       ██║   ██║   ██║    ██║╚██╔╝██║╚════██║██╔══██║
╚██████╔╝███████╗██║ ╚████║╚██████╔╝██║ ╚═╝ ██║███████╗███████║       ██║   ╚██████╔╝    ██║ ╚═╝ ██║███████║██║  ██║
 ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝╚══════╝       ╚═╝    ╚═════╝     ╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝


Input a table with a list of genome files > output a MSA file aligned with the Sibeliaz ready for phylogenetic analysis (+ RapidNJ tree included)

Usage example : 

'''python3 Genome_to_MSA.py -g Yersinia_genomes/Yersinia_genome.tab -o Yersinia_Genome_to_MSA_dir -nsites 10000 -K 15 -t 16 --Run_NJ_phylogeny yes'''

Available options : 

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
  -K K, --K K           For small datasets, like bacteria, we recommend k=15, and for mammalian-sized genomes k=25. The default is 25.
  -t THREADS, --threads THREADS
                        Number of threads to use.
  -Taxon TAXON, --Taxon TAXON
                        comma separated list of taxa for wich we want to force to include all sites in the MSA, in that case number of final sites >=
                        nsites. Name of Taxa should be the name of the Genome file without the full path and remove everythin after the '.'
  -NJ RUN_NJ_PHYLOGENY, --Run_NJ_phylogeny RUN_NJ_PHYLOGENY
                        Set to yes if you want a quick phylogeny made with rapidnj.
