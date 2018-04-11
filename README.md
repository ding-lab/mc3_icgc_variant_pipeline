## Setup

Require either [miniconda3] or [Anaconda3] and Python 3.5+.

Create a new conda virtual environment with name `venv`:

    conda env create -n venv -f environment.yml

Note that current `environment.yml` works only on Linux. Alternatively, re-create the env by:

    conda create -n venv r python                   \
        bedtools bedops                             \
        snakemake r-data.table r-stringr r-plyr     \
        r-lme4


Modify all the paths to data files in `config.yaml`.


## Usage

    snakemake --cores 4 -p  # launch Snakemake using 4 CPU cores
    fab clean               # remove all generated output files

    snakemake -n -p         # dry-run and print out the exact commands to be executed
    snakemake --summary     # print out the table of all related files
                            # and their status of whether being re-evaluated or not

    snakemake -n -p -k 	    # continue without dying if a files is not found 

    snakemake -n -p <rule>  # without "rule" name run up to a certain rule


## MATT's Attempts
    ### Problems with memory

    snakemake --cores 4 -p all 2> snakemake_bailey_115.log # This is my run a time 1:15 and it looks like Denali can only handle about 4 cores at any given time. This is going to take a bit longer than expected. 

    ### Working 
    snakemake --cores 4 -p merge_exome_wig_reduced_exome_mafs 2> bailey.snakemake.150.log 
    And the sort /tmp directory is too small. Rerun without sort. 
    

    
 
