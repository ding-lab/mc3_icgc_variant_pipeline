## Setup

Require either [miniconda3] or [Anaconda3] and Python 3.5+.

Create a new conda virtual environment with name `venv`:

    conda env create -n venv -f environment.yml

Note that current `environment.yml` works only on Linux. Alternatively, re-create the env by:

    conda create -n venv r python                   \
        bedtools                                    \
        snakemake r-data.table r-stringr r-plyr


Modify all the paths to data files in `config.yaml`.


## Usage

    snakemake --cores 4 -p  # launch Snakemake using 4 CPU cores
    fab clean               # remove all generated output files

    snakemake -n -p         # dry-run and print out the exact commands to be executed
    snakemake --summary     # print out the table of all related files
                            # and their status of whether being re-evaluated or not
