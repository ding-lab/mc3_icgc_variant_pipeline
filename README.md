## Setup

Require either [miniconda3] or [Anaconda3] and Python 3.5+.

Create a new conda virtual environment with name `venv`:

    conda create -n venv -f environment.yml


Check if all the paths in `config.yaml` exist. If not, modify them with the
correct paths.


## Usage

    snakemake --cores 4 -p  # launch Snakemake using 4 CPU cores
    fab clean               # remove all generated output files

    snakemake -n -p         # dry-run and print out the exact commands to be executed
    snakemake --summary     # print out the table of all related files
                            # and their status of whether being re-evaluated or not
