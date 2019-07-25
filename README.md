## Setup
We recommend [miniconda3] to install all the required R/Python dependencies.  

All the dependencies can be installed in a new conda environment using `environment.yml` (Linux only):

    conda env create -n venv -f environment.yml

Alternatively, re-create the env by:

    conda create -n venv r "python>=3.5"            \
        bedtools bedops                             \
        snakemake r-data.table r-stringr r-plyr     \
        r-lme4

Modify all the paths to data files in `config.yaml`.

[miniconda3]: https://docs.conda.io/en/latest/miniconda.html


## Pipeline execution 
The pipeline is managed by [Snakemake].  Here is the basic usage of Snakemake:

     # Dry-run and print out the exact commands to be executed
    snakemake -n -p <target rule name and/or output(s)>  

    # Launch Snakemake using 16 CPU cores
    snakemake -j 16 -p <target rule name and/or output(s)>

[Snakemake]: https://snakemake.readthedocs.io/en/stable/
