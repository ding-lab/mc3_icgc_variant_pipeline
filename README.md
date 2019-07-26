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

Since the full pipeline will take substantial time to complete, it is advised to run the snakemake by "stages".

    # Get all sequencing depth coverage files within region of interest
    snakemake all_exome_beds all_genome_beds

    # Step 1 reduce to shared samples
    snakemake reduce_exome_maf reduce_genome_maf get_orig_maf_header_exome get_orig_maf_header_genome

    # Step 2 reduce to calls with WES target region
    snakemake broadbed_exome broadbed_genome

    # Step 3 reduce to GENCODE v19 exon regions
    snakemake gafbed_exome gafbed_genome

    # Step 4: reduce the variants by exome wig files
    snakemake merge_exome_wig_reduced_exome_mafs merge_exome_wig_reduced_genome_mafs

    # STEP 5: Reduce the variants by genome wig files
    snakemake merge_self_exome merge_self_genome

    # Make SQLite database of variants
    snakemake make_sqlite_db

    # Generate variant matching outputs for downstream analysis
    snakemake make_analysis_file clean_analysis_file

[Snakemake]: https://snakemake.readthedocs.io/en/stable/


## Downstream analysis and figures generation
All the figures in the publication can be generated via Snakemake:

    snakemake all_figures
