configfile: 'config.yaml'

GENCODE_PROTEIN_CODING_BED = 'annotations/gencode.proteincoding.known.gene.bed'
EXOME_SAMPLE_IDS = []
GENOME_SAMPLE_IDS = []

with open(config['OVERLAPPED_SAMPLE_INFO']) as f:
    header = next(f)
    for line in f:
        exome_id, genome_id, *_ = line.rstrip().split(' ')
        EXOME_SAMPLE_IDS.append(exome_id)
        GENOME_SAMPLE_IDS.append(genome_id)

rule all:
    input:
        expand(
            'processed_data/GAFexon/{sample_id}.exome.broadgaf.maf',
            sample_id=EXOME_SAMPLE_IDS,
        ),
        expand(
            'processed_data/GAFexon/{sample_id}.genome.broadgaf.maf',
            sample_id=GENOME_SAMPLE_IDS,
        )


rule gencode_bed_protein_coding_only:
    # Create protein only region as BED file from GENCODE v19 annotation
    input: config['GENCODE_GTF']
    output: GENCODE_PROTEIN_CODING_BED
    shell:
        r"""
        zcat {input} |                                      \
        awk -F'[;\t]'                                       \
            '{{OFS="\t"}}{{                                 \
                if(                                         \
                    $3 == "exon" &&                         \
                    $12 == " gene_status \"KNOWN\"" &&      \
                    $11 == " gene_type \"protein_coding\""  \
                )                                           \
                    print $1,$4,$5,$13                      \
            }}' |                                           \
        sed -E -e 's/^chr([0-9]+|X|Y|M)/\1/' > {output}
        """


def map_maf_filepath(wildcards):
    maf_name = wildcards['name']
    # replace all `.` in the file name with `_` 
    config_entry = maf_name.replace('.', '_')
    maf_filepath = config['SAMPLES'][config_entry]
    return maf_filepath


rule temp_link_partially_processed_maf:
    """Link the intermediate MAFs by looking up their filepaths in config
    
    This rule use a customized function to look up in config entries under
    SAMPLES to find the corresponding file path, which usually relies in the
    original analysis folders.
    
    This rule should be removed once all processing steps have been ported to
    Snakemake.
    """
    input: map_maf_filepath
    output: 'processed_data/{name}.maf'
    shell:
        'ln -s {input} {output}'


rule maf_protein_coding_only:
    # Filter MAF file with only GENCODE protein coding regions
    input:
        maf='processed_data/{name}.maf',
        bed=GENCODE_PROTEIN_CODING_BED,
    output: 'processed_data/{name}.gaf4bed.exon.maf'
    shell:
        "bedtools intersect -a {input.maf} -b {input.bed} | "
        "sort -u -s -k 1,1V -k 2,2n -k 3,3n > {output}"


rule split_exome_maf:
    input: 'processed_data/exome.broadbed.gaf4bed.exon.maf'
    output: 
        expand(
            'processed_data/GAFexon/{sample_id}.exome.broadgaf.maf',
            sample_id=EXOME_SAMPLE_IDS
        )
    shell:
        r"""
        awk -F"\t" '{{
            print > "processed_data/GAFexon/"$12".exome.broadgaf.maf"
        }}' {input}
        """


rule split_genome_maf:
    input: 'processed_data/genome.broadbed.gaf4bed.exon.maf'
    output: 
        expand(
            'processed_data/GAFexon/{sample_id}.genome.broadgaf.maf',
            sample_id=GENOME_SAMPLE_IDS
        )
    shell:
        r"""
        awk -F"\t" '{{
            print > "processed_data/GAFexon/"$46".genome.broadgaf.maf"
        }}' {input}
        """
