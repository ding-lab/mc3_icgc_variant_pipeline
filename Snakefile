configfile: 'config.yaml'

from pathlib import Path

GENCODE_PROTEIN_CODING_BED = 'annotations/gencode.proteincoding.known.gene.bed'
EXOME_SAMPLE_IDS_TO_GENOME = {}
GENOME_SAMPLE_IDS_TO_EXOME = {}


with open(config['OVERLAPPED_SAMPLE_INFO']) as f:
    header = next(f)
    for line in f:
        exome_id, genome_id, *_ = line.rstrip().split(' ')
        EXOME_SAMPLE_IDS_TO_GENOME[exome_id] = genome_id
        GENOME_SAMPLE_IDS_TO_EXOME[genome_id] = exome_id


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
            }}' > {output}
        """
        #     }}' |                                           \
        # sed -E -e 's/^chr([0-9]+|X|Y|M)/\1/' > {output}
        


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

    For example, if input maf is ``exome.broadbed.maf``, it looks up the
    ``config['SAMPLES']['exome_broadbed']`` in the config file, and soft link
    the file under ``processed_data/``. Note that all dots in the filename are
    converted to underscore to reduce ambiguity.

    This rule should be removed once all processing steps have been ported to
    Snakemake.
    """
    input: map_maf_filepath
    output: 'processed_data/{name}.maf'
    wildcard_constraints:
        name=r"[^/]+"
    shell:
        'ln -s {input} {output}'


rule maf_protein_coding_only:
    # Filter MAF file with only GENCODE protein coding regions
    input:
        maf='processed_data/{seq_type}.maf',
        bed=GENCODE_PROTEIN_CODING_BED,
    output: 'processed_data/{seq_type}.gaf4bed.exon.maf'
    shell:
        "bedtools intersect -a {input.maf} -b {input.bed} | "
        "sort -u > {output}"
        # "sort -u -s -k 1,1V -k 2,2n -k 3,3n > {output}"


rule split_exome_maf:
    input: 'processed_data/exome.broadbed.gaf4bed.exon.maf'
    output:
        expand(
            'processed_data/GAFexon/{sample_exome_id}.exome.broadgaf.maf',
            sample_exome_id=EXOME_SAMPLE_IDS_TO_GENOME.keys()
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
            'processed_data/GAFexon/{sample_genome_id}.genome.broadgaf.maf',
            sample_genome_id=GENOME_SAMPLE_IDS_TO_EXOME.keys()
        )
    shell:
        r"""
        awk -F"\t" '{{
            print > "processed_data/GAFexon/"$46".genome.broadgaf.maf"
        }}' {input}
        """


rule reduce_exome_maf_low_coverage_wig:
    """
    Take split exome MAF and remove any variants not in whole exome wig
    (no sequencing coverage)

    Adapted from reduceGAFexonsMaf.py 
    """
    input:
        exome_maf='processed_data/GAFexon/{exome_id}.exome.broadgaf.maf',
        exome_wig_bed=str(Path(config['WIG_BED_DIR'], '{exome_id}.bed')),
    output:
        'processed_data/GAFexonReduced/{exome_id}.broadgaf.wig.exome.reduced.maf'
    shell:
        r"""
        bedtools intersect \
            -b {input.exome_wig_bed} \
            -a {input.exome_maf} \
            -wb | \
        awk -F"\t" '{{OFS = "\t"}} {{if($NF == 1) print $0}}' \
        > {output}
        """


def retrieve_exome_wig_bed(wildcards):
    # Given the genome sample id, return the paired exome sample id
    exome_id = GENOME_SAMPLE_IDS_TO_EXOME[wildcards.genome_id]
    return str(Path(
        config['WIG_BED_DIR'],
        '{exome_id}.bed'.format(exome_id=exome_id)
    ))


rule reduce_genome_maf_low_coverage_wig:
    # Take split genome MAF and remove any variants not in whole exome wig
    # (no sequencing coverage)
    input:
        genome_maf='processed_data/GAFexon/{genome_id}.genome.broadgaf.maf',
        exome_wig_bed=retrieve_exome_wig_bed,
    output:
        'processed_data/GAFexonReduced/{genome_id}.broadgaf.wig.genome.reduced.maf'
    shell:
        r"""
        bedtools intersect \
            -b {input.exome_wig_bed} \
            -a {input.genome_maf} \
            -wb | \
        awk -F"\t" '{{OFS = "\t"}} {{if($NF == 1) print $0}}' \
        > {output}
        """


rule merge_exome_wig_reduced_exome_mafs:
    # Concatenate all the reduced exome MAFs
    # readme.txt:38
    input:
        expand(
            'processed_data/GAFexonReduced/{exome_id}.broadgaf.wig.exome.reduced.maf',
            exome_id=EXOME_SAMPLE_IDS_TO_GENOME.keys(),
        )
    output:
        'processed_data/exome.broadbed.gafe.wigs.maf'
    shell:
        "cat {input} | sort -u > {output}"


rule merge_exome_wig_reduced_genome_mafs:
    # Concatenate all the reduced genome MAFs
    # readme.txt:39
    input:
        expand(
            'processed_data/GAFexonReduced/{genome_id}.broadgaf.wig.genome.reduced.maf',
            genome_id=GENOME_SAMPLE_IDS_TO_EXOME.keys(),
        )
    output:
        'processed_data/genome.broadbed.gafe.wigs.maf'
    shell:
        "cat {input} | sort -u > {output}"


# rule genome_wig_to_bed:
#     # Convert the genome sequencing coverage WIG as BED format
#     input: config['FAKE_ICGC46_WIG']
#     output: 'processed_data/fake_icgc46.wig.bed'
#     shell:
#         # New way: 
#         # r"wig2bed --zero-indexed < {input} | "
#         # # chr1 -> 1 naming conversion and remove non canonical chroms
#         # r"sed -E -e 's/^chr([0-9]+|X|Y|M)/\1/g' -e '/^chr/d'> {output}" 


rule reduce_maf_low_genome_wig:
    # Remove variants from either sequeing type not covered in the whole genome
    # wig e.g. no sequencing coverage in whole genome sequencing.
    #
    # Note that in whole genome sequencing we don't have the BAM(Wig) file for
    # each paired sample, so we rely on a "pesudo" joint wig file.
    # readme.txt:42 and 44 
    input:
        maf='processed_data/{seq_type}.broadbed.gafe.wigs.maf',
        # genome_bed='processed_data/fake_icgc46.wig.bed',
        genome_bed=config['FAKE_ICGC46_WIG_BED'],
    output: 'processed_data/{seq_type}.broadbed.gafe.wigs.rep.maf'
    shell:
        "bedtools intersect -a {input.maf} -b {input.genome_bed} | sort -u > {output}"


rule reduce_maf_low_genome_wig_uniq_snv:
    # Uniqe and sort the MAF
    input: 'processed_data/{seq_type}.broadbed.gafe.wigs.rep.maf'
    output: 'processed_data/{seq_type}.broadbed.gafe.wigs.rep.snv.uniq.maf'
    shell:
        r"""
        awk -F"\t" '{{if($2 == $3) print $0}}' {input} \
            | sort -u | sort -k1,1 -k2,2n \
            > {output}
        """


rule bdeolap_nchar:
    # Overlapp
    # readme.txt:66 and 68 
    input: 
        sample_mapping=config['OVERLAPPED_SAMPLE_INFO'],
        py_script='scripts/bedolap.{nchar}.py',
        exome_maf='processed_data/exome.broadbed.gafe.wigs.rep.snv.uniq.maf',
        genome_maf='processed_data/genome.broadbed.gafe.wigs.rep.snv.uniq.maf',
    output:  
        'output/e.g.mymerge.{nchar}.maf'        
    shell:
        'python {input.py_script} {input.sample_mapping} {input.exome_maf} {input.genome_maf} > {output}'


rule all:
    input:
        # expand(
        #     'processed_data/{seq_type}.broadbed.gafe.wigs.rep.snv.uniq.maf',
        #     seq_type=['exome', 'genome']
        # )
        expand(
            'output/e.g.mymerge.{nchar}.maf',
            nchar=['12char', '25char']
        )

