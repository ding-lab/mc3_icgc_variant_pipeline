configfile: 'config.yaml'

from pathlib import Path

GENCODE_PROTEIN_CODING_BED = 'annotations/gencode.proteincoding.known.gene.bed'

EXOME_MAF = config['EXOME_MAF']
GENOME_MAF = config['GENOME_MAF']

# Set up GENOME MAPPING FOR MAF
GENOME_ALIQUOT_TO_PCAWG_DONORID = {}
with open(config['GENOME_ALIQUOT_INFO']) as f: 
    header = next(f)
    for line in f: 
        _, __, submitter_donor_id, donor_id, _____, aliquot_id, *_ = line.rstrip().split('\t')
        GENOME_ALIQUOT_TO_PCAWG_DONORID[aliquot_id] = donor_id 


# Set up EXOME MAPPING FOR WIGS 
EXOME_ALIQUOT_TO_ID = {}
with open(config['EXOME_ALIQUOT_INFO']) as f:
    header = next(f)
    for line in f:
        _, __, aliquot_id, ____, _____, exome_id, *_ = line.rstrip().split('\t')
        EXOME_ALIQUOT_TO_ID[aliquot_id] = exome_id

# Set up the EXOME 2 GENOME MAPPING
EXOME_SAMPLE_IDS_TO_GENOME = {}
GENOME_SAMPLE_IDS_TO_EXOME = {}

with open(config['OVERLAPPED_SAMPLE_INFO']) as f:
    header = next(f)
    for line in f:
        pcawg_aliquot, mc3_aliquot, exome_id, genome_id = line.rstrip().split('\t')
 #       if exome_id in SUBSETS: 
        EXOME_SAMPLE_IDS_TO_GENOME[exome_id] = GENOME_ALIQUOT_TO_PCAWG_DONORID[pcawg_aliquot]
        GENOME_SAMPLE_IDS_TO_EXOME[GENOME_ALIQUOT_TO_PCAWG_DONORID[pcawg_aliquot]] = exome_id


# Find all the wigs and map them from their exome id.
EXOME_SAMPLE_IDS_TO_WIG = {}
for pth in Path(config['EXOME_WIG_DIR']).glob('*.coverage.wig.txt.gz'):
    aliquot_id_0, *_, aliquot_id_1, aliquot_id_2 = pth.name.split('.', 1)[0].split('_')
    exome_id = EXOME_ALIQUOT_TO_ID.get(aliquot_id_0) or  EXOME_ALIQUOT_TO_ID.get(aliquot_id_1) or EXOME_ALIQUOT_TO_ID.get(aliquot_id_2)
    if exome_id is not None:
        EXOME_SAMPLE_IDS_TO_WIG[exome_id] = str(pth)

# Find all the wigs and map them from genome id
GENOME_SAMPLE_IDS_TO_WIG = {}
for pth in Path(config['GENOME_WIG_DIR']).glob('*.coverage.gz'):
    aliquot_id = pth.name.split('.', 1)[0]
    if GENOME_ALIQUOT_TO_PCAWG_DONORID[aliquot_id] in GENOME_SAMPLE_IDS_TO_EXOME:
        exome_id = GENOME_SAMPLE_IDS_TO_EXOME[GENOME_ALIQUOT_TO_PCAWG_DONORID[aliquot_id]]
        if exome_id is not None:
            GENOME_SAMPLE_IDS_TO_WIG[exome_id] = str(pth)

#TWO SAMPLES ARE MISSING COVERAGE FILES 
#      1 TCGA-B5-A11I-01A-11D-A10M-09
#      1 TCGA-B6-A0RT-01A-21D-A099-09


#These rules take wigs and puts them into an exome regions only BED file format. 
rule convert_exome_wig_to_bed:
    input: 
        wig=lambda wildcards: EXOME_SAMPLE_IDS_TO_WIG[wildcards.exome_id],
        myroi=config['COMBINED_ROI_BED'],

    output: str(Path(config['EXOME_WIG_BED_DIR'], '{exome_id}.bed'))
    shell:
        "gunzip -c {input.wig} | wig2bed -d | bedtools intersect -a stdin -b {input.myroi} | python scripts/condenseBed2Wig.py > {output}"

# Temporary list of exome IDs without WIG
with open('./exome_ids_without_wig.list') as f:
    WIG_BLACKLIST = f.read().splitlines()


rule all_exome_beds:
    input: expand(str(Path(config['EXOME_WIG_BED_DIR'], '{exome_id}.bed')), exome_id=[x for x in EXOME_SAMPLE_IDS_TO_GENOME.keys() if x not in WIG_BLACKLIST])


rule convert_genome_wig_to_bed:
    input: 
        wig=lambda wildcards: GENOME_SAMPLE_IDS_TO_WIG[wildcards.exome_id],
        myroi=config['COMBINED_ROI_BED'],

    output: str(Path(config['GENOME_WIG_BED_DIR'], '{exome_id}.bed'))
    shell:
        "gunzip -c {input.wig} | tr -s \"\n\" | wig2bed -d | bedtools intersect -a stdin -b {input.myroi} | python scripts/condenseBed2Wig.py > {output}"

#This step missed a couple ids that I need to put back 
# gunzip -c /diskmnt/Datasets/PCAWG/Coverage/compressed_coverage_tracks/42de1441-6f3c-4b4d-b8b8-ad91e4b1dbe2.coverage.gz | tr -s '\n' | wig2bed -d | bedtools intersect -a stdin -b /diskmnt/Projects/ICGC_MC3/Reduction_Beds/combined.bitgt.gencode19.merged.bed | python scripts/condenseBed2Wig.py > ../Coverage_Reduction_Genome/TCGA-B5-A11G-01A-13D-A122-09.bed


# Temporary list of exome IDs without WIG
with open('./genome_ids_without.wig.list') as f:
    G_WIG_BLACKLIST = f.read().splitlines()


rule all_genome_beds:
    input: expand(str(Path(config['GENOME_WIG_BED_DIR'], '{exome_id}.bed')), exome_id=[x for x in GENOME_SAMPLE_IDS_TO_EXOME.values() if x not in G_WIG_BLACKLIST])




#Reduce to samples STEP 1 
rule reduce_exome_maf:
    input: EXOME_MAF
    output: 'processed_data/exome.maf'
    run:
        with open(output[0], "w") as e:
            with open(input[0], "r") as f:
                ehead = f.readline()[:-1].split("\t")
                hugo = ehead.pop(0)
                ehead.append(hugo)
                del ehead[:3]
                for line in f:
                    l = line.strip().split("\t")
                    eid = l[15]
                    gene = l.pop(0)
                    l.append(gene)
                    del l[:3]
                    if eid in EXOME_SAMPLE_IDS_TO_GENOME:
                        e.write("\t".join(l))
                        e.write("\n")


rule reduce_genome_maf:
    input: GENOME_MAF
    output: 'processed_data/genome.maf'
    run:
        with open(output[0],"w") as g:
            with open(input[0],"r") as f:
                ghead = f.readline()[:-1].split("\t")
                for line in f:
                    l = line[:-1].split("\t")
                    gid = l[42]
                    gene = l.pop(0)
                    l.append(gene)
                    if gid in GENOME_SAMPLE_IDS_TO_EXOME:
                        g.write("\t".join(l))
                        g.write("\n")


#Get the orginal header to pase on later 
rule get_orig_maf_header_exome:
    input: EXOME_MAF
    output: 'processed_data/exome.maf.header'
    run:
        with open(output[0], "w") as e:
            with open(input[0], "r") as f:
                ehead = f.readline()[:-1].split("\t")
                hugo = ehead.pop(0)
                ehead.append(hugo)
                del ehead[:3]
                e.write("\t".join(ehead))


rule get_orig_maf_header_genome:
    input: GENOME_MAF
    output: 'processed_data/genome.maf.header'
    run:
        with open(output[0],"w") as g:
            with open(input[0],"r") as f:
                ghead = f.readline()[:-1].split("\t")
                hugo = ghead.pop(0)
                ghead.append(hugo)
                g.write("\t".join(ghead))


#Reduce to broad target bed
rule broadbed_exome:
    input:
        maf='processed_data/exome.maf',
        bed='../Reduction_Beds/gaf_20111020+broad_wex_1.1_hg19.bed',
    output: 'processed_data/exome.broadbed.maf'
    shell:
        "bedtools intersect -a {input.maf} -b {input.bed} -wa | "
        "sort -u > {output}"


rule broadbed_genome:
    input:
        maf='processed_data/genome.maf',
        bed='../Reduction_Beds/gaf_20111020+broad_wex_1.1_hg19.bed',
    output: 'processed_data/genome.broadbed.maf'
    shell:
        "bedtools intersect -a {input.maf} -b {input.bed} -wa | "
        "sort -u > {output}"

#Reduce to the gencode 19 exon regions 
rule gafbed_exome:
    input:
        maf='processed_data/exome.broadbed.maf',
        bed='../Reduction_Beds/gencode.v19.basic.exome.bed',
    output: 'processed_data/exome.broadbed.gaf.maf'
    shell:
        "bedtools intersect -a {input.maf} -b {input.bed} -wa | "
        "sort -u > {output}"


rule gafbed_genome:
    input:
        maf='processed_data/genome.broadbed.maf',
        bed='../Reduction_Beds/gencode.v19.basic.exome.bed',
    output: 'processed_data/genome_broadbed.gaf.maf'
    shell:
        "bedtools intersect -a {input.maf} -b {input.bed} -wa | "
        "sort -u > {output}"

#Split up mutations by sample nice AWK command
rule split_exome_maf:
    input: rules.gafbed_exome.output 
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
    input: rules.gafbed_genome.output
    output:
        expand(
            'processed_data/GAFexon/{sample_genome_id}.genome.broadgaf.maf',
            sample_genome_id=GENOME_SAMPLE_IDS_TO_EXOME.keys()
        )
    shell:
        r"""
        awk -F"\t" '{{
            print > "processed_data/GAFexon/"$42".genome.broadgaf.maf"
        }}' {input}
        """

#These steps take reduce the variants by resiprocol wig files MC3 by PCAWG and PCAWG by MC3 .
rule reduce_exome_maf_low_coverage_wig:
    """
    Take split exome MAF and remove any variants not in whole exome wig
    (no sequencing coverage)

    Adapted from reduceGAFexonsMaf.py
    """
    input:
        exome_maf='processed_data/GAFexon/{exome_id}.exome.broadgaf.maf',
        genome_wig_bed=str(Path(config['GENOME_WIG_BED_DIR'], '{exome_id}.bed')),
    output:
        'processed_data/GAFexonReduced/{exome_id}.broadgaf.wig.exome.reduced.maf'
    shell:
        'python scripts/reduceMAFusingWIG.py {input.exome_maf} {input.genome_wig_bed} exome > {output}'


def retrieve_exome_wig_bed(wildcards):
    # Given the genome sample id, return the paired exome sample id
    exome_id = GENOME_SAMPLE_IDS_TO_EXOME[wildcards.genome_id]
    return str(Path(
        config['EXOME_WIG_BED_DIR'],
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
        'python scripts/reduceMAFusingWIG.py {input.genome_maf} {input.exome_wig_bed} genome > {output}'


# Concatenate all the reduced exome and genome MAFs separately after reciprocol
rule merge_exome_wig_reduced_exome_mafs:
    input:
        expand(
            'processed_data/GAFexonReduced/{exome_id}.broadgaf.wig.exome.reduced.maf', 
            exome_id=[x for x in GENOME_SAMPLE_IDS_TO_EXOME.values() if x not in G_WIG_BLACKLIST],
        )
    output:
        'processed_data/exome.broadbed.gafe.wigs.maf'
    shell:
        "cat {input} > {output}"


rule merge_exome_wig_reduced_genome_mafs:
    input:
        expand(
            'processed_data/GAFexonReduced/{genome_id}.broadgaf.wig.genome.reduced.maf',
            genome_id=GENOME_SAMPLE_IDS_TO_EXOME.keys(),
        )
    output:
        'processed_data/genome.broadbed.gafe.wigs.maf'
    shell:
        "cat {input} > {output}"

#I was able to find WIGS For 745 of the 747 samples. These next two steps are irrelevant
#rule genome_wig_to_bed:

#However I did notice an interesting situation where BROAD didn't make a call but other tools did and broad said that the variant wasn't covered. 

#This will be an interesting exercise. 

#So I'm going to do a self reduction of the Mutations based on coverage that should work everything out. We'll have to see the difference between the two groups to get things done. 
rule self_reduce_exome:
    input:
        exome_maf='processed_data/GAFexonReduced/{exome_id}.broadgaf.wig.exome.reduced.maf',
        exome_wig_bed=str(Path(config['EXOME_WIG_BED_DIR'], '{exome_id}.bed')),
        genome_reduced='processed_data/genome.broadbed.gafe.wigs.maf',
    output:
        'processed_data/GAFselfReduced/{exome_id}.broadgaf.wig.exome.reduced.self.maf'
    shell:
        'python scripts/reduceMAFusingWIG.py {input.exome_maf} {input.exome_wig_bed} exome > {output}'


def retrieve_genome_wig_bed(wildcards):
    e = GENOME_SAMPLE_IDS_TO_EXOME[wildcards.genome_id]
    out = str(Path(
        config['GENOME_WIG_BED_DIR'],
        '{e}.bed'.format(e=e)
    ))
    return out

rule self_reduce_genome:
    input:
        genome_maf='processed_data/GAFexonReduced/{genome_id}.broadgaf.wig.genome.reduced.maf',
        genome_wig_bed=retrieve_genome_wig_bed,
        exome_reduced='processed_data/exome.broadbed.gafe.wigs.maf'
    output:
        'processed_data/GAFselfReduced/{genome_id}.broadgaf.wig.genome.reduced.self.maf'   
    shell:
        'python scripts/reduceMAFusingWIG.py {input.genome_maf} {input.genome_wig_bed} genome > {output}'


rule merge_self_exome:
    input:
        expand(
            'processed_data/GAFselfReduced/{exome_id}.broadgaf.wig.exome.reduced.self.maf',
            exome_id=[x for x in GENOME_SAMPLE_IDS_TO_EXOME.values() if x not in G_WIG_BLACKLIST],
        )
    output:
        'processed_data/exome.broadbed.gafe.wigs.self.maf'
    shell:
        'cat {input} > {output}'



rule merge_self_genome:
    input: 
        expand(
            'processed_data/GAFselfReduced/{genome_id}.broadgaf.wig.genome.reduced.self.maf',
            genome_id=[x for x in GENOME_SAMPLE_IDS_TO_EXOME.keys() if GENOME_SAMPLE_IDS_TO_EXOME[x] not in G_WIG_BLACKLIST],
        )
    output:
        'processed_data/genome.broadbed.gafe.wigs.self.maf'
    shell:
        'cat {input} > {output}'


#The next two rules can be used for to generate sorted SNV only mafs if you want. 
rule reduce_maf_low_genome_wig_uniq_snv:
    # Uniqe and sort the MAF
    input: 'processed_data/{seq_type}.broadbed.gafe.wigs.self.maf'
    output: 'processed_data/{seq_type}.broadbed.gafe.wigs.snv.uniq.self.maf'
    shell:
        r"""
        awk -F"\t" '{{if($2 == $3) print $0}}' {input} \
            | sort -u -T /diskmnt/Projects/ICGC_MC3/tmp -S 1G | sort -k1,1 -k2,2n -T /diskmnt/Projects/ICGC_MC3/tmp -S 1G \
            > {output}
        """


rule bdeolap_nchar:
    # Overlapp
    # readme.txt:66 and 68
    # Not needed in current version (20180725)
    input:
        sample_mapping=config['OVERLAPPED_SAMPLE_INFO'],
        py_script='scripts/bedolap.{nchar}.py',
        exome_maf='processed_data/exome.broadbed.gafe.wigs.snv.uniq.self.maf',
        genome_maf='processed_data/genome.broadbed.gafe.wigs.snv.uniq.self.maf',
    output:
        'output/e.g.mymerge.{nchar}.maf'
    shell:
        'python {input.py_script} {input.sample_mapping} {input.exome_maf} {input.genome_maf} > {output}'


rule make_sqlite_db:
    input:
        exome_maf='processed_data/exome.broadbed.gafe.wigs.self.maf',
        exome_maf_header=rules.get_orig_maf_header_exome.output,
        genome_maf='processed_data/genome.broadbed.gafe.wigs.self.maf',
        genome_maf_header=rules.get_orig_maf_header_genome.output,
        genome_id_mapping=config['GENOME_ALIQUOT_INFO_FOR_SQLITE']
    params:
        exome_pipe='processed_data/exome.pipe',
        genome_pipe='processed_data/genome.pipe',
        sql='scripts/make_db.sql'
    output:
        'output/mc3_icgc.sqlite'
    shell:
        '''
        sed -e 's#\tSTRAND\t#\tSTRAND_VEP\t#g' < {input.exome_maf_header} > processed_data/exome.maf.for_db.header
        mkfifo {params.exome_pipe}
        cut -f 1-111 {input.exome_maf} > {params.exome_pipe} &
        
        cut -f 1-43 {input.genome_maf} > {params.genome_pipe} &
        sqlite3 {output} < {params.sql}

        rm {params.exome_pipe} {params.genome_pipe}
        '''
rule make_analysis_file:
    input: 
        db='output/mc3_icgc.sqlite',
        sql_script='scripts/gen_full.sql'
    output:
        'output/full.tsv'
    shell:
        '''
        sqlite3 {input.db} < {input.sql_script}
        
        mv full.tsv output
        '''


rule clean_analysis_file:
    input:
        full='output/full.tsv',
        Rscript='scripts/optimize_full.R'
    output:
        'output/full_cleaned.tsv'
    shell:
        '''
        Rscript --vanilla {input.Rscript} {input.full} {output}
        ''' 
    
rule all:
    input: rules.clean_analysis_file.output



rule upSetR_olap_figure:
    input: 
        prepScript='scripts/genCallerMatrix.py',
        makeUpsetR='scripts/make_upset.R',
        full='output/full_cleaned.tsv'
    output:
        datPrep='processed_data/data.4.upset.txt',
        upsetPDF='figures/upSetR.byCaller.pdf',
        upsetPDFall='figures/upSetR.byCaller.all.pdf'
    shell:
        '''
        python {input.prepScript} {input.full} > {output.datPrep}
       
        Rscript --quiet --vanilla {input.makeUpsetR} {output.datPrep} {output.upsetPDF} {output.upsetPDFall} 
        '''

rule landscape_sample_figure:
    input:
        full='output/full_cleaned.tsv',
        id_similar=config['ID_SIMILAR01'],
        makeLanR='scripts/make_lan.R'
    output:
        lanPDF='figures/landscape.sampleconcordance.pdf',
        lanNOTES='figures/landscape.sampleconcordance.notes.txt'
    shell: 
         '''
         Rscript --quiet --vanilla {input.makeLanR} {input.full} {input.id_similar} {output.lanPDF} {output.lanNOTES}
         '''

rule landscape_cancer_figure:
    input:
        prepScript='scripts/genCallerMatrix.cancer.py',
        makeCancerR='scripts/make_cancer.R',
        full='output/full_cleaned.tsv',
        idcancer=config['TCGA12_CANCER']
        
    output:
        datPrep='processed_data/data.4.cancerOLAP.txt',
        lanCancerPDF='figures/landscapteByCancertype.pdf'
       
    shell:
        '''
        python {input.prepScript} {input.full} > {output.datPrep}

        Rscript --quiet --vanilla {input.makeCancerR} {output.datPrep} {input.idcancer} {output.lanCancerPDF}
        '''

rule simulation_figure:
    input:
        full='output/full_cleaned.tsv',
        makeSimR='scripts/make_simulations.R',
        
    output:
        datPrep='processed_data/pp3.vaf.txt',
        simPDF='figures/WMeyerson_VAFpredictsRecoveryRate.pdf'

    shell:
        '''
        cat {input.full} | cut -f 1-2,112-113,6-9,117-120,12,37,38,123,148,149,152 > {output.datPrep}
        
        Rscript --vanilla {input.makeSimR} {output.datPrep} {input.full} {output.simPDF}
        '''


rule clonality_figure:
#So at this step just run this in the GenFigures Clonality section with the files that I need. Anywho. 
    input:
        makeClone='scripts/make_clonality.R',
        fullClone='/diskmnt/Projects/ICGC_MC3/Gen_Figures/Clonality_Figure/Processed_data/Full_Clonality.tsv'
    output:
        gClonePDF='figures/Clonality.ICGC.subclusters.v6.pdf',
        eClonePDF='figures/Clonality.TCGA.subclusters.v6.pdf'

    shell:
        '''
        Rscript --quiet --vanilla {input.makeClone} {input.fullClone} {output.gClonePDF} {output.eClonePDF} 
        '''

rule vaf_figure:
    input:
        full='output/full_cleaned.tsv',
        makeVAFbins='scripts/make_vafbins.R'
    output:
        vafBinPDF='figures/VAF.bins.pdf'
    shell:
        '''
        Rscript --quiet --vanilla {input.makeVAFbins} {input.full} {output.vafBinPDF}

        '''

#rule general_gene_figure:
#    input:
#        a
#    output:
#        b
#    shell:
#        '''
          
#        '''        
    


 
#rule all_figures:
#    input: 
#        rules.upSetR_olap_figure.output
#        rules.landscape_sample_figure.output
#        rules.landscape_cancer_figure.output
#        rules.simulation_figure.output
#        rules.figure_750_540_figure.output





















