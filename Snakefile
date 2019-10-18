configfile: 'config.yaml'

from pathlib import Path

# Path to the WES based MAF (TCGA) and WGS based MAF (ICGC)
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

# Set up the mapping between exome and genome IDs 
EXOME_SAMPLE_IDS_TO_GENOME = {}
GENOME_SAMPLE_IDS_TO_EXOME = {}

with open(config['OVERLAPPED_SAMPLE_INFO']) as f:
    header = next(f)
    for line in f:
        pcawg_aliquot, mc3_aliquot, exome_id, genome_id = line.rstrip().split('\t')
        EXOME_SAMPLE_IDS_TO_GENOME[exome_id] = GENOME_ALIQUOT_TO_PCAWG_DONORID[pcawg_aliquot]
        GENOME_SAMPLE_IDS_TO_EXOME[GENOME_ALIQUOT_TO_PCAWG_DONORID[pcawg_aliquot]] = exome_id

# Find all the wigs and map them from their exome id.
EXOME_SAMPLE_IDS_TO_WIG = {}
for pth in Path(config['EXOME_WIG_DIR']).glob('*.coverage.wig.txt.gz'):
    aliquot_id_0, *_, aliquot_id_1, aliquot_id_2 = pth.name.split('.', 1)[0].split('_')
    exome_id = EXOME_ALIQUOT_TO_ID.get(aliquot_id_0) or EXOME_ALIQUOT_TO_ID.get(aliquot_id_1) or EXOME_ALIQUOT_TO_ID.get(aliquot_id_2)
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

# Currently 2 samples miss their coverage files:
# TCGA-B5-A11I-01A-11D-A10M-09 and TCGA-B6-A0RT-01A-21D-A099-09

# For CADD plots
CADDVARS=["cHmmBivFlnk","cHmmEnhBiv","cHmmEnhG","cHmmEnh","cHmmHet","cHmmQuies","cHmmReprPC","cHmmReprPCWk","cHmmTssAFlnk","cHmmTssA","cHmmTssBiv","cHmmTxFlnk","cHmmTx","cHmmTxWk","cHmmZnfRpts","CpG","dnaHelT","dnaMGW","dnaProT","dnaRoll","EncExp","EncH3K27Ac","EncH3K4Me1","EncH3K4Me3","EncNucleo","EncOCCombPVal","EncOCctcfPVal","EncOCctcfSig","EncOCDNasePVal","EncOCDNaseSig","EncOCFairePVal","EncOCFaireSig","EncOCmycPVal","EncOCmycSig","EncOCpolIIPVal","EncOCpolIISig","ESP_AFR","ESP_AF","ESP_EUR","fitConsc","GC","GerpN","GerpRSpval","GerpRS","GerpS","mamPhCons","mamPhyloP","mapAbility20bp","mapAbility35bp","mirSVR_E","mirSVR_Score","motifDist","motifEScoreChng","PHRED","PolyPhenVal","priPhCons","priPhyloP","RawScore","relcDNApos","relCDSpos","relProts","scoreSegDup","SIFTval","TFBSPeaksMax","TG_AFR","TG_AF","TG_AMR","TG_ASN","TG_EUR","verPhCons","verPhyloP"]


#These rules take wigs and puts them into an exome regions only BED file format.
rule convert_exome_wig_to_bed:
    """Convert one WIG to BED and subset to be within ROI."""
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
    """Convert one exome WIG to BED and subset to be within ROI."""
    input:
        wig=lambda wildcards: GENOME_SAMPLE_IDS_TO_WIG[wildcards.exome_id],
        myroi=config['COMBINED_ROI_BED'],

    output: str(Path(config['GENOME_WIG_BED_DIR'], '{exome_id}.bed'))
    shell:
        "gunzip -c {input.wig} | tr -s \"\n\" | wig2bed -d | bedtools intersect -a stdin -b {input.myroi} | python scripts/condenseBed2Wig.py > {output}"

# This step missed a couple ids that I need to put back
# gunzip -c /diskmnt/Datasets/PCAWG/Coverage/compressed_coverage_tracks/42de1441-6f3c-4b4d-b8b8-ad91e4b1dbe2.coverage.gz | tr -s '\n' | wig2bed -d | bedtools intersect -a stdin -b /diskmnt/Projects/ICGC_MC3/Reduction_Beds/combined.bitgt.gencode19.merged.bed | python scripts/condenseBed2Wig.py > ../Coverage_Reduction_Genome/TCGA-B5-A11G-01A-13D-A122-09.bed


# Temporary list of exome IDs without WIG
with open('./genome_ids_without.wig.list') as f:
    G_WIG_BLACKLIST = f.read().splitlines()


rule all_genome_beds:
    input: expand(str(Path(config['GENOME_WIG_BED_DIR'], '{exome_id}.bed')), exome_id=[x for x in GENOME_SAMPLE_IDS_TO_EXOME.values() if x not in G_WIG_BLACKLIST])


# Reduce to samples STEP 1
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
        with open(output[0], "w") as g:
            with open(input[0], "r") as f:
                ghead = f.readline()[:-1].split("\t")
                for line in f:
                    l = line[:-1].split("\t")
                    gid = l[42]
                    gene = l.pop(0)
                    l.append(gene)
                    if gid in GENOME_SAMPLE_IDS_TO_EXOME:
                        g.write("\t".join(l))
                        g.write("\n")


# Get the original header to pase on later
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


# Reduce to broad target bed STEP 2
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


# Reduce to the gencode 19 exon regions STEP 3
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


# Split up mutations by sample nice AWK command
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

# STEP 4: reduce the variants by reciprocal wig files,
# e.g., MC3 variants by PCAWG wig and PCAWG variants by MC3 wig.
rule reduce_exome_maf_low_coverage_wig:
    """
    Take split exome MAF and remove any variants not in whole exome wig
    (no sequencing coverage)

    Adapted from reduceGAFexonsMaf.py
    """
    input:
        exome_maf='processed_data/GAFexon/{exome_id}.exome.broadgaf.maf',
        genome_wig_bed=str(Path(config['GENOME_WIG_BED_DIR'], '{exome_id}.bed'))
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


# Concatenate all the reduced exome and genome MAFs together 
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

# STEP 5: reduce the genome and exome findings by the genome wiggle files
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


#This step makes the SQLite data base
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
        db='output/mc3_icgc.sqlite'
    shell:
        '''
        sed -e 's#\tSTRAND\t#\tSTRAND_VEP\t#g' < {input.exome_maf_header} > processed_data/exome.maf.for_db.header
        mkfifo {params.exome_pipe}
        cut -f 1-111 {input.exome_maf} > {params.exome_pipe} &
        mkfifo {params.genome_pipe}
        cut -f 1-43 {input.genome_maf} > {params.genome_pipe} &
        sqlite3 {output} < {params.sql}

        rm -f {params.exome_pipe} {params.genome_pipe}
        '''

#This is the rule that does the merging
rule make_analysis_file:
    input:
        db=rules.make_sqlite_db.output['db'],
        sql_script='scripts/gen_full.sql'
    output:
        'output/full.tsv'
    shell:
        '''
        sqlite3 {input.db} < {input.sql_script}
        '''

#This script cleans up duplicated, and other merge products
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
    input: rules.clean_analysis_file.output, rules.make_analysis_file.output


#Below are the various rules and analysis algorithms used to generate figures

#Generates an upsetR plot
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

#This is the gut for generating the figure found in MAFit
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

#This is the bar chart that was analyzed by cancer type
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

#This is code provided by William Meyerson to run the simulations
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

#Generates the clonality figures
rule clonality_figure:
#So at this step just run this in the GenFigures Clonality section with the files that I need. Anywho.
    input:
        makeClone='scripts/make_clonality.R',
        fullClone='/diskmnt/Projects/ICGC_MC3/Gen_Figures/Clonality_Figure/Processed_data/Full_Clonality.tsv'
    output:
        gClonePDF='figures/clonality.ICGC.subclusters.pdf',
        eClonePDF='figures/clonality.TCGA.subclusters.pdf'

    shell:
        '''
        Rscript --quiet --vanilla {input.makeClone} {input.fullClone} {output.gClonePDF} {output.eClonePDF}
        '''

#Generates the VAF figures that are characterized by different VAF bins
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


# Produce a figure that looks into the MC3 controlled MAF to identity filters MC3 variants that were replicated by ICGC
rule inverse_controlled_processing:
    input:
        full='output/full_cleaned.tsv',
        mc3ControlMAF=config['MC3_CONTROLLED'],
        makeControlled='scripts/capture_mc3controlled_icgcuniq.py'
    output:
        MC3controlledInPCAWG='processed_data/mc3_controlled_maf.of.icgcUnique.txt',
        PCAWGuniqNOTinMC3='processed_data/unique_icgc.notin.mc3_controlled.txt',
        PCAWGuniqInMC3='processed_data/unique_icgc.in.mc3_controlled.txt'
    shell:
        '''
        python {input.makeControlled} {input.full} {input.mc3ControlMAF} {output.MC3controlledInPCAWG} {output.PCAWGuniqNOTinMC3} {output.PCAWGuniqInMC3}
        '''


#Here we looked into the effects of singles variant callers
rule single_caller_figure:
    input:
        plotR="scripts/make_singlecaller.R",
        mc3control="processed_data/mc3_controlled_maf.of.icgcUnique.txt"
    output:
        singlePDF="figures/single.caller.filter.pdf"
    shell:
        '''
        Rscript --quiet --vanilla {input.plotR} {input.mc3control} {output.singlePDF}
        '''


#MC3 produced many variant call filters. This is how we approached this issue
rule filter_figure:
    input:
        likertProc="scripts/likert_reformat.py",
        full='output/full_cleaned.tsv',
        makelikert="scripts/make_likerFilter.R",
    output:
        likert='processed_data/data.4.likert.txt',
        likertPlot='figures/MC3filtersInPCAWG.likert.pdf'
    shell:
        '''
        python {input.likertProc} {input.full} > {output.likert}

        Rscript --quiet --vanilla {input.makelikert} {output.likert} {output.likertPlot}
        '''


#This is processing step to add CADD annotations to the mutations
rule CADD_annotation_split:
    input:
        full='output/full_cleaned.tsv',
        splitfull='scripts/divide_full_forCADD.R'
    output:
        o0='processed_data/psudoVCF.SNPs.MC3_PCAWG.txt',
        o1='processed_data/vars_chr1-chr6.psuedo.vcf',
        o2='processed_data/vars_chr6-chr12.psuedo.vcf',
        o3='processed_data/vars_chr12-chrY.psuedo.vcf'
    shell:
        '''
        Rscript --quiet --vanilla {input.splitfull} {input.full} {output.o1} {output.o2} {output.o3} {output.o0}
        '''


#And applied using VEP
rule CADD_annotation_vep:
    input:
        'processed_data/{cadd}.psuedo.vcf',
    output:
        'processed_data/{cadd}.psuedo.annotated.txt'
    shell:
        '''
        docker run --rm -t -i -u `id -u`:`id -g` -v /diskmnt/Datasets/VEP:/home/vep/.vep -v /diskmnt/Datasets/CADD:/annotations/CADD -v $PWD:/data ensemblorg/ensembl-vep:release_91.3 vep --cache --offline --fork 4 --minimal --pick --assembly GRCh37  -i /data/{input} -o /data/{output} --plugin CADD,/annotations/CADD/whole_genome_SNVs_inclAnno.tsv.gz
        '''

#This step combines the two previous steps.
rule combine_CADD:
    input:
        tocat=expand('processed_data/{cadd}.psuedo.annotated.txt', cadd=[i for i in ["vars_chr1-chr6","vars_chr6-chr12","vars_chr12-chrY"]]),
        reform='scripts/reformat_anno.py'
    output:
        'processed_data/CADD.annotations.fullvars.txt'
    shell:
        '''
        python {input.reform} {input.tocat} > {output}
        '''

#And this produces figures looking into GC content using CADD
rule CADD_plot:
    input:
        full='output/full_cleaned.tsv',
        makeCADD='scripts/make_cadd.R',
        caddanno=rules.combine_CADD.output[0],
        unique_p='processed_data/unique_icgc.notin.mc3_controlled.txt',
    output:
        dynamic("figures/CADD/{caddvars}_PCAWG_unique.pdf")
    shell:
        '''
        Rscript --quiet --vanilla {input.makeCADD} {input.full} {input.unique_p} {input.caddanno} figures/CADD/
        '''

#This will bring the previous 3 rules into one figure
rule CADD_figure:
    input: dynamic("figures/CADD/{caddvars}_PCAWG_unique.pdf")


#The next 4 rules look at various aspects of incorporating CADD stats to these mutations
rule CADD_plot_depth:
    input:
        full='output/full_cleaned.tsv',
        makeCADD='scripts/make_cadd_depth.R',
        caddanno='processed_data/CADD.annotations.fullvars.txt',
        unique_p='processed_data/unique_icgc.notin.mc3_controlled.txt',
    output:
        dynamic("figures/CADD/{caddvars}_depth.pdf")
    shell:
        '''
        Rscript --quiet --vanilla {input.makeCADD} {input.full} {input.unique_p} {input.caddanno} figures/CADD/
        '''

rule CADD_depth:
    input: dynamic("figures/CADD/{caddvars}_depth.pdf")


rule CADD_plot_depth_sd:
    input:
        full='output/full_cleaned.tsv',
        makeCADD='scripts/make_cadd_sd.R',
        caddanno='processed_data/CADD.annotations.fullvars.txt',
        unique_p='processed_data/unique_icgc.notin.mc3_controlled.txt',
    output:
        dynamic("figures/CADD/{caddvars}_depth_sd.pdf")
    shell:
        '''
        Rscript --quiet --vanilla {input.makeCADD} {input.full} {input.unique_p} {input.caddanno} figures/CADD/
        '''

rule CADD_sd:
    input: dynamic("figures/CADD/{caddvars}_depth_sd.pdf")


#This is a figure not shown in the manuscript that is a general gene model to determine if there are systematic differences at the beginning, end or middle of a gene.
rule general_gene_annotate:
#Note that this takes a different python environment that Liang-Bo set up.
    input:
        gencode=config['GENCODE_DB'],
        full='output/full_cleaned.tsv',
        annotate="scripts/general_gene_annotate.py",
        enviroPython="/diskmnt/Projects/Users/lwang/miniconda3/envs/genemodel/bin/python"
    output:
        annolog="logs/genemodel.annotation.log",
        gen_gene="processed_data/general_gene_model_annotation.tsv.gz",
    shell:
        '''
        {input.enviroPython} {input.annotate} {input.gencode} {input.full} 2> {output.annolog} | gzip -c > {output.gen_gene}
        '''


#Grouped with rule above
rule general_gene_figure:
    input:
        makeGenGene='scripts/make_genGene.R',
        full='output/full_cleaned.tsv',
        gen_gene="processed_data/general_gene_model_annotation.tsv.gz",
        unique_p='processed_data/unique_icgc.notin.mc3_controlled.txt',
    output:
        gen_gene_unzip = 'processed_data/general_gene_model_annotation.tsv',
        ggProportionPDF='figures/Proportional_genemodel.pdf',
        ggCountPDF='figures/Count_genemodel.pdf'
    shell:
        '''
        gunzip {input.gen_gene}

        Rscript --quiet --vanilla {input.makeGenGene} {input.full} {output.gen_gene_unzip} {input.unique_p} {output.ggProportionPDF} {output.ggCountPDF}
        '''


#This is a cheap way to generate a sunburst plot using R
rule sunburst_figure:
    input:
        prepDat='scripts/maf_complement.py',
        g_broadgaf_maf='processed_data/genome_broadbed.gaf.maf',
        g_self_maf='processed_data/genome.broadbed.gafe.wigs.self.maf',
        pseudo=config['PSEUDO_GENES_LIST'],
        makeSunburst='scripts/make_sunburst.R',
        can299=config['CANCER_299']
    output:
        rem_by_covg='processed_data/removed.by.coverage.genome.maf',
        utr3_canPDF='figures/UTR3.Cancer.Histogram.pdf',
        utr5_canPDF='figures/UTR5.Cancer.Histogram.pdf',
        miss_canPDF='figures/MISS.Cancer.Histogram.pdf',
        sunburstPDF='figures/Coverage_sunburst.pdf'
    shell:
        '''
        python {input.prepDat} {input.g_self_maf} {input.g_broadgaf_maf} > {output.rem_by_covg}

        Rscript --quiet --vanilla {input.makeSunburst} {output.rem_by_covg} {input.pseudo} {input.can299} {output.utr3_canPDF} {output.utr5_canPDF} {output.miss_canPDF} {output.sunburstPDF}
        '''


# Next few rules look at how CADD was integrated into this manuscript.
rule CADD_covg_split:
    input:
        rem_by_covg='processed_data/removed.by.coverage.genome.maf',
        splitcovg='scripts/divide_covg_forCADD.R'
    output:
        o0='processed_data/covgVCF.SNPs.MC3_PCAWG.txt',
        o1='processed_data/vars_chr1-chr6.covg.vcf',
        o2='processed_data/vars_chr6-chr12.covg.vcf',
        o3='processed_data/vars_chr12-chrY.covg.vcf'
    shell:
        '''
        Rscript --quiet --vanilla {input.splitcovg} {input.rem_by_covg} {output.o1} {output.o2} {output.o3} {output.o0}
        '''

rule CADD_covg_annotation_vep:
    input:
        'processed_data/{cadd}.covg.vcf',
    output:
        'processed_data/{cadd}.covg.annotated.txt'
    shell:
        '''
        docker run --rm -t -i -u `id -u`:`id -g` -v /diskmnt/Datasets/VEP:/home/vep/.vep -v /diskmnt/Datasets/CADD:/annotations/CADD -v $PWD:/data ensemblorg/ensembl-vep:release_91.3 vep --cache --offline --fork 4 --minimal --pick --assembly GRCh37  -i /data/{input} -o /data/{output} --plugin CADD,/annotations/CADD/whole_genome_SNVs_inclAnno.tsv.gz
        '''

rule combine_covg_CADD:
    input:
        tocat=expand('processed_data/{cadd}.covg.annotated.txt', cadd=[i for i in ["vars_chr1-chr6","vars_chr6-chr12","vars_chr12-chrY"]]),
        reform='scripts/reformat_anno.py'
    output:
        'processed_data/CADD.annotations.vars.rem_by_covg.txt'
    shell:
        '''
        python {input.reform} {input.tocat} > {output}
        '''

rule covg_GC_depth_figures:
    input:
        plotr='scripts/make_gc_depth.R',
        full='output/full_cleaned.tsv',
        rmcovg='processed_data/removed.by.coverage.genome.maf',
        caddcovg='processed_data/CADD.annotations.vars.rem_by_covg.txt',
        pcawg_uniq='processed_data/unique_icgc.notin.mc3_controlled.txt',
        caddfull='processed_data/CADD.annotations.fullvars.txt'
    output:
        density_pdf='figures/denisty_GC.pdf',
        gc_depth_binned='figures/binned_gc_depth.pdf'
    shell:
        '''
        Rscript --quiet --vanilla {input.plotr} {input.rmcovg} {input.caddcovg} {input.full} {input.pcawg_uniq} {input.caddfull} {output.gc_depth_binned} {output.density_pdf}
        '''


rule covg_CADD_figures:
    input:
        plotr='scripts/make_gc_depth.all.R',
        full='output/full_cleaned.tsv',
        rmcovg='processed_data/removed.by.coverage.genome.maf',
        caddcovg='processed_data/CADD.annotations.vars.rem_by_covg.txt',
        pcawg_uniq='processed_data/unique_icgc.notin.mc3_controlled.txt',
        caddfull='processed_data/CADD.annotations.fullvars.txt'
    output:
        dynamic("figures/CADD/{caddvars}_depth_sd.pdf")
    shell:
        '''
        Rscript --quiet --vanilla {input.plotr} {input.rmcovg} {input.caddcovg} {input.full} {input.pcawg_uniq} {input.caddfull} {output.density_pdf} {output.gc_depth_binned}
        '''


rule CADD_covg:
    input: dynamic("figures/CADD/{caddvars}_covg.pdf")


#This is just a side script for converting bed files to wiggle files.
rule bed_wigs_genome:
    input:
        bed=str(Path(config['GENOME_WIG_BED_DIR'], '{exome_id}.bed')),
        bed2wig='scripts/bed2fixedwig.py'
    output:
        str(Path(config['GENOME_REDUCED_WIG'], '{exome_id}.wig'))
    shell:
        '''
        python {input.bed2wig} {input.bed} {output}

        '''


#We wanted to look at SMG using MuSiC to determine if we could pick up some non-coding hits. The next few rule prep input for that effort
rule gen_SMG_wigs:
    input:
        expand(str(Path(config['GENOME_REDUCED_WIG'], '{exome_id}.wig')), exome_id=[x for x in GENOME_SAMPLE_IDS_TO_EXOME.values() if x not in G_WIG_BLACKLIST])


rule subset_MAF_byCancer:
    #This script needs to split up the MAF and generate a wig_list
    input:
        maf='processed_data/genome_broadbed.gaf.maf',
        samplemap=config['OVERLAPPED_SAMPLE_INFO'],
        splitMAF='scripts/divide_genomeMAF_bycan.R'
    output:
        dynamic("processed_data/SMG/{cantype}.genome.reduced.maf")
    shell:
        '''
        Rscript --quiet --vanilla {input.splitMAF} {input.maf} {input.samplemap}
        '''


rule get_ROI:
    input:'/diskmnt/Software/bin/2020plus/data/intersection_snvbox_gafbitgt.unique.bed'

    output:'processed_data/SMG/input/2020plus_roi'

    shell:
        '''
        cut -f 1-4 {input} > {output}
        '''


#Look into mutation spectrum,
rule mutation_spectrum_figure:
    input:
        mc3reduced2exonsMAF='processed_data/exome.broadbed.gaf.maf',
        pcawgreduced2exonMAF='processed_data/genome_broadbed.gaf.maf',
        mkMutSpec='scripts/make_mutspec.R'
    output:
        mc3_mutspec='figures/mc3_mutspec.pdf',
        pcawg_mutspec='figures/pcawg_mutspec.pdf',
        mutspec_notes='processed_data/mutspecNotes.txt'
    shell:
        '''
        Rscript --quiet --vanilla {input.mkMutSpec} {input.mc3reduced2exonsMAF} {input.pcawgreduced2exonMAF} {output.mc3_mutspec} {output.pcawg_mutspec} {output.mutspec_notes}
        '''


#Looking in the single, tri-nucleotide, and indel differences in these data
rule snp_tnp_indel_figure:
    input:
        full='output/full_cleaned.tsv',
        dnpfig='scripts/make_onp_indel.R'
    output:
        snpbar='figures/snp_dnp_tnp_indel.pdf'
    shell:
        '''
        Rscript --quiet --vanilla {input.dnpfig} {input.full} {output.snpbar}
        '''


# This is just looking into exome only variants
rule lowVAFexome_keys:
    input:
        full='output/full_cleaned.tsv',
        tcgaUniq='scripts/make_tcgaOnly.R',
        can299=config['CANCER_299']
    output:
        tcgaOnly='processed_data/advantagesOfExome.txt'
    shell:
        '''
        Rscript --quiet --vanilla {input.tcgaUniq} {input.full} {input.can299} {output.tcgaOnly}
        '''

# This rule looks into Broad/Mutect and MUSE specific calls in the dataset. 
rule museonly_figure:
    input:
        full='output/full_cleaned.tsv',
        genMuse='scripts/make_museonly.R'
    output:
        likertmuse='figures/museOnly.pdf'
    shell: 
        '''
        Rscript --quiet --vanilla {input.genMuse} {input.full} {output.likertmuse}
        '''

#This rule looks into the correlation between MutPerBP and Concordance 
#This rule still needs to be run "by-hand"

rule mutpmb_concordance_figure:
    input:
        prepScript='scripts/genCallerMatrix.cancer.py',
        idcancer=config['TCGA12_CANCER'],
        full='output/full_cleaned.tsv',
        genCorr='scripts/make_corr.R',
        mutimmune='data/MUTperMB.immunepaper.txt'
    output:
        datPrep='processed_data/data.4.cancerOLAP.txt',
        mxcfig='figures/CancerXConcordance.pdf'

    shell:
        '''
        python {input.prepScript} {input.full} > {output.datPrep}

        Rscript --quiet --vanilla {input.genCorr} {output.datPrep} {input.idcancer} {input.mutimmine} {output.mxcfig}
        '''


rule barchart_overview_figure:
#Currently variation exists for these reasons. I need to track down these variants and add them the is algorithm. VAF, (X)Sublonal, Filtering, Caller, GC content, Cancertype,
#
    input:
        genAll='scripts/make_allbar.R',
        full='output/full_cleaned.tsv',
        fullClone='processed_data/Full_Clonality.tsv',
        cantype=config['TCGA12_CANCER']
        caddfull='processed_data/CADD.annotations.fullvars.txt'
    output:
        allfig='figures/allbar.v1.pdf'
    shell:
        '''
        Rscript --quiet --vanilla {input.genCorr} 
        '''



#This final bit of code, can be added to un-commented to the run many of these rules after the full_cleaned.tsv was generated.
rule all_figures:
    input:
        rules.upSetR_olap_figure.output,
        rules.landscape_sample_figure.output,
        rules.landscape_cancer_figure.output,
        rules.simulation_figure.output,
        rules.clonality_figure.output,
        rules.vaf_figure.output,
        rules.inverse_controlled_processing.output,
        rules.single_caller_figure.output,
        rules.filter_figure.output,
        rules.CADD_figure.output,
        rules.CADD_depth.output,
        rules.CADD_sd.output,
        # rules.general_gene_figure.output,
        rules.sunburst_figure.output,
        rules.covg_GC_depth_figures.output,
        rules.CADD_covg.output,
        rules.mutation_spectrum_figure.output,
        rules.snp_tnp_indel_figure.output,
        rules.museonly_figure.output,
        rules.mutpmb_concordance_figure.output,
        rules.barchart_overview_figure.output
