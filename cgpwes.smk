# Rules in Snakefile based on https://github.com/cancerit/dockstore-cgpwgs/blob/develop/scripts/analysisWGS.sh

WORK="/work/sduvarcall/"
INPUT_BAM = WORK+"G90-2019-VHL-case-projekt/bam_cgpwes"
OUTPUT = WORK+"G90-2019-VHL-case-projekt/cgpwes_out"

# SAMPLES, = glob_wildcards(INPUT+"{sample}_R1_001.fastq.gz")

# NORMAL = "G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27"
# TUMOR = "G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22"

# NORMAL = "G56-blod-pt1_truseq-nano-genome_HGL2LDSXX_S5"
# TUMOR = "G56-sampleA1_truseq-nano-genome_HGL2LDSXX_S1"

configfile: "/work/sduvarcall/G90-2019-VHL-case-projekt/samples_matched_pairwise_cgpwes.yaml"

REF_FILES = "/work/sduvarcall/G37-2016-BRCAX-Genomes/"
ref_build = "GRCh37"

# Resources - paths inside docker
if ref_build == "GRCh37":
    # GRCh37 paths
    ref = REF_FILES+"reference_files/core_ref_GRCh37d5"
    res = REF_FILES+"reference_files"
else:
    # GRCh38 paths
    ref = REF_FILES+"/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
    res = REF_FILES+"/reference_files_GRCh38"

# Global variables
if ref_build == "GRCh37":
    SPECIES = "Human"   # Specify according to bam file e.g. Human or human
    ASSEMBLY = "NCBI37" # Specify according to bam file e.g. GRCh38 or NCBI37
    CONTIG_EXCLUDE = "NC_007605,hs37d5,GL%" # GRCh37
else:
    SPECIES = "human"   # Specify according to bam file e.g. Human or human
    ASSEMBLY = "GRCh38" # Specify according to bam file e.g. GRCh38 or NCBI37
    CONTIG_EXCLUDE = "HLA%,chrUn%" # GRCh38

PROTOCOL = "WXS"
CAVESPLIT = 350000

# CONTIG_EXCLUDE = ""
if CONTIG_EXCLUDE == "":
    CAVE_CONTIG_EXCLUDE = ""
    PINDEL_CONTIG_EXCLUDE = ""
else:
    CAVE_CONTIG_EXCLUDE = "-x " + CONTIG_EXCLUDE
    PINDEL_CONTIG_EXCLUDE = "-e " + CONTIG_EXCLUDE

ASCAT_ADD_ARGS='' # Let ASCAT compute purity and ploidy
# ASCAT_ADD_ARGS='-pu ASCAT_PURITY -pi ASCAT_PLOIDY'
# ASCAT_ADD_ARGS='-pu 0.8 -pi 3'

## CNV tool to use; choose either ascat or battenberg
CNV_TOOL = "ascat"
# CNV_TOOL = "battenberg"

# Get sample names from config file
SAMPLES = [sample for sample in config]
# SAMPLES = ["G90-VHL521-10-CNS1"]
print(SAMPLES)


onstart:
    shell("mkdir -p " + OUTPUT)
    shell("mkdir -p " + OUTPUT + "/finished")
    # shell("rm -f " + OUTPUT+"/finished/TEST_SUCCESS_"+TUMOR+"_vs_"+NORMAL+".txt")


rule all:
    input:
        [expand(OUTPUT+"/finished/SUCCESS__"+PROTOCOL+"__"+CNV_TOOL+"__{tumor}__vs__{normal}.txt", tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in SAMPLES]
        # [expand(OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.battenberg", tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in config]
        # expand(OUTPUT+"/finished/TEST_SUCCESS_{tumor}_vs_{normal}.txt", tumor=TUMOR, normal=NORMAL)
        # expand(OUTPUT+"/finished/SUCCESS_{tumor}_vs_{normal}.txt", tumor=TUMOR, normal=NORMAL)
        # expand(OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.CaVEMan_setup", tumor=TUMOR, normal=NORMAL)


###############################################################################
#### Helper functions
###############################################################################
#def get_normalContamination(samplestats):
    #with open(samplestats) as f:
def get_normalContamination(wildcards):
    if CNV_TOOL == "ascat":
        with open(checkpoints.ascat.get(**wildcards).output[1], 'r') as f:
            for line in f:
                if line.startswith("rho"):
                    rho = line.strip().split(" ")
                if line.startswith("NormalContamination"):
                    normCont = line.strip().split(" ")
    else:
        with open(checkpoints.battenberg.get(**wildcards).output[1], 'r') as f:
            for line in f:
                if line.startswith("rho"):
                    rho = line.strip().split(" ")
                if line.startswith("NormalContamination"):
                    normCont = line.strip().split(" ")
    
    res = 1-float(rho[1])
    if res > 0 and res < 1:
        return(res)
    res = 1-float(normCont[1])
    if res > 0 and res < 1:
        return(res)
    return(0)


###############################################################################
#### Parallel Block 1
#### CaVEMan setup
#### Genotype Check
#### VerifyBam Normal
###############################################################################
rule caveman_setup:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam"
    output:
        res1_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.CaVEMan_setup"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 1] CaVEMan setup"; \
        TMP={OUTPUT}/{wildcards.tumor}/tmp; \
        mkdir -p {OUTPUT}/{wildcards.tumor}; \
        mkdir -p $TMP; \
    
        cp {res}/caveman/flag.vcf.config.{PROTOCOL}.ini {OUTPUT}/{wildcards.tumor}/flag.vcf.config.{PROTOCOL}.ini;
        perl -alne 'print join(qq{{\\t}},$F[0],0,$F[1],2);' < {ref}/genome.fa.fai | tee $TMP/norm.cn.bed > $TMP/tum.cn.bed;
        
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/genome.fa.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s {SPECIES} \
        -sa {ASSEMBLY} \
        -t {threads} \
        -st {PROTOCOL} \
        -tc {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed \
        -nc {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed \
        -td 5 -nd 2 \
        -tb {input.tumor} \
        -nb {input.normal} \
        -c {OUTPUT}/{wildcards.tumor}/flag.vcf.config.{PROTOCOL}.ini \
        -f {res}/caveman/flagging/flag.to.vcf.convert.ini \
        -e {CAVESPLIT} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/caveman_{CNV_TOOL} \
        {CAVE_CONTIG_EXCLUDE} \
        -p setup' \
        2>&1 | tee {output}
        """
        # 2>&1 | tee {output}
        # >& {output}


rule genotype_check:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam"
    output:
        res1_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.geno"
    shell:
        """
        echo -e "[Parallel block 1] Genotype Check"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=laran \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v compareBamGenotypes.pl \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/genotyped \
        -j {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/genotyped/result.json \
        -nb {input.normal} \
        -tb {input.tumor} \
        -s {res}/qcGenotype/general.tsv \
        -g {res}/qcGenotype/gender.tsv' \
        2>&1 | tee {output}
        """


rule verifybam_normal:
    input:
        normal=INPUT_BAM+"/{normal}.bam"
    output:
        res1_3=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{normal}.time.verify_WT"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 1] VerifyBam Normal"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v verifyBamHomChk.pl -d 25 \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.normal}/contamination \
        -b {input.normal} \
        -t {threads} \
        -j {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.normal}/contamination/result.json \
        -s {res}/qcGenotype/verifyBamID_snps.vcf.gz' \
        2>&1 | tee {output}
        """
        # mv {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}__{wildcards.tumor}.time.verify_WT {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}__{wildcards.normal}.time.verify_WT


###############################################################################
#### Parallel Block 2
#### CaVEMan split
###############################################################################
rule caveman_split:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam",
        res1_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.CaVEMan_setup",
        res1_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.geno",
        res1_3=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{normal}.time.verify_WT"
    output:
        res2_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.CaVEMan_split"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 2] CaVEMan split"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/genome.fa.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s {SPECIES} \
        -sa {ASSEMBLY} \
        -t {threads} \
        -st {PROTOCOL} \
        -tc {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed \
        -nc {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed \
        -td 5 -nd 2 \
        -tb {input.tumor} \
        -nb {input.normal} \
        -c {OUTPUT}/{wildcards.tumor}/flag.vcf.config.{PROTOCOL}.ini \
        -f {res}/caveman/flagging/flag.to.vcf.convert.ini \
        -e {CAVESPLIT} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/caveman_{CNV_TOOL} \
        {CAVE_CONTIG_EXCLUDE} \
        -p split' \
        2>&1 | tee {output}
        """


###############################################################################
#### Parallel Block 3
#### ASCAT
#### BRASS_input
#### BRASS_cover
###############################################################################
checkpoint ascat:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam"
    output:
        res3_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.ascat",
        samplestats=OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/ascat/{tumor}_vs_{normal}.samplestatistics.txt"
    # params:
    #     samplestats=OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/ascat/{tumor}.samplestatistics.txt"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 3] ASCAT"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v ascat.pl \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/ascat \
        -t {input.tumor} \
        -n {input.normal} \
        -sg {res}/ascat/SnpGcCorrections.tsv \
        -r {ref}/genome.fa \
        -q 20 \
        -g L \
        -l {res}/qcGenotype/gender.tsv \
        -rs {SPECIES} \
        -ra {ASSEMBLY} \
        -pr {PROTOCOL} \
        -pl ILLUMINA \
        -c {threads} \
        -force \
        {ASCAT_ADD_ARGS}' \
        2>&1 | tee {output.res3_1};

        cp {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/ascat/{wildcards.tumor}.samplestatistics.txt \
        {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/ascat/{wildcards.tumor}_vs_{wildcards.normal}.samplestatistics.txt
    
        cut -f 2,3,4,5 -d "," \
        {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
        | tr -s "," "\\t" > {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed;
        cut -f 2,3,4,7 -d "," \
        {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
        | tr -s "," "\\t" > {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed

        cp {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed {OUTPUT}/{wildcards.tumor}/tmp/norm_ascat.cn.bed
        cp {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed {OUTPUT}/{wildcards.tumor}/tmp/tum_ascat.cn.bed
        """
        # >& {output.res3_1};
        
        # -g XY \
        # -gc chrY \


rule BRASS_input:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam"
    output:
        res3_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.BRASS_input"
    threads: 2
    shell:
        """
        echo -e "[Parallel block 3] BRASS_input"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v brass.pl -j 4 -k 4 -c {threads} \
        -d {res}/brass/HiDepth.bed.gz \
        -f {res}/brass/brass_np.groups.gz \
        -g {ref}/genome.fa \
        -s {SPECIES} -as {ASSEMBLY} -pr {PROTOCOL} -pl ILLUMINA \
        -g_cache {res}/vagrent/vagrent.cache.gz \
        -vi {res}/brass/viral.genomic.fa.2bit \
        -mi {res}/brass/all_ncbi_bacteria \
        -b {res}/brass/500bp_windows.gc.bed.gz \
        -ct {res}/brass/CentTelo.tsv \
        -cb {res}/brass/cytoband.txt \
        -t {input.tumor} \
        -n {input.normal} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/brass_{CNV_TOOL} \
        -p input'
        2>&1 | tee {output}
        """


rule BRASS_cover:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam"
    output:
        res3_3=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.BRASS_cover"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 3] BRASS_cover"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v nice -n 10 brass.pl -j 4 -k 4 -c {threads} \
        -d {res}/brass/HiDepth.bed.gz \
        -f {res}/brass/brass_np.groups.gz \
        -g {ref}/genome.fa \
        -s {SPECIES} -as {ASSEMBLY} -pr {PROTOCOL} -pl ILLUMINA \
        -g_cache {res}/vagrent/vagrent.cache.gz \
        -vi {res}/brass/viral.genomic.fa.2bit \
        -mi {res}/brass/all_ncbi_bacteria \
        -b {res}/brass/500bp_windows.gc.bed.gz \
        -ct {res}/brass/CentTelo.tsv \
        -cb {res}/brass/cytoband.txt \
        -t {input.tumor} \
        -n {input.normal} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/brass_{CNV_TOOL} \
        -p cover' \
        2>&1 | tee {output}
        """


###############################################################################
#### Parallel Block 4
#### cgpPindel
#### CaVEMan
###############################################################################
rule pindel:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam"
    output:
        res4_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.cgpPindel",
        bed=OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/pindel/{tumor}_vs_{normal}.germline.bed",
    threads: 24
    shell:
        """
        echo -e "[Parallel block 4] cgpPindel"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v pindel.pl \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/pindel \
        -r {ref}/genome.fa \
        -t {input.tumor} \
        -n {input.normal} \
        -s {res}/pindel/simpleRepeats.bed.gz \
        -u {res}/pindel/pindel_np.gff3.gz \
        -f {res}/pindel/{PROTOCOL}_Rules.lst \
        -g {res}/vagrent/codingexon_regions.indel.bed.gz \
        -st {PROTOCOL} \
        -as {ASSEMBLY} \
        -sp {SPECIES} \
        {PINDEL_CONTIG_EXCLUDE} \
        -b {res}/pindel/HiDepth.bed.gz \
        -c {threads} \
        -sf {res}/pindel/softRules.lst' \
        2>&1 | tee {output.res4_1}
        """

rule caveman:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam",
        res2_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.CaVEMan_split",
        res3_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time."+CNV_TOOL,
    output:
        res4_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.caveman"
    params:
        #samplestats=OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/ascat/{tumor}.samplestatistics.txt"
        # NORM_CONTAM=lambda wildcards: get_normalContamination(OUTPUT+"/"+
        #     wildcards.tumor+"/{PROTOCOL}__"+wildcards.tumor+
        #     "/ascat/"+wildcards.tumor+".samplestatistics.txt")
        # NORM_CONTAM=lambda wildcards: get_normalContamination(expand(OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/ascat/{tumor}.samplestatistics.txt", tumor=wildcards.tumor))
        NORM_CONTAM=get_normalContamination
    threads: 24
    shell:
        """
        echo -e "[Parallel block 4] CaVEMan"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/genome.fa.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s {SPECIES} \
        -sa {ASSEMBLY} \
        -t {threads} \
        -st {PROTOCOL} \
        -tc {OUTPUT}/{wildcards.tumor}/tmp/tum_{CNV_TOOL}.cn.bed \
        -nc {OUTPUT}/{wildcards.tumor}/tmp/norm_{CNV_TOOL}.cn.bed \
        -td 5 -nd 2 \
        -tb {input.tumor} \
        -nb {input.normal} \
        -c {OUTPUT}/{wildcards.tumor}/flag.vcf.config.{PROTOCOL}.ini \
        -f {res}/caveman/flagging/flag.to.vcf.convert.ini \
        -e {CAVESPLIT} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/caveman_{CNV_TOOL} \
        {CAVE_CONTIG_EXCLUDE} \
        -k {params.NORM_CONTAM} \
        -no-flagging -noclean' \
        2>&1 | tee {output}
        """
        # NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){{print 1-$1;}}' {params.samplestats}`; \
        # -k ${{NORM_CONTAM}} \

###############################################################################
#### Parallel Block 5
#### BRASS
#### Pindel_annot
#### cgpFlagCaVEMan
###############################################################################
rule brass:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam",
        samplestats=OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/{CNV_TOOL}/{tumor}_vs_{normal}.samplestatistics.txt",
        res3_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.{CNV_TOOL}",
        res3_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.BRASS_input",
        res3_3=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.BRASS_cover"
    output:
        res5_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.BRASS"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 5] BRASS"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v brass.pl -j 4 -k 4 -c {threads} \
        -d {res}/brass/HiDepth.bed.gz \
        -f {res}/brass/brass_np.groups.gz \
        -g {ref}/genome.fa \
        -s {SPECIES} -as {ASSEMBLY} -pr {PROTOCOL} -pl ILLUMINA \
        -g_cache {res}/vagrent/vagrent.cache.gz \
        -vi {res}/brass/viral.genomic.fa.2bit \
        -mi {res}/brass/all_ncbi_bacteria \
        -b {res}/brass/500bp_windows.gc.bed.gz \
        -ct {res}/brass/CentTelo.tsv \
        -cb {res}/brass/cytoband.txt \
        -t {input.tumor} \
        -n {input.normal} \
        -ss {input.samplestats} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/brass_{CNV_TOOL}' \
        2>&1 | tee {output}
        """

rule Pindel_annot:
    input:
        res4_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.cgpPindel"
    output:
        res5_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.cgpPindel_annot"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 5] Pindel_annot"

        # ensure no annotated pindel
        rm -f {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.annot.vcf.gz*

        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v AnnotateVcf.pl -t \
        -c {res}/vagrent/vagrent.cache.gz \
        -i {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.flagged.vcf.gz \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.annot.vcf' \
        2>&1 | tee {output}
        """

rule cgpFlagCaVEMan:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam",
        bed=OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/pindel/{tumor}_vs_{normal}.germline.bed",
        res4_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.cgpPindel",
        res4_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.caveman",
    output:
        res5_3=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.cgpFlagCaVEMan"
    params:
        # samplestats=OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/ascat/{tumor}.samplestatistics.txt",
        # NORM_CONTAM=lambda wildcards: get_normalContamination(OUTPUT+"/"+
        #     wildcards.tumor+"/{PROTOCOL}__"+wildcards.tumor+
        #     "/ascat/"+wildcards.tumor+".samplestatistics.txt"),
        # GERMLINE_BED=lambda wildcards: expand(OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/pindel/{tumor}_vs_{normal}.germline.bed",
        #     tumor=wildcards.tumor, normal=wildcards.normal)
        NORM_CONTAM=get_normalContamination
    threads: 24
    shell:
        """
        echo -e "[Parallel block 5] CaVEMan flag"; \
        
        # need to sort and index pindel germline
        sort -k1,1 -k2,2n -k3,3n {input.bed} | bgzip -c > {input.bed}.gz;
        tabix -p bed {input.bed}.gz
        
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/genome.fa.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s {SPECIES} \
        -sa {ASSEMBLY} \
        -t {threads} \
        -st {PROTOCOL} \
        -tc {OUTPUT}/{wildcards.tumor}/tmp/tum_{CNV_TOOL}.cn.bed \
        -nc {OUTPUT}/{wildcards.tumor}/tmp/norm_{CNV_TOOL}.cn.bed \
        -td 5 -nd 2 \
        -tb {input.tumor} \
        -nb {input.normal} \
        -c {OUTPUT}/{wildcards.tumor}/flag.vcf.config.{PROTOCOL}.ini \
        -f {res}/caveman/flagging/flag.to.vcf.convert.ini \
        -e {CAVESPLIT} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/caveman_{CNV_TOOL} \
        {CAVE_CONTIG_EXCLUDE} \
        -k {params.NORM_CONTAM} \
        -in {input.bed}.gz \
        -p flag' \
        2>&1 | tee {output}
        """
        # NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){{print 1-$1;}}' {params.samplestats}`;
        # rm -f $GERMLINE_BED
        
        # sort -k1,1 -k2,2n -k3,3n {params.GERMLINE_BED} | bgzip -c > {params.GERMLINE_BED}.gz
        # tabix -p bed {params.GERMLINE_BED}.gz
        
        # -in {params.GERMLINE_BED}.gz \

###############################################################################
#### Parallel Block 6
#### CaVEMan_annot
#### VerifyBam tumor
###############################################################################
rule caveman_annot:
    input:
        res5_3=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.cgpFlagCaVEMan"
    output:
        res6_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.CaVEMan_annot"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 6] CaVEMan_annot";

        # annotate caveman
        rm -f {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/caveman_{CNV_TOOL}/{wildcards.tumor}_vs_{wildcards.normal}.annot.muts.vcf.gz*
    
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v AnnotateVcf.pl -t -c {res}/vagrent/vagrent.cache.gz \
        -i {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/caveman_{CNV_TOOL}/{wildcards.tumor}_vs_{wildcards.normal}.flagged.muts.vcf.gz \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/caveman_{CNV_TOOL}/{wildcards.tumor}_vs_{wildcards.normal}.annot.muts.vcf' \
        2>&1 | tee {output}
        """

rule verifyBam_tumor:
    input:
        tumor=INPUT_BAM+"/{tumor}.bam",
        res3_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.{CNV_TOOL}"
    output:
        res6_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.verify_MT"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 6] VerifyBam Tumour"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v verifyBamHomChk.pl -d 25 \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/contamination_{CNV_TOOL} \
        -b {input.tumor} \
        -t {threads} \
        -a {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/{CNV_TOOL}/{wildcards.tumor}.copynumber.caveman.csv \
        -j {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/contamination_{CNV_TOOL}/result.json \
        -s {res}/qcGenotype/verifyBamID_snps.vcf.gz' \
        2>&1 | tee {output}
        """
        # mv {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}__{wildcards.tumor}.time.verify_MT {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}__{wildcards.tumor}.time.verify_MT


###############################################################################
#### Execution of battenberg copynumber (Manual)
###############################################################################
checkpoint battenberg:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam"
    output:
        res3_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.battenberg",
        samplestats=OUTPUT+"/{tumor}/{PROTOCOL}__{tumor}/battenberg/{tumor}_vs_{normal}.samplestatistics.txt"
    threads: 24
    shell:
        """
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi  \
        --volume={WORK}:{WORK}  \
        cgpbattenberg353 \
        bash -c '/usr/bin/time -v battenberg.pl \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/battenberg \
        -r {ref}/genome.fa.fai \
        -u {res}/battenberg/battenberg_1000genomesloci2012_v3 \
        -e {res}/battenberg/battenberg_impute_v3/impute_info.txt \
        -c {res}/battenberg/probloci_270415.txt.gz \
        -ig {res}/battenberg/ignore_contigs.txt \
        -gc {res}/battenberg/battenberg_wgs_gc_correction_1000g_v3 \
        -ge L \
        -gl {res}/qcGenotype/gender.tsv \
        -tb {input.tumor} \
        -nb {input.normal} \
        -rs {SPECIES} \
        -ra {ASSEMBLY} \
        -nl 50 \
        -t {threads}' \
        2>&1 | tee {output.res3_1}
        
        cat {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/battenberg/{wildcards.tumor}_{{normal_contamination,cellularity_ploidy,rho_and_psi}}.txt | \
        awk -F'\t' '{{ \
        if (NR==1) printf("NormalContamination %s\\n",$1); \
        if (NR==3) printf("Ploidy %s\\n",$2); \
        if (NR==6) printf("rho %s\\npsi %s\\ngoodnessOfFit %s\\nGenderChr Y\\nGenderChrFound Y\\n",$2,$3,$5) \
        }}' > {output.samplestats}
        
        zcat {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/battenberg/{wildcards.tumor}_battenberg_cn.vcf.gz | \
        grep -v ^# | cut -f 1,2,8,10 | \
        awk '{{for(i=1;i<=NF;i++){{if (i==3) {{split($i,a,"=");printf("%s\\t",a[3])}} else if (i==4) {{split($i,a,":");printf("%s\\n", a[2])}} else {{printf("%s\\t",$i)}} }}}}' > \
        {OUTPUT}/{wildcards.tumor}/tmp/norm_battenberg.cn.bed
        
        zcat {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/battenberg/{wildcards.tumor}_battenberg_cn.vcf.gz | \
        grep -v ^# | cut -f 1,2,8,11 | \
        awk '{{for(i=1;i<=NF;i++){{if (i==3) {{split($i,a,"=");printf("%s\\t",a[3])}} else if (i==4) {{split($i,a,":");printf("%s\\n", a[2])}} else {{printf("%s\\t",$i)}} }}}}' > \
        {OUTPUT}/{wildcards.tumor}/tmp/tum_battenberg.cn.bed
    
        zcat {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/battenberg/{wildcards.tumor}_battenberg_cn.vcf.gz | \
        grep -v ^# | cut -f 1,2,8,10,11 | \
        awk '{{printf("%s,",NR);for(i=1;i<=NF;i++){{if (i==3) {{split($i,a,"=");printf("%s,",a[3])}} else if (i==4) {{split($i,a,":");printf("%s,%s,", a[2],a[3])}} else if (i==5) {{split($i,a,":");printf("%s,%s\\n", a[2],a[3])}} else {{printf("%s,",$i)}} }}}}' > \
        {OUTPUT}/{wildcards.tumor}/{PROTOCOL}__{wildcards.tumor}/battenberg/{wildcards.tumor}.copynumber.caveman.csv
        """
        # -ge L \
        # -gl {res}/qcGenotype/gender.tsv \


###############################################################################
#### Tying it all together with max possible paralization
###############################################################################
rule run_parallel:
    input:
        res1_3=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{normal}.time.verify_WT",
        res5_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.BRASS",
        res5_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.cgpPindel_annot",
        res6_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.CaVEMan_annot",
        res6_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{CNV_TOOL}_{tumor}_vs_{normal}.time.verify_MT",
        battenberg=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.battenberg"
    output:
        OUTPUT+"/finished/SUCCESS__{PROTOCOL}__{CNV_TOOL}__{tumor}__vs__{normal}.txt"
    shell:
        """
        echo 'Successfully completed entire cgp pipeline!' > {output}
        """

rule run_test:
    input:
        res1_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.CaVEMan_setup",
        # res1_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.geno",
        # res1_3=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{normal}.time.verify_WT",
        # res5_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.BRASS",
        # res5_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.cgpPindel_annot",
        # res6_1=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.CaVEMan_annot",
        # res6_2=OUTPUT+"/{tumor}/timings/{PROTOCOL}__{tumor}_vs_{normal}.time.verify_MT",
    output:
        OUTPUT+"/finished/TEST_SUCCESS_{tumor}_vs_{normal}.txt"
    shell:
        """
        echo 'Successfully completed entire cgp pipeline!' > {output}
        """
    