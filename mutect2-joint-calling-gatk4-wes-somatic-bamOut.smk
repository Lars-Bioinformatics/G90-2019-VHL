__title__ = "Pipeline for Somatic Joint Variant Calling with Mutect2 - VHL pipeline"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "14/05/2020"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "/work/sduvarcall/G90-2019-VHL-case-projekt/samples_somatic.yaml"

# Explicit paths for external input files
ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
gnomad = "/work/sduvarcall/knownSNPs/gnomad/af-only-gnomad.raw.sites.b37.vcf.gz"
interval = "/work/sduvarcall/Resources/GRCh37/MedExome_target_regions/MedExome_GRCh37_capture_targets.interval_list"
common_variants = "/work/sduvarcall/knownSNPs/gnomad/somatic-b37_small_exac_common_3.vcf"

# split interval files
target_regions = "/work/sduvarcall/Resources/GRCh37/MedExome_target_regions/split_regions/"
INTERVALS, = glob_wildcards(target_regions+"{interval}.interval_list")
INTERVALS = sorted(INTERVALS)
print(INTERVALS)

# list of normals for Panel of normals (PON)
NORMALS = list(set([config[sample]["normal"] for sample in config]))
print(NORMALS)

# Get sample names from config file
# SAMPLES = [sample for sample in config]
SAMPLES = ["G90-VHL521-10", "G90-VHL526-10"]
print(SAMPLES)


#########################################################
####                      Output                     ####
#########################################################
# log_file = "log_gatk4-wes-somatic.txt"
output_somatic = "mutect2_somatic_joint_calling_bamOut/"
# Panel of normals location
pon_location = output_somatic+"panel_of_normals/"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
# totim = time.time()
# timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
# mem = "-Xmx12g" # login nodes
# mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx50g"
mem = "-Xmx12g"

onstart:
    shell("mkdir -p "+output_somatic)
    shell("mkdir -p "+output_somatic+"split/")
    shell("mkdir -p "+output_somatic+"split_f1r2/")
    shell("mkdir -p "+output_somatic+"contamination/")
    shell("mkdir -p "+output_somatic+"bam_split/")
    shell("mkdir -p "+output_somatic+"realigned_bam/")
    shell("mkdir -p "+pon_location)
    shell("mkdir -p "+pon_location+"split")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
rule all:
    input:
        # pon_location+"pon.vcf.gz"
        expand(output_somatic+"{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=SAMPLES),


#########################################################
####       Create Somatic Panel of Normals           ####
#########################################################
# '''
# Run Mutect2 in tumor-only mode for each normal sample
# '''
# rule Mutect2_tumor_only_pon:
#     input:
#         bam="bam/{sample}_recal.bam",
#         intervals=target_regions+"{interval}.interval_list"
#     output:
#         vcf=pon_location+"split/{sample}_for_pon__{interval}__split.vcf.gz",
#         vcf_stats=pon_location+"split/{sample}_for_pon__{interval}__split.vcf.gz.stats",
#         # idx=pon_location+"split/{sample}_for_pon__{interval}__split.vcf.idx"
#     threads: 24
#     shell:
#         """
#         gatk --java-options {mem} Mutect2 \
#         -R {ref} \
#         -I {input.bam} \
#         -max-mnp-distance 0 \
#         --native-pair-hmm-threads {threads} \
#         -L {input.intervals} \
#         -O {output.vcf}
#         """
#
# rule merge_normal_vcf:
#     input:
#         vcf_subfile=expand(pon_location+"split/{{sample}}_for_pon__{interval}__split.vcf.gz", interval=INTERVALS)
#     output:
#         vcf=pon_location+"{sample}_for_pon.vcf.gz",
#         tbi=pon_location+"{sample}_for_pon.vcf.gz.tbi"
#     params:
#         vcf_subfile=expand("-I "+pon_location+"split/{{sample}}_for_pon__{interval}__split.vcf.gz", interval=INTERVALS)
#     shell:
#         """
#         gatk --java-options {mem} GatherVcfs \
#         {params.vcf_subfile} \
#         -O {output.vcf}
#
#         tabix -p vcf {output.vcf}
#         """
#
# rule merge_normal_vcf_stats:
#     input:
#         vcf_subfile=expand(pon_location+"split/{{sample}}_for_pon__{interval}__split.vcf.gz.stats", interval=INTERVALS)
#     output:
#         vcf_stats=pon_location+"{sample}_for_pon.vcf.gz.stats",
#     params:
#         vcf_subfile=expand("-stats "+pon_location+"split/{{sample}}_for_pon__{interval}__split.vcf.gz.stats", interval=INTERVALS)
#     shell:
#         """
#         gatk --java-options {mem} MergeMutectStats \
#         {params.vcf_subfile} \
#         -O {output.vcf_stats}
#         """
#
# rule GenomicsDB:
#     input:
#         vcf=expand(pon_location+"{sample}_for_pon.vcf.gz", sample=NORMALS),
#         vcf_stats=expand(pon_location+"{sample}_for_pon.vcf.gz.stats", sample=NORMALS)
#     output:
#         db=pon_location+"pon_db/callset.json"
#     params:
#         vcf=expand("-V "+pon_location+"{sample}_for_pon.vcf.gz", sample=NORMALS),
#         db=pon_location+"pon_db"
#     shell:
#         """
#         export LC_ALL="en_US.UTF-8" # Make sure 1000-sep is '.'
#
#         gatk --java-options {mem} GenomicsDBImport \
#         -R {ref} \
#         -L {interval} \
#         --genomicsdb-workspace-path {params.db} \
#         --max-num-intervals-to-import-in-parallel {threads} \
#         --merge-input-intervals true \
#         {params.vcf}
#         """
#
# '''
# Combine the normal calls using CreateSomaticPanelOfNormals.
# '''
# rule CreateSomaticPanelOfNormals:
#     input:
#         pon_location+"pon_db/callset.json"
#         # vcf=expand(pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
#     output:
#         pon=pon_location+"pon.vcf.gz"
#     params:
#         db=pon_location+"pon_db"
#         # vcf=expand("-V "+pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
#     shell:
#         """
#         gatk --java-options {mem} CreateSomaticPanelOfNormals \
#         -R {ref} \
#         -V gendb://{params.db} \
#         -O {output}
#         """
#         # {params.vcf} \


##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_matched:
    input:
        normal=lambda wildcards: expand("bam/{normal}_recal.bam", normal=config[wildcards.sample]["normal"]),
        tumors=lambda wildcards: expand("bam/{tumor}_recal.bam", tumor=config[wildcards.sample]["tumors"]),
        pon=pon_location+"pon.vcf.gz",
        intervals=target_regions+"{interval}.interval_list"
    output:
        vcf=output_somatic+"split/{sample}_somatic_mutect2__{interval}__split.vcf.gz",
        tbi=output_somatic+"split/{sample}_somatic_mutect2__{interval}__split.vcf.gz.tbi",
        vcf_stats=output_somatic+"split/{sample}_somatic_mutect2__{interval}__split.vcf.gz.stats",
        f1r2=output_somatic+"split_f1r2/{sample}_f1r2__{interval}__split.tar.gz",
        bam_subfile=output_somatic+"bam_split/{sample}__{interval}_mutect2.bam",
    params:
        tumors=lambda wildcards: expand("-I bam/{tumor}_recal.bam", tumor=config[wildcards.sample]["tumors"]),
        normal_name=lambda wildcards: expand("-normal {normal}", normal=config[wildcards.sample]["normal"].split("_",1)[0])
    threads: 24
    shell:
        """
        gatk --java-options {mem} Mutect2 \
        -R {ref} \
        {params.tumors} \
        -I {input.normal} \
        {params.normal_name} \
        -pon {input.pon} \
        --germline-resource {gnomad} \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --native-pair-hmm-threads {threads} \
        --f1r2-tar-gz {output.f1r2} \
        -bamout {output.bam_subfile} \
        -L {input.intervals} \
        -O {output.vcf}
        """

rule merge_mutect2_bam:
  input:
        bam_subfile=expand("bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
  output:
        bam=output_somatic+"realigned_bam/{sample}.bam"
  params:
        bam_subfile=expand("-I bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
  shell:
        """
        gatk --java-options {mem} MergeSamFiles \
        {params.bam_subfile} \
        --VALIDATION_STRINGENCY LENIENT \
        --USE_THREADING true \
        -O {output.bam} \
        --CREATE_INDEX true
        """

rule merge_somatic_vcf:
    input:
        vcf_subfile=expand(output_somatic+"split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
    output:
        vcf=output_somatic+"{sample}_somatic_mutect2.vcf.gz",
        idx=output_somatic+"{sample}_somatic_mutect2.vcf.gz.tbi",
    params:
        vcf_subfile=expand("-I "+output_somatic+"split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
    shell:
        """
        gatk --java-options {mem} GatherVcfs \
        {params.vcf_subfile} \
        -O {output.vcf}
    
        tabix -p vcf {output.vcf}
        """

rule merge_somatic_vcf_stats:
    input:
        vcf_subfile=expand(output_somatic+"split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
    output:
        vcf_stats=output_somatic+"{sample}_somatic_mutect2.vcf.gz.stats",
    params:
        vcf_subfile=expand("-stats "+output_somatic+"split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
    shell:
        """
        gatk --java-options {mem} MergeMutectStats \
        {params.vcf_subfile} \
        -O {output.vcf_stats}
        """

#########################################################
####          Learn Read Orientation Bias            ####
####                                                 ####
#### Note: Used to fix orientation bias artifacts    ####
####       from Formalin-Fixed Paraffin-Embedded     ####
####       (FFPE) samples - i.e. not needed for      ####
####       frozen tissue                             ####
#########################################################
rule learnReadOrientationModel:
    input:
        f1r2=expand(output_somatic+"split_f1r2/{{sample}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    output:
        f1r2_model=output_somatic+"{sample}_read-orientation-model.tar.gz"
    params:
        f1r2=expand("-I "+output_somatic+"split_f1r2/{{sample}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    shell:
        """
        gatk LearnReadOrientationModel \
        {params.f1r2} \
        -O {output}
        """

#########################################################
####           Create Contamination table            ####
#########################################################
rule GetPileupSummaries:
    input:
        bam="bam/{sample}_recal.bam"
    output:
        pileup=output_somatic+"contamination/{sample}_pileup.table"
    shell:
        """
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

rule CalculateContamination:
    input:
        normal=output_somatic+"contamination/{normal}_pileup.table",
        tumor=output_somatic+"contamination/{tumor}_pileup.table"
    output:
        contamination=output_somatic+"contamination/{tumor}_vs_{normal}_contamination.table"
    shell:
        """
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """

rule JoinContaminationTables:
    input:
        cont_tables=lambda wildcards: expand(output_somatic+"contamination/{tumor}_vs_{normal}_contamination.table", 
                    normal=config[wildcards.sample]["normal"], tumor=config[wildcards.sample]["tumors"])
    output:
        joined_table=output_somatic+"{sample}_merged_contamination.table"
    shell:
        """
        echo -e 'sample\tcontamination\terror' > {output};
        awk -F'\t' 'FNR == 2' {input} >> {output}
        """


#########################################################
####              Filter Mutect2 Calls               ####
#########################################################
rule FilterMutectCalls:
    input:
        vcf=output_somatic+"{sample}_somatic_mutect2.vcf.gz",
        idx=output_somatic+"{sample}_somatic_mutect2.vcf.gz.tbi",
        bam=output_somatic+"realigned_bam/{sample}.bam",
        stats=output_somatic+"{sample}_somatic_mutect2.vcf.gz.stats",
        contamination=output_somatic+"{sample}_merged_contamination.table",
        read_orientation=output_somatic+"{sample}_read-orientation-model.tar.gz"
    output:
        vcf=output_somatic+"{sample}_somatic_mutect2_filterFlag.vcf.gz",
        # tbi=output_somatic+"{sample}_somatic_mutect2_filterFlag.vcf.gz.tbi",
        vcf_pass=output_somatic+"{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
        # tbi_pass=output_somatic+"{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz.tbi"
    shell:
        """
        gatk --java-options {mem} FilterMutectCalls \
        -R {ref} \
        -V {input.vcf} \
        --contamination-table {input.contamination} \
        --orientation-bias-artifact-priors {input.read_orientation} \
        --stats {input.stats} \
        -L {interval} \
        -O {output.vcf}
        
        zcat {output.vcf} | egrep '^#|PASS' | bgzip -c > {output.vcf_pass}
        tabix -p vcf {output.vcf_pass}
        """
        # tabix -p vcf {output.vcf}
