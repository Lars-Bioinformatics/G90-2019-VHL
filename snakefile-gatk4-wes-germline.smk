__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "14/05/2020"
__version__ = "1.0"

import time, os, sys

configfile: "/work/sduvarcall/G90-2019-VHL-case-projekt/samples_germline.yaml"

# SAMPLES, = glob_wildcards("fastq/{sample}_R1.fastq.gz")
# SAMPLES = [sample.split("_")[0].rsplit("-",1)[0] for sample in SAMPLES]
# SAMPLES = ["G73-05B5_nimblegen-medexome_HJ5NKBGX9"]
# FAMNAME = os.getcwd().rsplit("/",1)[1]
# FAMNAME = "fam18"
# FAMNAME = [fam for fam in config]
# print(FAMNAME)

totim = time.time()
timeFormat = time.strftime("%Y%m%d-%H%M%S") # - python3
log_file = "log_gatk4-wes-germline_" + timeFormat + ".txt"

mem = "-Xmx12g" # login nodes
# mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx50g"

# Resources
ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
bed = "/work/sduvarcall/Resources/target_regions/MedExome_target_regions/MedExome_GRCh37_capture_targets.bed"
interval = "/work/sduvarcall/Resources/target_regions/MedExome_target_regions/MedExome_GRCh37_capture_targets.interval_list"
dbsnp = "/work/sduvarcall/knownSNPs/dbsnp_150.b37.vcf.gz"
mills_1000G = "/work/sduvarcall/knownSNPs/Mills_and_1000G_gold_standard.indels.b37.vcf"
cosmic = "/work/sduvarcall/cosmic/Cosmic-combined_v81_b37.vcf"
hapmap = "/work/sduvarcall/knownSNPs/hapmap_3.3.b37.vcf"
omni = "/work/sduvarcall/knownSNPs/1000G_omni2.5.b37.vcf"
phase1_1000G = "/work/sduvarcall/knownSNPs/1000G_phase1.indels.b37.vcf"


onstart:
    shell("echo 'Started execution of pipeline:' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
    shell("mkdir -p gvcf_files")

onsuccess:
    fiTime = 'Total time:', str((time.time()-totim) / 60), 'minutes'
    shell("echo 'Finished execution on' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
    shell("echo {fiTime} >> {log_file}")

onerror:
    fiTime = 'Total time:', str((time.time()-totim) / 60), 'minutes'
    shell("echo 'Finished execution on' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
    shell("echo {fiTime} >> {log_file}")
    shell("echo 'ERROR OCCURED, PLEASE REFER TO SLURM LOGFILE FOR DETAILS' >> {log_file}")


###############################################################################
### Rule all                                                                ###
###############################################################################
rule all:
    input:
        [expand("haplotypeCaller_germline/{fam_name}_variants.vcf.gz", fam_name=fam) for fam in config]


###############################################################################
### Call germline variants                                                  ###
###############################################################################
'''
HaplotypeCaller
'''
rule HaplotypeCaller:
    input:
        bam="bam/{sample}_recal.bam"
    output:
        gvcf="gvcf_files/{sample}_variants.g.vcf.gz"
    shell:
        """
        gatk --java-options {mem} HaplotypeCaller \
        -R={ref} \
        -I={input.bam} \
        -O={output.gvcf} \
        -ERC=GVCF \
        --dbsnp={dbsnp} \
        -L={interval} \
        -A Coverage \
        -A DepthPerAlleleBySample \
        -A BaseQualityRankSumTest
        """
        # Coverage - Total depth of coverage per sample and over all samples
        # DepthPerAlleleBySample - Unfiltered alternative allele depth (AD)
        # BaseQualityRankSumTest - Rank Sum Test of REF versus ALT base quality scores

'''
Create single multi-sample g.vcf file    
'''
rule CombineGVCFs:
    input:
        gvcfs=lambda wildcards: expand("gvcf_files/{sample}_variants.g.vcf.gz", sample=config[wildcards.fam_name])
    output:
        gvcf="gvcf_files/{fam_name}_variants_combined.g.vcf.gz"
    params:
        gvcfs=lambda wildcards: expand("-V=gvcf_files/{sample}_variants.g.vcf.gz", sample=config[wildcards.fam_name]),
    shell:
        """
        gatk --java-options {mem} CombineGVCFs \
        -R {ref} \
        {params.gvcfs} \
        -O {output.gvcf}
        """

'''
GenotypeGVCFs
'''
rule GenotypeGVCFs:
    input:
        gvcf="gvcf_files/{fam_name}_variants_combined.g.vcf.gz"
    output:
        vcf="haplotypeCaller_germline/{fam_name}_variants.vcf.gz"
    shell:
        """
        gatk --java-options {mem} GenotypeGVCFs \
        -R={ref} \
        -V={input.gvcf} \
        --dbsnp={dbsnp} \
        -L={interval} \
        -O={output} \
        """


# '''
# Recalibrate Variant Quality Score
# '''
# rule VariantRecalibrator:
#     input:
#         vcf="{fam_name}_variants.vcf.gz"
#     output:
#         recal="{fam_name}_variants.recal",
#         tranches="{fam_name}_variants.tranches",
#         rscript="{fam_name}_variants.plots.R"
#     shell:
#         """
#         gatk --java-options {mem} VariantRecalibrator \
#         -R {ref} \
#         -V {input} \
#         -O {output.recal} \
#         --resource hapmap,known=false,training=true,truth=true,prior=15.0:{hapmap} \
#         --resource omni,known=false,training=true,truth=false,prior=12.0:{omni} \
#         --resource 1000G,known=false,training=true,truth=false,prior=10.0:{phase1_1000G} \
#         --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{dbsnp} \
#         -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
#         -mode BOTH \
#         --tranches-file {output.tranches} \
#         --rscript-file {output.rscript}
#         """
#
#
# '''
# Apply Variant Recalibration
# '''
# rule ApplyVQSR:
#     input:
#         vcf="{fam_name}_variants.vcf.gz",
#         recal="{fam_name}_variants.recal",
#         tranches="{fam_name}_variants.tranches"
#     output:
#         vcf="{fam_name}_variants_recal_gatk.vcf.gz"
#     shell:
#         """
#         gatk --java-options {mem} ApplyVQSR \
#         -R {ref} \
#         -V {input.vcf} \
#         -O {output} \
#         --tranches-file {input.tranches} \
#         --recal-file {input.recal} \
#         --mode BOTH
#         """

###############################################################################
### Create resouces
###############################################################################
# rule bedToInterval:
#     input:
#         bed="SeqCap_EZ_Exome_v3_GRCh37_capture_targets.bed"
#     output:
#         interval="SeqCap_EZ_Exome_v3_GRCh37_capture_targets.interval_list"
#     shell:
#         """
#         gatk BedToIntervalList \
#         -I={input.bed} \
#         -O={output} \
#         -SD={ref}
#         """