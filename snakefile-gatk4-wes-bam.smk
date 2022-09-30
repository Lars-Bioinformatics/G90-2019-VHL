__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "14/05/2020"
__version__ = "1.0"

import time, os

SAMPLES, = glob_wildcards("fastq/{sample}_R1.fastq.gz")
# print(SAMPLES)
# FAMNAME = os.getcwd().rsplit("/",1)[1]
# FAMNAME = [fam for fam in config]
# print(FAMNAME)

totim = time.time()
timeFormat = time.strftime("%Y%m%d-%H%M%S") # - python3
log_file = "log_gatk4-wes-bam_" + timeFormat + ".txt"

# mem = "-Xmx12g" # login nodes
# mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
mem = "-Xmx50g"

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
    shell("mkdir -p Metrics; mkdir -p bam")

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
        # BAM file
        expand("bam/{sample}_recal.bam", sample=SAMPLES),
        # Metrics and Statistics
        expand("Metrics/{sample}_post_recalibration.grp", sample=SAMPLES),
        expand("Metrics/{sample}.HS_Metrics.txt", sample=SAMPLES),
        expand("Metrics/{sample}_unsorted.quality_distribution_metrics", sample=SAMPLES),
        expand("Metrics/{sample}.alignment_summary_metrics", sample=SAMPLES),
        expand("Metrics/{sample}_aggregation.alignment_summary_metrics", sample=SAMPLES)
        # Resources
        # "SeqCap_EZ_Exome_v3_GRCh37_capture_targets.interval_list"



###############################################################################
### Create recalibrated bam file                                            ###
###############################################################################
'''
Map reads to reference genome with bwa
'''
rule MapAndSort:
    input:
        f1 = "fastq/{sampleid}_{protocol}_{flowcell}_R1.fastq.gz",
        f2 = "fastq/{sampleid}_{protocol}_{flowcell}_R2.fastq.gz"
    output:
        bam = temp("{sampleid}_{protocol}_{flowcell}_sorted.bam"),
        bai = temp("{sampleid}_{protocol}_{flowcell}_sorted.bai")
    params:
        rgid = "{sampleid}_{protocol}_{flowcell}",
        rglb = "{protocol}",
        rgsm = "{sampleid}",
        rgpl = "Illumina",
        rgpu = "{flowcell}",
    threads: 24
    shell:
        """
        bwa mem -M -t {threads} \
        -R \"@RG\\tID:{params.rgid}\\tLB:{params.rglb}\\tSM:{params.rgsm}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\" \
        {ref} {input} | \
        gatk --java-options {mem} SortSam \
        --INPUT=/dev/stdin \
        --OUTPUT={output.bam} \
        --VALIDATION_STRINGENCY=LENIENT \
        --SORT_ORDER=coordinate \
        --CREATE_INDEX=TRUE
        """


'''
Remove duplicate reads
'''
rule MarkDuplicates:
    input:
        bam = "{sample}_sorted.bam",
        bai = "{sample}_sorted.bai"
    output:
        bam = temp("{sample}_dedup.bam"),
        bai = temp("{sample}_dedup.bai"),
        met = "Metrics/{sample}_duplicate_metrics.txt"
    shell:
        """
        gatk --java-options {mem} MarkDuplicates \
        --INPUT={input.bam} \
        --OUTPUT={output.bam} \
        --METRICS_FILE={output.met} \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        --CREATE_INDEX=true \
        """


'''
Obtain recalibration information
'''
rule BaseRecalibrator:
    input:
        bam = "{sample}_dedup.bam",
        bai = "{sample}_dedup.bai",
    output:
        table = "Metrics/{sample}_pre_recalibration.grp"
    shell:
        """
        gatk --java-options {mem} BaseRecalibrator \
        -R={ref} \
        -I={input.bam} \
        --known-sites={dbsnp} \
        --known-sites={mills_1000G} \
        --known-sites={phase1_1000G} \
        -L={interval} \
        -O={output} \
        """


'''
Apply Recalibration
'''
rule ApplyRecalibration:
    input:
        bam="{sample}_dedup.bam",
        bai="{sample}_dedup.bai",
        table = "Metrics/{sample}_pre_recalibration.grp"
    output:
        bam="bam/{sample}_recal.bam",
        bai="bam/{sample}_recal.bai"
    shell:
        """
        gatk --java-options {mem} ApplyBQSR \
        -R={ref} \
        -I={input.bam} \
        --bqsr-recal-file={input.table} \
        -L={interval} \
        -O={output.bam} \
        """

###############################################################################
### Collect Statistics                                                      ###
###############################################################################
'''
Collect recalibration metrics
'''
rule RecalibrationMetrics:
    input:
        bam = "bam/{sample}_recal.bam",
    output:
        table = "Metrics/{sample}_post_recalibration.grp"
    shell:
        """
        gatk --java-options {mem} BaseRecalibrator \
        -R={ref} \
        -I={input.bam} \
        --known-sites={dbsnp} \
        --known-sites={mills_1000G} \
        --known-sites={phase1_1000G} \
        -L={interval} \
        -O={output} \
        """

'''
Collect Hybrid Selection (HS) metrics
'''
rule collectHsMetrics:
    input:
        bam="bam/{sample}_recal.bam"
    output:
        "Metrics/{sample}.HS_Metrics.txt"
    shell:
        """
        gatk --java-options {mem} CollectHsMetrics \
        --INPUT={input.bam} \
        --REFERENCE_SEQUENCE={ref} \
        --OUTPUT={output} \
        --BAIT_INTERVALS={interval} \
        --TARGET_INTERVALS={interval} \
        """



rule CollectUnsortedReadgroupBamQualityMetrics:
    input:
        bam="bam/{sample}_recal.bam"
    output:
        "Metrics/{sample}_unsorted.quality_distribution_metrics"
        # And many more, not written here for simplicity
    params:
        out="Metrics/{sample}_unsorted"
    shell:
        """
        gatk --java-options {mem} CollectMultipleMetrics \
        --INPUT={input.bam} \
        --OUTPUT={params} \
        --ASSUME_SORTED=true \
        --PROGRAM="null" \
        --PROGRAM="CollectBaseDistributionByCycle" \
        --PROGRAM="CollectInsertSizeMetrics" \
        --PROGRAM="MeanQualityByCycle" \
        --PROGRAM="QualityScoreDistribution" \
        --METRIC_ACCUMULATION_LEVEL="null" \
        --METRIC_ACCUMULATION_LEVEL="ALL_READS" \
        """


rule CollectReadgroupBamQualityMetrics:
    input:
        bam="bam/{sample}_recal.bam"
    output:
        "Metrics/{sample}.alignment_summary_metrics"
        # And many more, not written here for simplicity
    params:
        out="Metrics/{sample}"
    shell:
        """
        gatk --java-options {mem} CollectMultipleMetrics \
        --INPUT={input.bam} \
        --REFERENCE_SEQUENCE={ref} \
        --OUTPUT={params} \
        --ASSUME_SORTED=true \
        --PROGRAM="null" \
        --PROGRAM="CollectAlignmentSummaryMetrics" \
        --PROGRAM="CollectGcBiasMetrics" \
        --METRIC_ACCUMULATION_LEVEL="null" \
        --METRIC_ACCUMULATION_LEVEL="READ_GROUP" \
        """
    
rule CollectAggregationMetrics:
    input:
        bam="bam/{sample}_recal.bam"
    output:
        "Metrics/{sample}_aggregation.alignment_summary_metrics"
        # And many more, not written here for simplicity
    params:
        out="Metrics/{sample}_aggregation"
    shell:
        """
        gatk --java-options {mem} CollectMultipleMetrics \
        --INPUT={input.bam} \
        --REFERENCE_SEQUENCE={ref} \
        --OUTPUT={params} \
        --ASSUME_SORTED=true \
        --PROGRAM="null" \
        --PROGRAM="CollectAlignmentSummaryMetrics" \
        --PROGRAM="CollectInsertSizeMetrics" \
        --PROGRAM="CollectGcBiasMetrics" \
        --PROGRAM="CollectSequencingArtifactMetrics" \
        --PROGRAM="QualityScoreDistribution" \
        --METRIC_ACCUMULATION_LEVEL="null" \
        --METRIC_ACCUMULATION_LEVEL="SAMPLE" \
        --METRIC_ACCUMULATION_LEVEL="LIBRARY"
        """


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