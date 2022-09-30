WORK="/work/sduvarcall/G90-2019-VHL-case-projekt/bam/"

SAMPLES, = glob_wildcards(WORK+"{sample}.bam")

# Resources - paths inside docker
ref_fai = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta.fai"

# onstart:
#     shell("mkdir -p " + OUTPUT)

rule all:
    input:
        expand(WORK+"{sample}.bam.bas", sample=SAMPLES)


rule bam_stats:
    input:
        bam=WORK+"{sample}.bam"
    output:
        bas=WORK+"{sample}.bam.bas"
    shell:
        """
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume=/work:/work \
        cgpwgs \
        bam_stats \
        -i {input} \
        -o {output} \
        -r {ref_fai} \
        -@ 12
        """