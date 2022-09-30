from snakemake.utils import R

WORK="/work/sduvarcall/G90-2019-VHL-case-projekt/"
INPUT_BAM = WORK+"bam_cgpwes"
OUTPUT = WORK+"sequenza_cgpwes"

# configfile: WORK+"samples_matched_pairwise.yaml"
configfile: WORK+"samples_matched_pairwise_cgpwes.yaml"

ref_build = "GRCh37"

# Resources - paths inside docker
if ref_build == "GRCh37":
    # GRCh37 paths
    # ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
    # wig = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.gc50Base.wig.gz"
    ref = "/work/sduvarcall/G37-2016-BRCAX-Genomes/reference_files/core_ref_GRCh37d5/genome.fa"
    wig = "/work/sduvarcall/G37-2016-BRCAX-Genomes/reference_files/core_ref_GRCh37d5/genome.gc50Base.wig.gz"
else:
    # GRCh38 paths
    ref = REF_FILES+"/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"


# Get sample names from config file
SAMPLES = [sample for sample in config]
# SAMPLES = ["G56-sampleA1_nimblegen-medexome_HGL2LDSXX"]
print(SAMPLES)

window = "50"

onstart:
    shell("mkdir -p "+OUTPUT)

rule all:
    input:
        [expand(OUTPUT+"/Sequenza_{tumor}_vs_{normal}_bin"+window+"/Sequenza_{tumor}_vs_{normal}_chromosome_view.pdf", 
            tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in SAMPLES],



###############################################################################
#### Sequenza
###############################################################################
rule gc_wiggle:
    input:
        ref = ref
    output:
        wig = wig
    shell:
        """
        sequenza-utils gc_wiggle --fasta {input.ref} -w 50 -o {output.wig}
        """

rule bam2seqz:
    input:
        # normal=INPUT_BAM+"/{normal}_recal.bam",
        # tumor=INPUT_BAM+"/{tumor}_recal.bam",
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam",
        wig = wig
    output:
        seqz=OUTPUT+"/{tumor}_vs_{normal}.seqz.gz"
    shell:
        """
        sequenza-utils bam2seqz -n {input.normal} -t {input.tumor} \
        -q 20 -N 20 -f "illumina" --fasta {ref} -gc {wig} -o {output}
        """

rule seqz_binning:
    input:
        seqz=OUTPUT+"/{tumor}_vs_{normal}.seqz.gz"
    output:
        seqz_bin=OUTPUT+"/{tumor}_vs_{normal}.seqz.bin{window}.gz"
    shell:
        """
        sequenza-utils seqz_binning --seqz {input.seqz} -w {window} -o {output.seqz_bin}
        """

rule sequenza_analysis:
    input:
        seqz_bin=OUTPUT+"/{tumor}_vs_{normal}.seqz.bin{window}.gz"
    output:
        res=OUTPUT+"/Sequenza_{tumor}_vs_{normal}_bin{window}/Sequenza_{tumor}_vs_{normal}_chromosome_view.pdf"
    run:
        R("""
        library(sequenza)
        seqzdata = sequenza.extract("{input.seqz_bin}", chromosome.list = c(1:22,"X","Y"))

        CP = sequenza.fit(seqzdata)
        
        id = "Sequenza_{wildcards.tumor}_vs_{wildcards.normal}"
        out.dir = paste0("{OUTPUT}/", id, "_bin{wildcards.window}")
        sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = id, out.dir=out.dir)        
        """)
