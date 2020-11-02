rule picard_reheading:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "reheaded_reads/{sample}.bam"
    log:
        "logs/picard/replace_rg/{sample}.log"
    params:
        "RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20"
    wrapper:
        "0.67.0/bio/picard/addorreplacereadgroups"

rule picard_mark_duplicates:
    input:
        "reheaded_reads/{sample}.bam"
    output:
        bam="dedup/{sample}.bam",
        metrics="dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.66.0/bio/picard/markduplicates"

rule create_dict:
    input:
        "data/genome.fa"
    output:
        "data/genome.dict"
    log:
        "logs/picard/create_dict.log"
    params:
        extra=""  # optional: extra arguments for picard.
    wrapper:
        "0.66.0/bio/picard/createsequencedictionary"

rule gatk_baserecalibrator:
    input:
        bam="dedup/{sample}.bam",
        ref="data/genome.fa",
        dic="data/genome.dict",
        known="data/db/dummysnp.vcf"  # optional known sites
    output:
        recal_table="recal/{sample}.grp"
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    params:
        extra="",  # optional
        java_opts="", # optional
    wrapper:
        "0.66.0/bio/gatk/baserecalibrator" 

rule gatk_applybqsr:
    input:
        bam="dedup/{sample}.bam",
        ref="data/genome.fa",
        dictio="data/genome.dict",
        recal_table="recal/{sample}.grp"
    output:
        bam="recal/{sample}.bam"
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log"
    params:
        extra="",  # optional
        java_opts="", # optional
    wrapper:
        "0.66.0/bio/gatk/applybqsr"