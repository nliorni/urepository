#this pack of rules is yet to be implemented in the workflow... they are  not even necessary in all cases
#DO NOT INCLUDE THIS RULE IF NOT NECESSARY, PLEASE, TRY TO GET ALL FILES BEFORE!!!

rule get_genome:
    output:
        "resources/genome.fasta"
    log:
        "logs/get-genome.log"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    wrapper:
        "0.59.2/bio/reference/ensembl-sequence"

checkpoint genome_faidx:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.fasta.fai"
    log:
        "logs/genome-faidx.log"
    wrapper:
        "0.59.2/bio/samtools/faidx"

rule genome_dict:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.dict"
    log:
        "logs/samtools/create_dict.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools dict {input} > {output} 2> {log} " ###NON AVEVO IL WRAPPER PER GATK??


rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/genome.fasta.fai"
    output:
        vcf="resources/variation.vcf.gz"
    log:
        "logs/get-known-variants.log"
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all"
    wrapper:
        "0.59.2/bio/reference/ensembl-variation"

rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz"
    output:
        "resources/variation.noiupac.vcf.gz"
    log:
        "logs/fix-iupac-alleles.log"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"

rule tabix_known_variants:
    input:
        "resources/variation.noiupac.vcf.gz"
    output:
        "resources/variation.noiupac.vcf.gz.tbi"
    log:
        "logs/tabix/variation.log"
    params:
        "-p vcf"
    wrapper:
        "0.59.2/bio/tabix"

rule bwa_index:
    input:
        "resources/genome.fasta"
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    resources:
        mem_mb=369000
    wrapper:
        "0.59.2/bio/bwa/index"