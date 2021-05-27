## Gatk HaplotypeCaller ##
use rule gatk_HaplotypeCaller as s_gatk_HaplotypeCaller with:
    output:
        gvcf=temp("results/s_calls/{sample}.g.vcf")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"

## Gatk GenotypeGvcfs ##
use rule gatk_GenotypeGvcfs as s_gatk_GenotypeGvcfs with:
    input: 
        gvcf="results/s_calls/{sample}.g.vcf",
        ref=config["reference"]
    output:
        vcf="results/s_genotyped/{sample}.vcf"
    log:
        "logs/gatk/genotype_gvcfs/{sample}.log"

#VARIANT FILTRATION
include: "single_filtering.smk"

## SnpEff Annotate ##
use rule snpeff_Annotate as s_snpeff_Annotate with:
    input:
        calls="results/merged/{sample}.vcf",
        db=config["snpeff"]
    output:
        calls="results/s_annotated/{sample}.vcf",
        stats="results/stats/snpeff_ann/{sample}.html",
        csvstats="results/stats/snpeff_ann/{sample}.csv"
    log:
        "logs/snpeff/{sample}.log"

## SnpSift Annotate ##
use rule snpsift_Annotate as s_snpsift_Annotate with:
    input: 
        call="results/s_annotated/{sample}.vcf",
        database=config["dbsnp"]
    output:
        call=temp("results/s_sift_annotated/{sample}.vcf")
    log:
        "logs/snpsift/annotate/{sample}.log"


## SnpSift VariantType ##
use rule snpsift_VarType as s_snpsift_VarType with:
    input: 
        vcf="results/s_sift_annotated/{sample}.vcf"
    output:
        vcf=temp("results/s_sift_vartyped/{sample}.vcf")
    log:
        "logs/snpsift/vartype/{sample}.log"

## SnpSift dbNSFP ##
use rule snpsift_dbNSFP as s_snpsift_dbNSFP with:
    input:
        call="results/s_sift_vartyped/{sample}.vcf",
        dbNSFP = config["dbnsfp"]
    output:
        call="results/s_sift_dbNSFP/{sample}.vcf"
    log:
        "logs/dbNSFP/{sample}.log"

## SnpSift ExtractFields ##
rule s_snpsift_ExtractFields:
    input:
        "results/s_sift_dbNSFP/{sample}.vcf"
    output:
        "results/csvfile/{sample}.csv"
    message:
        "Running SnpSift ExtractFields. Extracting fields of interest from the completly annotated vcf file '{input}' into '{output}'"
    shell:
        "python3 scripts/extractfields.py --input {input} > {output}"