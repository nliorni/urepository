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


## Annovar VCF to Human ##
use rule Annovar_vcf2human as s_Annovar_vcf2human with:
    input:
        vcf="results/merged/{sample}.vcf"
    output:
        human="results/annovar/{sample}.human"

## Annovar Process Human ##
use rule Annovar_processhuman as s_Annovar_processhuman with:
    input:
        human="results/annovar/{sample}.human"
    output:
        "results/annovar/{sample}_cleaned_single_freq_standard_monoallelic_annotated_final.tsv"

## dbNSFP Annotator ##
use rule dbNSFP_Annotator as s_dbNSFP_Annotator with:
    input:
        "results/annovar/{sample}_cleaned_single_freq_standard_monoallelic_annotated_final.tsv"
    output:
        "results/annovar/{sample}_cleaned_single_freq_standard_monoallelic_annotated_final_dbnsfp.xlsx"