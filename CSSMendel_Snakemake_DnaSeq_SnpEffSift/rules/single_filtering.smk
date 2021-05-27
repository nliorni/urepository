## Gatk SelectVariants (Snps) ##
use rule gatk_snps_SelectVariants as s_gatk_snps_SelectVariants with:
    input:
        vcf="results/s_genotyped/{sample}.vcf",
        ref=config["reference"]
    output:
        vcf="results/select/{sample}.snps.vcf"
    log:
        "logs/gatk/select/{sample}.snps.log"


## Gatk SelectVariants (Indels) ##
use rule gatk_indel_SelectVariants as s_gatk_indel_SelectVariants with:
    input:
        vcf="results/s_genotyped/{sample}.vcf",
        ref=config["reference"]
    output:
        vcf="results/select/{sample}.indel.vcf"
    log:
        "logs/gatk/select/{sample}.indel.log"

## Gatk VariantFiltration (Snps) ##   
use rule gatk_snps_VariantFiltration as s_gatk_snps_VariantFiltration with:
    input:
        vcf="results/select/{sample}.snps.vcf",
        ref=config["reference"]
    output:
        vcf="results/filtered/{sample}.snps.filtered.vcf"
    log:
        "logs/gatk/filter/{sample}.snps.log"
        
## Gatk VariantFiltration (Indels) ##
use rule gatk_indel_VariantFiltration as s_gatk_indel_VariantFiltration with:
    input:
        vcf="results/select/{sample}.indel.vcf",
        ref=config["reference"]
    output:
        vcf="results/filtered/{sample}.indel.filtered.vcf"
    log:
        "logs/gatk/filter/{sample}.indel.log"

## Picard MergeVcfs ##
rule s_picard_MergeVcfs:
    input:
        vcfs=["results/filtered/{sample}.snps.filtered.vcf", "results/filtered/{sample}.indel.filtered.vcf"]
    output:
        "results/merged/{sample}.vcf"
    log:
        "logs/picard/{sample}.mergevcfs.log"
    message:
        "Running Gatk MergeVcfs. Merging the filtered snps and indel variants {input.vcfs} into {output}"
    params:
        extra=""
    wrapper:
        "0.72.0/bio/picard/mergevcfs"
