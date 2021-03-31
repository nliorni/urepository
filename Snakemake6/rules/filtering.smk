## Gatk SelectVariants (Snps) ##
rule gatk_snps_SelectVariants:
    input:
        vcf="results/genotyped/all.vcf",
        ref=config["reference"],
    output:
        vcf="results/select/all.snps.vcf"
    log:
        "logs/gatk/select/snvs.log"
    message:
        "Running Gatk SelectVariants. Selecting snps variant from the vcf file {input.vcf} and storing them in {output.vcf}"
    params:
        extra=config["gatk_SelectVariants"]["snps"]["extra"],  # optional filter arguments, see GATK docs
        java_opts=""
    wrapper:
        "0.72.0/bio/gatk/selectvariants"
    
## Gatk SelectVariants (Indels) ##
rule gatk_indel_SelectVariants:
    input:
        vcf="results/genotyped/all.vcf",
        ref=config["reference"],
    output:
        vcf="results/select/all.indel.vcf"
    log:
        "logs/gatk/select/indel.log"
    message:
        "Running Gatk SelectVariants. Selecting indel variant from the vcf file {input.vcf} and storing them in {output.vcf}"
    params:
        extra=config["gatk_SelectVariants"]["indel"]["extra"],  # optional filter arguments, see GATK docs
        java_opts=""
    wrapper:
        "0.72.0/bio/gatk/selectvariants"

## Gatk VariantFiltration (Snps) ##   
rule gatk_snps_VariantFiltration:
    input:
        vcf="results/select/all.snps.vcf",
        ref=config["reference"],
    output:
        vcf="results/filtered/all.snps.filtered.vcf"
    log:
        "logs/gatk/filter/snps.log"
    message:
        "Running Gatk VariantFiltration. Filtering snps variants in {input.vcf} and storing the result in {output.vcf}"
    params:
        filters=config["gatk_FilterVariants"]["snps"]["filters"],
        extra="",  
        java_opts=""
    wrapper:
        "0.72.0/bio/gatk/variantfiltration"

## Gatk VariantFiltration (Indels) ##
rule gatk_indel_VariantFiltration:
    input:
        vcf="results/select/all.indel.vcf",
        ref=config["reference"],
    output:
        vcf="results/filtered/all.indel.filtered.vcf"
    log:
        "logs/gatk/filter/indel.log"
    message:
        "Running Gatk VariantFiltration. Filtering indel variants in {input.vcf} and storing the result in {output.vcf}"
    params:
        filters=config["gatk_FilterVariants"]["indels"]["filters"],
        extra="",  
        java_opts=""
    wrapper:
        "0.72.0/bio/gatk/variantfiltration"

## Picard MergeVcfs ##
rule picard_MergeVcfs:
    input:
        vcfs=["results/filtered/all.snps.filtered.vcf", "results/filtered/all.indel.filtered.vcf"]
    output:
        "results/merged/all.vcf"
    log:
        "logs/picard/mergevcfs.log"
    message:
        "Running Gatk MergeVcfs. Merging the filtered snps and indel variants {input.vcfs} into {output}"
    params:
        extra=""
    wrapper:
        "0.72.0/bio/picard/mergevcfs"

