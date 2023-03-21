## Gatk HaplotypeCaller ##
rule gatk_HaplotypeCaller:
    input:
        # single or list of bam files
        bam=outputDir+"results/{sample}/{sample}.bqsr.bam",
        bai=outputDir+"results/{sample}/{sample}.bqsr.bam.bai",
        ref=config["reference"],
        dic=config["dict"],
        regions=config["bed"],
        known=config["dbsnp"] # optional
    output:
        gvcf=outputDir+"results/{sample}/{sample}.g.vcf",
        idx=outputDir+"results/{sample}/{sample}.g.vcf.idx"
    message:
        "running gatk haplotypecaller for sample {wildcards.sample}"
    params:
        java_opts=config["gatk_HaplotypeCaller"]["params"]["java_opts"],
        extra=config["gatk_HaplotypeCaller"]["params"]["extra"]
    shell:
        "gatk --java-options {params.java_opts} HaplotypeCaller -I {input.bam} -R {input.ref} -L {input.regions} -ERC GVCF {params.extra} -O {output.gvcf} --dbsnp {input.known}"
    
 ## Gatk GenotypeGvcfs ##
rule gatk_GenotypeGvcfs:
    input:
        gvcf=outputDir+"results/{sample}/{sample}.g.vcf",  # combined gvcf over multiple samples calls/all.g.vcf
        ref=config["reference"]
    output:
        vcf=temp(outputDir+"results/{sample}/{sample}.vcf"),
        idx=temp(outputDir+"results/{sample}/{sample}.vcf.idx")
    message:
        "running gatk genotypegvcfs for sample {wildcards.sample}"
    params:
        extra=config["gatk_GenotypeGvcfs"]["params"]["extra"],  
        java_opts=config["gatk_GenotypeGvcfs"]["params"]["java_opts"], 
    threads: config["gatk_GenotypeGvcfs"]["threads"]
    shell:
        "gatk --java-options {params.java_opts} GenotypeGVCFs --variant {input.gvcf} --reference {input.ref} --output {output.vcf}"

## Gatk SelectVariants (Snps) ##
rule gatk_snps_SelectVariants:
    input:
        vcf=outputDir+"results/{sample}/{sample}.vcf",
        ref=config["reference"],
    output:
        vcf=temp(outputDir+"results/{sample}/{sample}.snps.vcf"),
        idx=temp(outputDir+"results/{sample}/{sample}.snps.vcf.idx")
    message:
        "running gatk selectvariants (snps) for sample {wildcards.sample}"
    params:
        extra=config["gatk_SelectVariants"]["snps"]["extra"],  # optional filter arguments, see GATK docs
        java_opts=config["gatk_SelectVariants"]["snps"]["java_opts"]
    shell:
        "gatk SelectVariants --variant {input.vcf} --reference {input.ref} {params.extra} --output {output.vcf}"

## Gatk SelectVariants (Indels) ##
rule gatk_indel_SelectVariants:
    input:
        vcf=outputDir+"results/{sample}/{sample}.vcf",
        ref=config["reference"],
    output:
        vcf=temp(outputDir+"results/{sample}/{sample}.indel.vcf"),
        idx=temp(outputDir+"results/{sample}/{sample}.indel.vcf.idx")
    message:
        "running gatk selectvariants (indels) for sample {wildcards.sample}"
    params:
        extra=config["gatk_SelectVariants"]["indel"]["extra"],  # optional filter arguments, see GATK docs
        java_opts=config["gatk_SelectVariants"]["indel"]["java_opts"]
    shell:
        "gatk SelectVariants --variant {input.vcf} --reference {input.ref} {params.extra} --output {output.vcf}"

## Gatk VariantFiltration (Snps) ##   
rule gatk_snps_VariantFiltration:
    input:
        vcf=outputDir+"results/{sample}/{sample}.snps.vcf",
        vcf_idx = outputDir+"results/{sample}/{sample}.snps.vcf.idx",
        ref=config["reference"],
    output:
        vcf=temp(outputDir+"results/{sample}/{sample}.snps.filtered.vcf"),
        idx=temp(outputDir+"results/{sample}/{sample}.snps.filtered.vcf.idx")
    message:
        "running gatk variantfiltrations (snps) for sample: {wildcards.sample}"
    params:
        #filters=config["gatk_FilterVariants"]["snps"]["filters"],
        extra=config["gatk_FilterVariants"]["snps"]["extra"],  
        java_opts=config["gatk_FilterVariants"]["snps"]["java_opts"]
    shell:
        "gatk VariantFiltration --variant {input.vcf} --reference {input.ref} --cluster-size 3 --cluster-window-size 10 --filter-name '100QUAL' --filter-expression 'QUAL < 100' --filter-name 'SnpWarning' --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --output {output.vcf}"

## Gatk VariantFiltration (Indels) ##
rule gatk_indel_VariantFiltration:
    input:
        vcf=outputDir+"results/{sample}/{sample}.indel.vcf",
        vcf_idx = outputDir+"results/{sample}/{sample}.indel.vcf.idx",
        ref=config["reference"],
    output:
        vcf=temp(outputDir+"results/{sample}/{sample}.indel.filtered.vcf"),
        idx=temp(outputDir+"results/{sample}/{sample}.indel.filtered.vcf.idx")
    message:
        "running gatk variantfiltrations (indels) for sample: {wildcards.sample}"
    params:
        #filters=config["gatk_FilterVariants"]["indels"]["filters"],
        extra=config["gatk_FilterVariants"]["indels"]["extra"],  
        java_opts=config["gatk_FilterVariants"]["indels"]["java_opts"]
    shell:
        "gatk VariantFiltration --variant {input.vcf} --reference {input.ref} --cluster-size 3 --cluster-window-size 10 --filter-name '100QUAL' --filter-expression 'QUAL < 100' --filter-name 'IndelWarning' --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --output {output.vcf}"


## add variant recalibration step?


## Picard MergeVcfs ##
rule picard_MergeVcfs:
    input:
        snps=outputDir+"results/{sample}/{sample}.snps.filtered.vcf",
        indel=outputDir+"results/{sample}/{sample}.indel.filtered.vcf"
    output:
        vcf=outputDir+"results/{sample}/{sample}.filter.vcf",
        idx=temp(outputDir+"results/{sample}/{sample}.filter.vcf.idx")
    message:
        "running picard mergevcfs for sample {wildcards.sample}"
    params:
        extra=config["picard_MergeVcfs"]["extra"]
    shell:
        "picard MergeVcfs {params.extra} --INPUT {input.snps} --INPUT {input.indel} --OUTPUT {output.vcf}"   