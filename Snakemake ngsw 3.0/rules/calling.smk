#bam quality control with KING algorithm??
###introdurre l'opzione -L per il file BED
#BASE CALLING WITH GATK
rule haplotype_caller:
    input:
        # single or list of bam files
        #bed=config["bed"]
        bam="sorted_reads/{sample}.sorted.bam",
        bai="sorted_reads/{sample}.sorted.bam.bai",
        ref=config["reference"],
        dic=config["dict"]
        # known="dbsnp.vcf"  # optional
    output:
        gvcf="calls/{sample}.g.vcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    message:
        "Calling variants in '{input.bam}' refering to '{input.ref}' using Gatk HaplotypeCaller... -Prepare to Wait Edition-"
    params:
        #extra=dd  # optional
        java_opts="-Xmx8G" #-L {input.bed}"
    wrapper:
        "0.66.0/bio/gatk/haplotypecaller"

#COMBINE GVCFs
rule combine_gvcfs:
    input:
        gvcfs=expand("calls/{sample}.g.vcf", sample=SAMPLES),
        ref=config["reference"]
    output:
        gvcf="combined/all.g.vcf"
    log:
        "logs/gatk/combinegvcfs.log"
    message:
        "Combining '{input.gvcfs}', with reference '{input.ref}', into '{output.gvcf}'"
    params:
        extra="",  # optional
        java_opts="",  # optional
    wrapper:
        "0.67.0/bio/gatk/combinegvcfs"


#GENOTYPE GVCFs
rule genotype_gvcfs:
    input:
        gvcf="combined/all.g.vcf",  # combined gvcf over multiple samples calls/all.g.vcf
        ref=config["reference"]
    output:
        vcf="genotyped/all.vcf",
    log:
        "logs/gatk/genotypegvcfs.log"
    message:
        "Genotyping the combined g.vcf '{input.gvcf}', with reference '{input.ref}', into '{output.vcf}' "
    params:
        extra="",  # optional
        java_opts="", # optional
    threads: 4
    wrapper:
        "0.67.0/bio/gatk/genotypegvcfs"

#CNV CALLING???

#Variant filtration???
