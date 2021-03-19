#BASE CALLING WITH GATK
rule s_gatk_HaplotypeCaller:
    input:
        # single or list of bam files
        bam="sorted_reads/{sample}.sorted.bam",
        bai="sorted_reads/{sample}.sorted.bam.bai",
        ref=config["reference"],
        dic=config["dict"],
        regions=config["bed"]
        #known="dbsnp.vcf"  # optional
    output:
        gvcf=temp("s_calls/{sample}.vcf")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    message:
        "Calling variants in '{input.bam}' refering to '{input.ref}' in regions '{input.regions}' using Gatk HaplotypeCaller... -Prepare to Wait Edition (if there's no bed)-"
    params:
        #extra=dd  # optional
        java_opts="-Xmx8G" #-L {input.bed}"
    shell:
        "gatk HaplotypeCaller -I {input.bam}  -R {input.ref} -L {input.regions} -O {output.gvcf}"
        #"0.66.0/bio/gatk/haplotypecaller"

# #GENOTYPE GVCFs
# rule s_gatk_GenotypeGvcf:
#     input:
#         gvcf="s_calls/{sample}.g.vcf",
#         ref=config["reference"]
#     output:
#         vcf="s_genotyped/{sample}.vcf"
#     log:
#         "logs/gatk/genotype_gvcfs/{sample}.log"
#     message:
#         "Genotyping the g.vcf '{input.gvcf}', with reference '{input.ref}', into '{output.vcf}' with Gatk GenotypeGVcfs "
#     params:
#         extra="",  # optional
#         java_opts="", # optional
#     threads: 4
#     wrapper:
#         "0.67.0/bio/gatk/genotypegvcfs"

#ANNOTATE VCF WITH SNPEFF
rule s_snpeff_Annotate:
    input:
        calls="s_calls/{sample}.vcf", # (vcf, bcf, or vcf.gz)
        db=config["snpeff"] # path to reference db downloaded with the snpeff download wrapper
    output:
        #multiext("snpeff/{sample}", ".vcf", ".html", ".csv")
        calls=temp("s_annotated/{sample}.vcf"),   # annotated calls (vcf, bcf, or vcf.gz)
        stats="stats/{sample}.html",  # summary statistics (in HTML), optional
        csvstats="stats/{sample}.csv" # summary statistics in CSV, optional
    log:
        "logs/snpeff/{sample}.log"
    message:
        "Annotating '{input.calls}' with '{input.db}' to generate '{output.calls}', '{output.stats}' and '{output.csvstats}' with SnpEff Annotate"
    params:
        extra="-Xmx8g"           # optional parameters (e.g., max memory 4g)
    wrapper:
        "0.66.0/bio/snpeff/annotate"

#SNPSIFT ANNOTATE
rule s_snpsift_Annotate:
    input:
        call="s_annotated/{sample}.vcf",
        database=config["dbsnp"]
    output:
        call=temp("s_sift_annotated/{sample}.vcf")
    message:
        "Further annotating '{input.call}' using '{input.database}' creating '{output.call}' with SnpSift Annotate"
    log:
        "logs/snpsift/annotate/{sample}.log"
    wrapper:
        "0.67.0/bio/snpsift/annotate"


#SNPSIFT VARTYPE
rule s_snpsift_VarType:
    input:
        vcf="s_sift_annotated/{sample}.vcf"
    output:
        vcf=temp("s_sift_vartyped/{sample}.vcf")
    message:
        "Further annotating '{input.vcf}' creating '{output.vcf}' with SnpSift VarType"
    log:
        "logs/snpsift/vartype/{sample}.log"
    wrapper:
        "0.67.0/bio/snpsift/varType"

#dbNSFP
rule s_snpsift_dbNSFP:
   input:
       call = "s_sift_vartyped/{sample}.vcf",
       dbNSFP = config["dbnsfp"]
   output:
       call = "s_sift_dbNSFP/{sample}.vcf"
   message:
       "Further annotating '{input.call}' using '{input.dbNSFP}', creating '{output.call}' with SnpSift dbNSFP"
   log:
       "logs/dbNSFP/{sample}.log"
   wrapper:
       "0.67.0/bio/snpsift/dbnsfp"

#EXTRACT FIELDS

rule s_snpsift_ExtractFields:
    input:
        "s_sift_dbNSFP/{sample}.vcf"
    output:
        "csvfile/{sample}.csv"
    message:
        "Extracting fields of interest from the completly annotated vcf file {input} into {output}"
    shell:
        "python3 scripts/extractfields.py --input {input} > {output}"