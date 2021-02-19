#BASE CALLING WITH GATK
rule gatk_HaplotypeCaller:
    input:
        # single or list of bam files
        bam="sorted_reads/{sample}.sorted.bam",
        bai="sorted_reads/{sample}.sorted.bam.bai",
        ref=config["reference"],
        dic=config["dict"],
        regions=config["bed"]
        #known="dbsnp.vcf"  # optional
    output:
        gvcf=temp("calls/{sample}.g.vcf")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    message:
        "Calling variants in '{input.bam}' refering to '{input.ref}' in regions '{input.regions}' using Gatk HaplotypeCaller... -Prepare to Wait Editionn (if there's no bed)-"
    params:
        #extra=dd  # optional
        java_opts="-Xmx8G" #-L {input.bed}"
    shell:
        "gatk HaplotypeCaller -I {input.bam}  -R {input.ref} -L {input.regions} -ERC GVCF -O {output.gvcf}"
        #"0.66.0/bio/gatk/haplotypecaller"

#COMBINE GVCFs
rule gatk_CombineGvcfs:
    input:
        gvcfs=expand("calls/{sample}.g.vcf", sample=SAMPLES),
        ref=config["reference"]
    output:
        gvcf=temp("combined/all.g.vcf")
    log:
        "logs/gatk/combinegvcfs.log"
    message:
        "Combining '{input.gvcfs}', with reference '{input.ref}', into '{output.gvcf}' with Gatk CombineGVcfs"
    params:
        extra="",  # optional
        java_opts="",  # optional
    wrapper:
        "0.67.0/bio/gatk/combinegvcfs"


#GENOTYPE GVCFs
rule gatk_GenotypeGvcfs:
    input:
        gvcf="combined/all.g.vcf",  # combined gvcf over multiple samples calls/all.g.vcf
        ref=config["reference"]
    output:
        vcf=temp("genotyped/all.vcf")
    log:
        "logs/gatk/genotypegvcfs.log"
    message:
        "Genotyping the combined g.vcf '{input.gvcf}', with reference '{input.ref}', into '{output.vcf}' with Gatk GenotypeGVcfs "
    params:
        extra="",  # optional
        java_opts="", # optional
    threads: 4
    wrapper:
        "0.67.0/bio/gatk/genotypegvcfs"


#ANNOTATE VCF WITH SNPEFF
rule snpeff_Annotate:
    input:
        calls="genotyped/all.vcf", # (vcf, bcf, or vcf.gz)
        db=config["snpeff"] # path to reference db downloaded with the snpeff download wrapper
    output:
        #multiext("snpeff/{sample}", ".vcf", ".html", ".csv")
        calls="annotated/all.vcf",   # annotated calls (vcf, bcf, or vcf.gz)
        stats="annotated/all.html",  # summary statistics (in HTML), optional
        csvstats="annotated/all.csv" # summary statistics in CSV, optional
    log:
        "logs/snpeff/all.log"
    message:
        "Annotating '{input.calls}' with '{input.db}' to generate '{output.calls}', '{output.stats}' and '{output.csvstats}' with SnpEff Annotate"
    params:
        extra="-Xmx8g"           # optional parameters (e.g., max memory 4g)
    wrapper:
        "0.66.0/bio/snpeff/annotate"

#SNPSIFT ANNOTATE
rule snpsift_Annotate:
    input:
        call="annotated/all.vcf",
        database=config["dbsnp"]
    output:
        call=temp("sift_annotated/all.vcf")
    message:
        "Further annotating '{input.call}' using '{input.database}' creating '{output.call}' with SnpSift Annotate"
    log:
        "logs/snpsift/annotate/all.log"
    wrapper:
        "0.67.0/bio/snpsift/annotate"


#SNPSIFT VARTYPE
rule snpsift_VarType:
    input:
        vcf="sift_annotated/all.vcf"
    output:
        vcf=temp("sift_vartyped/all.vcf")
    message:
        "Further annotating '{input.vcf}' creating '{output.vcf}' with SnpSift VarType"
    log:
        "logs/snpsift/vartype/all.log"
    wrapper:
        "0.67.0/bio/snpsift/varType"

#dbNSFP
rule snpsift_dbNSFP:
   input:
       call = "sift_vartyped/all.vcf",
       dbNSFP = config["dbnsfp"]
   output:
       call = "sift_dbNSFP/all.vcf"
   message:
       "Further annotating '{input.call}' using '{input.dbNSFP}', creating '{output.call}' with SnpSift dbNSFP"
   log:
       "logs/dbNSFP/all.log"
   wrapper:
       "0.67.0/bio/snpsift/dbnsfp"