## Gatk HaplotypeCaller ##
rule gatk_HaplotypeCaller:
    input:
        # single or list of bam files
        bam="results/sorted_bam/{sample}.sorted.bam",
        bai="results/sorted_bam/{sample}.sorted.bam.bai",
        ref=config["reference"],
        dic=config["dict"],
        regions=config["bed"]
        #known="dbsnp.vcf"  # optional
    output:
        gvcf=temp("results/calls/{sample}.g.vcf")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    message:
        "Running Gatk HaplotypeCaller. Calling variants in '{input.bam}' refering to '{input.ref}' in regions '{input.regions}'."
    params:
        java_opts=config["gatk_HaplotypeCaller"]["params"]["java_opts"],
        extra=config["gatk_HaplotypeCaller"]["params"]["extra"]
    shell:
        "gatk HaplotypeCaller --java-options {params.java_opts} -I {input.bam}  -R {input.ref} -L {input.regions} -ERC GVCF {params.extra} -O {output.gvcf}"
        #"0.66.0/bio/gatk/haplotypecaller"

## Gatk CombineGvcfs ##
rule gatk_CombineGvcfs:
    input:
        gvcfs=expand("results/calls/{sample}.g.vcf", sample=config["samples"]),
        ref=config["reference"]
    output:
        gvcf=temp("results/combined/all.g.vcf")
    log:
        "logs/gatk/combinegvcfs/all.log"
    message:
        "Running Gatk CombineGVcfs. Combining '{input.gvcfs}', with reference '{input.ref}', into '{output.gvcf}'."
    params:
        extra=config["gatk_CombineGvcfs"]["params"]["extra"],  
        java_opts=config["gatk_CombineGvcfs"]["params"]["java_opts"],  
    wrapper:
        "0.67.0/bio/gatk/combinegvcfs"


## Gatk GenotypeGvcfs ##
rule gatk_GenotypeGvcfs:
    input:
        gvcf="results/combined/all.g.vcf",  # combined gvcf over multiple samples calls/all.g.vcf
        ref=config["reference"]
    output:
        vcf=temp("results/genotyped/all.vcf")
    log:
        "logs/gatk/genotypegvcfs/all.log"
    message:
        "Running Gatk GenotypeGVcfs. Genotyping the g.vcf '{input.gvcf}', with reference '{input.ref}', into '{output.vcf}'."
    params:
        extra=config["gatk_GenotypeGvcfs"]["params"]["extra"],  
        java_opts=config["gatk_GenotypeGvcfs"]["params"]["java_opts"], 
    threads: 10
    wrapper:
        "0.67.0/bio/gatk/genotypegvcfs"

#VARIANT FILTRATION
include: "filtering.smk"

## SnpEff Annotate ##
rule snpeff_Annotate:
    input:
        calls="results/merged/all.vcf", # (vcf, bcf, or vcf.gz)
        db=config["snpeff"] # path to reference db downloaded with the snpeff download wrapper
    output:
        calls="results/annotated/all.vcf",   # annotated calls (vcf, bcf, or vcf.gz)
        stats="results/annotated/all.html",  # summary statistics (in HTML), optional
        csvstats="results/annotated/all.csv" # summary statistics in CSV, optional
    log:
        "logs/snpeff/annotate/all.log"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=100000
    message:
        "Running SnpEff Annotate. Annotating '{input.calls}' with '{input.db}' to generate '{output.calls}', '{output.stats}' and '{output.csvstats}'."
    wrapper:
        "0.73.0/bio/snpeff/annotate"

## SnpSift Annotate ##
rule snpsift_Annotate:
    input:
        call="results/annotated/all.vcf",
        database=config["dbsnp"]
    output:
        call="results/sift_annotated/all.vcf"
    message:
        "Running SnpSift Annotate. Further annotating '{input.call}' using '{input.database}' creating '{output.call}'."
    params:
        java_opts=config["snpsift_Annotate"]["params"]["java_opts"],
        extra=config["snpsift_Annotate"]["params"]["extra"],
        java_path=config["java_path"],
        snpsift_path=config["snpsift_path"]
    log:
        "logs/snpsift/annotate/all.log"
    shell:
        "{params.java_path} -jar {params.java_opts} {params.snpsift_path} annotate {params.extra} {input.database} {input.call} > {output.call}"


## SnpSift VariantType ##
rule snpsift_VarType:
    input:
        vcf="results/sift_annotated/all.vcf"
    output:
        vcf=temp("results/sift_vartyped/all.vcf")
    message:
        "Running SnpSift VarType. Further annotating '{input.vcf}' creating '{output.vcf}'."
    params:
        java_opts=config["snpsift_VarType"]["params"]["java_opts"],
        extra=config["snpsift_VarType"]["params"]["extra"],
        java_path=config["java_path"],
        snpsift_path=config["snpsift_path"]
    log:
        "logs/snpsift/vartype/all.log"
    shell:
        "{params.java_path} -jar {params.snpsift_path} VarType {input.vcf} > {output.vcf}"

## SnpSift dbNSFP ##
rule snpsift_dbNSFP:
   input:
       call = "results/sift_vartyped/all.vcf",
       dbNSFP = config["dbnsfp"]
   output:
       call = "results/sift_dbNSFP/all.vcf"
   message:
       "Running SnpSift dbNSFP. Further annotating '{input.call}' using '{input.dbNSFP}', creating '{output.call}'."
   params:
        java_opts=config["snpsift_dbNSFP"]["params"]["java_opts"],
        extra=config["snpsift_dbNSFP"]["params"]["extra"],
        java_path=config["java_path"],
        snpsift_path=config["snpsift_path"]
   log:
       "logs/snpsift/dbNSFP/all.log"
   shell:
       "{params.java_path} -jar {params.java_opts} {params.snpsift_path} dbnsfp  {params.extra} -db {input.dbNSFP} {input.call} > {output.call}"


## SnpSift ExtractFields ##
rule snpsift_ExtractFields:
    input:
        "results/sift_dbNSFP/all.vcf"
    output:
        "results/csvfile/all.csv"
    log:
        "logs/snpsift/ExtractFields/all.log"
    message:
        "Running SnpSift ExtractFields. Extracting fields of interest from the completly annotated vcf file '{input}' into '{output}'"
    shell:
        "python3 scripts/extractfields.py --input {input} > {output}"

## SnpSift ExtractFields ##
# rule snpsift_ExtractFields_AG:
#     input:
#         "results/sift_dbNSFP/all.vcf"
#     output:
#         "results/csvfile/all.csv"
#     log:
#         "logs/snpsift/ExtractFields/all.log"
#     params: path=config["extractfields"]["path"]
#     message:
#         "Running SnpSift ExtractFields. Extracting fields of interest from the completly annotated vcf file '{input}' into '{output}'"
#     shell:
#         "python3 scripts/extractfields.py --path {params.path} --input {input} > {output}"

## Check Gene Lists ##
# rule genelists_check:
#     input:
#         csv="csvfile/all.csv",
#         genelists=config["genelists"]
#     output:
#         "genelistedcsv/all.csv"
#     message:
#         "Checking if the genes in the specified lists are present in the performed annotation"
#     shell:
#         "python3 genelists.py --input {input.csv} --genelists {input.genelists} > {output}"
