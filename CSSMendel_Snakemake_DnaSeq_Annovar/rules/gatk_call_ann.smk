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

## Annovar VCF to Human ##
rule Annovar_vcf2human:
    input: 
        vcf="results/merged/all.vcf"
    output:
        human="results/annovar/all.human"
    message:
        "Running Annovar Vcf2Human"
    shell:
        "python3 /software/VCF2Human/vcf2human.py -v {input.vcf} -u {output.human} "


## Annovar Process Human ##
rule Annovar_processhuman:
    input:
        human="results/annovar/all.human"
    output:
        "results/annovar/all_standard_monoallelic_annotated_final.tsv"
    params:
        ip="192.168.10.19",
        port="27018",
        annovar="/software/annovar/table_annovar.pl",
        humandb="/software/annovar/humandb",
        extra=" -sw 10 --complete False --genome-version hg19"
    message: 
        "Running Annovar ProcessHuman"
    shell:
        "python3 scripts/processhuman_NL.py {input} -a standard --annovarcmd {params.annovar} --humandb {params.humandb} --dbsnp_ver 151 --ip {params.ip} --port {params.port} {params.extra} -d results/annovar"

## dbNSFP Annotator ##
rule dbNSFP_Annotator:
    input: 
        "results/annovar/all_standard_monoallelic_annotated_final.tsv"
    output:
        "results/annovar/all_standard_monoallelic_annotated_final_dbnsfp.xlsx"
    params:
        ip="192.168.10.19",
        port="27018"
    message:
        "Running dbNSFP Annotator."
    shell:
        "python3 /software/dbNSFP_Annotator/dbnsfp_annotator_MT.py -v 40 -i {params.ip} -p {params.port} -r all -u {input}"

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
