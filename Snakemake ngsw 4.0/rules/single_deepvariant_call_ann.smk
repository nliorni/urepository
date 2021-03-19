# #DEEPVARIANT
rule s_deepvariant_calling:
   input:
       bam="sorted_reads/{sample}.sorted.bam",
       bai="sorted_reads/{sample}.sorted.bam.bai",
       ref=config["reference"],
       dic=config["dict"],
       regions=config["bed"]
   output:
       vcf="s_deepcalls/{sample}.vcf"    #try with g.vcf
   params:
       model="wgs",   # {wgs, wes}
       extra="--regions chr1:43389920-43450402" #prova --regions {regions}, regions=config["regions"]
   log:
       "logs/deepvariant/{sample}/stdout.log"
   message:
        "Trying out DeepVariant: deepcalling {input.bam} -based on {input.ref}- into {output.vcf}"
   wrapper:
       "0.68.0/bio/deepvariant"

#ANNOTATE VCF WITH SNPEFF
rule s_dv_snpeff_Annotate:
    input:
        calls="s_deepcalls/{sample}.vcf", # (vcf, bcf, or vcf.gz)
        db=config["snpeff"] # path to reference db downloaded with the snpeff download wrapper
    output:
        #multiext("snpeff/{sample}", ".vcf", ".html", ".csv")
        calls=temp("s_dv_annotated/{sample}.vcf"),   # annotated calls (vcf, bcf, or vcf.gz)
        stats="s_dv_annotated/{sample}.html",  # summary statistics (in HTML), optional
        csvstats="s_dv_annotated/{sample}.csv" # summary statistics in CSV, optional
    log:
        "logs/snpeff/{sample}.log"
    message:
        "Annotating '{input.calls}' with '{input.db}' to generate '{output.calls}', '{output.stats}' and '{output.csvstats}' with SnpEff Annotate"
    params:
        extra="-Xmx8g"           # optional parameters (e.g., max memory 4g)
    wrapper:
        "0.66.0/bio/snpeff/annotate"

#SNPSIFT ANNOTATE
rule s_dv_snpsift_Annotate:
    input:
        call="s_dv_annotated/{sample}.vcf",
        database=config["dbsnp"]
    output:
        call=temp("s_dv_sift_annotated/{sample}.vcf")
    message:
        "Further annotating '{input.call}' using '{input.database}' creating '{output.call}' with SnpSift Annotate"
    log:
        "logs/snpsift/annotate/{sample}.log"
    wrapper:
        "0.67.0/bio/snpsift/annotate"


#SNPSIFT VARTYPE
rule s_dv_snpsift_VarType:
    input:
        vcf="s_dv_sift_annotated/{sample}.vcf"
    output:
        vcf=temp("s_dv_sift_vartyped/{sample}.vcf")
    message:
        "Further annotating '{input.vcf}' creating '{output.vcf}' with SnpSift VarType"
    log:
        "logs/snpsift/vartype/{sample}.log"
    wrapper:
        "0.67.0/bio/snpsift/varType"

#dbNSFP
rule s_dv_snpsift_dbNSFP:
   input:
       call = "s_dv_sift_vartyped/{sample}.vcf",
       dbNSFP = config["dbnsfp"]
   output:
       call = "s_dv_sift_dbNSFP/{sample}.vcf"
   message:
       "Further annotating '{input.call}' using '{input.dbNSFP}', creating '{output.call}' with SnpSift dbNSFP"
   log:
       "logs/dbNSFP/{sample}.log"
   wrapper:
       "0.67.0/bio/snpsift/dbnsfp"

#EXTRACT FIELDS

rule s_dv_snpsift_ExtractFields:
    input:
        "s_dv_sift_dbNSFP/{sample}.vcf"
    output:
        "dv_csvfile/{sample}.csv"
    message:
        "Extracting fields of interest from the completly annotated vcf file {input} into {output}"
    shell:
        "python3 scripts/extractfields.py --input {input} > {output}"
