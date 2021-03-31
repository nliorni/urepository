## Google DeepVariant ##
use rule deepvariant_calling as s_deepvariant_calling with:
    output:
        vcf="results/s_deepcalls/{sample}.vcf" 

## SnpEff Annotate ##
use rule snpeff_Annotate as s_dv_snpeff_Annotate with:
    input:
        calls="results/s_deepcalls/{sample}.vcf",
        db=config["snpeff"]
    output:
        calls="results/s_dv_annotated/{sample}.vcf",
        stats="results/s_dv_annotated/{sample}.html",
        csvstats="results/s_dv_annotated/{sample}.csv"
    log:
        "logs/snpeff/{sample}.log"

## SnpSift Annotate ##
use rule snpsift_Annotate as s_dv_snpsift_Annotate with:
    input: 
        call="results/s_dv_annotated/{sample}.vcf",
        database=config["dbsnp"]
    output:
        call="results/s_dv_sift_annotated/{sample}.vcf"
    log:
        "logs/snpsift/annotate/{sample}.log"


## SnpSift VariantType ##
use rule snpsift_VarType as s_dv_snpsift_VarType with:
    input: 
        vcf="results/s_dv_sift_annotated/{sample}.vcf"
    output:
        vcf="results/s_dv_sift_vartyped/{sample}.vcf"
    log:
        "logs/snpsift/vartype/{sample}.log"

## SnpSift dbNSFP ##
use rule snpsift_dbNSFP as s_dv_snpsift_dbNSFP with:
    input:
        call="results/s_dv_sift_vartyped/{sample}.vcf",
        dbNSFP = config["dbnsfp"]
    output:
        call="results/s_dv_sift_dbNSFP/{sample}.vcf"
    log:
        "logs/dbNSFP/{sample}.log"

## SnpSift ExtractFields ##
rule s_dv_snpsift_ExtractFields:
    input:
        "results/s_dv_sift_dbNSFP/{sample}.vcf"
    output:
        "results/dv_csvfile/{sample}.csv"
    message:
        "Running SnpSift ExtractFields. Extracting fields of interest from the completly annotated vcf file '{input}' into '{output}'"
    shell:
        "python3 scripts/extractfields.py --input {input} > {output}"
