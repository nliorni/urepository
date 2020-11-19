#SNPEFF DOWNLOAD --> VERIFICARE CHE FUNZIONI, A CASA!!!!!!!!!
#rule snpeff_download:
#    output:
#        # wildcard {reference} may be anything listed in `snpeff databases`
#        directory("resources/snpeff/{reference}")
#    log:
#        "logs/snpeff/download/{reference}.log"
#    params:
#        reference="{reference}"
#    wrapper:
#        "0.66.0/bio/snpeff/download"

#ANNOTATE VCF WITH SNPEFF
rule snpeff_annotation:
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
        "Annotating '{input.calls}' with '{input.db}' to generate '{output.calls}', '{output.stats}' and '{output.csvstats}' with SNPeff Annotate"
    params:
        extra="-Xmx8g"           # optional parameters (e.g., max memory 4g)
    wrapper:
        "0.66.0/bio/snpeff/annotate"

#SNPSIFT ANNOTATE
rule snpsift_annotate:
    input:
        call="annotated/all.vcf",
        database=config["snpdb2"]
    output:
        call="sift_annotated/all.vcf"
    message:
        "Further annotating '{input.call}' using '{input.database}' creating '{output.call}'"
    log:
        "logs/snpsift/annotate/all.log"
    wrapper:
        "0.67.0/bio/snpsift/annotate"


#SNPSIFT VARTYPE
rule snpsift_vartype:
    input:
        vcf="sift_annotated/all.vcf"
    output:
        vcf="sift_vartyped/all.vcf"
    message:
        "Further annotating '{input.vcf}' creating '{output.vcf}'"
    log:
        "logs/snpsift/vartype/all.log"
    wrapper:
        "0.67.0/bio/snpsift/varType"

#GATK VARIANTSTOTABLE
rule gatk_variants2table:
    input:
        "sift_vartyped/all.vcf"
    output:
        "sift_vartyped/all_annotated.csv"
    log:
        "logs/v2t/all.log"
    message:
        "Converting '{input}' into '{output}' with Gatk VariantsToTable"
    shell:
        "gatk VariantsToTable -V {input} -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F ANN -O {output}"
