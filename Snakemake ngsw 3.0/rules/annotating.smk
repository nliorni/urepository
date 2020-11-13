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
        db="resources/snpeff/hg38" # path to reference db downloaded with the snpeff download wrapper
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



###MAKE ALL THE NOT STRICTLY NEEDED FILES TEMPORARY !!!!

rule gatk_variants2table:
    input:
        "annotated/all.vcf"
    output:
        "annotated/all_annotated.csv"
    log:
        "logs/v2t/all.log"
    message:
        "Converting {input} into {output} with Gatk VariantsToTable"
    shell:
        "gatk VariantsToTable -V {input} -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F ANN -O {output}"
