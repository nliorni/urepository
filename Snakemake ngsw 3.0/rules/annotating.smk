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

rule snpeff_annotation:
    input:
        calls="calls/{sample}.g.vcf", # (vcf, bcf, or vcf.gz)
        db="resources/snpeff/hg38" # path to reference db downloaded with the snpeff download wrapper
    output:
        #multiext("snpeff/{sample}", ".vcf", ".html", ".csv")
        calls="snpeff/{sample}.vcf",   # annotated calls (vcf, bcf, or vcf.gz)
        stats="snpeff/{sample}.html",  # summary statistics (in HTML), optional
        csvstats="snpeff/{sample}.csv" # summary statistics in CSV, optional
    log:
        "logs/snpeff/{sample}.log"
    params:
        extra="-Xmx8g"           # optional parameters (e.g., max memory 4g)
    wrapper:
        "0.66.0/bio/snpeff/annotate"
        
    #CNV variant calling???

###MAKE ALL THE NOT STRICTLY NEEDED FILES TEMPORARY !!!!