## VEP Annotate ##
rule VEP_annotate:
    input:
        calls=get_input_vcf(),
        #calls="results/{sample}/{sample}.filter.vcf",  # .vcf, .vcf.gz or .bcf
        cache=config["vep_annotate"]["cache"], # can be omitted if fasta and gff are specified
        plugins="/software/dataseq/databases/vep/plugins",
    output:
        calls=temp(outputDir+"results/{sample}/{sample}.vep.vcf"),  # .vcf, .vcf.gz or .bcf
        stats=outputDir+"results/{sample}/vep_stats/{sample}.html"
    params:
        plugins=config["vep_annotate"]["plugins"],
        extra=config["vep_annotate"]["extra"] # optional: extra arguments
    log:
        outputDir+"logs/vep/{sample}.vepanno.log"
    message: 
        "annotating with VEP for sample {wildcards.sample}"
    threads: 20
    wrapper:
        "v1.15.0/bio/vep/annotate"

## SnpSift VariantType ##
rule snpsift_vep_VarType:
    input:
        vcf=outputDir+"results/{sample}/{sample}.vep.vcf"
    output:
        vcf=temp(outputDir+"results/{sample}/{sample}.vep.final.vcf")
    message:
        "running snpsift vartype for sample {wildcards.sample}"
    params:
        java_opts=config["snpsift_VarType"]["params"]["java_opts"],
        extra=config["snpsift_VarType"]["params"]["extra"],
        snpsift_path=config["snpsift_path"]
    log:
        outputDir+"logs/snpsift/vartype/{sample}.log"
    shell:
        "java -jar {params.snpsift_path} VarType {input.vcf} > {output.vcf}"


## SnpSift ExtractFields - 1st parsing step ##
rule snpsift_vep_ExtractFields:
    input:
        outputDir+"results/{sample}/{sample}.vep.final.vcf"
    output:
        temp(outputDir+"results/{sample}/{sample}.vep.ssef.tsv")
    log:
        outputDir+"logs/snpsift/ExtractFields/{sample}.log"
    params:
        fieldsfile=config["vcf_fieldsfile"]
    message:
        "running snpsift extractfields (parsing vcf fields) fields for sample {wildcards.sample}"
    shell:
        "python3 workflow/scripts/extractfields.py --input {input} --fieldsfile {params.fieldsfile} > {output}"    

## VEP ExtractFields - 2nd parsing step ##
rule vep_extractFields:
    input:
        tsv=outputDir+"results/{sample}/{sample}.vep.ssef.tsv",
        vcf=outputDir+"results/{sample}/{sample}.vep.vcf"
    output:
        temp(outputDir+"results/{sample}/{sample}.vep.tsv")
    params:
        fieldsfile=config["vep_fieldsfile"]
    message:
        "running VEP extractfields (parsing CSQ fields) for sample {wildcards.sample}"
    shell:
        "python3 workflow/scripts/vep_extractfields.py --input {input.vcf} --fieldsfile {params.fieldsfile} --tsv {input.tsv}"


rule final_parsing_genelist:
    input:
        outputDir+"results/{sample}/{sample}.vep.tsv"
    output:
        outputDir+"results/{sample}/{sample}.xlsx"
    params:
        genelists = config["genelists"] if "genelists" in config and config["genelists"] else []
    message:
        "final operations for sample {wildcards.sample}: formatting columns, checking predictors and adding gene lists"
    shell:
        "python3 workflow/scripts/final_parsing.py -i {input} -g {params.genelists} -o {output}"
