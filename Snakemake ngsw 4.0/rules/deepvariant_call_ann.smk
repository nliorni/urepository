#DEEPVARIANT
rule deepvariant_calling:
   input:
       bam="sorted_reads/{sample}.sorted.bam",
       bai="sorted_reads/{sample}.sorted.bam.bai",
       ref=config["reference"],
       dic=config["dict"],
       regions=config["bed"]
   output:
       vcf="deepcalls/{sample}.g.vcf"    #try with g.vcf
   params:
       model="wgs",   # {wgs, wes}
       extra="--regions chr1:43389920-43450402" #prova --regions {regions}, regions=config["regions"]
   log:
       "logs/deepvariant/{sample}/stdout.log"
   message:
        "Trying out DeepVariant: deepcalling {input.bam} -based on {input.ref}- into {output.vcf}"
   wrapper:
       "0.68.0/bio/deepvariant"


#PICARD RENAME SAMPLES   
rule picard_RenameSampleInVcf:
    input:
        sample="deepcalls/{sample}.g.vcf"
    output:
        "deepnamed/{sample}_renamed.g.vcf"
    message:
        "Renaming the vcf file {input.sample} to make it compatible with glnexus"
    shell:
        "picard RenameSampleInVcf INPUT={input.sample} NEW_SAMPLE_NAME={wildcards.sample} OUTPUT={output}"

#GLNEXUS
rule combine_GLNexus:
    input:
        gvcfs=expand("deepnamed/{sample}_renamed.g.vcf", sample=config["samples"])
    output:
        "gl_combined/all.bcf"
    message: 
        "Combining the gvcf files {input.gvcfs} obtained from DeepVariant calling into {output} using GLNexus"
    shell:
        "scripts/glnexus_cli --config DeepVariant {input.gvcfs} > {output}"

#BCFTOOLS CONVERTING 
rule bcftools_View:
    input:
        "gl_combined/all.bcf"
    output:
        "gl_combined/all.vcf"
    message:
        "Converting {input} into {output} using bcftools View"
    shell:
        "bcftools view {input} > {output} "

#ANNOTATE VCF WITH SNPEFF
rule dv_snpeff_Annotate:
    input:
        calls="gl_combined/all.vcf", # (vcf, bcf, or vcf.gz)
        db=config["snpeff"] # path to reference db downloaded with the snpeff download wrapper
    output:
        #multiext("snpeff/{sample}", ".vcf", ".html", ".csv")
        calls="dv_annotated/all.vcf",   # annotated calls (vcf, bcf, or vcf.gz)
        stats="dv_annotated/all.html",  # summary statistics (in HTML), optional
        csvstats="dv_annotated/all.csv" # summary statistics in CSV, optional
    log:
        "logs/snpeff/all.log"
    message:
        "Annotating '{input.calls}' with '{input.db}' to generate '{output.calls}', '{output.stats}' and '{output.csvstats}' with SnpEff Annotate"
    params:
        extra="-Xmx8g"           # optional parameters (e.g., max memory 4g)
    wrapper:
        "0.66.0/bio/snpeff/annotate"

#SNPSIFT ANNOTATE
rule dv_snpsift_Annotate:
    input:
        call="dv_annotated/all.vcf",
        database=config["dbsnp"]
    output:
        call=temp("dv_sift_annotated/all.vcf")
    message:
        "Further annotating '{input.call}' using '{input.database}' creating '{output.call}' with SnpSift Annotate"
    log:
        "logs/snpsift/annotate/all.log"
    wrapper:
        "0.67.0/bio/snpsift/annotate"


#SNPSIFT VARTYPE
rule dv_snpsift_VarType:
    input:
        vcf="dv_sift_annotated/all.vcf"
    output:
        vcf=temp("dv_sift_vartyped/all.vcf")
    message:
        "Further annotating '{input.vcf}' creating '{output.vcf}' with SnpSift VarType"
    log:
        "logs/snpsift/vartype/all.log"
    wrapper:
        "0.67.0/bio/snpsift/varType"

#dbNSFP
rule dv_snpsift_dbNSFP:
   input:
       call = "dv_sift_vartyped/all.vcf",
       dbNSFP = config["dbnsfp"]
   output:
       call = "dv_sift_dbNSFP/all.vcf"
   message:
       "Further annotating '{input.call}' using '{input.dbNSFP}', creating '{output.call}' with SnpSift dbNSFP"
   log:
       "logs/dbNSFP/all.log"
   wrapper:
       "0.67.0/bio/snpsift/dbnsfp"

#EXTRACT FIELDS
rule dv_snpsift_ExtractFields:
    input:
        "dv_sift_dbNSFP/all.vcf"
    output:
        "dv_csvfile/all.csv"
    message:
        "Extracting fields of interest from the completly annotated vcf file {input} into {output}"
    shell:
        "python3 scripts/extractfields.py --input {input} > {output}"

