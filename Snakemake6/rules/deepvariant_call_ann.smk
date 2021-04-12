## Google Deepvariant ##
rule deepvariant_calling:
   input:
       bam="results/sorted_bam/{sample}.sorted.bam",
       bai="results/sorted_bam/{sample}.sorted.bam.bai",
       ref=config["reference"],
       dic=config["dict"],
       regions=config["bed"]
   output:
       vcf="results/deepcalls/{sample}.g.vcf"    
   params:
       model=config["deepvariant"]["model"],  
       extra=config["deepvariant"]["regions"]  
   log:
       "logs/deepvariant/{sample}/stdout.log"
   message:
        "Running Google DeepVariant. Deepcalling '{input.bam}' -based on '{input.ref}' - into '{output.vcf}'"
   wrapper:
       "0.68.0/bio/deepvariant"


## Picard RenameSampleInVcf ##  
rule picard_RenameSampleInVcf:
    input:
        sample="results/deepcalls/{sample}.g.vcf"
    output:
        "results/deepnamed/{sample}_renamed.g.vcf"
    message:
        "Running Picard RenameSampleInVcf. Renaming the vcf file '{input.sample}' to make it compatible with Glnexus"
    shell:
        "picard RenameSampleInVcf INPUT={input.sample} NEW_SAMPLE_NAME={wildcards.sample} OUTPUT={output}"

## GLNexus ##
rule combine_GLNexus:
    input:
        gvcfs=expand("results/deepnamed/{sample}_renamed.g.vcf", sample=config["samples"])
    output:
        "results/gl_combined/all.bcf"
    message: 
        "Running GLNexus. Combining the gvcf files '{input.gvcfs}' obtained from DeepVariant calling into '{output}' using GLNexus"
    shell:
        "scripts/glnexus_cli --config DeepVariant '{input.gvcfs}' > '{output}'"

## Bcftools View ##
rule bcftools_View:
    input:
        "results/gl_combined/all.bcf"
    output:
        "results/vcf/all.vcf"
    message:
        "Running bcftools View. Converting '{input}' into '{output}'."
    shell:
        "bcftools view {input} > {output} "


## SnpEff Annotate ##
rule dv_snpeff_Annotate:
    input:
        calls="results/vcf/all.vcf", # (vcf, bcf, or vcf.gz)
        db=config["snpeff"] # path to reference db downloaded with the snpeff download wrapper
    output:
        #multiext("snpeff/{sample}", ".vcf", ".html", ".csv")
        calls="results/dv_annotated/all.vcf",   # annotated calls (vcf, bcf, or vcf.gz)
        stats="results/dv_annotated/all.html",  # summary statistics (in HTML), optional
        csvstats="results/dv_annotated/all.csv" # summary statistics in CSV, optional
    log:
        "logs/snpeff/all.log"
    message:
        "Running SnpEff Annotate. Annotating '{input.calls}' with '{input.db}' to generate '{output.calls}', '{output.stats}' and '{output.csvstats}'."   
    wrapper:
        "0.66.0/bio/snpeff/annotate"

## SnpSift Annotate ##
rule dv_snpsift_Annotate:
    input:
        call="results/dv_annotated/all.vcf",
        database=config["dbsnp"]
    output:
        call=temp("results/dv_sift_annotated/all.vcf")
    message:
        "Running SnpSift Annotate. Further annotating '{input.call}' using '{input.database}' creating '{output.call}'."
    log:
        "logs/snpsift/annotate/all.log"
    wrapper:
        "0.67.0/bio/snpsift/annotate"


## SnpSift VariantType ##
rule dv_snpsift_VarType:
    input:
        vcf="results/dv_sift_annotated/all.vcf"
    output:
        vcf=temp("results/dv_sift_vartyped/all.vcf")
    message:
        "Running SnpSift VarType. Further annotating '{input.vcf}' creating '{output.vcf}'."
    log:
        "logs/snpsift/vartype/all.log"
    wrapper:
        "0.67.0/bio/snpsift/varType"

## SnpSift dbNSFP ##
rule dv_snpsift_dbNSFP:
   input:
       call = "results/dv_sift_vartyped/all.vcf",
       dbNSFP = config["dbnsfp"]
   output:
       call = "results/dv_sift_dbNSFP/all.vcf"
   message:
       "Running SnpSift dbNSFP. Further annotating '{input.call}' using '{input.dbNSFP}', creating '{output.call}'."
   log:
       "logs/dbNSFP/all.log"
   wrapper:
       "0.67.0/bio/snpsift/dbnsfp"

## SnpSift ExtractFields ##
rule dv_snpsift_ExtractFields:
    input:
        "results/dv_sift_dbNSFP/all.vcf"
    output:
        "results/dv_csvfile/all.csv"
    message:
        "Running SnpSift ExtractFields. Extracting fields of interest from the completly annotated vcf file '{input}' into '{output}'"
    shell:
        "python3 scripts/extractfields.py --input {input} > {output}"
