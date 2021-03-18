# #DEEPVARIANT
rule Deepvariant:
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
       extra="--regions chr1:43389920-43450402"
   log:
       "logs/deepvariant/{sample}/stdout.log"
   message:
        "Trying out DeepVariant: deepcalling {input.bam} -based on {input.ref}- into {output.vcf}"
   wrapper:
       "0.68.0/bio/deepvariant"


#PENSA A BCFTOOLS CONCAT AL POSTO DI GLNEXUS!!!!! FORSE FUNZIONA!   https://github.com/dnanexus-rnd/DeepVariant-GLnexus-WDL



#al posto di picard rename samples è possibile e consigliato usare bcftools reheader
#OPPURE you can also use the --sample_name flag when running DeepVariant to manually specify the sample name

# Note: The BAM files should provide unique names for each sample in their `SM`
# header tag, which is usually derived from a command-line flag to the read
# aligner. If your BAM files don't have unique `SM` tags (and if it's not feasible
# to adjust the alignment pipeline), add the `--sample_name=XYZ` flag to
# `run_deepvariant` to override the sample name written into the gVCF file header.



#PICARD RENAME SAMPLES   #bisogna trovare un modo per chiamare differentmente i due campioni
rule picard_RenameSampleInVcf:
    input:
        sample="deepcalls/{sample}.g.vcf"
    output:
        "deepnamed/{sample}_renamed.g.vcf"
    message:
        "Renaming the vcf file {input.sample} to make it compatible with glnexus"
    shell:
        "picard RenameSampleInVcf INPUT={input.sample} NEW_SAMPLE_NAME={wildcard.sample} OUTPUT={output}"


#OPPURE

# rule bcftools_reheader:
#     input:
#         vcf="a.bcf",
#         ## new header, can be omitted if "samples" is set
#         header="header.txt",
#         ## file containing new sample names, can be omitted if "header" is set
#         samples="samples.tsv"
#     output:
#         "a.reheader.bcf"
#     params:
#         extra="",  # optional parameters for bcftools reheader
#         view_extra="-O b"  # add output format for internal bcftools view call
#     wrapper:
#         "0.72.0/bio/bcftools/reheader"


#GLNEXUS
# rule combine_GlNexus:
#     input:
#         vcf=expand("deepnamed/{sample}_renamed.vcf", sample=config["samples"]),
#         ref=config["reference"]
#     output:
#         "deepcombined/all.vcf"
#     log:
#         "logs/glnexus/allg.log"
#     message:
#         "Combining gvcfs {input.vcf} into {output} obtained from a DeepVariant calling"
#     shell:
#         "scripts/GLnexus/glnexus_cli {input.vcf} > {output}"
