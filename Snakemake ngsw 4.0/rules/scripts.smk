#SAMTOOLS DEPTH
rule samtools_Depth:
    input:
        bams="recal/{sample}.bam",
        bed=config["bed"], # optional
    output:
        "stats/{sample}_samtools_stats.csv"
    log:
        "logs/samtools/depth/{sample}.log"
    message:
        "Creating '{output}' from '{input.bams}' and '{input.bed}' with Samtools Depth"
    params:
        # optional bed file passed to -b
        extra="" # optional additional parameters as string
    wrapper:
        "0.67.0/bio/samtools/depth"

#MULTIQC
#rule multiqc:
#    input:
#        expand("samtools_stats/{sample}.txt", sample=["a", "b"])
#    output:
#        "qc/multiqc.html"
#    params:
#        ""  # Optional: extra parameters for multiqc.
#    log:
#        "logs/multiqc.log"
#    wrapper:
#        "0.67.0/bio/multiqc"


#BEDTOOLS COVERAGE
rule bedtools_Coverage:
   input:
       bed=config["bed"],
       ref=config["filtered"],
       bam="recal/{sample}.bam"
   output:
       "stats/{sample}_bedtools_stats.csv"
   log:
       "logs/bedtools/coverage/{sample}.log"
   message:
        "Creating '{output}' from '{input.bed}', '{input.ref}' and '{input.bam}' with Bedtools Coverage"
   shell:
       "bedtools coverage -sorted -d -a {input.bed} -g {input.ref} -b {input.bam} > {output}"

#MEQC
rule Meqc:
    input:
        samstats="stats/{sample}_samtools_stats.csv",
        bedstats="stats/{sample}_bedtools_stats.csv"
    log:
        "logs/meqc/{sample}.log"
    output:
        "stats/{sample}_meqc_stats.txt"
    message:
        "Creating '{output}' from '{input.samstats}' and '{input.bedstats}' with the homemade script MeQC "
    shell:
        "scripts/MEQC.R 1,2,3,5,10,20,30,50,100,500,1000,5000 {output} {input.samstats} {input.bedstats}"

#COVERAGEBED
rule bedtools_CoverageBed:
    input:
        a=config["bed"],
        b="recal/{sample}.bam"
    output:
        "stats/{sample}.cov"
    log:
        "logs/coveragebed/{sample}.log"
    params:
        extra=""  # optional parameters
    threads: 8
    message:
        "Creating '{output}', from '{input.a}' and '{input.b}' with Bedtools CoverageBed"
    wrapper:
        "0.67.0/bio/bedtools/coveragebed"

#iCOVERAGE
rule iCoverage:
    input: 
        "stats/{sample}.cov"
    output:
        "stats/{sample}_iCoverage_output_file.xlsx"
    log:
        "logs/icoverage/{sample}.log"
    message:
        "Creating {output} from {input} with iCoverage homemade script"
    shell:
        "python3 scripts/icoverage.py -c {input} -t 30 -o {output}"
