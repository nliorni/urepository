#SAMTOOLS DEPTH
rule samtools_depth:
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

#BEDTOOLS COVERAGE
rule bedtools_coverage:
   input:
       bed=config["bed"],
       ref=config["filtered"],
       bam="recal/{sample}.bam"
   output:
       "stats/{sample}_bedtools_stats.csv"
   log:
       "logs/bedtools/coverage/{sample}.log"
   message:
        "Creating '{output}' with Bedtools Coverage, from '{input.bed}', '{input.ref}' and '{input.bam}'"
   shell:
       "bedtools coverage -sorted -d -a {input.bed} -g {input.ref} -b {input.bam} > {output}"

#MEQC
rule meqc:
    input:
        samstats="stats/{sample}_samtools_stats.csv",
        bedstats="stats/{sample}_bedtools_stats.csv"
    log:
        "logs/meqc/{sample}.log"
    output:
        "meqc/{sample}_meqc_stats.txt"
    shell:
        "scripts/MEQC.R 1,2,3,5,10,20,30,50,100,500,1000,5000 {output} {input.samstats} {input.bedstats}"

#COVERAGEBED
rule coverageBed:
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
        "Trying out coverageBed, creating '{output}', from '{input.a}' and '{input.b}'"
    wrapper:
        "0.67.0/bio/bedtools/coveragebed"

#iCOVERAGE
#rule icoverage:
#    input: 
#        "stats/{sample}.cov"
#    output:
#        "icov/{sample}_output_file.xlsx"
#    log:
#        "logs/icoverage/{sample}.log"
#    message:
#        "Creating {output} from {input} using iCoverage script"
#    shell:
#        "/scripts/icoverage.py -c {input} -t 30 -o {output}"
