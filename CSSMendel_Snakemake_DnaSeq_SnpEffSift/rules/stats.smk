## Samtools Depth ##
rule samtools_Depth:
    input:
        bams="results/bqsr_bam/{sample}.bam",
        bed=config["bed"], # optional
    output:
        "results/stats/{sample}_samtools_stats.csv"
    log:
        "logs/samtools/depth/{sample}.log"
    message:
        "Running Samtools Depth. Creating '{output}' from '{input.bams}' and '{input.bed}'."
    params:
        # optional bed file passed to -b
        extra="" # optional additional parameters as string
    wrapper:
        "0.67.0/bio/samtools/depth"


## Bedtools Coverage ##
rule bedtools_Coverage:
   input:
       bed=config["bed"],
       ref=config["filtered"],
       bam="results/bqsr_bam/{sample}.bam"
   output:
       "results/stats/{sample}_bedtools_stats.csv"
   log:
       "logs/bedtools/coverage/{sample}.log"
   message:
        "Running Bedtools Coverage. Creating '{output}' from '{input.bed}', '{input.ref}' and '{input.bam}'."
   shell:
       "bedtools coverage -sorted -d -a {input.bed} -g {input.ref} -b {input.bam} > {output}"

## MeQC ##
rule Meqc:
    input:
        samstats="results/stats/{sample}_samtools_stats.csv",
        bedstats="results/stats/{sample}_bedtools_stats.csv"
    log:
        "logs/meqc/{sample}.log"
    output:
        "results/stats/{sample}_meqc_stats.txt"
    message:
        "Running MeQC. Creating '{output}' from '{input.samstats}' and '{input.bedstats}'."
    shell:
        "scripts/MEQC.R 1,2,3,5,10,20,30,50,100,500,1000,5000 {output} {input.samstats} {input.bedstats}"

## Bedtools CoverageBed ##
rule bedtools_CoverageBed:
    input:
        a=config["bed"],
        b="results/bqsr_bam/{sample}.bam"
    output:
        "results/stats/{sample}.cov"
    log:
        "logs/coveragebed/{sample}.log"
    params:
        extra=""  # optional parameters
    threads: 8
    message:
        "Running Bedtools CoverageBed. Creating '{output}', from '{input.a}' and '{input.b}' with."
    wrapper:
        "0.67.0/bio/bedtools/coveragebed"

## iCoverage ##
rule iCoverage:
    input: 
        "results/stats/{sample}.cov"
    output:
        "results/stats/{sample}_iCoverage_output_file.xlsx"
    log:
        "logs/icoverage/{sample}.log"
    message:
        "Running iCoverage. Creating '{output}' from '{input}'."
    shell:
        "python3 scripts/icoverage.py -c {input} -t 30 -o {output}"
