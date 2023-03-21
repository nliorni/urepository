## Samtools Depth ##
rule samtools_depth:
    input:
        bams=outputDir+"results/{sample}/{sample}.bqsr.bam",
        bed=config["bed"], # optional
    output:
        temp(outputDir+"results/{sample}/qc/{sample}_samtools_stats.csv")
    message:
        "running samtools depth for sample {wildcards.sample}"
    params:
        extra=config["samtools_Depth"]["extra"] # optional additional parameters as string
    shell:
        "samtools depth -a {input.bams} -b {input.bed} > {output}"

## Bedtools Coverage ##
rule bedtools_coverage:
   input:
       bed=config["bed"],
       ref=config["genomeFile"],
       bam=outputDir+"results/{sample}/{sample}.bqsr.bam"
   output:
       temp(outputDir+"results/{sample}/qc/{sample}_bedtools_stats.csv")
   message:
        "running bedtools coverage for sample {wildcards.sample}"
   shell:
       "bedtools coverage -sorted -d -a {input.bed} -g {input.ref} -b {input.bam} > {output}"

    
## meQC ##
rule meQC:
    input:
        samstats=outputDir+"results/{sample}/qc/{sample}_samtools_stats.csv",
        bedstats=outputDir+"results/{sample}/qc/{sample}_bedtools_stats.csv"
    output:
        outputDir+"results/{sample}/qc/{sample}_meqc_stats.txt"
    params:
        ranges=config["MeQC"]["ranges"]
    message:
        "running meqc for sample {wildcards.sample}"
    shell:
        "workflow/scripts/MEQC.R {params.ranges} {output} {input.samstats} {input.bedstats}"

## Bedtools CoverageBed ##
rule bedtools_coveragebed:
    input:
        bed=config["bed"],
        bam=outputDir+"results/{sample}/{sample}.bqsr.bam",
	    genomeFile=config["genomeFile"]
    output:
        temp(outputDir+"results/{sample}/qc/{sample}.cov")
    params:
        extra=config["bedtools_CoverageBed"]["extra"]  
    message:
        "running bedtools coveragebed for sample {wildcards.sample}"
    shell:
        "/software/bedtools_2.27/bin/coverageBed -a {input.bed} -b {input.bam} -d -sorted -g {input.genomeFile} > {output}"

## iCoverage ##
rule iCoverage:
    input: 
        outputDir+"results/{sample}/qc/{sample}.cov"
    output:
        outputDir+"results/{sample}/qc/{sample}_iCoverage_output_file.xlsx"
    params: 
        threshold=config["iCoverage"]["threshold"]
    message:
        "running icoverage for sample {wildcards.sample}"
    shell:
        "python3 workflow/scripts/icoverage.py -c {input} -t {params.threshold} -o {output}"

rule iCoverage_csv:
    input:
        outputDir+"results/{sample}/qc/{sample}_iCoverage_output_file.xlsx"
    output:
        temp(outputDir+"results/{sample}/qc/{sample}_iCoverage_output_file.csv")
    params:
        ""
    message:
        "get the icoverage csv for {wildcards.sample}"
    shell:
        "python3 workflow/scripts/to_csv.py --input {input} --output {output}"