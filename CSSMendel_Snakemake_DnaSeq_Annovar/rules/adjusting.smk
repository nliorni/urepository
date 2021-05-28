## Picard AddOrReplaceReadGroups ##
rule picard_AddOrReplaceGroups:
    input:
        "results/bam/{sample}.bam"
    output:
        temp("results/reheaded_bam/{sample}.bam")
    log:
        "logs/picard/replace_rg/{sample}.log"
    message:
        "Running Picard AddOrReplaceGroups. Reheading '{input}' into '{output}'."
    params:
        config["picard_AddOrReplaceGroups"]["params"]
    resources:
        mem_mb=8192
    wrapper:
        "0.67.0/bio/picard/addorreplacereadgroups"

## Picard ReorderSam ##
rule picard_ReorderSam:
    input:
        bam="results/reheaded_bam/{sample}.bam",
        ref=config["reference"]
    output:
        temp("results/reordered_bam/{sample}.bam")
    params:
        java_opts="",
        extra=""
    message:
        "Running Picard ReorderSam."
    shell:
        "picard ReorderSam {params.extra} -I {input.bam} -O {output} -SD {input.ref}"

## Samtools View  ##
rule samtools_View:
    input:
        "results/reordered_bam/{sample}.bam"
    output:
        temp("results/processed_bam/{sample}.bam")
    message:
        "Running Samtools View. Cleaning BAM file according to chromosomes in BED file."
    log:
        "logs/samtools/view/{sample}.log"
    params:
        config["samtools_View"]["params"]
    wrapper:
        "0.72.0/bio/samtools/view"

## Picard MarkDuplicates ##
rule picard_MarkDuplicates:
    input:
        "results/processed_bam/{sample}.bam"
    output:
        bam=temp("results/md_bam/{sample}.bam"),
        metrics="results/stats/{sample}.md.metrics.txt"
    log:
        "logs/picard/dedup/{sample}_dedup.log"
    message:
        "Running Picard MarkDuplicates. Marking [and removing] duplicates from '{input}' into '{output.bam}' and writing metrics into '{output.metrics}'."
    params:
        config["picard_MarkDuplicates"]["params"]
    resources:
        mem_mb=8192
    wrapper:
        "0.66.0/bio/picard/markduplicates"


## Gatk BaseRecalibrator ##
rule gatk_BaseRecalibrator:
    input:
        bam="results/md_bam/{sample}.bam",
        ref=config["reference"],
        dic=config["dict"],
        known=config["dbsnp"]  # optional known sites
    output:
        recal_table="results/recal_table/{sample}.grp"
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    message:
        "Running Gatk BaseRecalibrator. Recalibration of '{input.bam}' based on '{input.known}' creating '{output}'."
    params:
        extra=config["gatk_BaseRecalibrator"]["params"]["extra"],  
        java_opts=config["gatk_BaseRecalibrator"]["params"]["java_opts"]
    resources:
        mem_mb=8192 
    wrapper:
        "0.66.0/bio/gatk/baserecalibrator"

## Gatk ApplyBQSR ##
rule gatk_ApplyBQSR:
    input:
        bam="results/md_bam/{sample}.bam",
        ref=config["reference"],
        dictio=config["dict"],
        recal_table="results/recal_table/{sample}.grp"
    output:
        bam=temp("results/bqsr_bam/{sample}.bam")
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log"
    message:   
        "Running Gatk ApplyBQSR. Applying the '{input.recal_table}' to '{input.bam}' outputing '{output}'."
    params:
        extra=config["gatk_ApplyBQSR"]["params"]["extra"],  
        java_opts=config["gatk_ApplyBQSR"]["params"]["java_opts"]
    resources:
        mem_mb=8192
    wrapper:
        "0.66.0/bio/gatk/applybqsr"


