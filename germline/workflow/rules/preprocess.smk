## Picard AddOrReplaceReadGroups ##
rule picard_addorreplacereadgroups:
    input:
        get_input_bam()
    output:
        temp(outputDir+"results/{sample}/{sample}.aorg.bam")
    message:
        "running picard addorreplacereadgroups for sample {wildcards.sample}"
    params:
        RGLB = "lib1",
        RGPL = "ILLUMINA",
        RGPU = "unit1",
        RGSM = lambda wildcards: f"{[wildcards.sample]}" # change to wildcards.sample (otherwise cohort mode won't work)
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} RGLB={params.RGLB} RGPL={params.RGPL} RGPU={params.RGPU} RGSM={params.RGSM}"

## Samtools View  ##
rule samtools_view:
    input:
        bam=outputDir+"results/{sample}/{sample}.aorg.bam",
        bed=config["bed"]
    output:
        temp(outputDir+"results/{sample}/{sample}.sview.bam")
    message:
        "running samtools view for sample {wildcards.sample}"
    shell:
        "samtools view -b -h -L {input.bed} {input.bam} > {output}"

## Picard MarkDuplicates ##
rule picard_markduplicates:
    input:
        outputDir+"results/{sample}/{sample}.sview.bam"
    output:
        bam=temp(outputDir+"results/{sample}/{sample}.dedup.bam"),
        metrics=outputDir+"results/{sample}/qc/{sample}.md.metrics.txt"
    message:
        "running picard markduplicates for sample {wildcards.sample}"
    params:
        config["picard_MarkDuplicates"]["params"]
    resources:
        mem_mb=config["picard_MarkDuplicates"]["mem_mb"]
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.metrics} {params}"

## Gatk BaseRecalibrator ##
rule gatk_baserecalibrator:
    input:
        bam=outputDir+"results/{sample}/{sample}.dedup.bam",
        ref=config["reference"],
        dic=config["dict"],
        known=config["dbsnp"]  # optional known sites
    output:
        recal_table=temp(outputDir+"results/{sample}/{sample}.grp")
    message:
        "running gatk baserecalibrator for sample {wildcards.sample}"
    params:
        extra=config["gatk_BaseRecalibrator"]["params"]["extra"],  
        java_opts=config["gatk_BaseRecalibrator"]["params"]["java_opts"]
    resources:
        mem_mb=config["gatk_BaseRecalibrator"]["mem_mb"] 
    shell:
        "gatk BaseRecalibrator --input {input.bam} --reference {input.ref} --known-sites {input.known} --output {output.recal_table}"

## Gatk ApplyBQSR ## 
rule gatk_applyBQSR:
    input:
        bam=outputDir+"results/{sample}/{sample}.dedup.bam",
        ref=config["reference"],
        dictio=config["dict"],
        recal_table=outputDir+"results/{sample}/{sample}.grp"
    output:
        bam=outputDir+"results/{sample}/{sample}.bqsr.bam",
        bai=temp(outputDir+"results/{sample}/{sample}.bqsr.bai")
    message:   
        "running gatk applybqsr for sample {wildcards.sample}"
    params:
        extra=config["gatk_ApplyBQSR"]["params"]["extra"],  
        java_opts=config["gatk_ApplyBQSR"]["params"]["java_opts"]
    resources:
        mem_mb=config["gatk_ApplyBQSR"]["mem_mb"]
    shell:
        "gatk ApplyBQSR --input {input.bam} --bqsr-recal-file {input.recal_table} --reference {input.ref} {params.extra} --output {output.bam}"

rule samtools_flagstat:
    input:
        bam = outputDir+"results/{sample}/{sample}.bqsr.bam",
        bai = outputDir+"results/{sample}/{sample}.bqsr.bai"
    output:
        stats = outputDir+"results/{sample}/qc/flagstat_{sample}.txt"
    params:
        ""
    message:
        "retrieving stats from the BWA-aligned BAM for {wildcards.sample}"
    shell:
        "samtools flagstat {input.bam} > {output.stats}"

## Samtools Index ##
rule samtools_index:
    input:
        outputDir+"results/{sample}/{sample}.bqsr.bam"
    output:
        outputDir+"results/{sample}/{sample}.bqsr.bam.bai",
    message:
        "running samtools index for sample {wildcards.sample}"
    params:
        config["samtools_Index"]["params"] 
    shell:
        "samtools index {params} {input}"