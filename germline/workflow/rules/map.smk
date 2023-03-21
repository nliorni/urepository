rule bwa_mem_sort:
    input:
        reads = lambda wildcards: expand(f"{config['samples'][wildcards.sample]}_{{unit}}.{{extension}}", unit=config["unit"], extension=config["extension"])
    output:
        bam = temp(outputDir+"results/{sample}/{sample}.bam")
    params:
        reference = config["reference"],
        sort_criteria = " ", # -n if sort by name
        threads = 20
    message:
        "running bwa mem for sample {wildcards.sample}"
    shell:
        "bwa mem -t {params.threads} {params.reference} {input.reads} | samtools sort {params.sort_criteria} -o {output.bam}"


# rule samtools_flagstat:
#     input:
#         bam = outputDir+"results/{sample}/{sample}.bam"
#     output:
#         stats = outputDir+"results/{sample}/qc/flagstat_{sample}.txt"
#     params:
#         ""
#     message:
#         "retrieving stats from the BWA-aligned BAM for {wildcards.sample}"
#     shell:
#         "samtools flagstat {input.bam} > {output.stats}"