rule picard_sortsam:
    input:
        "recal/{sample}.bam"
    output:
        "sorted_reads/{sample}.sorted.bam"
    log:
        "logs/picard/sort_sam/{sample}.log"
    params:
        sort_order="coordinate",
        extra="VALIDATION_STRINGENCY=LENIENT" # optional: Extra arguments for picard.
    wrapper:
        "0.66.0/bio/picard/sortsam"

rule samtools_index:
    input:
        "sorted_reads/{sample}.sorted.bam"
    output:
        "sorted_reads/{sample}.sorted.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.66.0/bio/samtools/index"