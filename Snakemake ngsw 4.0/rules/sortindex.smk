#SORTING WITH PICARD
rule picard_Sortsam:
    input:
        "recal/{sample}.bam"
    output:
        temp("sorted_reads/{sample}.sorted.bam")
    log:
        "logs/picard/sort_sam/{sample}.log"
    message:
        "Sorting reads in '{input}' into '{output}' using Picard Sortsam"
    params:
        sort_order="coordinate",
        extra="VALIDATION_STRINGENCY=LENIENT" # optional: Extra arguments for picard.
    wrapper:
        "0.66.0/bio/picard/sortsam"

#INDEXING WITH SAMTOOLS
rule samtools_Index:
    input:
        "sorted_reads/{sample}.sorted.bam"
    output:
        temp("sorted_reads/{sample}.sorted.bam.bai")
    message:
        "Indexing the '{input}' into '{output}' using Samtools Index"
    params:
        "" # optional params string
    wrapper:
        "0.66.0/bio/samtools/index"