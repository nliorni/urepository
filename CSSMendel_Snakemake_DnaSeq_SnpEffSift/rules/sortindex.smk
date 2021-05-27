## Picard SortSam ##
rule picard_Sortsam:
    input:
        "results/bqsr_bam/{sample}.bam"
    output:
        "results/sorted_bam/{sample}.sorted.bam"
    log:
        "logs/picard/sort_sam/{sample}.log"
    message:
        "Sorting reads in '{input}' into '{output}' using Picard Sortsam"
    params:
        sort_order=config["picard_Sortsam"]["params"]["sort_order"],
        extra=config["picard_Sortsam"]["params"]["extra"]
    wrapper:
        "0.66.0/bio/picard/sortsam"
# considera la rimozione di questa regola: bwa mem usa subito anche samtools sort. E' davvrero necessario? #

## Samtools Index ##
rule samtools_Index:
    input:
        "results/sorted_bam/{sample}.sorted.bam"
    output:
        temp("results/sorted_bam/{sample}.sorted.bam.bai")
    log:
        "logs/samtools/index/{sample}.log"
    message:
        "Indexing the '{input}' into '{output}' using Samtools Index"
    params:
        config["samtools_Index"]["params"] 
    wrapper:
        "0.66.0/bio/samtools/index"