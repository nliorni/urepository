#MAPPING WITH BWA
rule bwa_map:
    input:
        reads=["data/samples/{sample}_1.fastq", "data/samples/{sample}_2.fastq"]
    output:
        "mapped_reads/{sample}.bam"
    params:
        index="data/genome.fa",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    log:
        "logs/bwa/map/{sample}.log"
    wrapper:
        "0.66.0/bio/bwa/mem"