#MAPPING WITH BWA
rule bwa_map:
    input:
        reads=["data/samples/{sample}_1.fastq", "data/samples/{sample}_2.fastq"]
    output:
        "mapped_reads/{sample}.bam"
    params:
        index=config["reference"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    message:
        "Mapping '{input}' on reference genome with BWA mem"
    log:
        "logs/bwa/map/{sample}.log"
    threads: 4
    wrapper:
        "0.66.0/bio/bwa/mem"
