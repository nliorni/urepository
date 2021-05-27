## Bwa Mem ##
rule bwa_Mem:
    input:
        reads=lambda wildcards: expand(f"{config['samples'][wildcards.sample]}_{{unit}}.fastq.gz", unit=["L001_R1_001","L001_R2_001"]) 
    output:
        "results/bam/{sample}.bam"
    params:
        index=config["reference"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort=config["bwaMem"]["params"]["sort"],             # Can be 'none', 'samtools' or 'picard'.
        sort_order=config["bwaMem"]["params"]["sort_order"],  # Can be 'queryname' or 'coordinate'.
        sort_extra=config["bwaMem"]["params"]["sort_extra"]            # Extra args for samtools/picard.
    message:
        "Running BWA mem. Mapping '{input}' on reference genome '{params.index}' producing {output}."
    threads: 20
    log:
        "logs/bwa/map/{sample}.log"
    wrapper:
        "0.66.0/bio/bwa/mem"