rule gatk_Mutect2:
    input:
        fasta=config["reference"],
        map="sorted_reads/{sample}.sorted.bam"
    output:
        vcf = "soma_calls/{sample}.vcf"
    message:
        "Testing Mutect2 with {wildcards.sample}"
    resources:
        mem_mb=1024
    log:
        "logs/mutect_{sample}.log"
    wrapper:
         "v0.69.0/bio/gatk/mutect"


#add mitochondrial java option