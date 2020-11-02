#bam quality control with KING algorithm??
###introdurre l'opzione -L per il file BED
rule haplotype_caller:
    input:
        # single or list of bam files
        #bed=data/bed/placeholder.bed
        bam="sorted_reads/{sample}.sorted.bam",
        bai="sorted_reads/{sample}.sorted.bam.bai",
        ref="data/genome.fa",
        dic="data/genome.dict"
        # known="dbsnp.vcf"  # optional
    output:
        gvcf="calls/{sample}.g.vcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra="",  # optional
        java_opts="-Xmx8G", # optional #-L {input.bed}?????
    wrapper:
        "0.66.0/bio/gatk/haplotypecaller"


#should I merge the vcfs from sample A and sample B? individual VCF files can be merged using BCFtools or
#or similar package
#should I do that now or after the snpeff?
#should I filter this to remove artifacts?? (Integrative Genome Viewer?)

#is this the rule for converting correctly  a g.vcf file to a vcf???
#rule genotype_gvcfs:
#    input:
#        gvcf="calls/{sample}.g.vcf",  # combined gvcf over multiple samples calls/all.g.vcf
#        ref="genome.fa"
#    output:
#        vcf="calls/{sample}.vcf",
#    log:
#        "logs/gatk/genotypegvcfs.log"
#    params:
#        extra="",  # optional
#        java_opts="", # optional
#    wrapper:
#        "0.67.0/bio/gatk/genotypegvcfs"
