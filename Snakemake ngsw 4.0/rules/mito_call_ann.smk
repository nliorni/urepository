# gatk Mutect2 \
#   			-R reference.fa \
#  			-L chrM \
#  			--mitochondria-mode \
#   			-I mitochondria.bam \
#   			-O mitochondria.vcf.gz


rule Mitocall_Mutect2:
    input:
        ref=config["reference"],
        bam="mitochondria.bam",
        chr="chrM"
    output:
        vcf="mitochondria.vcf.gz"
    shell:
        "gatk Mutect2 -R {input.ref} -L {input.chr} --mitochondria-mode -I {input.bam} -O {output.vcf}"