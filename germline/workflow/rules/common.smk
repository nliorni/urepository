def get_input_bam():
    workflow = config["workflow"]
    if workflow=="StartWithBam":
        inputCmd=lambda wildcards: expand(f"{config['samples'][wildcards.sample]}.bam")
    else:
        inputCmd=outputDir+"results/{sample}/{sample}.bam"
    return inputCmd

def get_input_vcf():
    workflow = config["workflow"]
    if workflow=="StartWithVcf":
        inputCmd=lambda wildcards: expand(f"{config['samples'][wildcards.sample]}.vcf")
    else:
        inputCmd=outputDir+"results/{sample}/{sample}.filter.vcf"
    return inputCmd


