rule fastqc:
    input:
        reads = lambda wildcards: expand(f"{config['samples'][wildcards.sample]}_{{unit}}.{{extension}}", unit=config["unit"], extension=config["extension"])
    output:
        html=outputDir+"results/{sample}/qc/{sample}_{unit}_fastqc.html",
        zip = outputDir+"results/{sample}/qc/{sample}_{unit}_fastqc.zip"
    message:
        "running fastqc for sample {wildcards.sample}"
    params:
        outdir = outputDir+"results/{sample}/qc/",
        threads = 20
    shell:
        "fastqc -t {params.threads} --outdir {params.outdir} {input.reads}"

rule multiqc:
    input:
        expand([outputDir+"results/{sample}/{sample}.xlsx" ,outputDir + "results/{sample}/qc/{sample}_{unit}_fastqc.zip", outputDir + "results/{sample}/vep_stats/{sample}.html", outputDir + "results/{sample}/qc/{sample}.md.metrics.txt", outputDir + "results/{sample}/qc/{sample}_meqc_stats.txt", outputDir+"results/{sample}/qc/{sample}_iCoverage_output_file.csv", outputDir+"results/{sample}/qc/flagstat_{sample}.txt"],  sample = config["samples"], unit = config["unit"])
    output:
        outputDir + "results/multiqc/multiqc_report.html"
    message:
        "running multiqc"
    # conda: 
    #     "multiqc_less"
    params:
        extra = "",
        outdir = outputDir+"results/multiqc/"
    shell:
        "multiqc {params.extra} --force -o {params.outdir} {input}"
