## Optional outputs ##
fastqc_output=[]
if config["fastqc"]==True:
    fastqc_output.append("results/qc/fastqc/{sample}_{num}.html")

meqc_output=[]
if config["meqc"]==True:
    meqc_output.append("results/stats/{sample}_meqc_stats.txt")

icov_output=[]
if config["icov"]==True:
    icov_output.append("results/stats/{sample}_iCoverage_output_file.xlsx")


## Target rules ##  ## Workflow options ##
rule cohort:
    input:
        expand(fastqc_output, sample=config["samples"], num=['1', '2']), #cambia con unit=["R1", "R2"]
        expand(meqc_output, sample=config["samples"]),
        expand(icov_output, sample=config["samples"]),
        "results/csvfile/all.csv"
    message:
        "Your gatk workflow annotating variant in multiple combined samples is complete"


rule single:
    input:
        expand(fastqc_output, sample=config["samples"], num=['1', '2']),
        expand(meqc_output, sample=config["samples"]),
        expand(icov_output, sample=config["samples"]),
        expand("results/csvfile/{sample}.csv", sample=config["samples"])
    message:
        "Your gatk workflow annotating variant on samples independently is complete"
        



## Modules ##
include: "rules/mapping.smk"
include: "rules/qc.smk"
include: "rules/adjusting.smk"
include: "rules/sortindex.smk"
include: "rules/stats.smk"
include: "rules/gatk_call_ann.smk"
include: "rules/single_gatk_call_ann.smk"




